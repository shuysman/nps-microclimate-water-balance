#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
May 2024
Modified Water Balance Model to run on 1m LiDAR data
Proof of Concept for downscaled Water Balance runs at smaller scales.  Used to examin planting microsites for
Summer 2024 Whitebark Pine planting with Nancy Bockino in Shoshone NF

December 2019
Gridded Cloud Water Balance Model Version 1.5: Retroactively removing Senay et al. NDVI method of calculating AET BUT retaining Jennings et al. 2018 T50 coefficients for estimating snow.
This version also removes the IGRID veg layer precip correction. The relevant lines have just been commented out. 
--Modified from the "season test" code of November 2019.
This version reads the reprojected tif files from the MACA GCMs instead of daymet netCDF files.
"""
import datetime
import itertools
import multiprocessing
import glob
import os
import pickle
import sys
import time
from pathlib import Path
from statistics import mean

import netCDF4
import numpy as np
import pandas as pd
import xarray as xr
import rioxarray
import utm
from osgeo import gdal

np.seterr(invalid = 'ignore') #Don't let nans in rasters raise a warning.
#np.seterr(over = 'ignore') #These have been checked and screen valid exceptions.
np.seterr(divide = 'ignore') #could comment out for debugging if wish. 

def est_snow():
    swe_m = melt_one_day()
    out_swe = accum_one_day(swe_m)
    return out_swe
        
def melt_one_day(): # Melt factor used is an average of table 1 in Hock 2003. Journal of Hydrology 282:104-115. Math based on equation 2 of this publication.
    global melt_factor, low_temperature_differences, accumswe
    melt_delta = low_temperature_differences * melt_factor
    melt_delta = np.where(melt_delta < 0,0,melt_delta) # When Tavg < low_thresh_temperatures, melt_delta < 0, therefore this takes care of condition Td <= T0 in Hock 2003.
    end_swe = accumswe - melt_delta
    end_swe = np.where(end_swe < 0,0,end_swe)
    return end_swe
        
def accum_one_day(start_swe): # Dingman et al. 2002., Lutz et al. 2010.
    global precip_fraction, tmean, low_temperature_differences, low_thresh_temperatures, high_thresh_temperatures, precip
    rain = precip_fraction * low_temperature_differences
    rain = np.where(tmean <= low_thresh_temperatures, 0,rain)
    rain = np.where(tmean > high_thresh_temperatures, 1.0, rain)
    snow_increment = (1.0 - rain) * precip
    #snow_increment = np.where(snow_increment < 0,0,snow_increment)
    end_swe = start_swe + snow_increment
    return end_swe

def namefixlen(v):
    s = str(v)
    while len(s) < 3:
        s = '0' + s
    return s
 
def get_tmean(tmin,tmax):
    tmax = tmax - 273.15 # MACA GCM data provide temperatures in degrees Kelvin
    tmin = tmin - 273.15
    tmean = (tmax + tmin) /2
    return tmean,tmax,tmin
    
def leapyearlist():
    #Generate list of leap years 
    leapyearlist=list(range(1804,2300,4))
    leapyearlist.remove(1900);leapyearlist.remove(2100);leapyearlist.remove(2200)
    return leapyearlist   

def onlyabsneg(arr):
    ret = np.where(arr > 0,0,arr)
    ret = np.absolute(ret)
    return ret

def nonneg(arr):
    return np.where(arr <0,0,arr)

def hamon_evapotrans(): #daytime_length in hours, temperatures in degrees C.
    global daylength, tmean
    es = calc_Es() # Lutz et al 2010
    bottom = tmean+273.3
    bracket = es/bottom
    pet = 29.8*daylength * bracket # Days = 1 because daily timestep so not included.
    return pet #mm / day

def calc_inverse_rel_distance(doy): #CHECKS OK
    doy = float(doy)
    bracket_term = ((2*np.pi)/365) * doy
    retval = (.033 * np.cos(bracket_term)) + 1
    return retval #radians

def calc_solar_declination(doy): #CHECKS OK
    doy = float(doy)
    bracket_term = (((2*np.pi)/365) * doy) - 1.39
    retval = .409 * np.sin(bracket_term)
    return retval # radians

def calc_sunset_hour_angle(latitude,solar_declination): # CHECKS OK
    #Both inputs need to be in radians
    bracket_term = np.tan(latitude) * -1 * np.tan(solar_declination) # All trig function use radians
    bracket_term = np.where(bracket_term < -1,-1,bracket_term)
    bracket_term = np.where(bracket_term > 1, 1, bracket_term) #arccos is undefined outside the range [-1,1] and high latitudes in the dataset lead to too large values of np.tan(latitude)
    retval = np.arccos(bracket_term)
    return retval #radians

def calc_Rso(Ra): #CHECKS OK
    global elevation #meters
    # Equation 37 in ch3 of FAO doc.
    # Ra as calc from equation 21.
    # Elevation in meters
    correction_term = (.00002 * elevation) + .75
    Rso = Ra * correction_term
    return Rso
    
def calculate_Ra(inverse_rel_distance,sunset_hour_angle,latitude,solar_declination): #CHECKS OK
    #Ra = extra-terrestrial radiation in MJ/m2/day. Equation 21 of Ch3 FAO doc. Radiation hitting top of Earth's atmosphere.
    # Latitude (j) in radians. Positive for northern hemisphere and negative for southern hemisphere.
    # Sunset hour angle (Ws) in radians. Solar declination (d) in radians.
    Gsc = .0820 # Solar constant : MJ/m^2/day
    first_term = (float((24*60))/np.pi) * Gsc #37.586
    first_half_bracket_term = (np.sin(latitude) * np.sin(solar_declination)) * sunset_hour_angle
    second_half_bracket_term = (np.cos(latitude) * np.cos(solar_declination) * np.sin(sunset_hour_angle))
    full_bracket_term = first_half_bracket_term + second_half_bracket_term
    Ra = first_term * full_bracket_term * inverse_rel_distance
    return Ra

def calc_Rnl(Ea,Rs,Rso): #CHECKS OK
    #Rnl = net long-wave radiation in MJ/m2/day. Equation 39 Ch3 FAO doc.
    #Ea = actual vapor pressure (Kpa). Rs = Solar radiation in MJ/m2/day as calc by equation 35.
    #Rso = clear sky radiation in MJ/m2/day as calculated by equation 37.
    global tmax, tmin
    relative_shortwave_radiation = Rs/Rso
    #if relative_shortwave_radiation > 1 : relative_shortwave_radiation = 1
    relative_shortwave_radiation = np.where(relative_shortwave_radiation > 1, 1, relative_shortwave_radiation)
    tmax_k = tmax + 273.16;tmin_k = tmin + 273.16 # Convert to degrees Kelvin
    tmax_k = tmax_k**4;tmin_k=tmin_k**4
    left_bracket_term = (tmax_k+tmin_k)/2;left_bracket_term = left_bracket_term * .000000004903 # Stefan-Boltzmann constant
    middle_term = 0.34 - (.14 * np.sqrt(Ea))
    right_term = (1.35*relative_shortwave_radiation) - 0.35
    Rnl = left_bracket_term * middle_term * right_term
    return Rnl

def calc_Rn(Rnl,Rs): # CHECKS OK
    # Equation 40 of ch3 FAO doc.
    # Using estimate of Rns from equation 38 :> Rns = (1-a)*Rs. a is default set to 0.25 unless a regression has been done.
    Rns = 0.75* Rs
    Rn = Rns - Rnl
    return Rn

def calc_todays_PET(Et_type): 
    #print('Calculating PET')
    global tmax,tmin,tmean,precip,accumswe,year,day_index,daylength,latitude
    if Et_type == 'oudin': 
        inverse_rel_distance = calc_inverse_rel_distance(day_index + 1) # day_of_year = day_index + 1
        solar_declination = calc_solar_declination(day_index + 1)
        sunset_hour_angle = calc_sunset_hour_angle(latitude,solar_declination)
        Ra = calculate_Ra(inverse_rel_distance,sunset_hour_angle,latitude,solar_declination)
        Et = calc_oudin_pet(Ra)
    if Et_type == 'hamon': Et = hamon_evapotrans()
    elif Et_type == 'penman_montieth':
        inverse_rel_distance = calc_inverse_rel_distance(day_index + 1) # day_of_year = day_index + 1
        solar_declination = calc_solar_declination(day_index + 1)
        sunset_hour_angle = calc_sunset_hour_angle(latitude,solar_declination)
        Ra = calculate_Ra(inverse_rel_distance,sunset_hour_angle,latitude,solar_declination)
        srad = get_param_one_day(year, 'srad',day_index) # Incidient short wave radiation flux density, taken as an average over the daylight period of the day in MJ/m2/day.
        Rs= (srad * daylength) / 1000000 # Convert to daily total radiation in MJ/m2/day. See ref here: https://daac.ornl.gov/DAYMET/guides/Daymet_V3_CFMosaics.html
        # Rs is now equivalent to Rs in FAO procedures for calculating Penman Montieth: http://www.fao.org/docrep/X0490E/x0490e07.htm#meteorological%20factors%20determining%20et
        G = 0 # Soil heat flux is assumed to be zero for daily calculations
        Rso = calc_Rso(Ra) # Clear sky solar radiation
        
        Es_minus_Ea,Ea = calc_Es_minus_Ea() 
        Rnl = calc_Rnl(Ea,Rs,Rso) #The ratio of Rs/Rso constrained to 1.
        
        D = calc_D()
        Rn = calc_Rn(Rnl,Rs)
        Et = calc_Penman_Montieth_Eto(Rn,G,Es_minus_Ea,D)
    Et = np.where(accumswe > 0,0,Et) # No PET can happen when snow on ground.
    return Et

def calc_gamma():#CHECKS OK
    # Gamma is the psychrometric constant used for Penman_Montieth Eto
    global atmospheric_pressure
    gamma = 0.665 * .001 * atmospheric_pressure
    return gamma

def calc_atmospheric_pressure(): # CHECKS OK
    global elevation
    tr = .0065*elevation
    top = 293-tr
    inside = top/293
    right = inside**5.26
    atmos_pressure = 101.3*right
    return atmos_pressure

def calc_D(): 
    #Slope of the vapor pressure curve
    # t is AVERAGE air temperature in degrees C
    # Equation 13 in FAO doc  
    global tmean
    tk = tmean + 237.3
    bracket_term = (17.27*tmean)/(tk)
    right_term = 0.6108 * np.exp(bracket_term)
    top_term = 4098*right_term
    bottom_term = (tk)**2
    D = top_term/bottom_term
    return D


def calc_Es(): #CHECKS OK
    #Equation 12 in FAO doc
    # This is just saturation vapor pressure for Tmean
    global tmax, tmin
    e_tmax = calc_saturation_vapor_pressure(tmax)
    e_tmin = calc_saturation_vapor_pressure(tmin)
    Es = (e_tmax + e_tmin)/2
    return Es #kPa

def calc_Es_minus_Ea(): 
    # Will use humidity data or estimate Ea by assuming that dewpoint is two degrees below daily Tmin.
    global tmax, tmin
    Es = calc_Es()
    tmin_adjusted = tmin - 2 # Suggested for arid regions by FAO. In the absence of humidity data, assume dewpoint is 2 degrees below Tmin.
    Ea = calc_saturation_vapor_pressure(tmin_adjusted)
    Es_minus_Ea = Es - Ea
    return Es_minus_Ea,Ea

def calc_saturation_vapor_pressure(t): #CHECKS OK
    # This is Es for a given t
    # t = degrees C, equation 11 in FAO doc
    bracket_term = (17.27*t)/(t+237.3)
    right_term = np.exp(bracket_term)
    saturation_vapor_pressure = 0.6108 * right_term
    return saturation_vapor_pressure #kPa


def calc_Penman_Montieth_Eto(Rn,G,Es_minus_Ea,D,u2 = 2): #CHECKS OK WITH EXAMPLE 17 CHAPTER 4.
    #Rn = net radiation (MJ/m2/day), G = soil heat flux density (MJ/m2/day),Tavg in degrees C
    #Es_minus_Ea = vapor pressure deficit in KPa, D = slope of vapor pressure curve
    #gamma = psychrometric constant (Kpa / degrees C)
    #u2 = wind speed at 2m height, set to 2 m/s by default
    global gamma,tmean
    topleft_term = .408*D*(Rn-G)
    corrected_tavg = tmean + 273
    topmiddle_term = (900/corrected_tavg)*gamma
    topright_term = u2*Es_minus_Ea
    overall_top_term = (topmiddle_term*topright_term) + topleft_term
    wind_correction = 0.34*u2
    bottomright_term = (1 + wind_correction) * gamma
    overall_bottom_term = D + bottomright_term
    Eto = overall_top_term / overall_bottom_term
    Eto = np.where(Eto < 0,0,Eto)
    return Eto
    
def calc_oudin_pet(Ra): # Equation 3 in Oudin et al. 2010
    global tmean 
    # Since this is daily timestep, mdays = 1 and not included here.
    top = Ra * (tmean + 5) * .408
    bottom = 100
    PET = top/bottom
    PET = np.where(tmean > -5, PET, 0)
    return PET # mm / d

#These two funcs were used to make the heat load layer. Not needed for model runs. 
def fold_aspect(aspect): # Used to make the folded aspect array for heat load. Can be deleted or saved for reference.
    inner = np.absolute(aspect - 225)
    ret = np.absolute(180 - inner)
    ret = ret * (np.pi/180) # convert to radians
    return ret

def calc_heat_load(lat,slope,aspect): # S1 appendix to Lutz et al. 2010 # Used to make the heat load layer. Saved for reference.
    first = 0.339 + (0.808 * np.cos(lat) * np.cos(slope))
    print('first',np.nanmin(first),np.nanmax(first))
    second = 0.196 * (np.sin(lat) * np.sin(slope))
    print('second',np.nanmin(second),np.nanmax(second))
    third = 0.482 * (np.cos(aspect) * np.sin(slope)) # Folded aspect used, converted to radians for np funcs
    print('third',np.nanmin(third),np.nanmax(third))
    HL = first - second - third
    print('HL',np.nanmin(HL),np.nanmax(HL))
    return HL

def heat_load_adjust_pet():
    global heat_load, PET
    PET_HL = PET * heat_load
    return PET_HL

def calc_remove_fraction():
    # When pet = w, exponent = 0, therefore multiplier = 0, therefore amount removed from soil = 0
    # When w = 0, multiplier becomes 1 - (1/e^(pet/soil_whc))
    # When w < pet, multiplier becomes  1 - (1/e^(reduced_pet/soil_whc)), where reduced pet = pet - w
    # When pet = soil_whc, still can remove only a frac = 1 - (1/e)
    global soil_whc, soil_water, PET_adjusted, w 
    exponent = (nonneg((PET_adjusted - w)) / soil_whc) * -1 #B2
    exponent = np.where(np.isinf(exponent) == True, 0, exponent) # Some soil_whc cells indicate zero holding capacity. Division by zero in prev line creates inf, must be screened.
    multiplier = (1.0 - np.e**exponent) #B3
    decrements = soil_water * multiplier #B4
    return decrements

def calc_aet_and_runoff_and_soilw(): # soil_water and associated vars are 2D arrays. Does not use NDVI correction.
    # w = water input to soil = rain + melt
    global soil_water, PET_adjusted, w,year,day_index
    soilw_inputs = np.where(w > PET_adjusted, w - PET_adjusted, 0) # If PET >= w, assume no water input to soil. If w > PET then water input to soil = "excess", where excess = w-aet, but aet = pet when w > PET, so excess = w - pet. See equation 12 in Lutz et al. appendix 
    remove_fraction = calc_remove_fraction() #Delta soil in Lutz appendix eq 13.
    soilw_remove_amounts = np.where(w <= PET_adjusted, remove_fraction, 0) #If w > PET then no water removed from soil. If w <= PET, then Lutz eq 13 (calc_remove_fraction) used to determined amount removed.
    soil_water = soil_water + soilw_inputs #Since we are working with 2D arrays, some locations will have positive inputs, some will be zero.
    #if np.nanmin(soilw_remove_amounts) < 0 : raise Exception('{y}  {d} Negative soil water remove amounts detected. Terminating program'.format(y = year, d = day_index))
    #if np.nanmin(soilw_inputs) < 0 : raise Exception('{y}  {d} Negative soil water remove amounts detected. Terminating program'.format(y = year, d = day_index))
    soil_water = soil_water - soilw_remove_amounts #Since we are working with 2D arrays, some locations will have positive remove amounts, some will be zero.
    soil_water = nonneg(soil_water)
    runoff = nonneg(soil_water - soil_whc) # The nonneg function modifier ensures that only locations with positive runoff values pass through.
    soil_water = np.where(soil_water > soil_whc, soil_whc,soil_water) #Bring all locations with soil_water > whc back to whc. Previous line has gotten rid of that excess as runoff.
    total_evap_avail = w + remove_fraction #B4
    aet = np.where(PET_adjusted >= total_evap_avail,total_evap_avail,PET_adjusted) #From Lutz appendix: "AET equals the smaller of PET or (delta_soil + w)"; total_evap_avail = delta_soil + w
    aet = np.where(np.isinf(aet) == True, 0, aet)
    aet = np.where(np.isnan(PET_adjusted) == True, np.nan,aet)
    soil_water = np.where(np.isnan(PET_adjusted) == True, np.nan, soil_water)
    return aet,runoff

'''
# This is the version2 function which in corporates Senay et al. NDVI correction.
def calc_aet_and_runoff_and_soilw(): # soil_water and associated vars are 2D arrays
    # w = water input to soil = rain + melt
    global soil_water, PET_adjusted, w,year,day_index,rain, melt, soil_whc
    #------ Code below is from older version of model that used Lutz et al. methods--------------
    # remove_fraction = calc_remove_fraction() #Delta soil in Lutz appendix eq 13.
    # soilw_remove_amounts = np.where(w <= PET_adjusted, remove_fraction, 0) #If w > PET then no water removed from soil. If w <= PET, then Lutz eq 13 (calc_remove_fraction) used to determined amount removed.
    # soilw_inputs = np.where(w > PET_adjusted, w - PET_adjusted, 0) # If PET >= w, assume no water input to soil. If w > PET then water input to soil = "excess", where excess = w-aet, but aet = pet when w > PET, so excess = w - pet. See equation 12 in Lutz et al. appendix 
    # soil_water = soil_water + soilw_inputs #Since we are working with 2D arrays, some locations will have positive inputs, some will be zero.
    # soil_water = soil_water - soilw_remove_amounts #Since we are working with 2D arrays, some locations will have positive remove amounts, some will be zero.
    # ----End old model code -----
    # --------------------------------------------------------------------------------------------------------------------------
    # Lines below are from Senay et al. Use NDVI to calculate AET
    SWi = soil_water + rain + melt # This is different from Senay et al. They use SWi = previous_soil_water + effective_precipitation. This is previous soil water + fraction of precip falling as rain + melt from snow.
    # Notes we are no longer using the calc_remove_fraction() to apply the exponential decay of soil water removal.
    varA = 1.25
    varB = 0.2
    ndvi_name = 'ndvi' + namefixlen(day_index+1) + '.npy.npz'
    ndvi_array = np.load(ndvi_name)
    ndviF = ndvi_array['ndvi']
    etasw1A = (varA * ndviF + varB) * PET_adjusted
    etasw1B = (varA * ndviF) * PET_adjusted
    etasw1 = np.where(ndviF > 0.4, etasw1A, etasw1B)
    etasw2 = (SWi / (0.5 * soil_whc)) * etasw1
    etasw3 = np.where(SWi > (0.5 *soil_whc), etasw1, etasw2)
    etasw4 = np.where(etasw3 > SWi, SWi, etasw3)
    aet = np.where(etasw4 > soil_whc, soil_whc, etasw4)
    etasw1A = 0 # reset to save memory
    etasw1b = 0
    etasw1 = 0
    etasw2 = 0
    etasw3 = 0
    etasw4 = 0
    # End Senay et al. NDVI Code
    # -----------------------------------------------------
    aet = np.where(np.isinf(aet) == True, 0, aet) # Fix unrealistic values
    aet = np.where(np.isnan(PET_adjusted) == True, np.nan,aet) # Fix unrealistic values
    soil_water_temp = SWi - aet # per Senay et al. pers comm.
    soil_water = np.where(SWi > soil_whc, soil_whc - aet, soil_water_temp)
    soil_water = nonneg(soil_water)
    runoff = nonneg(SWi - soil_whc) # This has changed from previous version of model. Was runoff = soil_water - soil_whc
    soil_water = np.where(soil_water > soil_whc, soil_whc,soil_water) #Bring all locations with soil_water > whc back to whc. Previous lines have gotten rid of that excess as runoff.
    #total_evap_avail = w + remove_fraction #B4 # From old model with Lutz equations.
    #aet = np.where(PET_adjusted >= total_evap_avail,total_evap_avail,PET_adjusted) #From Lutz appendix: "AET equals the smaller of PET or (delta_soil + w)"; total_evap_avail = delta_soil + w # From old model
    soil_water = np.where(np.isnan(PET_adjusted) == True, np.nan, soil_water)
    return aet,runoff
'''

def chunks(l, n):
    #Yield successive n-sized chunks from l.
    for i in range(0, len(l), n):
        yield l[i:i + n]

def mp_write_daily(): # Subsets to just CONUS and multiplies the output_mult_factor if applicable
    global file_handles,day_index,var_dict, output_params, npz_cores,jobs, output_mult_factor, model, scenario
    chunk_size = npz_cores
    param_chunks = chunks(output_params,chunk_size)
    var_dict = {'PET':PET_adjusted,'AET':AET,'runoff':runoff,'Deficit':deficit,'rain':rain,'water_input_to_soil':w,'melt':melt,
               'days_snow':daily_snow,'agdd':agdd,'accum_precip':accum_precip,'accumswe':accumswe,'soil_water':soil_water}  
    for chunk in param_chunks:
        jobs = []
        for f in chunk:
            #output_chunk = var_dict[f][subsetting_indices[0][0]:subsetting_indices[0][1],subsetting_indices[1][0]:subsetting_indices[1][1]] # Subsetting spatially here to just CONUS
            output_chunk = var_dict[f]
            p = multiprocessing.Process(target = write_daily_to_npz,args=(model,scenario,f,output_chunk,output_mult_factor))
            p.start()
            jobs.append(p)
        for j in jobs:
            j.join()
            
def write_daily_to_npz(model,scenario,param,data, output_mult_factor):
    filename = output_data_path + model + '_' + scenario + '_' + param + '_' + str(current_date.year) + '_' + current_date.strftime("%j")  + '.npz'
    data = data * output_mult_factor
    if output_mult_factor > 1:
        data = np.round(data,0)
        data = np.where(np.isnan(data) == True, -9999, data)
    np.savez_compressed(filename, param = data)
  

def launch_new_collation(year, variable):
    print(f"Collating {year} {variable}")
    collate_pool.apply_async(collate_into_netcdf, args = (year, variable))
    
def collate_into_netcdf(year, variable):
    global model, scenario, reference_ds
    npz_files = glob.glob(output_data_path + f"{model}_{scenario}_{variable}_{year}_*.npz")
    npz_files.sort()

    signif_digits = 1
   
    encoding = {
        'compression': 'zlib',
        #'shuffle': True,
        'complevel': 3,
        #'fletcher32': False,
        #'contiguous': False,
        #'chunksizes': chunksizes
        '_FillValue': -9999.0
    }

    dataarrays = []
    
    for file in npz_files:
        ## Ex CanESM2_rcp45_AET_2006_006.npz
        yday = file.split("_")[-1].split(".")[0]
        file_date = datetime.datetime.strptime(f"{year}-{yday}", "%Y-%j")

        arr = np.load(file)['param']
        arr = np.around(arr, decimals = signif_digits)
        if all_ints: 
            arr = arr.astype('i2')
            
        da = xr.DataArray(data = arr,
                          dims = ("y", "x"),
                          coords = dict(
                              x = ("x", reference_ds.x.data),
                              y = ("y", reference_ds.y.data),
                              time = ((), file_date),
                          ),
                          name = variable,
                          attrs = dict(
                              description=variable,
                              units = "mm")
                          )
        da.encoding = encoding
        dataarrays.append(da)

    combined = xr.concat(dataarrays, dim = "time") \
                 .rio.write_crs(str(reference_ds.rio.crs), inplace = True) \
                 .rio.write_grid_mapping(inplace = True) \
                     .rio.write_coordinate_system(inplace = True)

    print("Writing: " + output_data_path + f"{model}_{scenario}_{year}_{variable}.nc")
    
    with xr.set_options(keep_attrs = True):
        combined.to_netcdf(path = output_data_path + f"{model}_{scenario}_{year}_{variable}.nc",
                           mode = "w",
                           engine = "netcdf4",
                           format = "NETCDF4",
                           )

    for f in npz_files:
        os.remove(f)
       
def pixel2coord(x, y, raster):
    """Returns global coordinates from pixel x, y coords"""
    # GDAL affine transform parameters, According to gdal documentation xoff/yoff are image left corner, a/e are pixel wight/height and b/d is rotation and is zero if image is north up. 
    xoff, a, b, yoff, d, e = raster.GetGeoTransform()
    
    xp = a * x + b * y + xoff
    yp = d * x + e * y + yoff
    return(xp, yp)

def raster2xy(raster):
    """Returns xy for each pixel in units of original crs, i.e., lat long or meters"""
    rows, colms = raster.ReadAsArray().shape
    x, y = np.fromfunction(pixel2coord, shape = (rows, colms), raster = raster)
    return (x, y)

def make_correction_factors(elevation_error):
    ## Average of N/S lapse rates in Tercek 2021 paper
    high_lapse_rates = (
        ### Average N/S slopes then divide by 1000 to convert C/km to C/m
        mean([-5.98, -5.91]) / 1000,
        mean([-6.49, -7.89]) / 1000,
        mean([-7.33, -10.52]) / 1000,
        mean([-9.28, -11.29]) / 1000,
        mean([-9.03, -10.78]) / 1000,
        mean([-8.28, -10.02]) / 1000,
        mean([-7.82, -9.65]) / 1000,
        mean([-7.32, -8.40]) / 1000,
        mean([-7.71, -8.99]) / 1000,
        mean([-8.47, -9.20]) / 1000,
        mean([-6.67, -7.97]) / 1000,
        mean([-6.24, -6.33]) / 1000
    )

    low_lapse_rates = (
        mean([-1.62, -0.92]) / 1000,
        mean([-2.03, -2.46]) / 1000,
        mean([-2.78, -2.65]) / 1000,
        mean([-4.08, -4.52]) / 1000,
        mean([-3.50, -3.23]) / 1000,
        mean([-1.67, -0.70]) / 1000,
        mean([-0.34, -0.32]) / 1000,
        mean([-0.68, 0.41]) / 1000,
        mean([-0.63, -0.37]) / 1000,
        mean([-2.42, -3.31]) / 1000,
        mean([-2.62, -2.95]) / 1000,
        mean([-2.68, -2.51]) / 1000
    )

    high_correction_factors = list(map(lambda lapse: lapse * elevation_error, high_lapse_rates))
    low_correction_factors = list(map(lambda lapse: lapse * elevation_error, low_lapse_rates))

    return (high_correction_factors, low_correction_factors)
  
if __name__ == '__main__':
    model = sys.argv[1]
    scenario = sys.argv[2]
    site = sys.argv[3]
    
    input_data_path = f"{os.environ['HOME']}/out/{site}/daily-split/" ## daily-split directory with daily T/P layers
    output_data_path = f"{os.environ['HOME']}/out/{site}/wb/" ## Output Directory for daily water balance npz
    site_data_path = f"../data/input/{site}/"

    Path(output_data_path).mkdir(parents = True, exist_ok = True)
    
    print(f"Model: {model}, Scenario: {scenario}")
    print(f"Input dir: {input_data_path}")
    print(f"Output dir: {output_data_path}")

    site_data = pd.read_csv("sites.csv")
    site_data = site_data[site_data['site'] == site]
    ### Get the elevation the metdata cell containing the site uses.
    ### This will be compared with the "real" elevations from the DEM
    ### to get an elevation error for lapse rate corrections
    metdata_elev = site_data.values[0][3] ### ugly, need to learn how to lookup in pandas df or make a dict from the csv
    
    web = True
    npz_cores = 1
    multiprocessing.set_start_method('fork')  ## Fork is POSIX only, will not run on Windows
    collate_cores = 4 # This can be raised once the model loops finish.

    melt_factor = 4.0
    precip_fraction = 0.167
    if web == False: textfilenamelist = 'canesm2_list.txt'
    else: textfilenamelist = 'null'
    start_time = time.time()  
    PET_type = 'oudin'
    agdd_base = 5.5 # Degrees C
    output_mult_factor = 10.0 # See def mp_write_daily(). For smaller file sizes, this multiplies all output values by 10 and stores them as 16 bit integer.
    all_ints = True # If True, all output netcdfs will be int16. If False, then categories in def launch_new_collation() are used.

    elev = gdal.Open(site_data_path + "dem/dem_nad83.tif")
    elevation = elev.ReadAsArray()
    elevation_error = elevation - metdata_elev ## How much higher or lower input DEM is compared to metdata elevations

    high_correction_factors, low_correction_factors = make_correction_factors(elevation_error)
    
    slope = gdal.Open(site_data_path + "dem/slope_nad83.tif").ReadAsArray()
    aspect = gdal.Open(site_data_path + "dem/aspect_nad83.tif").ReadAsArray()

    reference_ds = rioxarray.open_rasterio(site_data_path + "dem/dem_nad83.tif")
    
    ## Slope in radians, aspect in degrees
    ## Aspect is used to run old_aspect() which takes aspects
    ## in degrees and returns folded aspect in radians
    ## 
    ## calc_heat_load() takes latitude, slope, and
    ## folded aspect in radians
    ##
    ## Heat load should be from 0.03 - 1.11
    aspect = np.where(aspect < 0, np.nan, aspect)
    aspect_folded = fold_aspect(aspect)
    slope = np.where(slope < 0, np.nan, slope)
    slope = np.where(slope > 60, np.nan, slope) ## Slope is 0-60 for eq3 in McCune 2002
    slope = slope * (np.pi / 180)

    height = elevation.shape[0]
    width = elevation.shape[1]
    
    new_x, new_y = raster2xy(elev)
    ## Determine lat and lons from utm values
    new_lats, new_lons = utm.to_latlon(new_x, new_y, 12, 'N')
    if (np.nanmax(new_lats) > 60) or (np.nanmin(new_lats) < 30) : raise Exception("Latitude outside of range for heatload function (30-60 degrees). Terminating.")
    
    ##output_params = ['soil_water','PET','AET','Deficit','runoff','agdd','accumswe', 'rain']
    output_params = ['AET', 'Deficit']
    output_units = {'PET':'mm','AET':'mm','Deficit':'mm','accumswe':'mm','melt':'mm','days_snow':'mm','rain':'mm','water_input_to_soil':'mm','runoff':'mm','agdd':'C','accum_precip':'mm'}       
    accum_precip = np.zeros((height,width)).astype(np.float32)
    agdd = np.zeros((height,width)).astype(np.float32)
    last_accumswe = np.zeros((height,width)).astype(np.float32)

    latitude = (np.pi / 180) * new_lats
    if (np.nanmax(latitude) > 1.5) or (np.nanmin(latitude) < 0.1) : raise Exception('Latitude is not in radians or wrong file has been used. Terminating.')
    heat_load = calc_heat_load(latitude, slope, aspect_folded)
    if (np.nanmax(heat_load) > 1.11) or (np.nanmin(heat_load) < 0.03) : raise Exception("Heat load exceeds specified values. Terminating.")
    np.savez_compressed(f"{site}-{model}-{scenario}-heatload.npz", heat_load)
    
    #Global vars used here to save duplication in memory associated with passing to funcs
    ##elevation = np.load(input_data_path + 'etopo1_aligned_array.npy') # In meters
    ##multiplier_mask = np.load(input_data_path + 'multiplier_mask.npz')['arr_0']
    multiplier_mask = np.full((height, width), 1).astype(np.float32)
    if (np.nanmax(elevation) > 5700) or (np.nanmin(elevation) < 0): raise Exception('Do you have the wrong elevation file? Terminating.')

    soil_whc = gdal.Open(site_data_path + "soil/soil_whc_025.tif").ReadAsArray()
    if (np.nanmax(soil_whc) < 10) : raise Exception("Soil values appear to be in cm, convert to mm before use. Terminating") ### Soil values are between around 10-100mm for 25cm depths.  Need to adjust this if deeper soil depths are used
    soil_water = np.copy(soil_whc) # Initialize soil values at full.

    ##snow_thresh_file = np.load(input_data_path + 'jennings_t50_coefficients.npz') # Jennings, K. et al. 2018. Spatial variation of the rain-snow temperature threshold across the northern hemisphere. Nature Communications 9: 1148. DOI: 10.1038/s41467-018-03629-7
    ##snow_thresh_temperatures = snow_thresh_file['t50']
    snow_thresh_temperatures = gdal.Open(site_data_path + "jennings_t50_coefficients.tif").ReadAsArray()
    low_thresh_temperatures = snow_thresh_temperatures - 3.0 # The Jennings coefficients are T50, i.e. where precip is half snow, half rain
    high_thresh_temperatures = snow_thresh_temperatures + 3.0 # This sets up a 6 degree span, which corresponds to the 1/6 = 0.167 precip fraction.
    
    # if PET_type == 'penman_montieth':
    #     #slope = np.load('aligned_slope_in_radians.npy') # radians needed for numpy trig functions
    #     #aspect = np.load('aligned_folded_aspect_in_radians.npy')
    #     atmospheric_pressure = calc_atmospheric_pressure() # Used to calculate gamma. 
    #     gamma = calc_gamma()
    
    jobs = []
    
    log_file = open(f'{site}-{model}-{scenario}-logfile.csv','w')

    if model == "historical":
        climate_data = pd.read_csv(site_data_path + "gridmet_1979_2023.csv")
        year = 1979
        pr_data = climate_data.pr
        tmax_data = climate_data.tmmx
        tmin_data = climate_data.tmmn
    else:
        climate_data = pd.read_csv(site_data_path + "macav2metdata_2006_2099.csv")
        year = 2006
        pr_data = climate_data.get(f"pr_{model}_r1i1p1_{scenario}")
        tmax_data = climate_data.get(f"tasmax_{model}_r1i1p1_{scenario}")
        tmin_data = climate_data.get(f"tasmin_{model}_r1i1p1_{scenario}")
        
    dates = climate_data.date

    accumswe = np.zeros((height, width)) # Set snow to zero everywhere for start of run. Note first year "Spin up" will not have accurate snow values so must consider only data beginning in start_year + 1.

    day_index = 0

    collate_pool = multiprocessing.Pool(processes = collate_cores)
    
    for index in range(0, len(pr_data)):
                
        precip = np.full((height, width), pr_data[index])
        tmax = np.full((height, width), tmax_data[index])
        tmin = np.full((height, width), tmin_data[index])

        last_year = year
        current_date = datetime.datetime.strptime(dates[index],  "%Y-%m-%dT%H:%M:%SZ").date()
        year = current_date.year
        month = current_date.month

        if year > last_year:
            day_index = 0
            print(f"{last_year} complete")
            for variable in output_params:
                    launch_new_collation(last_year, variable)
                    
        else:
            day_index += 1
        
        var_dict = {}
        print(f"Calculating {str(current_date)}, Day Index: {day_index}, Year: {year}")
        log_file.write(str(current_date))

        ## Lapse Rate Corrections
        tmax = tmax + high_correction_factors[month - 1]
        tmin = tmin + low_correction_factors[month - 1] 

        tmean,tmax,tmin = get_tmean(tmin, tmax)
            
        low_temperature_differences = tmean - low_thresh_temperatures
        
        accum_precip = accum_precip + precip
        gdd = nonneg(tmean - agdd_base)
        agdd = agdd + gdd
        
        accumswe = est_snow()

        PET = calc_todays_PET(PET_type)
        PET_adjusted = heat_load_adjust_pet()        
        snow_diff = accumswe - last_accumswe
        melt = onlyabsneg(snow_diff)
        daily_snow = nonneg(snow_diff)
        rain = nonneg(precip - daily_snow)
        w = rain + melt 
        
        AET,runoff = calc_aet_and_runoff_and_soilw() 
        deficit = nonneg(PET_adjusted - AET)
        last_accumswe = accumswe
        
        mp_write_daily()


    ### Collate last year, loop ends at max_year - 1
    for variable in output_params:
        print("Collating final year")
        launch_new_collation(year, variable)

    collate_pool.close()
    collate_pool.join()

    log_file.close()
            
    print('Done!!')
