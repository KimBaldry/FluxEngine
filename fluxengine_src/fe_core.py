#!/usr/bin/env python2

# ofluxghg-flux-calc.py IGA working version - now using GIT at OCEANFLUX-SHARED/workspace/i.g.c.ashton_exeter.ac.uk/FE_IGA_new
#test
# utility to load various input netcdf datasets and determine
# air-sea flux of CO2 based on parameterisations in the ESA STSE OceanFlux Greenhouse Gases Technical Specification (TS) 


# History
# date, decsription, author, email, institution
# v0 10/2012  Basic structure and output, Jamie Shutler, jams@pml.ac.uk, Plymouth Marine Laboratory.
# v1 24/10/2012  improved error handling and some method defs, Jamie Shutler, jams@pml.ac.uk, Plymouth Marine Laboratory.
# v3 03/12/2012  additional functionality for uncertainty analyses, improved error handling, jams@pml.ac.uk, Plymouth Marine Laboratory.
# v4 06/08/2013 approaching final form and almost compliant with TS, jams@pml.ac.uk, Plymouth Marine Laboratory
# v5 06/08/2014 version used for final runs and FluxEngine publication
#Further changes marked by #IGA within the text. Version control within a GIT repo

# ToDo (wishList)
# Add uncertainty values as an option through configuration file
# Incorporate rain noise through configuration file
# Pressure units often vary - account for this through config?
# Line 1736 (re.match) returns an error if runParams.pco2_prod is 'pco2'
# 

 # netcdf bits
import pandas as pd
import sys
from math import log, exp, pow, isnan;
from numpy import size, flipud, mean, zeros, nonzero, array, resize, ma, arange, dtype, ones, meshgrid, where;
from numpy import any as npany;
from random import normalvariate
import logging;
from os import path;

from datalayer_edit import DataLayer, DataLayerMetaData;
from settings import Settings;
from debug_tools import calc_mean; #calculate mean ignoring missing values.

from datetime import timedelta, datetime;

 # debug mode switches
DEBUG = False
DEBUG_PRODUCTS = True
VERIFICATION_RUNS = False # forces flux calculations to use Takahashi SST_t as the SST dataset for the pco2 data
DEBUG_LOGGING = True; # TH: Added to make debugging easier. Sets up a logger object and writes a to the specified path (or cwd)
#workingDirectory = os.getcwd(); #Only used to make a full path for logs, every other path is suppled as an absolute filepath.
                                #Should really be passed as an argument.


 # missing value set for intermediate data sets and output dataset
missing_value = -999.0
fill_value = -999.0
missing_value_int = -999
fill_value_int = -999

# # valid data ranges for testing output products and flagging violations
# # orignally set in TS, some have been modified based on the actual input datasets
# # these are used in the NetCDF metadata and to check the valid data ranges to create the
# # OFA11 process indicator layer
#OF_min = -0.75
#OF_max = 0.6
#OK1_min = 0.0
#OK1_max = 100.0
#OK3_min = 0.0
#OK3_max = 100.0
#OKD_min = 0.0
#OKD_max = 100.0
#OKR_min = 0.0
#OKR_max = 100.0
#OKB1_min = 0.0
#OKB1_max = 100.0
#OKS1_min = 15.0
#OKS1_max = 40.0
#OKT1_min = -1.8
#OKT1_max = 30.5
#OFWR_min = -0.75
#OFWR_max = 0.6
#
# # these thresholds are in g-C m^-3 (rather than ppm as in the TS)
#OSC1_min = 0.095
#OSC1_max = 0.27
#OIC1_min = 0.11
#OIC1_max = 0.22
#
# # these thresholds are in uatm
#OBPC_min = 200
#OBPC_max = 800

# function name definition for stderr and stdout messages
function = "(ofluxghg-flux-calc, main)"



#Stores a set of parameters required to run the FluxEngine, allowing easy centralised access.
class RunParameters:
    def set_parameters(self, parameterDict):
        for key in vars(self).keys(): #remove all previously stored values
            delattr(self, key);
        for key in parameterDict.keys(): #Create instance variables for every key:value pair
            setattr(self, key, parameterDict[key]);


#
# method definitions
#
# writing the final netcdf output
### change to writing a .txt file K.B.
def write_txt(fluxEngineObject, verbose=False):
    
    
    dataLayers = fluxEngineObject.data;
    runParams = fluxEngineObject.runParams;

    #outputChunk = int(runParams.run_count % runParams.output_temporal_chunking);
    #fluxEngineObject.logger.debug("Writing netCDF output with output_chunk = %d", outputChunk);
    
#if outputChunk == 0: #Create a new file and write output to it.
    n = fluxEngineObject.n;
    latitudeData = fluxEngineObject.latitude_data;
    longitudeData = fluxEngineObject.longitude_data;

    function="write_txt";

    data = {'lat':latitudeData, 'lon':longitudeData} 
    # Create DataFrame 
    fulldf = pd.DataFrame(data) 
    metadf = pd.DataFrame(columns=['variable','attribute','value'])
    #data layers
    #
    for dataLayerName in dataLayers:
        try:
            if verbose:
                print "Writing datalayer '"+dataLayerName+"' to .tsv file as "+dataLayers[dataLayerName].netCDFName;

            data = dataLayers[dataLayerName].fdata; #fdata is usually a view by sometimes a copy so it has to be done this way. There is probably a better way to do this.
            fulldf[dataLayers[dataLayerName].netCDFName] = data;

        except AttributeError as e:
            print "%s:No netCDFName or data attribute found in DataLayer '%s'." % (function, dataLayerName);
            raise e;
        except ValueError as e:
            print type(e), e.args;
            print "%s: Cannot resize datalayer '%s'" % (function, dataLayers[dataLayerName].name);
            raise e;
# Need to write attributes as a header, or as another file

        #variable.missing_value = missing_value;
        scale_factor = 1.0;
        add_offset = 0.0;
        metadf = metadf.append({'variable': dataLayers[dataLayerName].name, 'attribute': 'missing_value' , 'value': str(missing_value)}, ignore_index=True)
        metadf = metadf.append({'variable': dataLayers[dataLayerName].name, 'attribute': 'scale_factor' , 'value': str(scale_factor)}, ignore_index=True)
        metadf = metadf.append({'variable': dataLayers[dataLayerName].name, 'attribute': 'add_offset' , 'value': str(add_offset)}, ignore_index=True)
        try:
            if dataLayers[dataLayerName].units != None:
                #variable.units = dataLayers[dataLayerName].units;
                metadf = metadf.append({'variable': dataLayerName, 'attribute': 'units' , 'value': dataLayers[dataLayerName].units}, ignore_index=True)
        except AttributeError:
            print "%s: No units found for datalayer named '%s'." % (function, dataLayerName);

        try:
            if dataLayers[dataLayerName].minBound != None:
                #variable.valid_min = dataLayers[dataLayerName].minBound;
                metadf = metadf.append({'variable': dataLayerName, 'attribute': 'valid_min' , 'value': dataLayers[dataLayerName].minBound}, ignore_index=True)
        except AttributeError:
            print "%s: No minBound found for datalayer named '%s'." % (function, dataLayerName);

        try:
            if dataLayers[dataLayerName].maxBound != None:
                #variable.valid_max = dataLayers[dataLayerName].maxBound;
                metadf = metadf.append({'variable': dataLayerName, 'attribute': 'valid_max' , 'value': dataLayers[dataLayerName].maxBound}, ignore_index=True)
        except AttributeError:
            print "%s: No maxBound found for datalayer named '%s'." % (function, dataLayerName);

        try:
            if dataLayers[dataLayerName].standardName != None:
                #variable.standard_name = dataLayers[dataLayerName].standardName;
                metadf = metadf.append({'variable': dataLayerName, 'attribute': 'standard_name' , 'value': dataLayers[dataLayerName].standardName}, ignore_index=True)
        except AttributeError:
            print "%s: No standardName found for datalayer named '%s'." % (function, dataLayerName);

        try:
            if dataLayers[dataLayerName].longName != None:
                #variable.long_name = dataLayers[dataLayerName].longName;
                metadf = metadf.append({'variable': dataLayerName, 'attribute': 'long_name' , 'value': dataLayers[dataLayerName].longName}, ignore_index=True)
        except AttributeError:
            print "%s: No longName found for datalayer named '%s'." % (function, dataLayerName);
        
    #set some global attributes
    metadf = metadf.append({'variable': 'Conventions', 'attribute': '' , 'value': 'CF-1.6'}, ignore_index=True)
    metadf = metadf.append({'variable': 'Institution', 'attribute': '' , 'value': 'Originally developed by the partners of the ESA OceanFlux GHG and OceanFlux GHG Evolution projects. Now continued by the CarbonLab team at the University of Exeter. Modified by Kimberlee Baldry at the University of Tasmania'}, ignore_index=True)
    metadf = metadf.append({'variable': 'Contact', 'attribute': '' , 'value': 'email: j.d.shutler@exeter.ac.uk and kimberlee.baldry@utas.edu.au'}, ignore_index=True)


    #Output all the parameters used in this run.
    for paramName in vars(runParams).keys():
        paramValue = getattr(runParams, paramName);
        if paramValue is not None:
            if type(paramValue) is bool: #netCDF does not support bool types.
                paramValue = int(paramValue);
            elif isinstance(paramValue, timedelta): #netCDF does not support object instances.
                paramValue = str(paramValue);
            elif paramValue == None: #netCDF does not support None type.
                paramValue = "None";
                
#                 setattr(ncfile, paramName, paramValue);
        
#         if int(fluxEngineObject.runParams.output_temporal_chunking) != 1:
#             setattr(ncfile, "start_year", fluxEngineObject.runParams.year);
#             setattr(ncfile, "start_month", fluxEngineObject.runParams.month);
#             setattr(ncfile, "start_day", fluxEngineObject.runParams.day);
#             setattr(ncfile, "start_hour", fluxEngineObject.runParams.hour);
#             setattr(ncfile, "start_minute", fluxEngineObject.runParams.minute);
#             setattr(ncfile, "start_second", fluxEngineObject.runParams.second);
     
           
    fulldf.to_csv(runParams.output_path,sep = '\t',index = False)
    front, end = runParams.output_path.split('.')
    metadf.to_csv(front+"_meta.tsv",sep = '\t',index = False)
    return 0;


#Calculates solubility of distilled water (i.e. assuming 0 salinity) given global temperature.
#overwrites solubilityDistilled with the calculated value.
#TODO: no need to pass nx, ny into all these functions.
#TODO: rain_wet_deposition_switch isn't used!
def calculate_solubility_distilled(solubilityDistilled, salinity, rain_wet_deposition_switch, sstskin, deltaT, n):
    #First create a 0 salinity dataset
    salDistil = array([missing_value] * len(salinity))
    for i in arange(n):
        if (salinity[i] != missing_value):
            salDistil[i] = 0.0
        else:
            salDistil[i] = missing_value
    
    #Next calculate the solubility using the zero salinity 'distilled water' dataset
    solubilityDistilled = solubility(sstskin, salDistil, deltaT, n, True);
    return solubilityDistilled;

#whitecapping data using relationship from TS and parameters from table 1 of Goddijn-Murphy et al., 2010, equation r1
#Takes two DataLayers as input. whitecap is written in place.
def calculate_whitecapping(windu10, whitecap):
    if (not isinstance(windu10, DataLayer)) or (not isinstance(whitecap, DataLayer)):
        raise ValueError("ofluxghg_flux_calc.calculate_whitecapping: Invalid arguments. Arguments must be DataLayer type.");
    if windu10.n != whitecap.n:
        raise ValueError("ofluxghg_flux_calc.calculate_whitecapping: Invalid arguments. windu10 and whitecap dimensions do not match (%d, %d vs %d, %d)." % (windu10.n, whitecap.n));
    
    for i in arange(windu10.n):
        if (windu10.fdata[i] != missing_value):
            whitecap.fdata[i] = 0.00159 * pow(windu10.fdata[i], 2.7)
        else:
            whitecap.fdata[i] = missing_value
    return whitecap.fdata;


def add_noise(data, rmse_value, mean_value, no_elements):
# randomly adding noise to data array, based log normal distribution
# rmse value is used as the standard deviation of the noise function

   # intialise the random number generator
  for i in arange(no_elements):
     if ( (data[i] != missing_value) and (data[i] != 0.0) ):
      orig = data[i]
      value = log(data[i])
      stddev = float(rmse_value/orig)
      noise = normalvariate(0,stddev) # determines the random noise value based on the input data and the standard deviation of the uncertainty
      value = value + noise
      data[i] = exp(value)
  return data

######Ians alternative version (thanks PLand) - Untested#############

# def add_noise(data, err_value, mean_value, no_elements):
#  # randomly adding noise to data array, based log normal distribution
#  # rmse value is used as the standard deviation of the noise function

#     if type(err_value) is not float:#If not a float, its an array of spatially varying error variances (rain)
#         sizetup = data.shape
#         data = (data + np.random.normal(size = sizetup) * err_value).clip(0)
#     else:
#         # intialise the random number generator
#         for i in arange(no_elements):
#             if ( (data[i] != missing_value) and (data[i] != 0.0) ):         
#                 orig = data[i]
#                 value = log(data[i])
#                 stddev = float(err_value/orig)
#                 noise = normalvariate(0,stddev) # determines the random noise value based on the rain input data and the standard deviation of the uncertainty
#                 value = value + noise
#                 data[i] = exp(value)
#         return data
######Ians altered version############END

def add_bias_k_biology_wind(data, bias_k_value, biology_data, biology_value, wind_data, wind_value, no_elements, percent_switch):
 # adding bias offset to data array based on biology and wind conditions

   for i in arange(no_elements):
      if ( (data[i] != missing_value) and (biology_data[i] != missing_value) and (wind_data[i] != missing_value) ):         
         if ( (biology_data[i] >= biology_value) and (wind_data[i] <= wind_value) ):
            if percent_switch == 1:
               data[i] -= (data[i]*(bias_k_value/100.0))
            else:
               data[i] += bias_k_value
   return data
   
def add_bias(data, bias_value, no_elements):
 # adding bias offset to data array

   for i in arange(no_elements):
      if ( (data[i] != missing_value) ):         
         orig = data[i]
         data[i] = orig + bias_value
   return data

def add_sst_rain_bias(data, bias_value, rain_intensity, rain_data, wind_speed, wind_data, no_elements):
 # adding bias offset to data array based on rain intensity and wind speed

   for i in arange(no_elements):
      if ( (data[i] != missing_value) ):
         if ((rain_data[i] >= rain_intensity) and (wind_data[i] <= wind_speed)):        
            orig = data[i]
            data[i] = orig + bias_value
   return data


#determine the schmidt number
#based on Schmid relationship from Wanninkhof1992 - Relationship between wind speed and gas exchange over the ocean, JGR Oceans
def schmidt_Wanninkhof1992(sstC_fdata, n, gas):
#calculating the schmidt data

   sc_fdata = array([missing_value] * n)
   if 'o2' in gas.lower():
       for i in arange(n):
          if (sstC_fdata[i] != missing_value):
             sc_fdata[i] = 1953.4 - (128.0 * sstC_fdata[i]) + (3.9918 * (sstC_fdata[i] * sstC_fdata[i])) - (0.050091 * (sstC_fdata[i] * sstC_fdata[i] * sstC_fdata[i]))#IGA-O2
          else:
        # assigning invalid values
             sc_fdata[i] = missing_value

   if 'n2o' in gas.lower():
       for i in arange(n):
          if (sstC_fdata[i] != missing_value):            
             sc_fdata[i] = 2301.1 - (151.1 * sstC_fdata[i]) + (4.7364 * (sstC_fdata[i] * sstC_fdata[i])) - (0.059431 * (sstC_fdata[i] * sstC_fdata[i] * sstC_fdata[i]))#IGA-N2O
          else:
        # assigning invalid values
             sc_fdata[i] = missing_value

   if 'ch4' in gas.lower():
       for i in arange(n):
          if (sstC_fdata[i] != missing_value):
              sc_fdata[i] = 2039.2 - (120.31 * sstC_fdata[i]) + (3.4209 * (sstC_fdata[i] * sstC_fdata[i])) - (0.040437 * (sstC_fdata[i] * sstC_fdata[i] * sstC_fdata[i]))#IGA-CH4
          else:
        # assigning invalid values
             sc_fdata[i] = missing_value
   if 'co2' in gas.lower():
       for i in arange(n):
          if (sstC_fdata[i] != missing_value):
              # relationship is only valid for temperatures <=30.0 oC
              sc_fdata[i] = 2073.1 - (125.62 * sstC_fdata[i]) + (3.6276 * (sstC_fdata[i] * sstC_fdata[i])) - (0.043219 * (sstC_fdata[i] * sstC_fdata[i] * sstC_fdata[i]))#IGA-CO2
          else:
        # assigning invalid values
             sc_fdata[i] = missing_value
   return sc_fdata

#based on Schmid relationship from Wanninkhof2014 - Relationship between wind speed and gas exchange over the ocean revisited, Limnology and Oceanography
def schmidt_Wanninkhof2014(sstC_fdata, n, gas):
    sc_fdata = array([missing_value] * n)
    if 'o2' in gas.lower():
        for i in arange(n):
            if (sstC_fdata[i] != missing_value):
                sc_fdata[i] = 1920.4 + (-135.6 * sstC_fdata[i]) + (5.2122 * sstC_fdata[i]**2) + (0.10939 * sstC_fdata[i]**3) + (0.00093777 * sstC_fdata[i]**4);
            else:
            # assigning invalid values
                sc_fdata[i] = missing_value

    if 'n2o' in gas.lower():
        for i in arange(n):
            if (sstC_fdata[i] != missing_value):            
                sc_fdata[i] = 2356.2 + (-166.38 * sstC_fdata[i]) + (6.3952 * sstC_fdata[i]**2) + (-0.13422 *sstC_fdata[i]**3) + (0.0011506 * sstC_fdata[i]**4);
            else:
            # assigning invalid values
                sc_fdata[i] = missing_value

    if 'ch4' in gas.lower():
        for i in arange(n):
            if (sstC_fdata[i] != missing_value):
                sc_fdata[i] = 2101.2 + (-131.54 * sstC_fdata[i]) + (4.4931 * sstC_fdata[i]**2) + (-0.08676 * sstC_fdata[i]**3) + (0.00070663 * sstC_fdata[i]**4);
            else:
            # assigning invalid values
                sc_fdata[i] = missing_value
    if 'co2' in gas.lower():
        for i in arange(n):
            if (sstC_fdata[i] != missing_value):
                # relationship is only valid for temperatures <=30.0 oC
                sc_fdata[i] = 2116.8 + (-136.25 * sstC_fdata[i]) + (4.7353 * sstC_fdata[i]**2) + (-0.092307 * sstC_fdata[i]**3) + (0.0007555 * sstC_fdata[i]**4);
            else:
            # assigning invalid values
                sc_fdata[i] = missing_value
    return sc_fdata
    

#solubility calculation equation from Table A2 of Wanninkkhof, JGR, 1992
def solubility_Wanninkhof1992(sstK, sal, deltaT, n, flux_calc, gas):
    #solubility calculation
    #equation from Table A2 of Wanninkkhof, JGR, 1992
    sol = array([missing_value] * n)
    if gas == 'co2':
        for i in arange(n):
            if ( (sstK[i] != missing_value) and (sal[i] != missing_value) and (sstK[i] > 0.0) ):
                sol[i] = -60.2409 + ( 93.4517*(100.0 / sstK[i]) ) + (23.3585 * (log(sstK[i]/100.0))) + (sal[i] * (0.023517 + ( (-0.023656)*(sstK[i]/100.0)) + (0.0047036*( (sstK[i]/100.0)*(sstK[i]/100.0) ) ) ) );
                sol[i] = exp(sol[i])
                #runParams.flux_calc is a switch to remove Delta_T-Sb component - ie selects use of RAPID or EQUILIBRIUM flux models from Woolf et al., 2016
                if flux_calc != 2:
                    deltaT[i] = 0.0
                sol[i] = sol[i] * (1 - (0.015*deltaT[i]))
            else:
                sol[i] = missing_value
    elif gas == 'o2':
        for i in arange(n):
            if ( (sstK[i] != missing_value) and (sal[i] != missing_value) and (sstK[i] > 0.0) ):
                sol[i] = -58.3877 + ( 85.8079*(100.0 / sstK[i]) ) + (23.8439 * (log(sstK[i]/100.0))) + (sal[i] * (-0.034892 + ( (0.015568)*(sstK[i]/100.0)) + (-0.0019387*( (sstK[i]/100.0)*(sstK[i]/100.0) ) ) ) );
                sol[i] = exp(sol[i])
                #runParams.flux_calc is a switch to remove Delta_T-Sb component - ie selects use of RAPID or EQUILIBRIUM flux models from Woolf et al., 2016
                if flux_calc != 2:
                    deltaT[i] = 0.0
                sol[i] = sol[i] * (1 - (0.015*deltaT[i]))
            else:
                sol[i] = missing_value
        
    elif gas == 'n2o':
        for i in arange(n):
            if ( (sstK[i] != missing_value) and (sal[i] != missing_value) and (sstK[i] > 0.0) ):
                sol[i] = -64.8539 + ( 100.2520*(100.0 / sstK[i]) ) + (25.2049 * (log(sstK[i]/100.0))) + (sal[i] * (-0.062544 + ( (0.035337)*(sstK[i]/100.0)) + (-0.0054699*( (sstK[i]/100.0)*(sstK[i]/100.0) ) ) ) );
                sol[i] = exp(sol[i])
                #runParams.flux_calc is a switch to remove Delta_T-Sb component - ie selects use of RAPID or EQUILIBRIUM flux models from Woolf et al., 2016
                if flux_calc != 2:
                    deltaT[i] = 0.0
                sol[i] = sol[i] * (1 - (0.015*deltaT[i]))
            else:
                sol[i] = missing_value
        
    elif gas == 'ch4':
        for i in arange(n):
            if ( (sstK[i] != missing_value) and (sal[i] != missing_value) and (sstK[i] > 0.0) ):
                sol[i] = -68.8862 + ( 101.4956*(100.0 / sstK[i]) ) + (28.7314 * (log(sstK[i]/100.0))) + (sal[i] * (-0.076146 + ( (0.043970)*(sstK[i]/100.0)) + (-0.0068672*( (sstK[i]/100.0)*(sstK[i]/100.0) ) ) ) );
                sol[i] = exp(sol[i])
                #runParams.flux_calc is a switch to remove Delta_T-Sb component - ie selects use of RAPID or EQUILIBRIUM flux models from Woolf et al., 2016
                if flux_calc != 2:
                    deltaT[i] = 0.0
                sol[i] = sol[i] * (1 - (0.015*deltaT[i]))
            else:
                sol[i] = missing_value
    
    return sol

#solubility calculation equation Table 2 of Wanninkhof, Rik. "Relationship between wind speed and gas exchange over the ocean revisited." Limnology and Oceanography: Methods 12.6 (2014): 351-362.
def solubility_Wanninkhof2014(sstK, sal, deltaT, n, flux_calc, gas):
    #solubility calculation
    #equation from Table A2 of Wanninkkhof, JGR, 1992
    sol = array([missing_value] * n)
    if gas == 'co2':
        for i in arange(n):
            if ( (sstK[i] != missing_value) and (sal[i] != missing_value) and (sstK[i] > 0.0) ):
                sol[i] = -58.0931 + ( 90.5069*(100.0 / sstK[i]) ) + (22.2940 * (log(sstK[i]/100.0))) + (sal[i] * (0.027766 + ( (-0.025888)*(sstK[i]/100.0)) + (0.0050578*( (sstK[i]/100.0)*(sstK[i]/100.0) ) ) ) );
                sol[i] = exp(sol[i])
                #runParams.flux_calc is a switch to remove Delta_T-Sb component - ie selects use of RAPID or EQUILIBRIUM flux models from Woolf et al., 2016
                if flux_calc != 2:
                    deltaT[i] = 0.0
                sol[i] = sol[i] * (1 - (0.015*deltaT[i]))
            else:
                sol[i] = missing_value
    elif gas == 'o2':
        for i in arange(n):
            if ( (sstK[i] != missing_value) and (sal[i] != missing_value) and (sstK[i] > 0.0) ):
                sol[i] = -58.3877 + ( 85.8079*(100.0 / sstK[i]) ) + (23.8439 * (log(sstK[i]/100.0))) + (sal[i] * (-0.034892 + ( (0.015568)*(sstK[i]/100.0)) + (-0.0019387*( (sstK[i]/100.0)*(sstK[i]/100.0) ) ) ) );
                sol[i] = exp(sol[i])
                #runParams.flux_calc is a switch to remove Delta_T-Sb component - ie selects use of RAPID or EQUILIBRIUM flux models from Woolf et al., 2016
                if flux_calc != 2:
                    deltaT[i] = 0.0
                sol[i] = sol[i] * (1 - (0.015*deltaT[i]))
            else:
                sol[i] = missing_value
        
    elif gas == 'n2o':
        for i in arange(n):
            if ( (sstK[i] != missing_value) and (sal[i] != missing_value) and (sstK[i] > 0.0) ):
                sol[i] = -62.7062 + ( 97.3066*(100.0 / sstK[i]) ) + (24.1406 * (log(sstK[i]/100.0))) + (sal[i] * (-0.058420 + ( (0.033193)*(sstK[i]/100.0)) + (-0.0051313*( (sstK[i]/100.0)*(sstK[i]/100.0) ) ) ) );
                sol[i] = exp(sol[i])
                #runParams.flux_calc is a switch to remove Delta_T-Sb component - ie selects use of RAPID or EQUILIBRIUM flux models from Woolf et al., 2016
                if flux_calc != 2:
                    deltaT[i] = 0.0
                sol[i] = sol[i] * (1 - (0.015*deltaT[i]))
            else:
                sol[i] = missing_value
        
    elif gas == 'ch4':
        for i in arange(n):
            if ( (sstK[i] != missing_value) and (sal[i] != missing_value) and (sstK[i] > 0.0) ):
                sol[i] = -68.8862 + ( 101.4956*(100.0 / sstK[i]) ) + (28.7314 * (log(sstK[i]/100.0))) + (sal[i] * (-0.076146 + ( (0.043970)*(sstK[i]/100.0)) + (-0.0068672*( (sstK[i]/100.0)*(sstK[i]/100.0) ) ) ) );
                sol[i] = exp(sol[i])
                #runParams.flux_calc is a switch to remove Delta_T-Sb component - ie selects use of RAPID or EQUILIBRIUM flux models from Woolf et al., 2016
                if flux_calc != 2:
                    deltaT[i] = 0.0
                sol[i] = sol[i] * (1 - (0.015*deltaT[i]))
            else:
                sol[i] = missing_value
    
    return sol

#Calculate the mass boundary layer concentration (ie concentration in the water)
#Each argument should be supplied as fdata matrix (flattened matrix)
def calculate_concw(concFactor, foundationSolubility, fCO2water, concw):
    for i in range(len(concw)):
        if ( (foundationSolubility[i] != missing_value) and (fCO2water[i] != missing_value)):
            concw[i] = foundationSolubility[i] * concFactor * fCO2water[i];
        else:
            concw[i] = missing_value;

#Calculate the interfacial concentration (i.e. at the interface between the ocean and the atmosphere)
#Each argument should be supplied as fdata matrix (flattened matrix)
def calculate_conca(concFactor, skinSolubility, fCO2air, conca):
    for i in range(len(conca)):
        if ( (skinSolubility[i] != missing_value) and (fCO2air[i] != missing_value)):
            conca[i] = skinSolubility[i] * concFactor * fCO2air[i];
        else:
            conca[i] = missing_value;


#Copies missing_values from 'master' to 'derived' e.g. for making mean and stddev datasets consistent.
#master and derived should be DataLayers
def copy_missing_values(master, derived, missingValue=DataLayer.missing_value):
    if master.n != derived.n:
        raise ValueError("copy_missing_values: master and derived DataLayers do not have the same dimensions: (%d, %d) versus (%d, %d)." % (master.n, derived.n));
    for i in arange(master.n):
        if (master.fdata[i] == missingValue):
            derived.fdata[i] = missing_value;

#Creates / adds to the failed_quality_fdata layer for each element of 'data' which exceeds the specified min/max range.
def check_output_dataset(datalayer, failed_quality_fdata):
    #function = "(check_output, main)"
    #checking the contents of a dataset
    for i in arange(len(datalayer.fdata)):
        if (datalayer.fdata[i] != DataLayer.missing_value):
            if datalayer.maxBound != None and datalayer.fdata[i] > datalayer.maxBound:
                #print "\n%s Dataset %s fails OceanFluxGHG TS table 7 valid limits, (defined min/max are: %lf/%lf, found %lf at grid point %d (%lf) exiting" % (function, name, min_range, max_range, data[i], i, i/nx)
                #print "\n%s Coincident data values sstskinC_fdata:%lf sstfndC_fdata:%lf windu10_fdata:%lf sig_wv_ht_fdata:%lf solfnd_fdata:%lf solskin_fdata:%lf k_fdata:%lf concw_fdata:%lf conca_fdata:%lf sal_fdata:%lf" % (function, sstskinC_fdata[i], sstfndC_fdata[i], windu10_fdata[i], sig_wv_ht_fdata[i], solfnd_fdata[i], solskin_fdata[i], k_fdata[i], concw_fdata[i], conca_fdata[i], sal_fdata[i])
                if failed_quality_fdata[i] == DataLayer.missing_value:
                    failed_quality_fdata[i] = 1; # first entry so need to initiase
                else:
                    failed_quality_fdata[i] += 1;
            elif datalayer.minBound != None and datalayer.fdata[i] < datalayer.minBound:
                #print "\n%s Dataset %s fails OceanFluxGHG TS table 7 valid limits, (defined min/max are: %lf/%lf, found %lf at grid point %d (%lf) exiting" % (function, name, min_range, max_range, data[i], i, i/nx)
                #print "\n%s Coincident data values sstskinC_fdata:%lf sstfndC_fdata:%lf windu10_fdata:%lf sig_wv_ht_fdata:%lf solfnd_fdata:%lf solskin_fdata:%lf k_fdata:%lf concw_fdata:%lf conca_fdata:%lf sal_fdata:%lf" % (function, sstskinC_fdata[i], sstfndC_fdata[i], windu10_fdata[i], sig_wv_ht_fdata[i], solfnd_fdata[i], solskin_fdata[i], k_fdata[i], concw_fdata[i], conca_fdata[i], sal_fdata[i])
                if failed_quality_fdata[i] == DataLayer.missing_value:
                    failed_quality_fdata[i] = 1; # first entry so need to initiase
                else:
                    failed_quality_fdata[i] += 1;




#Returns false if dataLayer dimensions do not match the reference dimensions.
def check_dimensions(dataLayer, ref_n, DEBUG=False):
   function = "(check_dimensions, main)"
   if len(dataLayer.data) == ref_n:
      if DEBUG:
         print "\n%s Input data (%s) have identical dimensions to reference values (%s, %s) "% (function, dataLayer.name,len(dataLayer.data))
         return True;
   else:
      print "\n%s Input data ('%s') dimensions are non-identical to reference (new: %s, %s is not equal to: %s, %s)." % (function, dataLayer.name,len(dataLayer.data), ref_n)
      return False;


class FluxEngine:
    def __init__(self, parameterDict):
        self.runParams = RunParameters();
        self.runParams.set_parameters(parameterDict);
        
        #self.rootDirectory = path.abspath(path.dirname(inspect.stack()[0][1]));
        #TODO: Shouldn't used src_home, this is going to be removed from config file soon. Use rootDirectory instead.
        self.defaultSettings = Settings(path.join(path.dirname(__file__), "settings.xml")); #Load default settings metadata for datalayers.
        self.data = {};
        self.kParameterisationFunctors = [] #List of functor objects which encapsulates the k rate calculation.
                                            #The k rate calculation can be extended in a modular fashion.
                                            #Each functor is called in order.
        self.processIndicatorFunctors = []; #List of functor objects which calculate process indicator layers.
                                            #Each functor is called in order, hence functors later in the list can use outputs from those earlier as inputs.
        
        self._load_lon_lat_time();
    
    #Adds a data layer by reading a netCDF file. Optionally adds stddev and count data from the same file.
    #preprocessing is an optional list of functions to transform the data before storing.
    #TODO: transposeData should now be added as a preprocessing function
    def add_data_layer(self, name, infile, prod, stddevProd=None, countProd=None, transposeData=False, preprocessing=None):
        function = "(fe_core, FluxEngine.add_data_layer)";
        
        try:
            self._add_single_data_layer(name, infile, prod, transposeData, preprocessing=preprocessing);
            
        except KeyError as e:
            print "\n%s: Data variable '%s' is missing from %s input (%s)" % (function, prod, name, infile);
            print e, e.args;
            raise e;
           # return False;
#        except ValueError as e: #E.g. incorrect number of dimensions
#            print "\n%s: %s" % (function, e.args);
#            raise e;
#            #return False;
       
            
        #Datalayer was successfully added.
        return True;
    
    #Adds a single datalayer and acquires metadata
    #TODO: transposeData is now be handled by preprocessing, so this cna be removed...
    def _add_single_data_layer(self, name, infile, prod, transposeData=False, preprocessing=None):
        #function = "(ofluxghg_flux_calc, FluxEngine._add_single_data_layer)";
        
        metaData = self._extract_data_layer_meta_data(name);

        dl = DataLayer.create_from_file(name, infile, prod, metaData, transposeData=transposeData, preprocessing=preprocessing);
        self.data[name] = dl;
    
    #Creates a DataLayer which is filled (by default) with DataLayer.missing_value.
    #The DataLayer will be added to self.data and is accessable by it's 'name', e.g. self.data["new_datalayer"].
    #If no dimensions (nx, ny) are specified, then the default self.nx and self.ny dimensions are used.
    #If there is known metadata for the datalayer (e.g. min, max, units etc.) these will be loaded from settings.xml
    #Metadata can be overwritten in the config file using the the datalayer name and the xml attribute from the settings.xml file, e.g.:
    #       datalayername_units = C m^2s^-1
    #       datalayername_maxBound = 100.0
    def add_empty_data_layer(self, name, n=None, fillValue=DataLayer.missing_value):
        if n==None: n=self.n;
        
        
        metaData = self._extract_data_layer_meta_data(name);        
        dl = DataLayer.create_empty_datalayer(name, n, metaData, fillValue=fillValue);
        
        self.data[name] = dl;
    
    #Returns a DataLayerMetaData instance for the named DataLayer.
    #First it will search for default values specified in the settings xml file
    #Failing that the default default values are used
    #Finally any metadata values which are specified in the config file are used to overwrite default values
    def _extract_data_layer_meta_data(self, name):
        #Extract default metadata
        if name in self.defaultSettings.allDataLayers: #If default values are specified in the settings.xml file
            metaData = self.defaultSettings.allDataLayers[name];
        else: #Otherwise create metadata with 'default' default values
            metaData = DataLayerMetaData(name);
        
        #overwrite default metadata with any found in the run parameters that correspond to the specified name.
        for attribute in vars(metaData):
            if attribute != name:
                if name+"_"+attribute in vars(self.runParams):
                    if attribute in ["temporalChunking", "temporalSkipInterval"]: #Integers
                        setattr(metaData, attribute, int(getattr(self.runParams, name+"_"+attribute))); #TODO: should be handled in setup_tools::create_run_parameters really.
                    else: #floats
                        setattr(metaData, attribute, getattr(self.runParams, name+"_"+attribute));
        
        return metaData;
    
    #Defines the k parameterisation to be used. Requires an callable object which:
    #   implements input_names() and output_names()
    #   and contains standard_name and long_name string attributes.
    def add_k_parameterisation_component(self, kObject):
        function = "(FluxEngine.add_k_parameterisation_component)";
        if kObject != None:
            self.kParameterisationFunctors.append(kObject);
        else:
            raise ValueError("%s: Trying to add 'None' k parameterisation component." % function);
    
    def add_process_indicator_functor(self, piObject):
        function = "(FluxEngine.add_process_indicator_functor)";
        if piObject != None:
            self.processIndicatorFunctors.append(piObject);
        else:
            raise ValueError("%s: Trying to add 'None' process indicator layer component." % function);

    def run(self):
        status = self._check_datalayers(); #Check for consistency of datalayers (e.g. dimensions all match one another).
        if status == True:
            return self._run_fluxengine(self.runParams);
        else:
            print "Not all datalayers are consistent (check dimensions of data).";
            return status;
 
    #Reads longitude, latitude, time and dimension sizes (nx, ny).
    def _load_lon_lat_time(self):
        function = "(ofluxghg_flux_calc, FluxEngine._load_lon_lat_time)";
        
        #Find the path of the relevent datalayer.
        try:
            axesDatalayerInfile = getattr(self.runParams, self.runParams.axes_data_layer+"_infile");
        except ValueError as e:
            print "Couldn't find file path for axes_data_layer (%s). Check that this is correctly defined in the config file (it must match the name of another datalayer)." % self.runParams.axes_data_layer;
            print e.args;
        
        #Read longitude, latitude and time data from the sstskin infile, so open this file.
        try:
            dataset = pd.read_table(axesDatalayerInfile, sep="\t");
        except IOError as e:
            print "\n%s: axes_data_layer (%s) inputfile %s does not exist" % (function, self.runParams.axes_data_layer, axesDatalayerInfile);
            print type(e), e.args;
        
        #Read lat and lat
        try:
            self.latitude_data = dataset[self.runParams.latitude_prod];
            self.longitude_data = dataset[self.runParams.longitude_prod];
            self.year_data = dataset[self.runParams.year_prod];
            if self.latitude_data[0]<0: #IGA - it is a vector that is in opposite orientation to 'taka'
                self.latitude_data = flipud(self.latitude_data);
        except KeyError as e:
            raise ValueError ("%s: Couldn't find longitude (%s) and/or latitude (%s) variables in %s. Have you set longitude_prod and latitude_prod correctly in your configuration file?" % (function, self.runParams.longitude_prod, self.runParams.latitude_prod, axesDatalayerInfile));

                  
        #set time (since 1st Jan 1970)
#         try:
#             #curDatetime = datetime(self.runParams.year, self.runParams.month, self.runParams.day, self.runParams.hour, self.runParams.minute, self.runParams.second);
#             #self.time_data = (curDatetime - datetime(1970, 1, 1)).total_seconds();
#             #self.time_data = dataset.variables[self.runParams.time_prod][:];
#         except KeyError as e:
#             raise ValueError("%s: Couldn't find time (%s%) variables in %s. Have you set time_prod correctly in your configuration file?" % (function, self.runParams.time_prod, self.runParams.sstskin_infile));

        #set dimention
        self.n = len(self.latitude_data)
    #Ran after each datalayer is read in. Checks for consistency between datalayers.
    #Rescales data layers which can be rescaled.
    #returns true if all data layers are successfully validated.
    def _check_datalayers(self):
        #Check that dimensions match
        for key in self.data:
            if check_dimensions(self.data[key], self.n,  DEBUG) == False:
                print e.args;
                return False;
        return True;
    
    #Applies the mask to all data layers
    def _apply_mask(self):
        if "mask" in self.data:
            toIgnore = where(self.data["mask"].fdata == 0);
            
            #Apply mask to each data layer
            for dataLayerName in self.data:
                if dataLayerName != "mask":
                    self.data[dataLayerName].fdata[toIgnore] = self.data[dataLayerName].missing_value = True;


    def _run_fluxengine(self, runParams):
        function = "(ofluxghg_flux_calc, FluxEngine._run_fluxengine_)";
        
        #Set up logging object
        try:
            self.logger = logging.getLogger('FluxEngine_debug_log');
            #hdlr = logging.FileHandler(os.path.join(workingDirectory, runParams.LOG_PATH), filemode='w')
            hdlr = logging.FileHandler(runParams.LOG_PATH);
            formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s');
            hdlr.setFormatter(formatter);
            self.logger.addHandler(hdlr);
            if DEBUG_LOGGING == True:
                self.logger.setLevel(logging.DEBUG);
            else:
                self.logger.setLevel(logging.DEBUG);
        except:
            print "\n%s Couldn't initialise logger at path %s" % (function, runParams.LOG_PATH);
        
        #TODO: Replace directly with self.nx, self.ny, no need to use local variables here.
        n = self.n;
        
        #Apply mask to all data layers, if applicable.
        self._apply_mask(); #Checks for existance of mask internally.
        
        self.add_empty_data_layer("windu10_moment2");
        self.add_empty_data_layer("windu10_moment3");
        for i in arange(self.n):
          if self.data["windu10"].fdata[i] != missing_value:
             self.data["windu10_moment2"].fdata[i] = self.data["windu10"].fdata[i]*self.data["windu10"].fdata[i]
             self.data["windu10_moment3"].fdata[i] = self.data["windu10"].fdata[i]*self.data["windu10"].fdata[i]*self.data["windu10"].fdata[i]
     
        if runParams.TAKAHASHI_DRIVER == True:
           #need to generate the moment2 and moment3 data
           self.add_empty_data_layer("windu10_moment2");
           self.add_empty_data_layer("windu10_moment3");
           for i in arange(self.n):
              if self.data["windu10"].fdata[i] != missing_value:
                 self.data["windu10_moment2"].fdata[i] = self.data["windu10"].fdata[i]*self.data["windu10"].fdata[i]
                 self.data["windu10_moment3"].fdata[i] = self.data["windu10"].fdata[i]*self.data["windu10"].fdata[i]*self.data["windu10"].fdata[i]
        
        #SOCATv4 - using input foundation temperature as the SST temp-------------START
        #If there is no pco2_sst data, we need to get it / generate it.
        if "pco2_sst" not in self.data: #SOCATv4
            try:
                print "No pco2_sst data was supplied."
                self.add_empty_data_layer("pco2_sst");
                self.data["pco2_sst"].fdata = self.data["sstfnd"].fdata-273.15; #copy/convert sstfnd
            except (IOError, KeyError, ValueError) as e:
                print "pco2_sst data not available and could read sstfnd so cannot proceed.";
                print type(e), "\n"+e.args;
                return 1;   
        
        #some specific pco2 conditions
        #TODO: This shouldn't be randomly here.
        if "pco2_sw" in self.data:
            self.data["pco2_sw"].fdata[abs(self.data["pco2_sw"].fdata)<0.1]=DataLayer.missing_value#IGA_SOCATv4

            
        #initialising some data structures for the calculations
        self.add_empty_data_layer("scskin");
        self.add_empty_data_layer("scfnd");
        self.add_empty_data_layer("solubility_skin");
        self.add_empty_data_layer("solubility_fnd");
        self.add_empty_data_layer("pH2O");
        
        #apply additive saline skin value, but why selectively apply it?
        #TODO: Why are these hard-coded values. Should be using minBound and maxBound?
        self.add_empty_data_layer("salinity_skin");
        for i in arange(self.n):
            if (self.data["salinity"].fdata[i] >= 0.0) and (self.data["salinity"].fdata[i] <= 50.0):
                if ( (self.data["salinity"].fdata[i] + runParams.saline_skin_value) <= 50.0):
                    self.data["salinity_skin"].fdata[i] = self.data["salinity"].fdata[i] + runParams.saline_skin_value
            else:
                self.data["salinity"].fdata[i] = missing_value
                self.data["salinity_skin"].fdata[i] = missing_value
        if (runParams.saline_skin_value != 0.0):
            print "%s Using the saline skin model (%lf psu added to skin salinities)" % (function, runParams.saline_skin_value)
        
        #conversion of rain data from mm day-1 to mm hr^-1
        if "rain" in self.data:
            for i in arange(self.n):   
               if (self.data["rain"].fdata[i] != DataLayer.missing_value):
                   self.data["rain"].fdata[i] /= 24.0;
        
        
        #interpreting fnd_data option
        ####Derive sstskin as necessary.
        #Two possible datasets: sstskin and sstfnd, from different parts of the water column.
        #One option: sst gradients
        #If using only one dataset then copy that over the other. Otherwise keep both.
        #runParams.cool_skin_difference is the assumed temperature difference between the skin and foundation layer
        
        #using sstfnd, so copy sstfnd into sstskin
        #sstskin = sstfnd
        if runParams.sst_gradients_switch == 0 and runParams.use_sstskin_switch == 0 and runParams.use_sstfnd_switch == 1:
           print "%s SST gradient handling is off, using SSTfnd data selection in configuration file for all components of the flux calculation (this will ignore any SSTskin data in configuration file)." % (function)
           #actually copy sstfnd data into the sstskin dataset to make sure
           for i in arange(n):
               if self.data["sstfnd"].fdata[i] != missing_value:
                   self.data["sstskin"].fdata[i] = self.data["sstfnd"].fdata[i]
                   
               else:
                   self.data["sstskin"].fdata[i] = missing_value;
        
        
        #IGA added for the case where only foundation is provided and gradients are on------------------------------
        #must estimate sstskin (= sstfnd - runParams.cool_skin_difference)
        elif runParams.sst_gradients_switch == 1 and runParams.use_sstskin_switch == 0 and runParams.use_sstfnd_switch == 1:
            print "%s Using SSTfnd data selection with correction for skin temperature (SSTskin = SSTfnd - %d)(ignoring SSTskin data in configuration file)." % (function, runParams.cool_skin_difference)
            #actually copy sstfnd data into the sstskin dataset to make sure
            if "sstskin" not in self.data: #Must add the sstskin layer first!
                self.add_empty_data_layer("sstskin");
            for i in arange(n): #sstdkin = sstfnd - runParams.cool_skin_difference
                if self.data["sstfnd"].fdata[i] != missing_value:
                    self.data["sstskin"].fdata[i] = self.data["sstfnd"].fdata[i]-runParams.cool_skin_difference
                else:
                    self.data["sstskin"].fdata[i] = missing_value;
        
        #Using sstskin, so calculate it from sstfnd.
        elif runParams.sst_gradients_switch == 0 and runParams.use_sstskin_switch == 1 and runParams.use_sstfnd_switch == 0:
           print "%s SST gradient handling is off, using SSTskin to derive SSTfnd (SSTfnd = SSTskin + %d) for flux calculation (ignoring SSTfnd data in configuration file)." % (function, runParams.cool_skin_difference)
           #setting sstfnd_ data fields to skin values
           for i in arange(n):
              if self.data["sstskin"].fdata[i] != missing_value:
                 
                 self.data["sstfnd"].fdata[i] = self.data["sstskin"].fdata[i] + runParams.cool_skin_difference
                 self.data["sstskin"].fdata[i] = self.data["sstskin"].fdata[i] + runParams.cool_skin_difference

              else:
                 self.data["sstfnd"].fdata[i] = missing_value
                 self.data["sstskin"].fdata[i] = missing_value
                 
        elif runParams.sst_gradients_switch == 1 and runParams.use_sstskin_switch == 1 and runParams.use_sstfnd_switch == 0:
           print "%s SST gradient handling is on, using SSTskin and SSTfnd = SSTskin + %d for flux calculation (ignoring SSTfnd data in configuration file)." % (function, runParams.cool_skin_difference)
            #setting sstfnd_ data fields to skin values  
           for i in arange(n):
              if self.data["sstskin"].fdata[i] != missing_value:
                 self.data["sstfnd"].fdata[i] = self.data["sstskin"].fdata[i] + runParams.cool_skin_difference
              else:
                 self.data["sstfnd"].fdata[i] = missing_value
        elif runParams.sst_gradients_switch == 0 and runParams.use_sstskin_switch == 1 and runParams.use_sstfnd_switch == 1:
           print "%s SST gradient handling is off, using SSTfnd and SSTskin from the configuration file." % (function)        
        elif runParams.sst_gradients_switch == 1 and runParams.use_sstskin_switch == 1 and runParams.use_sstfnd_switch == 1:
           print "%s SST gradient handling is on, using SSTfnd and SSTskin from the configuration file." % (function)
        else:
           print "\n%s sst_gradients_switch (%d), use_sstskin_switch (%d) and use_sstfnd_switch (%d) combination in configuration not recognised, exiting." % (function, runParams.sst_gradients_switch, runParams.use_sstskin_switch, runParams.use_sstfnd_switch)
           return 1;
       

        #quality filtering and conversion of SST datasets
        #Calculate sstskinC and sstfndC
        self.add_empty_data_layer("sstskinC");
        self.add_empty_data_layer("sstfndC");
        for i in arange(n):
            if self.data["sstskin"].fdata[i] != missing_value:
                self.data["sstskinC"].fdata[i] = self.data["sstskin"].fdata[i] - 273.15;
                self.data["sstfndC"].fdata[i] = self.data["sstfnd"].fdata[i] - 273.15;
        
         # ability to randomly perturb the input datasets
         # needed for the ensemble analyses
         # stddev of noise is using published RMSE for each dataset
         # all datasets are considered to have bias=0, hence mean of noise=0
        if (runParams.random_noise_windu10_switch == 1):
           add_noise(self.data["windu10"].fdata, 0.44, 0.0, n)
           add_noise(self.data["windu10_moment2"].fdata, 0.44, 0.0, n)
           add_noise(self.data["windu10_moment3"].fdata, 0.44, 0.0, n)
           print "%s Adding random noise to windu10_mean, windu10_moment2 and windu10_moment3 (mean 0.0, stddev 0.44 ms^-1 - assuming using ESA GlobWave data)" % (function)
        
        if (runParams.random_noise_sstskin_switch == 1):
           add_noise(self.data["sstskin"].fdata, 0.14, 0.0, n)
           print "%s Adding random noise to sstskin (mean 0.0, stddev 0.14 ^oC - assuming using ESA CCI ARC data)" % (function)
        
        if (runParams.random_noise_sstfnd_switch == 1):
           add_noise(self.data["sstfnd"].fdata, 0.6, 0.0, n)
           print "%s Adding random noise to sstfnd (mean 0.0, stddev 0.6 ^oC - assuming using OSTIA data)" % (function)
        
        if (runParams.random_noise_pco2_switch == 1):
           print "/n%s Shape of pco2 data",self.data["pco2_sw"].fdata.shape
           add_noise(self.data["pco2_sw"].fdata, 6, 0.0, n)
           print "%s Adding random noise to pco2/fco2 data (mean 0.0, stddev 6 uatm - Using Candyfloss data, value provided by J Shutler)" % (function)
           #print "%s Adding random noise to pco2/fco2 data (mean 0.0, stddev 2.0 uatm - assuming using SOCAT flag A and flag B data)" % (function)
        
        #Ians rain noise test
        #if (random_noise_rain == 1):
        # self.data["rain"].data = reshape(self.data["rain"].data,self.data["rain"].data.shape[0]*self.data["rain"].data.shape[1])
        # rain_err_data = reshape(rain_err_data,rain_err_data.shape[0]*rain_err_data.shape[1])
        # add_noise(self.data["rain"].data, rain_err_data, 0.0, rain_nx*rain_ny)
        # print "%s Adding random noise to rain data using variance field" % (function)
        
        #Ians rain noise test end
        
         # bias values to be added here
        if (runParams.bias_windu10_switch == 1):
           add_bias(self.data["windu10"].fdata, runParams.bias_windu10_value, n)
            # makes no sense to add bias to second and third order moments, as any bias in the system 
            # would only impact on the mean (which is the 1st order moment)
           print "%s Adding bias to windu10_mean (not to second and third order moments) (value %lf ms^-1)" % (function, runParams.bias_windu10_value)
        
        if (runParams.bias_sstskin_switch == 1):
           add_bias(self.data["sstskin"].fdata, runParams.bias_sstskin_value, n)
           print "%s Adding bias noise to sstskin (value %lf ^oC)" % (function, runParams.bias_sstskin_value)
        
        if (runParams.bias_sstfnd_switch == 1):
           add_bias(self.data["sstfnd"].fdata, runParams.bias_sstfnd_value, n)
           print "%s Adding bias noise to sstfnd (value %lf ^oC)" % (function, runParams.bias_sstfnd_value)
        
        if (runParams.bias_pco2_switch == 1):
           add_bias(self.data["pco2_sw"].fdata, runParams.bias_pco2_value, n)
           print "%s Adding bias noise to pco2w (value %lf uatm)" % (function, runParams.bias_pco2_value)
        
        
         # bias based on change in sstskin due to rain
        if (runParams.bias_sstskin_due_rain_switch == 1):
           add_sst_rain_bias(self.data["sstskin"].fdata, runParams.bias_sstskin_due_rain_value, runParams.bias_sstskin_due_rain_intensity, self.data["rain"].fdata, runParams.bias_sstskin_due_rain_wind, self.data["windu10"].fdata, n)
        
         # quality filtering of wind and Hs data
        for i in arange(n):
           if (self.data["windu10"].fdata[i] != missing_value):
               # valid range taken from GlobWave Product User Guide Phase 3 (PUG3) doucment
              if (self.data["windu10"].fdata[i] > 50.0) or (self.data["windu10"].fdata[i] < 0.0):
                 self.data["windu10"].fdata[i] = missing_value
           if "sig_wv_ht" in self.data and self.data["sig_wv_ht"].fdata[i] != missing_value:
               # valid range taken from GlobWave Product User Guide Phase 3 (PUG3) doucment
              if (self.data["sig_wv_ht"].fdata[i] > 20.0) or (self.data["sig_wv_ht"].fdata[i] < 0.0):
                 self.data["sig_wv_ht"].fdata[i] = missing_value
           #if (self.data["sigma0"].fdata[i] != missing_value):
           #   if (self.data["sigma0"].fdata[i] > ??.0):
            #     self.data["sigma0"].fdata[i] = missing_value
        
        if (runParams.pco2_data_selection != 0 and VERIFICATION_RUNS != True):
           # signifies that we're using SOCAT data or in-situ data
           #print "%s Using the SOCAT data " % (function)
           #if self.data["pco2_sst"].fdata[i] != missing_value:
           for i in arange(n):
              if isnan(self.data["pco2_sst"].fdata[i]) != True:# and self.data["pco2_sst"].fdata[i] > 0.0 ): #SOCATv4_IGA
                 if self.data["pco2_sst"].fdata[i] > 260:
                   self.data["pco2_sst"].fdata[i] = self.data["pco2_sst"].fdata[i] - 273.15#IGA - If statement added because in-situ SST data may not be in K!
              
                 if self.data["pco2_sst"].fdata[i] > 30.5 or self.data["pco2_sst"].fdata[i] < -1.8:
                    self.data["pco2_sst"].fdata[i] = missing_value
              else:
                 self.data["pco2_sst"].fdata[i] = missing_value
        
        # quality control/contrain all SST data
        # check all SST data are within -1.8 - 30.5^oC (or 271.35 - 303.65K)
        for i in arange(n):
           if ((self.data["sstskin"].fdata[i] != missing_value) and (self.data["sstskinC"].fdata[i] != missing_value) and (self.data["sstfnd"].fdata[i] != missing_value) and (self.data["sstfndC"].fdata[i] != missing_value)):
              if ( (self.data["sstskinC"].fdata[i] > 30.5) or (self.data["sstfndC"].fdata[i] > 30.5) or (self.data["sstskinC"].fdata[i] < -1.8) or (self.data["sstfndC"].fdata[i] < -1.8)):
                 self.data["sstfnd"].fdata[i] = missing_value
                 self.data["sstskin"].fdata[i] = missing_value
                 self.data["sstfndC"].fdata[i] = missing_value
                 self.data["sstskinC"].fdata[i] = missing_value
                 
        # ensure that the different SST data all cover the same spatial regions
        # convert all missing values into standard value, rather than variations that seem to exist in some of these data
        #note this is an intersect operation (compared to the above)
        for i in arange(n):
           if ((self.data["sstfnd"].fdata[i] == self.data["sstfnd"].fillValue) or (self.data["sstskin"].fdata[i] == self.data["sstskin"].fillValue)):
              self.data["sstfnd"].fdata[i] = missing_value
              self.data["sstskin"].fdata[i] = missing_value
        
        # ---------------------------------------------------------------------------------------------------------------------------------------
        # -                                                                                                                                     -
        # converting pressure data from Pascals to millibar #IGA Need to formalise pressure data units - this converts from Pa.-------------START
        if runParams.TAKAHASHI_DRIVER != True:
            if npany(self.data["pressure"].fdata > 10000.0):
                for i in arange(n):
                    if (self.data["pressure"].fdata[i] != missing_value):
                        self.data["pressure"].fdata[i] = self.data["pressure"].fdata[i] * 0.01
                print "Converted pressure data to mbar from Pa"
        # ------------------------------------------------------IGA_SOCATv4 - Was removed as using pressure data in mbar.-------------END
        # -                                                                                                                                     -
        # ---------------------------------------------------------------------------------------------------------------------------------------
        
        #########################################################
        # calculations for components of the flux calculation
        #########################################################
         
        # pCO2/fCO2 extrapolation (if turned on) from reference year # edited K.B. also added [i] t
        pco2_increment = (self.year_data - runParams.pco2_reference_year) * runParams.pco2_annual_correction;
        pco2_increment_air = pco2_increment;
            
        DeltaT_fdata = array([missing_value] * n)
        
        if runParams.flux_calc == 1:
           print "%s Using the RAPID model (from Woolf et al., 2016)" % (function)
        elif runParams.flux_calc == 2:
           print "%s Using the EQUILIBRIUM model (from Woolf et al., 2016)" % (function)
           for i in arange(n):
              if ( (self.data["sstfndC"].fdata[i] != missing_value) & (self.data["sstskinC"].fdata[i] != missing_value) ):
                 DeltaT_fdata[i] = self.data["sstfndC"].fdata[i] - self.data["sstskinC"].fdata[i]
              else:
                 DeltaT_fdata[i] = missing_value
        elif runParams.flux_calc == 3:
            print "%s Using the BULK model (ie Flux = k*S*delta_pCO2)" % (function)
            print "%s Note: this assumes that the SSTskin dataset is the only temperature dataset, and so this is what will be used to calculate k and the solubility" % (function)
        else:
           print "\n%s runParams.flux_calc from configuration not recognised, exiting." % (function)
           return 1;
        
        #Calculating the schmidt number at the skin and fnd
        if runParams.schmidt_parameterisation == "schmidt_Wanninkhof2014":
            self.data["scskin"].fdata = schmidt_Wanninkhof2014(self.data["sstskinC"].fdata, n, runParams.GAS)
            self.data["scfnd"].fdata = schmidt_Wanninkhof2014(self.data["sstfndC"].fdata, n, runParams.GAS)
        elif runParams.schmidt_parameterisation == "schmidt_Wanninkhof1992":
            self.data["scskin"].fdata = schmidt_Wanninkhof1992(self.data["sstskinC"].fdata, n, runParams.GAS)
            self.data["scfnd"].fdata = schmidt_Wanninkhof1992(self.data["sstfndC"].fdata, n, runParams.GAS)
        else:
            raise ValueError("Unrecognised schmidt/solubility parameterisation selected: "+runParams.schmidtParameterisation);
        
        #Calculate solubility
        if runParams.schmidt_parameterisation == "schmidt_Wanninkhof2014":
            #calculating the skin solubility, using skin sst and salinity
            self.data["solubility_skin"].fdata = solubility_Wanninkhof2014(self.data["sstskin"].fdata, self.data["salinity_skin"].fdata, DeltaT_fdata, n, True, runParams.GAS.lower());
            #calculating the interfacial solubility
            self.data["solubility_fnd"].fdata = solubility_Wanninkhof2014(self.data["sstfnd"].fdata, self.data["salinity"].fdata, DeltaT_fdata, n, runParams.flux_calc, runParams.GAS.lower());
        elif runParams.schmidt_parameterisation == "schmidt_Wanninkhof1992":
            #calculating the skin solubility, using skin sst and salinity
            self.data["solubility_skin"].fdata = solubility_Wanninkhof1992(self.data["sstskin"].fdata, self.data["salinity_skin"].fdata, DeltaT_fdata, n, True, runParams.GAS.lower());
            #calculating the interfacial solubility
            self.data["solubility_fnd"].fdata = solubility_Wanninkhof1992(self.data["sstfnd"].fdata, self.data["salinity"].fdata, DeltaT_fdata, n, runParams.flux_calc, runParams.GAS.lower());
        else:
            raise ValueError("Unrecognised schmidt/solubility parameterisation selected: "+runParams.schmidtParameterisation);
    
         # calculate pCO2 data using mean sea level pressure data
         # equation 26, Kettle et al, 2009, ACP
         # pCO2_air = X[CO2] ( P(t) - pH2O(t) )
         # where X[CO2] is from Takahashi climatology
         # pH2O(t) = 1013.25 exp (24.4543 - 67.4509 (100/Tk(t) - 4.8489 ln (Tk(t) / 100) - 0.000544 S)
         # S = salinity, Tk = temperature in Kelvin
         # pCO2_water = pCO2_water_tak exp (0.0423 (T_foundation - T_tak)
         # T_tak = Takahashi temperature
        
        sys.stdout.flush();

        
        ###############################
        # Calculate:                  #
        #   pH2O                      #
        ###############################
        for i in arange(n):
            #if ( (self.data["salinity_skin"].fdata[i] != missing_value) and (self.data["sstskin"].fdata[i] != missing_value) and (self.data["pressure"].fdata[i] != missing_value) and (self.data["vco2_air"].fdata[i] != missing_value) and (self.data["sstfnd"].fdata[i] != missing_value) and (self.data["pco2_sst"].fdata[i] != missing_value) and (self.data["pco2_sw"].fdata[i] != missing_value) and (self.data["sstskin"].fdata[i] !=0.0) ):
            if (self.data["salinity_skin"].fdata[i] != missing_value) and (self.data["sstskin"].fdata[i] != missing_value):
                self.data["pH2O"].fdata[i] = 1013.25 * exp(24.4543 - (67.4509 * (100.0/self.data["sstskin"].fdata[i])) - (4.8489 * log(self.data["sstskin"].fdata[i]/100.0)) - 0.000544 * self.data["salinity_skin"].fdata[i])
            else:
                self.data["pH2O"].fdata[i] = missing_value
                    
        #######################################################
        # Calculate corrected values for pco2_sw and pco2_air #
        # Only do this is pco2 data is provided               #
        #######################################################  
        # this may be needed when using SMOS salinity data
        # To-DO: awaiting info from Lonneke and David before implementing fully
        pCO2_salinity_term = 0;
        #dSalinity = salinity_rmse
        #if (salinity_option == 1):
        #  # need to determine dS/S using original salinity (its been modified above to add the salinity_rmse)
        #  # so (self.data["salinity"].fdata[i] - salinity_rmse) ensures that we are dealing with the original value of salinity
          #  # using Ys=1 as a global correction following Sarmiento and Gruber, 2006)
        #pCO2_salinity_term = 1.0*(dSalinity/(self.data["salinity"].fdata[i] - salinity_rmse) ) #TH: Commented this out, but should not be removed as may be used in the future.
        #else:
        # pCO2_salinity_term = 0.0
        
        #if pCO2 data in sea water is provided calculated corrected values.
        if "pco2_sw" in self.data: #Only calculate partial pressure data is available
            self.add_empty_data_layer("pco2_sw_cor");
            for i in range(len(self.data["pco2_sw"].fdata)):
                #If statement below to maintain a consistent calculation with previous versions. Perhaps not needed but would invalidate reference data otherwise.
                if ( (self.data["salinity_skin"].fdata[i] != missing_value) and (self.data["sstskin"].fdata[i] != missing_value) and (self.data["pressure"].fdata[i] != missing_value) and (self.data["sstfnd"].fdata[i] != missing_value) and (self.data["pco2_sst"].fdata[i] != missing_value) and (self.data["pco2_sw"].fdata[i] != missing_value) and (self.data["sstskin"].fdata[i] !=0.0) ):
                    if runParams.GAS == 'CO2' and runParams.pco2_data_selection != 3:
                        # correction to different years, correction is data and year specific.
                        # note for 2010, correction for SOCAT isn't strictly required. However the contents of the exponential will collapse
                        # to 1 (with some rounding error expected), so effectively no correction will be applied
                        self.data["pco2_sw_cor"].fdata[i] = pco2_increment[i] + (self.data["pco2_sw"].fdata[i] *exp( (0.0423*(self.data["sstfndC"].fdata[i] - self.data["pco2_sst"].fdata[i])) - (0.0000435*((self.data["sstfndC"].fdata[i]*self.data["sstfndC"].fdata[i]) - (self.data["pco2_sst"].fdata[i]*self.data["pco2_sst"].fdata[i]) )) + pCO2_salinity_term) );
                    else:
                        self.data["pco2_sw_cor"].fdata[i] = self.data["pco2_sw"].fdata[i];
        
        
        #if interface/air CO2 data is provided as molar fraction and there is no partial pressure then calculate partial pressure.
        if ("vco2_air" in self.data) and ("pco2_air" not in self.data):
            # vco2 in ppm * 1000000 = atm
            # result /1000 to then convert from atm to uatm
            # hence * 0.001 factor
            self.add_empty_data_layer("pco2_air");
            for i in range(len(self.data["vco2_air"].fdata)):
                #If statement below to maintain a consistent calculation with previous versions. Perhaps not needed but would invalidate reference data otherwise.
                if ( (self.data["salinity_skin"].fdata[i] != missing_value) and (self.data["sstskin"].fdata[i] != missing_value) and (self.data["pressure"].fdata[i] != missing_value) and (self.data["vco2_air"].fdata[i] != missing_value) and (self.data["sstfnd"].fdata[i] != missing_value) and (self.data["pco2_sst"].fdata[i] != missing_value) and (self.data["pco2_sw"].fdata[i] != missing_value) and (self.data["sstskin"].fdata[i] !=0.0) ):
                    if runParams.GAS == 'CO2' and runParams.ATMGAS == 'V':
                        self.data["pco2_air"].fdata[i] = (self.data["vco2_air"].fdata[i] * (self.data["pressure"].fdata[i] - self.data["pH2O"].fdata[i])) / 1013.25;
                    else:
                        self.data["pco2_air"].fdata[i] = self.data["vco2_air"].fdata[i]
        
        #Now calculate corrected values for pCO2 at the interface/air
        ###Converts from ppm to microatm THc
        self.add_empty_data_layer("pco2_air_cor");
        #If statement added below to maintain a consistent calculation with previous versions. Perhaps not needed but would invalidate reference data otherwise.
        if runParams.TAKAHASHI_DRIVER == False: #Different for takahashi run to maintain compatability with verification run. This will be updated when verification runs are updated
            for i in range(len(self.data["pco2_air"].fdata)):
                if ( (self.data["salinity_skin"].fdata[i] != missing_value) and (self.data["sstskin"].fdata[i] != missing_value) and (self.data["pressure"].fdata[i] != missing_value) and (self.data["sstfnd"].fdata[i] != missing_value) and (self.data["pco2_sst"].fdata[i] != missing_value) and (self.data["pco2_sw"].fdata[i] != missing_value) and (self.data["sstskin"].fdata[i] !=0.0) ):
                    if runParams.GAS == 'CO2' and runParams.ATMGAS == 'V':
                        ##THtodo: 1e-6 can be removed...
                        self.data["pco2_air_cor"].fdata[i] = self.data["pco2_air"].fdata[i] + pco2_increment_air[i];
                    else:
                        self.data["pco2_air_cor"].fdata[i] = self.data["pco2_air"].fdata[i]
        else: #runParams.TAKAHASHI_DRIVER==True #Added to maintain compatability with takahashi verification. Should be removed when the verification runs are updated.
            for i in range(len(self.data["vco2_air"].fdata)):
                #If statement below to maintain a consistent calculation with previous versions. Perhaps not needed but would invalidate reference data otherwise.
                if ( (self.data["salinity_skin"].fdata[i] != missing_value) and (self.data["sstskin"].fdata[i] != missing_value) and (self.data["pressure"].fdata[i] != missing_value) and (self.data["vco2_air"].fdata[i] != missing_value) and (self.data["sstfnd"].fdata[i] != missing_value) and (self.data["pco2_sst"].fdata[i] != missing_value) and (self.data["pco2_sw"].fdata[i] != missing_value) and (self.data["sstskin"].fdata[i] !=0.0) ):
                    if runParams.GAS == 'CO2' and runParams.ATMGAS == 'V':
                        #THtodo: 1e-6 can be removed...
                        self.data["pco2_air_cor"].fdata[i] = (self.data["vco2_air"].fdata[i] * 1e-6 * (self.data["pressure"].fdata[i] - self.data["pH2O"].fdata[i]) / (1e-6 * 1013.25)) + (pco2_increment_air[i])
                    else:
                        self.data["pco2_air_cor"].fdata[i] = self.data["vco2_air"].fdata[i]

                    
            
            #SOCAT, so: conversion of pCO2 to fCO2 from McGillis and Wanninkhof 2006, Marine chemistry with correction from Weiss 1974 (as the equation in 2006 paper has a set of brackets missing)
            #runParams.pco2_data_selection ==2 signifies SOCAT fCO2 data, so converting pCO2_air_cor_fdata to fCO2_air_cor_fdata      
            if runParams.pco2_data_selection == 2 or runParams.pco2_data_selection == 4 or runParams.pco2_data_selection == 45:
                b11_fdata = array([missing_value] * n);
                d12_fdata = array([missing_value] * n);
                for i in range(len(self.data["pco2_air_cor"].fdata)):
                    #If statement below to maintain a consistent calculation with previous versions. Perhaps not needed but would invalidate reference data otherwise.
                    if ( (self.data["salinity_skin"].fdata[i] != missing_value) and (self.data["sstskin"].fdata[i] != missing_value) and (self.data["pressure"].fdata[i] != missing_value) and (self.data["vco2_air"].fdata[i] != missing_value) and (self.data["sstfnd"].fdata[i] != missing_value) and (self.data["pco2_sst"].fdata[i] != missing_value) and (self.data["pco2_sw"].fdata[i] != missing_value) and (self.data["sstskin"].fdata[i] !=0.0) ):
                        # conversion of pCO2 to fCO2 from McGillis and Wanninkhof 2006, Marine chemistry with correction from Weiss 1974 (as the equation in 2006 paper has a set of brackets missing)
                        b11_fdata[i] = -1636.75 + (12.0408*self.data["sstskin"].fdata[i]) - (0.0327957*self.data["sstskin"].fdata[i]*self.data["sstskin"].fdata[i]) + (3.16528e-5 * self.data["sstskin"].fdata[i]*self.data["sstskin"].fdata[i]*self.data["sstskin"].fdata[i])
                        d12_fdata[i] = 57.7 - (0.118*self.data["sstskin"].fdata[i])
                         # gas constant
                        R = 82.0578 # in [cm^3 atm/(mol K)]
                        # 1/0.987 = 1.0131712 - conversion between bar and atm, so *1013.25 is the conversion from millibar to atm.
                        # the combination of the B11 and d12 terms are in cm^3/mol and so these cancel with the P/RT term (in mol/cm^3) so the whole of the exp term is dimensionless
                        self.data["pco2_air_cor"].fdata[i] = self.data["pco2_air_cor"].fdata[i] * exp((b11_fdata[i] + (2*d12_fdata[i]) ) * 1e-6 * ((self.data["pressure"].fdata[i] * 1013.25)/(R*self.data["sstskin"].fdata[i]) ))

        
        ######################################
        # Takahashi verification information #  ##TODO: CHECK IF THIS CAN BE REMOVED?
        ######################################
        if runParams.TAKAHASHI_DRIVER: #Assumes CO2 data input is not suppled as concentrations
            # debuggin differences in pH20 values
            pCO2a_diff_fdata = array([missing_value] * n)
            dpCO2_diff_fdata = array([missing_value] * n)
            for i in arange(n):
                #Additional pCO2 outputs for Takahashi verification
                if self.data["pco2_air"].fdata[i] != missing_value:
                    pCO2a_diff_fdata[i] = self.data["pco2_air_cor"].fdata[i] - self.data["pco2_air"].fdata[i]
                    dpCO2_diff_fdata[i] = (self.data["pco2_sw_cor"].fdata[i] - self.data["pco2_air_cor"].fdata[i]) - (self.data["pco2_sw_cor"].fdata[i] - self.data["pco2_air"].fdata[i])
            
            pH2O_takahashi_fdata = array([missing_value] * n)
            humidity_fdata = array([missing_value] * n)
            pH2O_diff_fdata = array([missing_value] * n)
            for i in arange(n):
                #Additional humidity outputs for Takahashi verification
                if self.data["pco2_air"].fdata[i] != missing_value and self.data["pH2O"].fdata[i] != missing_value and self.data["pressure"].fdata[i] != missing_value and self.data["vco2_air"].fdata[i] != missing_value:
                    pH2O_takahashi_fdata[i] = self.data["pressure"].fdata[i] -  (self.data["pco2_air"].fdata[i] *1e-6 * 1013.25) / (self.data["vco2_air"].fdata[i] * 1e-6)
                    humidity_fdata[i] = pH2O_takahashi_fdata[i]/self.data["pH2O"].fdata[i];     
                    pH2O_diff_fdata[i] = ((humidity_fdata[i])-1.0) * 100.0;    
                else:
                    pH2O_takahashi_fdata[i] = missing_value;
                    humidity_fdata[i] = missing_value;
                    pH2O_diff_fdata[i] = missing_value;



        #######################################
        # Calculating gas transfer velocity k #
        #######################################
        for kParameterisationFunctor in self.kParameterisationFunctors:
            #Check each input exists #TODO: This should go in the pre-run checks!
            for inputDataName in kParameterisationFunctor.input_names():
                if inputDataName not in self.data: #This additional check isn't really needed as it is done in the functor and in the driver script.
                    raise KeyError("Selected kParameterisation ("+kParameterisationFunctor.name+") requires input data layers which have not been provided. Required the following DataLayers:\n"+str(kParameterisationFunctor.input_names())+"\nmissing:\n"+str(inputDataName));
            
            #Before running it is necessary to create any non-existing output layers that are required by the k calculation
            for outputDataName in kParameterisationFunctor.output_names():
                if outputDataName not in self.data:
                    self.add_empty_data_layer(outputDataName);
            
            #Call the functor to calculate k
            kParameterisationOutput = kParameterisationFunctor(self.data);
            if kParameterisationOutput == False:
                raise RuntimeError("%s: k parameterisation component (%s) returned False indicating k has not been calculated successfully." % (function, kParameterisationFunctor.name));
        

         # ability to investigate bias on k due to surface biology/slicks
         # assumes that bias values are realistic and that they won't cause the k_fdata to become unrealistic
         # also assumes that self.data["biology"].fdata exist (ie it won't if they indicator layers are off)
        if (runParams.bias_k_switch == 1):
           add_bias_k_biology_wind(self.data["k"].fdata, runParams.bias_k_value, self.data["biology"].fdata, runParams.bias_k_biology_value, self.data["windu10"].fdata, runParams.bias_k_wind_value, n, runParams.bias_k_percent_switch)
           if runParams.bias_k_percent_switch==0:
              print "\n%s Adding bias to chosen k (k_fdata) parameterisation data (bias value of %lf ms^-1 added, where biology (biology fdata) is > %lf mg m^-3 and wind speed (windu10) is < %lf m s^-1)" % (function, runParams.bias_k_value, runParams.bias_k_biology_value, runParams.bias_k_wind_value)
           else:
              print "\n%s Adding bias to chosen k (k_fdata) parameterisation data (bias value of - %lf percent being used, where biology (biology fdata) is > %lf mg m^-3 and wind speed (windu10) is < %lf m s^-1)" % (function, runParams.bias_k_value, runParams.bias_k_biology_value, runParams.bias_k_wind_value)


         ############################
         # actual flux calculation  #
         ############################
         # determine flux based on kHO6_fdata derivation
         # CO2 flux = k * s * Delta_PCO2
         # s = salinity, Delta_PCO2 calculated above
          # solubility conversion between mol kg^-1 atm^-1 to g-C m^-3 uatm^-1
          # mol -> grams for CO2 x 12
          # kg=liter -> m^-3 = x 1000
          # atm -> uatm = /1000000
          # result = x12.0108/1000
          # provides solubility in g-C m^-3 uatm^-1
          
          # multiplying solubility (g-C m^-3 uatm^-1) by pCO2 (uatm) = concentration in g-C m^-3
          
          # kg to m^-3 = 1 liter = 1/1000 m^-3
          # 1 liter = 1 kg
          
          # k conversion from m/s to to m/day
          # k are already in cm/h)
          # x 24 (converts from hours to days)
          # x 10^-2 (converts from cm to metres)
          # 24 * 10**{-2} = *24/100
        
          # this wil give a daily flux in g-C m^-2 day^-1
          # expanding equations
        
        #k_factor = 36.0 * 24.0 / 100.0 # used if using 10^-4 ms units for k, whereas we are now using cm/h for k
        k_factor = 24.0 / 100.0
        
        #Somewhere to store the net flux
        self.add_empty_data_layer("FH06");
        #Update FH06 (OF/gas flux data layer) description to reflect the gas transfer velocity calculation used in the calculation.
        self.data["FH06"].longName = self.data["FH06"].longName % self.data["k"].name;

        #If using rain wet deposition, calculate the gas solubility in distilled water
        if runParams.rain_wet_deposition_switch:
            self.add_empty_data_layer("FKo07");
            self.add_empty_data_layer("solubility_distilled");
            calculate_solubility_distilled(self.data["solubility_distilled"].fdata, self.data["salinity"].fdata,
                                           runParams.rain_wet_deposition_switch, self.data["sstskin"], DeltaT_fdata, self.n);
        
        if ((runParams.kb_asymmetry != 1.0) and (runParams.k_parameterisation == 3)):
           print "%s kb asymetry has been enabled (runParams.kb_asymmetry:%lf and runParams.k_parameterisation:%d)" % (function, runParams.kb_asymmetry, runParams.k_parameterisation)
           if runParams.flux_calc == 3:
               raise ValueError("kb_asymmetry is not supported for the 'bulk' calculation. Try using 'rapid' or 'equilibrium' instead.");

        ###################################################
        # If concentration data are not provided as input #
        #    calculate them from corrected pco2 data      #
        ###################################################
        concFactor = (12.0108/1000.0);
        if "concw" not in self.data:
            self.add_empty_data_layer("concw");
            if runParams.flux_calc == 3: #Bulk calculation, so should use the same solubility as conca
                calculate_concw(concFactor, self.data["solubility_skin"].fdata, self.data["pco2_sw_cor"].fdata, self.data["concw"].fdata); #calculate concw
            else: #Not bulk, so use the solubility at foundation layer
                calculate_concw(concFactor, self.data["solubility_fnd"].fdata, self.data["pco2_sw_cor"].fdata, self.data["concw"].fdata); #calculate concw
        if "conca" not in self.data:
            self.add_empty_data_layer("conca");
            calculate_conca(concFactor, self.data["solubility_skin"].fdata, self.data["pco2_air_cor"].fdata, self.data["conca"].fdata); #calculate conca
        
        
        ##############################
        # Main flux calculation loop # #assume corrected pco2 data at the moment #################################
        ##############################
        for i in arange(n):
            if ( (self.data["k"].fdata[i] != missing_value) and (self.data["concw"].fdata[i] != missing_value) and (self.data["conca"].fdata[i] != missing_value) ):
                #flux calculation
                #TODO: This is partially K-parameterisation dependent. Need to decouple this...
                
                #Rapid and equilibrium flux calculations
                if runParams.flux_calc == 1 or runParams.flux_calc == 2:
                    if ((runParams.kb_asymmetry != 1.0) and (runParams.k_parameterisation == 3)):
                        kd_component = self.data["kd"].fdata[i] * k_factor
                        kb_component = self.data["kb"].fdata[i] * k_factor
                        self.data["FH06"].fdata[i] = (kd_component * (self.data["concw"].fdata[i] - self.data["conca"].fdata[i])) + (kb_component * (self.data["concw"].fdata[i] - (runParams.kb_asymmetry *self.data["conca"].fdata[i]) ) )
                    else:
                        self.data["FH06"].fdata[i] = (self.data["k"].fdata[i] * k_factor) * (self.data["concw"].fdata[i] - self.data["conca"].fdata[i])

                # using simplified flux calculation with no separation of temperature between airside and waterside CO2
                # assumes that the skin temperature dataset is the only temperature dataset
                #Bulk model (F = k*solubility*(pCO2_water - pCO2_air)
                elif runParams.flux_calc == 3:
                    #Note that that concw is calculated differently if flux_calc ==3 vs !=3, uses skin solubility (same as conca) if flux_calc==3.
                    #The below calculation is therefore not identical to RAPID, and is the correct BULK forumulation.
                    self.data["FH06"].fdata[i] = self.data["k"].fdata[i] * k_factor * (self.data["concw"].fdata[i] - self.data["conca"].fdata[i]);
                
                #Some unexpected flux calculation was specified. Should never get this far as fe_setup.py will check this.
                else:
                    raise ValueError("Unrecognised flux calculation. Currently supported options are 'rapid' (1), 'equilibrium' (2) and 'bulk' (3). Recieved: "+str(runParams.flux_calc));
                
                
                # calculating and adding in flux component for wet deposition due to rain
                if runParams.rain_wet_deposition_switch:
                    if "pco2_air_cor" not in self.data: #Check that pCO2 / vCO2 data were supplied
                        raise ValueError("Cannot use wet deposition (rain_wet_deposition_switch) without specifying pCO2 or vCO2 data.");
                    else:
                        # relationship from Komori et al.m 2007
                        # need solubility of CO2 in fresh water
                        # solubility calculated using sstkin
                        # 24.0/1000.0 = conversion from mm hr^-1 to m day^-1
                        # flux is then in g C m^-2 day^-1
                        # flux is always negative, ie going into the ocean
                        self.data["FKo07"].fdata[i] = -(self.data["rain"].fdata[i] * (24.0/1000.0)) * (concFactor * self.data["solubility_distilled"].fdata[i]) * self.data["pco2_air_cor"].fdata[i]
                        if self.data["FH06"].fdata[i] != missing_value:
                            self.data["FH06"].fdata[i] += self.data["FKo07"].fdata[i]
                        else:
                            self.data["FH06"].fdata[i] = self.data["FKo07"].fdata[i]
            else:
                self.data["FH06"].fdata[i] = missing_value
        
        
        #Adding verification data outputs at same units as T09
        if runParams.TAKAHASHI_DRIVER == True:
          solskin_takadata = array([missing_value] * n)
          FH06_takadata = array([missing_value] * n)
          for i in arange(n):
           if self.data["solubility_skin"].fdata[i] != missing_value:
              solskin_takadata[i] = self.data["solubility_skin"].fdata[i]*1000 #from 'mol kg-1 atm-1' to 'mmol kg-1 atm-1'
           if self.data["FH06"].fdata[i] != missing_value:
              FH06_takadata[i] = self.data["FH06"].fdata[i]/30.5 #from flux per day to flux per month
        else:
          solskin_takadata = array([missing_value] * n)
          FH06_takadata = array([missing_value] * n)
        
        
        #
        # quality control of datasets, following TS (OceanFluxGHG_TS_D2-9_v1.8-signed.pdf) table 7
        #
        #calculate takahashi style DpCO2 and check range
        if "pco2_sw_cor" in self.data and "pco2_air_cor" in self.data:
            self.add_empty_data_layer("dpco2_cor");
            for i in arange(self.n):
               if ( (self.data["pco2_sw_cor"].fdata[i] != missing_value) and (self.data["pco2_air_cor"].fdata[i] != missing_value) ):
                  self.data["dpco2_cor"].fdata[i] = (self.data["pco2_sw_cor"].fdata[i] - self.data["pco2_air_cor"].fdata[i]) 
               else:
                  self.data["dpco2_cor"].fdata[i] = missing_value

        self.add_empty_data_layer("dpconc_cor");        
        for i in arange(self.n):
           if ( (self.data["concw"].fdata[i] != missing_value) and (self.data["conca"].fdata[i] != missing_value) ):
              self.data["dpconc_cor"].fdata[i] = (self.data["concw"].fdata[i] - self.data["conca"].fdata[i]) 
           else:
              self.data["dpconc_cor"].fdata[i] = missing_value

        #Calculate the total number of quality violations per grid location
        self.add_empty_data_layer("failed_quality")
        #self.d = {};
        for outputVar in ["FH06", "kt", "k", "kd", "kb", "salinity", "sstskinC", "concw", "conca", "dpco2_cor", "sstfndC", "pco2_sw_cor"]:
            if outputVar in self.data: #Only check outputs that we've actually created...
                check_output_dataset(self.data[outputVar], self.data["failed_quality"].fdata);
                #self.d[outputVar] = self.data["failed_quality"].fdata.copy();
                #self.d[outputVar].shape = self.data["failed_quality"].data.shape;

        #whitecapping data using relationship from TS and parameters from table 1 of Goddijn-Murphy et al., 2010, equation r1
        self.add_empty_data_layer("whitecap");
        self.data["whitecap"].fdata = calculate_whitecapping(self.data["windu10"], self.data["whitecap"]);
        
        
        #runParams.pco2_data_selection == 1 signifies that we're using SOCAT data 
        #this dataset is filled with missing_values prior to output, otherwise non-socat data will appear in the SFUG field in the netcdf
        #which will confuse the user
        #enabled for TAKAHASHI_DRIVER to enable checking
        if runParams.TAKAHASHI_DRIVER != True:
           if (self.runParams.pco2_data_selection == 0):
              self.data["pco2_sw"].fdata = array([missing_value] * self.n)
              self.data["pco2_sw_stddev"].fdata = array([missing_value] * self.n)
        
        #
        #procesing indictor attribute layers
        #
        #temp if clause using the -l flag for now #TODO: remove this, config should just specify required process indicator layers
        if self.runParams.use_sstfnd_switch == True:
            for piFunctor in self.processIndicatorFunctors:
                try:
                    for inputDataLayer in piFunctor.input_names(): #Check all required inputs exist.
                        if inputDataLayer not in self.data:
                            raise ValueError("%s: Missing input DataLayer (%s) for process indicator functor %s." % (function, inputDataLayer, piFunctor.name));
    
                    for outputDataLayer in piFunctor.output_names(): #Add any output DataLayers that don't already exist
                        if outputDataLayer not in self.data:
                            self.add_empty_data_layer(outputDataLayer);
                    #Execute the process indicator layer functor.
                    piFunctor(self.data);
                except ValueError as e:
                    print e.args;
                    print "Exiting...";
                    return 1;
        
        #
        # write out results
        #        
        #Temp refactoring:
        #TODO: Remove this and test.
#        for name in ["krain", "kt", "kb", "kd"]:
#            if name not in self.data:
#                self.add_empty_data_layer(name);
        
        
        #write out the final ouput to netcdf
        write_txt(self);        
        print "%s SUCCESS writing file %s" % (function, runParams.output_path)
#        
#        #Finally, close the logger
        handlers = self.logger.handlers[:];
        for handler in handlers:
            handler.close();
        
        sys.stdout.flush();
        return 0;

