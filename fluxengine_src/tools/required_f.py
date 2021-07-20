import os
import numpy
import numpy.lib.recfunctions
import datetime
import argparse
import netCDF4
import glob
import multiprocessing
import pandas as pd;

#Returns standard column names, dtypes and columns as string names or string column numbers   
def construct_column_info(year_col, month_col, day_col, hour_col, minute_col, second_col, longitude_col, latitude_col, \
                                   salinity_col, salinity_sub_col, SST_C_col, Tequ_col, air_pressure_col, air_pressure_sub_col, air_pressure_equ_col, \
                                   fCO2_col, expocode_col, socatversion, notsocatformat):
    stndColNames = ["expocode", "year", "month", "day", "hour", "minute", "second", "longitude", "latitude", "salinity", "sst", "T_equ", \
                    "air_pressure", "air_pressure_equ", "salinity_sub", "air_pressure_sub", "fCO2", "fCO2_qc_flag"];
    colDTypes = ['U24', '<i8', '<i8', '<i8', '<i8', '<i8', '<i8', '<f8', '<f8', '<f8', '<f8', '<f8', '<f8', '<f8', '<f8', '<f8', '<f8', '<i8'];
    #colDTypes = [str, int, int, int, int, int, int, float, float, float, float, float, float, float, float, float, float, int];
    
    if notsocatformat: #Specify the columns as given by the command line parameters.
        colIdentifiers = [expocode_col, year_col, month_col, day_col, hour_col, minute_col, second_col, longitude_col, latitude_col, \
                                   salinity_col, SST_C_col, Tequ_col, air_pressure_col, air_pressure_equ_col, salinity_sub_col, air_pressure_sub_col, \
                                   fCO2_col, None];
    else: #using socat so use socatversion to determine correct columns
        if socatversion == 2:
            colIdentifiers = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '10', '11', '12', '13', '14', '15', '16', '20', '22'];
#            return ['SOCAT_DOI',
#             'QC_ID',
#             'yr',
#             'mon',
#             'day',
#             'hh',
#             'mm',
#             'ss',
#             'latitude [dec.deg.N]',
#             'sample_depth [m]',
#             'sal',
#             'SST [deg.C]',
#             'Tequ [deg.C]',
#             'PPPP [hPa]',
#             'Pequ [hPa]',
#             'd2l [km]',
#             'fCO2rec [uatm]'];
        elif socatversion in [3, 4, 5, 6]:
            colIdentifiers = ['0', '3', '4', '5', '6', '7', '8', '9', '11', '12', '13', '14', '15', '16', '17', '18', '22', '24'];
            #colIdentifiers = ['0', '4', '5', '6', '7', '8', '9', '10', '11', '13', '18', '14', '15', '16', '19', '17', '23', '25'];
#            return ['yr',
#             'mon',
#             'day',
#             'hh',
#             'mm',
#             'ss',
#             'longitude [dec.deg.E]',
#             'latitude [dec.deg.N]',
#             'sal',
#             'SST [deg.C]',
#             'Tequ [deg.C]',
#             'PPPP [hPa]',
#             'Pequ [hPa]',
#             'WOA_SSS',
#             'NCEP_SLP [hPa]',
#             'fCO2rec [uatm]',
#             'fCO2rec_flag']
        else:
            raise ValueError("No value columns could be generated. Only socat version 2, 3, 4, 5, and 6 are supported.");
        
    columnInfo = [];
    for i in range(len(colIdentifiers)):
        columnInfo.append( (stndColNames[i], colDTypes[i], colIdentifiers[i]) );
    return columnInfo;


def convert_column_id_to_index(header, columnIdentifier):
    if columnIdentifier == None:
        return None;
    else:
        try:
            if int(columnIdentifier) <= len(header):
                index = int(columnIdentifier);
            else:
                raise IndexError("Index (%d) is out of range for header length of %d." % (int(columnIdentifier), len(header)));
        except ValueError:
            try:
                index = header.index(columnIdentifier);
            except ValueError:
                raise ValueError("column index '%s' could not be determined in header %s" % (columnIdentifier, str(header)));
    return index;

def ReadInData(inputfile, columnInfo, socatversion, delimiter='\t'):
    """
    Reads in the data from the SOCAT v2 or v3 files
       inputfile - the input SOCAT ascii csv file
       delimiter - the delimiter in the file (default to tabs)
    """
    with open(inputfile) as FILE:
       linestoskip=0;
       #Find how many lines we want to skip in the input SOCAT file by searching for 2 column headings
       if socatversion != None:
           for preline in FILE:
               if 'Expocode' not in preline or 'yr' not in preline:
                   linestoskip+=1
               else:
                   break;
           #Extract the header
           header = [s.strip() for s in preline.strip().split(delimiter)];
       else: #Using insitu, assume header is first line.
           header = [s.strip() for s in FILE.readline().split(delimiter)];
    

    columnInfoToExtract = [info for info in columnInfo if info[2] != None]; #These will be extracted from the datafile
    
    #Convert columns into indices
    indicesToExtract = [convert_column_id_to_index(header, info[2]) for info in columnInfoToExtract];
    namesOfExtracted = [info[0] for info in columnInfoToExtract];
    order = numpy.argsort(indicesToExtract); #pandas ignores the column order so we need to rearrange the column names accordingly
    namesOfExtracted = [namesOfExtracted[i] for i in order];
    
    print namesOfExtracted;
    
    #dtypesOfExtracted = [info[1] for info in columnInfoToExtract];

    #Read in the columns we want into a data array
    print "Reading in data from SOCAT file: %s"%inputfile
    print "This can take a while for large datasets.\n";

    data = pd.read_table(inputfile, skiprows=linestoskip+1, sep=delimiter, engine='c', usecols=indicesToExtract, names=namesOfExtracted, low_memory=False);#, dtype=dtypesOfExtracted);
    
    #Now insert additional columns filled with nan
    colNamesToInsert = [info[0] for info in columnInfo if info[2] == None];
    toInsert = numpy.full((len(data), len(colNamesToInsert)), numpy.nan);
    toInsert = pd.DataFrame(toInsert, columns=colNamesToInsert);
    data = data.join(toInsert);
    
    #Reorder columns
    orderedColNames = [info[0] for info in columnInfo];
    data = data[orderedColNames];

    #data.columnInfo = columnInfo; #attach metadata to the data
    data = data.to_records(index=False);
    if "fCO2_qc_flag" in colNamesToInsert: #nan values must be float type, fCO2_qc_flag is usually int, so have to split here.
        data = data.astype([('expocode', 'S24'), ('year','<i8'), ('month','<i8'), ('day','<i8'), ('hour','<i8'), ('minute','<i8'), ('second','<i8'),
                     ('longitude','<f8'), ('latitude','<f8'), ('salinity', '<f8'),
                    ('SST', '<f8'), ('T_equ', '<f8'), ('air_pressure', '<f8'), ('air_pressure_equ', '<f8'),
                    ('salinity_sub', '<f8'), ('air_pressure_sub', '<f8'), ('fCO2', '<f8'), ('fCO2_qc_flag', '<f8')]);
    else: #fCO2_qc_flag was provided by the data file, so set type as int.
        data = data.astype([('expocode', 'S24'), ('year','<i8'), ('month','<i8'), ('day','<i8'), ('hour','<i8'), ('minute','<i8'), ('second','<i8'),
                         ('longitude','<f8'), ('latitude','<f8'), ('salinity', '<f8'),
                        ('SST', '<f8'), ('T_equ', '<f8'), ('air_pressure', '<f8'), ('air_pressure_equ', '<f8'),
                        ('salinity_sub', '<f8'), ('air_pressure_sub', '<f8'), ('fCO2', '<f8'), ('fCO2_qc_flag', '<i8')]);
    
    return data;
