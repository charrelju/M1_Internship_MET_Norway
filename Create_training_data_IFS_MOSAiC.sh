#!/bin/bash -f
#$ -N Training_data_IFS_MOSAiC
#$ -l h_rt=00:15:00
#$ -S /bin/bash
#$ -pe shmem-1 1
#$ -l h_rss=4G,mem_free=4G
#$ -q research-r8.q
#$ -t 1-12
##$ -j y
##$ -m ba
#$ -o /home/justinec/Documents/OUT/OUT_$JOB_NAME.$JOB_ID_$TASK_ID
#$ -e /home/justinec/Documents/ERR/ERR_$JOB_NAME.$JOB_ID_$TASK_ID
##$ -R y
##$ -r y

source /modules/rhel8/conda/install/etc/profile.d/conda.sh
conda activate /lustre/storeB/users/justinec/myconda

echo "Got $NSLOTS slots for job $SGE_TASK_ID."

cat > "/home/justinec/Documents/PROG/Training_data_IFS_MOSAiC_""$SGE_TASK_ID"".py" << EOF

################################################################################################

#!/usr/bin/env python
# coding: utf-8

# This script collocate the data of the MOSAiC campaign's radiosoundings and of the IFS's operational forecasts for one radiosounding. 
# 
# It can plot temperature, relative humidity, u, v, wind direction and wind speed of both. 
# 
# This scrip is used in a bash script to reproduce figures and collocated data for multiple radiosoundings.

# ### IMPORT

# In[11]:


import numpy as np
from netCDF4 import Dataset as ncfile
from netCDF4 import num2date
import matplotlib.dates as dates
import datetime
import os
import itertools
from itertools import chain
import glob


# ### SELECT LEAD TIME AND LIST ALL IFS FILES

# In[12]:


lead_time = [0,6,12,18,24,30,36,42,48,54,60,66]
task_ID = lead_time[$SGE_TASK_ID-1] #$SGE_TASK_ID
#print('lead time :', task_ID)
task_ID_str = str(task_ID)
if len(task_ID_str) < 2 :
    task_ID_str = '0'+task_ID_str


# In[13]:


paths = []
path_input = '/lustre/storeB/users/justinec/master_internship/data/IFS/op_10m_wind_T2m_RH_PL/'
for year in range(2019, 2021):
    yearstr = str(year)
    for month in range(1,13):
        monthstr = "{:02d}".format(month)
        p = path_input + yearstr + '/' + monthstr + '/'
        paths.append(p)


# In[14]:


dataset = []
for path, dirs, files in chain.from_iterable(os.walk(path) for path in paths):
        nh_files = [s for s in files if "ECMWF_op_PL_" in s]
        sorted_files=sorted(nh_files)
        for i in range(0,len(sorted_files)):
                dataset.append(sorted(glob.glob(path+sorted_files[i])))


# ### FUNCTIONS

# In[15]:


"""Function that convert times of a dataset depemding on its unit and calendar 
type to a datetime then to Matplotlib dates.
---exemple---
IFS :  'hours since 2020-1-1 00:00:00' become datetime of type 2020-01-01 00:00:00 then 18262
MOSAiC : 'seconds since 2020-01-01T05:44:27.397Z' become datetime of type 2020-01-01 05:44:27.397000 then 18262.23920598
"""

def convert_time(dataset) :
    dataset_time = dataset.variables['time']
    t_unit = dataset_time.units
    t_cal = dataset_time.calendar
    dataset_nctime=[]; 
    dataset_nctime.append(num2date(dataset_time,units = t_unit,calendar = t_cal, only_use_python_datetimes=True, only_use_cftime_datetimes=False)) #datetime of type 2019-10-01 00:00:00
    #dates.date2num --> convert datetime objects to Matplotlib dates (better for compare ERA5 and MOSAiC times and for figures)
    #np.squeeze --> pass Matplotlib dates in column instead of in line
    dataset_nctimenum = np.squeeze(dates.date2num(dataset_nctime))
    return dataset_nctime,dataset_nctimenum,t_unit
# -----------------------------------------------------------------------------------------

# -----------------------------------------------------------------------------------------
"""Function that finds the indexes of the closest values between IFS and MOSAiC.
   Useful for:
- keep the latitudes, longitudes and times of IFS corresponding to MOSAiC
- keep the MOSAiC pressure levels corresponding to IFS"""
def corresponding_index(ifs,mosaic) :
    diff = np.absolute(ifs-mosaic)
    indx = diff.argmin()
    return indx
# -----------------------------------------------------------------------------------------

# -----------------------------------------------------------------------------------------
"""Function that finds the values of t, rh, u and v corresponding to the desired pressure levels.
It add a value when the pressure level isn't available in IFS.
This value is the mean between the two closest surrounding pressure levels.
It delete a value when the pressure level in IFS isn't in the desire pressure levels."""
def interpolate_IFS_MOSAiC() :
    ifs_temp = [] ; ifs_rh = [] ; ifs_u = [] ; ifs_v = [] ; err=0
    for i in range(0, len(indx_lat)) :
        
        if ifs_pres_list_cut[i] not in desired_levels :
            del(ifs_pres_list_cut[i])
            err=err+1
            ind = ifs_pres_list_cut.index(desired_levels[i])+err
            t = ifs.variables['t'][indx_time][ind][indx_lat[i]][indx_lon[i]]
            ifs_temp.append(t)
            rh = ifs.variables['r'][indx_time][ind][indx_lat[i]][indx_lon[i]]
            ifs_rh.append(rh)
            u = ifs.variables['u'][indx_time][ind][indx_lat[i]][indx_lon[i]]
            ifs_u.append(u)
            v = ifs.variables['v'][indx_time][ind][indx_lat[i]][indx_lon[i]]
            ifs_v.append(v)
        
        elif desired_levels[i] in ifs_pres_list_cut :
            ind = ifs_pres_list_cut.index(desired_levels[i])+err
            t = ifs.variables['t'][indx_time][ind][indx_lat[i]][indx_lon[i]]
            ifs_temp.append(t)
            rh = ifs.variables['r'][indx_time][ind][indx_lat[i]][indx_lon[i]]
            ifs_rh.append(rh)
            u = ifs.variables['u'][indx_time][ind][indx_lat[i]][indx_lon[i]]
            ifs_u.append(u)
            v = ifs.variables['v'][indx_time][ind][indx_lat[i]][indx_lon[i]]
            ifs_v.append(v)

        else :
            ifs_pres_list_cut.insert(i,desired_levels[i])
            ind = ifs_pres_list_cut.index(desired_levels[i])+err
            t = ifs.variables['t'][indx_time][ind][indx_lat[i+1]][indx_lon[i+1]]
            ifs_temp.insert(i,(np.mean([ifs_temp[i-1],t])))
            rh = ifs.variables['r'][indx_time][ind][indx_lat[i+1]][indx_lon[i+1]]
            ifs_rh.insert(i,(np.mean([ifs_rh[i-1],rh])))
            u = ifs.variables['u'][indx_time][ind][indx_lat[i+1]][indx_lon[i+1]]
            ifs_u.insert(i,(np.mean([ifs_u[i-1],u])))
            v = ifs.variables['v'][indx_time][ind][indx_lat[i+1]][indx_lon[i+1]]
            ifs_v.insert(i,(np.mean([ifs_v[i-1],v])))
            err=err-1
    
    return ifs_temp, ifs_rh, ifs_u, ifs_v
# -----------------------------------------------------------------------------------------

# -----------------------------------------------------------------------------------------
"""Function that plot the profiles of temperature, relative humidity, zonal wind u, meridional wind v, wind direction and wind speed of ERA5 and MOSAiC"""
def plotprofiles_IFS_MOSAiC() :
    
    import matplotlib.gridspec as gridspec
    from matplotlib import pyplot as plt

    parameters = {'axes.labelsize':20, 'axes.titlesize':20, 'legend.fontsize':16, 'xtick.labelsize':20, 'ytick.labelsize':20, 
              'legend.title_fontsize':20, 'figure.titlesize':25}
    plt.rcParams.update(parameters)

    fig, axs = plt.subplots(1,6,figsize=(30,15))
    axs[0].plot(mosaic_temp,desired_levels, 'C0')
    axs[0].plot(ifs_temp,desired_levels, 'r')
    axs[0].invert_yaxis()
    axs[0].set_xlabel ('T [K]')
    axs[0].set_ylabel ('Pressure [hPa]')

    axs[1].plot(mosaic_rh, desired_levels, 'C0')
    axs[1].plot(ifs_rh,desired_levels, 'r')
    axs[1].invert_yaxis()
    axs[1].set_xlabel ('RH [%]')

    axs[2].plot(mosaic_u,desired_levels, 'C0')
    axs[2].plot(ifs_u,desired_levels, 'r')
    axs[2].invert_yaxis()
    axs[2].set_xlabel (r'u [$m.s^{-1}$]')

    axs[3].plot(mosaic_v,desired_levels, 'C0')
    axs[3].plot(ifs_v,desired_levels, 'r')
    axs[3].invert_yaxis()
    axs[3].set_xlabel (r'v [$m.s^{-1}$]')

    axs[4].plot(mosaic_wdir,desired_levels, 'C0')
    axs[4].plot(ifs_wdir,desired_levels, 'r')
    axs[4].set_xlabel (r'Wind direction [degree]')
    axs[4].invert_yaxis()

    axs[5].plot(mosaic_wspeed,desired_levels, 'C0', label='MOSAiC')
    axs[5].plot(ifs_wspeed,desired_levels, 'r', label='IFS')
    axs[5].invert_yaxis()
    axs[5].set_xlabel (r'Wind speed [$m.s^{-1}$]')
    
    ifs_start_time = 'IFS operational forecast start time : ' + start_time[6:8] + '/' + start_time[4:6] + '/' + start_time[0:4] + ' ' + start_time[8:10] + ':' + start_time[10:12] + '\n'
    #ifs_forecast_date = str(int(ifs_nctime[0][indx_time].strftime('%Y%m%d%H%M%S')))
    ifs_forecast = 'IFS operational forecast time : ' + day + '/' + month + '/' + year + ' ' + hour + ':' + forecast_time[10:12] + '\n'
    mosaic_obs = 'MOSAiC observation time : ' + day + '/' + month + '/' + year + ' ' + hour + ':00:00'
    fig.legend()
    fig.suptitle(ifs_start_time+ifs_forecast+mosaic_obs, size=16)
    fig.tight_layout()
    fig.subplots_adjust(top=0.90)
    path_fig = '/lustre/storeB/users/justinec/master_internship/figures/ifs_mosaic_profiles/'
    plt.savefig(path_fig+'/fig_profiles_tropo_'+start_time+'_'+task_ID_str+'.png')
# -----------------------------------------------------------------------------------------

# -----------------------------------------------------------------------------------------    
def ncfile_IFS_MOSAiC() :
    # Create a ncfile in a path
    path_output = '/lustre/storeB/users/justinec/master_internship/data/IFS_MOSAiC/'
    output_filename = path_output + '/collocated_IFS_MOSAiC_' + start_time + '_' + task_ID_str + '.nc'
    output_netcdf = ncfile(output_filename, 'w', format = 'NETCDF4')

    # Dimension
    pres = output_netcdf.createDimension('pres',len(desired_levels))

    # Variables
    pres = output_netcdf.createVariable('pres', 'd', ('pres'))
    temp_ifs = output_netcdf.createVariable('temp_ifs', 'd', ('pres'))
    temp_mosaic = output_netcdf.createVariable('temp_mosaic', 'd', ('pres'))
    rh_ifs = output_netcdf.createVariable('rh_ifs', 'd', ('pres'))
    rh_mosaic = output_netcdf.createVariable('rh_mosaic', 'd', ('pres'))
    wspeed_ifs = output_netcdf.createVariable('wspeed_ifs', 'd', ('pres'))
    wspeed_mosaic = output_netcdf.createVariable('wspeed_mosaic', 'd', ('pres'))
    wdir_ifs = output_netcdf.createVariable('wdir_ifs', 'd', ('pres'))
    wdir_mosaic = output_netcdf.createVariable('wdir_mosaic', 'd', ('pres'))
    lat_mosaic = output_netcdf.createVariable('lat_mosaic', 'd', ('pres'))
    lon_mosaic = output_netcdf.createVariable('lon_mosaic', 'd', ('pres'))
    time_mosaic = output_netcdf.createVariable('time_mosaic', 'd', ('pres'))
    time_forecast = output_netcdf.createVariable('forecast_time', 'd')
    time_start = output_netcdf.createVariable('start_time', 'd')
    
    # Information
    pres.units = 'millibars'
    pres.long_name = 'Pressure level'
    pres.standard_name = 'pressure_level'
    temp_ifs.units = 'K'
    temp_ifs.long_name = 'Temperature'
    temp_ifs.standard_name = 'air_temperature'
    temp_mosaic.units = 'K'
    temp_mosaic.long_name = 'Temperature'
    temp_mosaic.standard_name = 'air_temperature'
    rh_ifs.units = 'percent'
    rh_ifs.long_name = 'Relative Humidity'
    rh_ifs.standard_name = 'relative_humidity'
    rh_mosaic.units = 'percent'
    rh_mosaic.long_name = 'Relative Humidity'
    rh_mosaic.standard_name = 'relative_humidity'
    wspeed_ifs.units = 'm s-1'
    wspeed_ifs.long_name = 'Wind speed'
    wspeed_ifs.standard_name = 'wind_speed'
    wspeed_mosaic.units = 'm s-1'
    wspeed_mosaic.long_name = 'Wind speed'
    wspeed_mosaic.standard_name = 'wind_speed'
    wdir_ifs.units = 'degree'
    wdir_ifs.long_name = 'Wind direction'
    wdir_ifs.standard_name = 'wind_from_direction'
    wdir_ifs.comment = 'Wind direction with 0°:north, 90°:east, 180°:south, 270°:west'
    wdir_mosaic.units = 'degree'
    wdir_mosaic.long_name = 'Wind direction'
    wdir_mosaic.standard_name = 'wind_from_direction'
    wdir_mosaic.comment = 'Wind direction with 0°:north, 90°:east, 180°:south, 270°:west'
    lat_mosaic.units = 'degree_North'
    lat_mosaic.long_name = 'Latitude'
    lat_mosaic.standard_name = 'latitude'
    lon_mosaic.units = 'degree_East'
    lon_mosaic.long_name = 'Longitude'
    lon_mosaic.standard_name = 'longitude'
    time_mosaic.units = t_unit
    time_mosaic.long_name = 'Time'
    time_mosaic.standard_name = 'time'
    time_forecast.units = 'yyyymmddhhmm'
    time_forecast.long_name = 'Forecast time'
    time_forecast.standard_name = 'forecast_time'
    time_start.units = 'yyyymmddhhmm'
    time_start.long_name = 'Forecast start time'
    time_start.standard_name = 'start_time'
    
    # Assignment of variables
    pres[:] = desired_levels
    temp_ifs[:] = ifs_temp
    temp_mosaic[:] = mosaic_temp
    rh_ifs[:] = ifs_rh
    rh_mosaic[:] = mosaic_rh
    wspeed_ifs[:] = ifs_wspeed
    wspeed_mosaic[:] = mosaic_wspeed
    wdir_ifs[:] = ifs_wdir
    wdir_mosaic[:] = mosaic_wdir
    lat_mosaic[:] = mosaic_lat_collocated
    lon_mosaic[:] = mosaic_lon_collocated
    time_mosaic[:] = mosaic_time_collocated
    time_forecast[:] = forecast_time
    time_start[:] = start_time

    output_netcdf.close() #close the netcdf file


# ### BROWSES ALL IFS FILES, COLLOCATE IFS-MOSAiC FOR THE SELECTED LEAD TIME, PLOT FIGURES AND CREATE NCFILE FOR EACH TIME

# In[17]:


for i in range(0,len(dataset)) :

# ----- READ DATA
    start_time = str(dataset[i])[101:113]                                  #start time of IFS operational forescast 
    start_datetime = datetime.datetime.strptime(start_time, '%Y%m%d%H%M')
    forecast_time =  start_datetime + datetime.timedelta(hours = task_ID)
    forecast_time = forecast_time.strftime('%Y%m%d%H%M')                   #time of IFS operational forecast and of MOSAiC observations

    year=forecast_time[0:4] ; month=forecast_time[4:6] ; day=forecast_time[6:8] ; hour=forecast_time[8:10]
    
    ifs_link  = dataset[i][0]

    ppidir_mosaic = '/lustre/storeB/users/maltem/Arctic/MOSAiC/radiosondes/'+year+'/'+month+'/'
    mosaic_link  = ppidir_mosaic + 'PST-RS-01_2_RS41-GDP_001_'+year+month+day+'T'+hour+'0000_1-000-001.nc'
    if os.path.isfile(mosaic_link) == False:
        mosaic_link  = ppidir_mosaic + 'PST-RS-01_2_RS41-GDP_001_'+year+month+day+'T'+hour+'0000_1-000-002.nc'
    if os.path.isfile(mosaic_link) == False:
        mosaic_link  = ppidir_mosaic + 'PST-RS-01_2_RS41-GDP_001_'+year+month+day+'T'+hour+'0000_1-000-003.nc'
    try:
        ifs = ncfile(ifs_link,'r')          #dataset of radiosoundings of IFS
        mosaic = ncfile(mosaic_link,'r')    #dataset of radiosoundings of MOSAiC
    except FileNotFoundError as e:
        print(f"FileNotFoundError successfully handled\n"f"{e}")
        continue #as MOSAiC file is missing for this date, we pass to the next date
        
# ----- RADIOSONDE OPERATIONAL FORESCAST IFS -----
    ifs_lat = ifs.variables['lat'][:]
    ifs_lon = ifs.variables['lon'][:]
    ifs_pres = ifs.variables['plev'][:]*10e-3  #Pa in hPA
    ifs_time = ifs.variables['time'][:]

# ----- RADIOSONDE MOSAiC -----
    mosaic_lat = mosaic.variables['lat'][:]
    mosaic_lon = mosaic.variables['lon'][:]
    mosaic_pres = mosaic.variables['press'][:]
    mosaic_time = mosaic.variables['time'][:]

# ----- CONVERT TIME -----
    ifs_nctime,ifs_nctimenum, t_unit = convert_time(ifs)
    mosaic_nctime, mosaic_nctimenum, t_unit = convert_time(mosaic)

# ----- CORRESPONDING INDEXES OF LATITUDE, LONGITUDE, TIME AND OF PRESSURE LEVELS BETWEEN IFS AND MOSAIC -----
    desired_levels = list(range(300, 1025, 50)) #from 300 hPa to 1000 hPa by 50 hPa
    #desired_levels=np.flip(ifs_pres) #if we want the same levels as IFS simply decomment this line
    desired_levels = np.flip(desired_levels)

    # Considering a single MOSAiC time because radiosounding lasts 1h30. It's short enough for all radiosounding times to be close to the same IFS time [00:00, 6:00, 12:00, 18:00]
    indx_time = corresponding_index(ifs_nctimenum,mosaic_nctimenum[0])
    #print('ifs forecast time :', ifs_nctime[0][indx_time])
    #print('mosaic time :', mosaic_link[87:102])
    
    # Indexes of MOSAiC pressure levels closest to the desired pressure levels
    indx_level = []
    for i in range(0,len(desired_levels)) :
        indx_level.append(corresponding_index(mosaic_pres,desired_levels[i]))

    # Latitude and longitude of IFS closest to MOSAiC

    #considering the latitude and longitude of MOSAiC of each desired pressure level
    indx_lat = []
    for i in range(0,len(mosaic_lat[indx_level])) :
        indx_lat.append(corresponding_index(ifs_lat,mosaic_lat[indx_level][i]))
    indx_lon = []
    for i in range(0,len(mosaic_lon[indx_level])) :
        indx_lon.append(corresponding_index(ifs_lon,mosaic_lon[indx_level][i]))

    mosaic_lat_collocated = mosaic_lat[indx_level]    #keep mosaic's latitudes of each desired pressure level
    mosaic_lon_collocated = mosaic_lon[indx_level]    #keep mosaic's longitudes of each desired pressure level
    mosaic_time_collocated = mosaic_time[indx_level]  #keep mosaic's times of each desired pressure level

    ifs_lat_collocated = ifs_lat[indx_lat] #latitude closest to mosaic
    ifs_lon_collocated = ifs_lon[indx_lon] #longitude closest to mosaic

# ----- CORRESPONDING TEMPERATURE, RELATIVE HUMIDITY, ZONAL WIND U, MERIDIONAL WIND V, WIND DIRECTION AND WIND SPEED OF IFS AND MOSAIC -----
    
    # MOSAiC
    mosaic_temp = mosaic.variables['temp'][indx_level]
    mosaic_rh = mosaic.variables['rh'][indx_level]
    mosaic_u = mosaic.variables['wzon'][indx_level]
    mosaic_v = mosaic.variables['wmeri'][indx_level]
    mosaic_wdir = mosaic.variables['wdir'][indx_level]
    mosaic_wspeed = mosaic.variables['wspeed'][indx_level]

    # IFS
    #considering latitude and longitude of each pressure level
    ifs_pres_list=ifs_pres.tolist()
    ifs_pres_list_cut = ifs_pres_list[:ifs_pres_list.index(desired_levels[-1])+1]
    ifs_temp, ifs_rh, ifs_u, ifs_v = interpolate_IFS_MOSAiC()
    ifs_wdir = []
    ifs_wspeed = []
    for i in range(0, len(indx_lat)) :
        ifs_wdir.append(((180/np.pi) * np.arctan2(ifs_u[i], ifs_v[i]) + 180) % 360)
        ifs_wspeed.append(np.sqrt(ifs_u[i]**2+ifs_v[i]**2))
    
# ----- PLOT
    #plotprofiles_IFS_MOSAiC() #comment this line to not plot figure

# ----- CREATE COLLOCATED NCFILE
    ncfile_IFS_MOSAiC()


##################################################################################################

EOF
python3 "/home/justinec/Documents/PROG/Training_data_IFS_MOSAiC_""$SGE_TASK_ID"".py"



