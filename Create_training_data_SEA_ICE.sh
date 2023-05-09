#!/bin/bash -f
#$ -N Training_data_SEA_ICE
#$ -l h_rt=00:15:00
#$ -S /bin/bash
#$ -pe shmem-1 1
#$ -l h_rss=4G,mem_free=4G
#$ -q research-r8.q
#$ -t 1-120
##$ -j y
##$ -m ba
#$ -o /home/justinec/Documents/OUT/OUT_$JOB_NAME.$JOB_ID_$TASK_ID
#$ -e /home/justinec/Documents/ERR/ERR_$JOB_NAME.$JOB_ID_$TASK_ID
##$ -R y
##$ -r y

source /modules/rhel8/conda/install/etc/profile.d/conda.sh
conda activate /lustre/storeB/users/justinec/myconda

echo "Got $NSLOTS slots for job $SGE_TASK_ID."

cat > "/home/justinec/Documents/PROG/Training_data_SEA_ICE_""$SGE_TASK_ID"".py" << EOF

################################################################################################


#!/usr/bin/env python
# coding: utf-8

# This script collocate the data of the MOSAiC campaign's radiosoundings and of the OSI SAF's ice concentration for one radiosounding. 
# 
# This scrip is used in a bash script to collocate data for multiple radiosoundings.

# ### IMPORT

# In[5]:


from netCDF4 import Dataset as ncfile
import datetime
import os
from netCDF4 import num2date
import numpy as np
import matplotlib.dates as dates
import matplotlib.pyplot as plt
import pyproj
import glob


# ### READ DATA

# In[98]:


date_min = '20191101000000'
date_max = '20191201180000'
task_ID = $SGE_TASK_ID


# In[99]:


"""Function that create a list of dates and take the date corresponding to the index giving by the bash script."""
def task_date(date_min, date_max, task_ID):
    current_date = datetime.datetime.strptime(date_min, '%Y%m%d%H%M%S')
    end_date = datetime.datetime.strptime(date_max, '%Y%m%d%H%M%S')
    list_date = []
    while current_date <= end_date:
        list_date.append(current_date.strftime('%Y%m%d%H%M%S'))
        current_date = current_date + datetime.timedelta(hours = 6)
    date_task = list_date[task_ID - 1]
    return(date_task)


# In[100]:


date_task = task_date(date_min,date_max,task_ID) #str
datetask = datetime.datetime.strptime(date_task, '%Y%m%d%H%M%S') #datetime


# In[101]:


print('date task :', datetask)


# In[102]:


datemin = datetask - datetime.timedelta(days = 3)
datemax = datetask + datetime.timedelta(days = 3)
date_min = datemin.strftime('%Y%m%d')
date_max = datemax.strftime('%Y%m%d')


# In[103]:


year=date_task[0:4] ; month=date_task[4:6] ; day=date_task[6:8] ; hour=date_task[8:10]

# ------ EUMETSAT OSISAF ------
ppidir_osi = '/lustre/storeB/project/copernicus/osisaf/data/reprocessed/ice/conc/v3p0/'
osi_link  = ppidir_osi + year + '/' + month + '/ice_conc_nh_ease2-250_cdr-v3p0_'+year+month+day+'1200.nc'
osi = ncfile(osi_link,'r')        #dataset of EUMETSAT OSISAF

# ------ MOSAIC ------
ppidir_mosaic = '/lustre/storeB/users/maltem/Arctic/MOSAiC/radiosondes/'+year+'/'+month+'/'
mosaic_link  = ppidir_mosaic + 'PST-RS-01_2_RS41-GDP_001_'+year+month+day+'T'+hour+'0000_1-000-001.nc'
if os.path.isfile(mosaic_link) == False:
    mosaic_link  = ppidir_mosaic + 'PST-RS-01_2_RS41-GDP_001_'+year+month+day+'T'+hour+'0000_1-000-002.nc'
if os.path.isfile(mosaic_link) == False:
    mosaic_link  = ppidir_mosaic + 'PST-RS-01_2_RS41-GDP_001_'+year+month+day+'T'+hour+'0000_1-000-003.nc'
mosaic = ncfile(mosaic_link,'r')  #dataset of radiosoundings of MOSAiC

# ------ CRYOSAT-SMOS ------
ppidir_cryosat = '/lustre/storeB/users/yuriib/thredds/carra-sice-verif/products/CryoSat-2_SMOS/'
datemin = datetask - datetime.timedelta(days = 3)
datemax = datetask + datetime.timedelta(days = 3)
date_min = datemin.strftime('%Y%m%d')
date_max = datemax.strftime('%Y%m%d')
cryosat_link = ppidir_cryosat + year + '/' + month + '/W_XX-ESA,SMOS_CS2,NH_25KM_EASE2_'+date_min+'_'+date_max+'_r_v204_01_l4sit.nc'
if os.path.isfile(cryosat_link) == False:
    diff = []
    three_days = datetime.timedelta(days = 3)
    liste = sorted(glob.glob(ppidir_cryosat + year + '/' + month + '/W_XX-ESA,SMOS_CS2,NH_25KM_EASE2_*.nc'))
    for i in range(0,len(liste)) :
        datemin = datetime.datetime.strptime(liste[i][117:125], '%Y%m%d')
        middle = datemin + three_days
        if datetask > middle :
            diff.append(datetask - middle)
        if datetask < middle :
            diff.append(middle - datetask)
    print('diff :', diff)
    if np.min(diff) < three_days :
        print('diff min :', np.min(diff))
        ind = np.argmin(diff)
        print('ind :', ind)
        cryosat_link = liste[ind]
        cryosat = ncfile(cryosat_link,'r')  #dataset of CRYOSAT-SMOS
    else :
        print('no cryosat-smos data available for', datetask)
cryosat = ncfile(cryosat_link,'r')  #dataset of CRYOSAT-SMOS
print(cryosat_link)


# In[104]:


# ----- EUMETSAT OSISAF -----
osi_time = osi.variables['time'][:]


# In[105]:


# ----- RADIOSOUNDINGS MOSAiC -----
mosaic_lat = mosaic.variables['lat'][:]
mosaic_lon = mosaic.variables['lon'][:]
mosaic_pres = mosaic.variables['press'][:]
mosaic_time = mosaic.variables['time'][:]


# In[106]:


# ----- CRYOSAT-SMOS -----
cryosat_time = cryosat.variables['time'][:]


# ### CONVERT TIME

# In[107]:


"""Function that convert times of a dataset depemding on its unit and calendar 
type to a datetime then to Matplotlib dates.
---exemple---
ERA5 :  'hours since 1900-01-01 00:00:00.0' become datetime of type 2019-10-01 00:00:00 then 18170
MOSAiC : 'seconds since 2019-10-31T22:55:09.757Z' become datetime of type 2019-10-31 22:55:09.757000 then 18200.954974039352
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


# In[108]:


osi_nctime,osi_nctimenum, t_unit = convert_time(osi)
mosaic_nctime, mosaic_nctimenum, t_unit = convert_time(mosaic)
cryosat_nctime, cryosat_nctimenum, t_unit = convert_time(cryosat)


# In[109]:


print('osisaf : {} \n mosaic : {} \n cryosat-smos : {}'.format(osi_nctime[0], mosaic_nctime[0], cryosat_nctime[0]))


# ### CORRESPONDING INDEXES OF TIME AND OF PRESSURE LEVELS BETWEEN OSISAF AND MOSAIC

# In[110]:


"""Function that finds the indexes of the closest values between OSISAF and MOSAiC.
   Useful for:
- keep the latitudes, longitudes and times of OSISAF corresponding to MOSAiC
- keep the MOSAiC pressure levels corresponding to desired levels"""
def corresponding_index(osi,mosaic) :
    diff = np.absolute(osi-mosaic)
    indx = diff.argmin()
    return indx


# In[111]:


desired_levels = list(range(300, 1025, 25)) #from 300 hPa to 1000 hPa by 25 hPa
#desired_levels=era5_pres #if we want the same levels as ERA5 simply decomment this line


# In[112]:


# Considering a single MOSAiC time because radiosounding lasts 1h30. It's short enough for all radiosounding times to be close to the same OSISAF time [12:00]
indx_time = corresponding_index(osi_nctimenum,mosaic_nctimenum[0])
    
# Indexes of MOSAiC pressure levels closest to the desired pressure levels
indx_level = []
for i in range(0,len(desired_levels)) :
    indx_level.append(corresponding_index(mosaic_pres,desired_levels[i]))


# In[113]:


mosaic_lat_collocated = mosaic_lat[indx_level]    #keep mosaic's latitudes of each desired pressure level
mosaic_lon_collocated = mosaic_lon[indx_level]    #keep mosaic's longitudes of each desired pressure level
mosaic_time_collocated = mosaic_time[indx_level]  #keep mosaic's times of each desired pressure level


# ### CONVERT LATITUDE AND LONGITUDE OF MOSAIC IN COORDINATES OF PROJECTION OF OSISAF

# In[114]:


proj_osi = "+proj=laea +lon_0=0 +datum=WGS84 +ellps=WGS84 +lat_0=90.0"
crs_osi = pyproj.CRS.from_proj4(proj_osi)
proj_mosaic = "+proj=latlon"
crs_mosaic = pyproj.CRS.from_proj4(proj_mosaic)
transform_mosaic_to_osi = pyproj.Transformer.from_crs(crs_mosaic, crs_osi, always_xy = True)
xx_osiproj, yy_osiproj = transform_mosaic_to_osi.transform(mosaic_lon_collocated[-1], mosaic_lat_collocated[-1])
xx_osiproj = xx_osiproj*1e-3
yy_osiproj = yy_osiproj*1e-3


# In[115]:


osi_xc = osi.variables['xc'][:]
osi_yc = osi.variables['yc'][:]

cryosat_xc = cryosat.variables['xc'][:]
cryosat_yc = cryosat.variables['yc'][:]


# ### CORRESPONDING INDEXES OF COORDINATES OF PROJECTION XC AND YC BETWEEN OSISAF AND MOSAIC

# In[116]:


osi_indx_xc = corresponding_index(osi_xc,xx_osiproj)
osi_indx_yc = corresponding_index(osi_yc,yy_osiproj)


# In[117]:


cryosat_indx_xc = corresponding_index(cryosat_xc,xx_osiproj)
cryosat_indx_yc = corresponding_index(cryosat_yc,yy_osiproj)


# ### SEA ICE CONCENTRATION FROM OSISAF AND CRYOSAT-SMOS, SEA ICE THICKNESS AND TYPE FROM CRYOSAT-SMOS

# In[120]:


#si for sea ice
si_conc_osi = osi.variables['ice_conc'][0][osi_indx_yc][osi_indx_xc]
si_conc_cryosat = cryosat.variables['sea_ice_concentration'][0][cryosat_indx_yc][cryosat_indx_xc]
si_thickness = cryosat.variables['analysis_sea_ice_thickness'][0][cryosat_indx_yc][cryosat_indx_xc]
si_type = cryosat.variables['sea_ice_type'][0][cryosat_indx_yc][cryosat_indx_xc]


# ### SAVE DATA IN A .TXT FILE

# In[121]:


path_file = '/lustre/storeB/users/justinec/master_internship/data/SEA_ICE/'


# In[122]:


Output_file = path_file + 'file.txt'
if os.path.isfile(Output_file) == False:
    output_header = 'date'+'\t'+'sea_ice_concentration_osi[%]'+'\t'+'sea_ice_concentration_cryosat-smos[%]'+'\t'+'sea_ice_thickness[m]'+'\t'+'sea_ice_type[2 young, 3 old]'
    #
    Output = open(Output_file, 'w')
    Output.write(output_header + '\n')
    Output.close()
#
output_str = date_task + '\t' + str(si_conc_osi) + '\t' + str(si_conc_cryosat) + '\t' + str(si_thickness) + '\t' + str(si_type)

#
Output = open(Output_file, 'a')
Output.write(output_str + '\n')
Output.close()


##################################################################################################

EOF
python3 "/home/justinec/Documents/PROG/Training_data_SEA_ICE_""$SGE_TASK_ID"".py"