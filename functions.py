# File of functions used in Create_training_data_IFS_MOSAiC.ipynb

import numpy as np
from netCDF4 import Dataset as ncfile
from netCDF4 import num2date
import matplotlib.dates as dates
import datetime

# -----------------------------------------------------------------------------------------

# -----------------------------------------------------------------------------------------
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