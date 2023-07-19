import os
import sys
import numpy as np
from netCDF4 import Dataset
from scipy.spatial.distance import *
function_path='/lustre/storeB/users/cyrilp/Python_functions/'
sys.path.insert(0, function_path)
from position_sea_ice_edge_including_coastlines import *
##########################################
# Constants
##########################################
file_T4 = '/lustre/storeB/project/copernicus/ARC-MFC/ARC-METNO-ARC-TOPAZ4_2_PHYS-FOR/arctic/mersea-class1/20210122_dm-metno-MODEL-topaz4-ARC-b20210114-fv02.0.nc'
path_output = '/lustre/storeB/project/copernicus/svalnav/Data/TOPAZ4/'
##########################################
# Dataset
##########################################
nc_T4 = Dataset(file_T4, 'r')
x_T4 = nc_T4.variables['x'][:]
y_T4 = nc_T4.variables['y'][:]
lat_T4 = nc_T4.variables['latitude'][:,:]
lon_T4 = nc_T4.variables['longitude'][:,:]
temperature_T4 = nc_T4.variables['temperature'][0,0,:,:] 
#
LSM_T4 = np.zeros(np.shape(lat_T4))
LSM_T4[temperature_T4.mask == False] = 1
##########################################
# Distance to land
##########################################
xx_T4, yy_T4 = np.meshgrid(x_T4 * 100 * 1000, y_T4 * 100 * 1000)
xx_T4_flat = np.ndarray.flatten(xx_T4)
yy_T4_flat = np.ndarray.flatten(yy_T4)
Coord_T4 = np.array([xx_T4_flat, yy_T4_flat]).T
#
coastlines_T4_flat = np.ndarray.flatten(position_sea_ice_edge_including_coastlines(LSM_T4))
Coord_coastlines_T4 = np.array([xx_T4_flat[coastlines_T4_flat == 1], yy_T4_flat[coastlines_T4_flat == 1]]).T
#
D_mat = cdist(Coord_T4, Coord_coastlines_T4, metric = 'euclidean')
Dist_coastlines = np.nanmin(D_mat, axis = 1)
Dist_coastlines_mat = Dist_coastlines.reshape(np.shape(temperature_T4))
Dist_coastlines_mat[LSM_T4 == 0] = np.nan
##########################################
# Output netCDF file
##########################################
output_filename = path_output + 'TOPAZ4_land_sea_mask.nc'
output_netcdf = Dataset(output_filename, 'w', format = 'NETCDF4')
#
x = output_netcdf.createDimension('x',len(x_T4))
y = output_netcdf.createDimension('y',len(y_T4))
#
x = output_netcdf.createVariable('x', 'd', ('x'))
y = output_netcdf.createVariable('y', 'd', ('y'))
latitude = output_netcdf.createVariable('latitude', 'd', ('y','x'))
longitude = output_netcdf.createVariable('longitude', 'd', ('y','x'))
LSM = output_netcdf.createVariable('LSM', 'd', ('y','x'))
distance_to_land = output_netcdf.createVariable('distance_to_land', 'd', ('y','x'))
coastlines = output_netcdf.createVariable('coastlines', 'd', ('y','x'))
#
x.units = '100 km'
x.standard_name = 'projection_x_coordinate'
y.units = '100 km'
y.standard_name = 'projection_y_coordinate'
latitude.units = 'degrees_north'
latitude.standard_name = 'latitude'
longitude.units = 'degrees_east'
longitude.standard_name = 'longitude'
LSM.units = 'mask'
LSM.standard_name = 'Land sea mask (1: ocean, 0: land)'
distance_to_land.units = 'meters'
distance_to_land.standard_name = 'distance to land'
coastlines.standard_name = 'coastlines'
#
x[:] = x_T4
y[:] = y_T4
longitude[:,:] = lon_T4
latitude[:,:] = lat_T4
LSM[:,:] = LSM_T4
distance_to_land[:,:] = Dist_coastlines_mat
coastlines[:,:] = coastlines_T4
#
output_netcdf.grid_projection = "+proj=stere +lon_0=-45. +lat_ts=90. +lat_0=90. +a=6378273. +b=6378273. ellps=sphere"
output_netcdf.close()
