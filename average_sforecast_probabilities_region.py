'''
Averages seasonal forecast probabilities over a region
bounded by a shapefile (in this case, Kenya)

Forecast probabilities are saved in a text file.

Scripts that generate data required for this (in run order):
1) 'save_subset_sforecast.py'
2) 'save_sforecast_outcomes.py'
3) 'save_sforecast_probabilities.py'

M. Young 7/05/2020
'''
from __future__ import division
import glob
import numpy as np
from netCDF4 import Dataset
from netCDF4 import MFDataset
import time as tt
from datetime import datetime, timedelta
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from mpl_toolkits.basemap import Basemap,shiftgrid
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits import basemap
from scipy.stats.stats import pearsonr
import netCDF4 as nc4
import calendar
import os.path
import os
import itertools as itt
import sys
import shapely
import geopandas as gpd

execfile('date_str.py')


dir_main = '/gws/nopw/j04/ncas_climate_vol1/users/vboult/seasonal_forecasting/'
dir_in = dir_main+'c3s_hindcasts/region_subsets/forecast_probabilities/'
dir_out = dir_main+''
shapefile = dir_main+'/GIS/Kenya.shp'


years = np.arange(1993,2016+1,1)
c3s_nmembers = [25,7,25,30]
center_ids = ['ecmwf_sys5','ukmo_sys13','meteofrance_sys6','dwd_sys2']
tercile_names = ['Below normal','Normal','Above normal']
tercile_short_name = ['Dry','Normal','Wet']
# define percentile categories to analyse
percentiles = [100./3.,200./3.]
p_name = 'tercile' # give categories an appropriate name


model=0
region='Kenya'
month_start=12
month_end=12
month_lead=0

p_file = dir_in+center_ids[model]+'_total_rainfall_'+p_name+'_probabilities_'+region+'_s'+mon_string(month_start)+'_e'+mon_string(month_end)+'_lead'+str(month_lead)+'.nc'
nc_fid = Dataset(p_file, 'r')
c3s_lat = np.flip(np.array(nc_fid.variables['latitude'][:]))  # extract/copy the data
c3s_lon = np.array(nc_fid.variables['longitude'][:])
p_all = np.array(nc_fid.variables['p'][:]) #  probabilities
nc_fid.close()

mask_shape = gpd.read_file(shapefile)
[lon2d,lat2d] = np.meshgrid(c3s_lon,c3s_lat)
# lon2 = lon2d.reshape(-1)
# lat2 = lat2d.reshape(-1)
lon2=lon2d.flatten()
lat2=lat2d.flatten()
# Create empty mask
mask = []
# Fill mask with true / false
for lat, lon in zip(lat2, lon2):
  this_point = shapely.geometry.Point(lon, lat)
  res = mask_shape.contains(this_point)
  mask.append(res.values[0])

# Reshape mask
mask = np.array(mask).reshape(lon2d.shape)

# p_mask_test = np.ma.masked_array(p_all[0,0,:,:],mask=mask)
mask_final = np.copy(mask)*1.0
mask_final[mask_final==0]=np.nan
# mask_final[mask_final==True]=1
p_masked=p_all*mask_final

p_region = np.nanmean(p_masked,axis=(2,3))

# save as text file
f_out = dir_in+center_ids[model]+'_total_rainfall_'+p_name+'_probabilities_avg_'+region+'_s'+mon_string(month_start)+'_e'+mon_string(month_end)+'_lead'+str(month_lead)+'.txt'
headers=' '.join(tercile_short_name)
headers='Year '+headers
fm_list = ['%i','%10.5f','%10.5f','%10.5f']
np.savetxt(f_out,np.hstack((np.expand_dims(years,axis=1),p_region)),fmt=fm_list,delimiter=' ',header=headers)

#
# plt.subplot(1,2,1)
# mymap = Basemap(projection='cyl',resolution='l',\
#         llcrnrlat=np.min(c3s_lat),urcrnrlat=np.max(c3s_lat),\
#         llcrnrlon=np.min(c3s_lon),urcrnrlon=np.max(c3s_lon))
# # mymap.drawcoastlines(linewidth=1)
# # mymap.drawcountries(linewidth=1)
# mymap.readshapefile(shapefile.replace('.shp',''),'kenya')
# x, y = mymap(*np.meshgrid(c3s_lon,c3s_lat))
# mymap.pcolormesh(x,y,p_all[0,0,:,:])
# # mymap.colorbar()
# plt.subplot(1,2,2)
# mymap = Basemap(projection='cyl',resolution='l',\
#         llcrnrlat=np.min(c3s_lat),urcrnrlat=np.max(c3s_lat),\
#         llcrnrlon=np.min(c3s_lon),urcrnrlon=np.max(c3s_lon))
# # mymap.drawcoastlines(linewidth=1,color='grey')
# # mymap.drawcountries(linewidth=1)
# mymap.readshapefile(shapefile.replace('.shp',''),'kenya')
#
# x, y = mymap(*np.meshgrid(c3s_lon,c3s_lat))
# mymap.pcolormesh(x,y,p_mask_test)
# # mymap.colorbar()
# plt.show()
