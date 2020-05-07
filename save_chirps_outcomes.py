'''
Saves CHIRPS  outcomes (1/0) for certain
event (e.g. lower tercile)

Saves as netcdf file

Submit this to lotus using the shell script
'submit_save_chirps_outcomes.sh'

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
from sklearn.metrics import roc_curve, roc_auc_score
import itertools as itt
import sys
execfile('date_str.py')


def mask_percentiles(data,percentile,arg):
  '''
  Function: 'mask_percentiles':
  Function takes in a 2d spatial forecast array and 2d spatial array of corresponding percentiles
  and masks the forecast array according to input arguments, 'arg':
  'between': mask the forecast between 2 percentiles
  'below': mask the forecast below a percentile
  'above': mask the forecast above a percentile
  '''
  data_mask= np.copy(data)

  if arg == 'between': # find between 2 percentiles
    # the method of masking below doesn't work
    #data_mask[(data_mask >= percentile[0,:,:]) & (data_mask <= percentile[1,:,:])] = np.nan
    idx = np.where((data_mask >= percentile[0,:,:]) & (data_mask <= percentile[1,:,:]))
    data_mask[idx] = np.nan
    data_mask[np.isnan(data_mask)==False] = 0
    data_mask[np.isnan(data_mask)==True] = 1

  elif arg == 'below': # find data below percentile
    data_mask[data_mask < percentile] = np.nan
    data_mask[data_mask >= percentile] = 0
    data_mask[np.isnan(data_mask)==True] = 1

  elif arg == 'above': # find data above percentile
    data_mask[data_mask > percentile] = np.nan
    data_mask[data_mask <= percentile] = 0
    data_mask[np.isnan(data_mask)==True] = 1

  return data_mask


dir_main = '/gws/nopw/j04/ncas_climate_vol1/users/vboult/seasonal_forecasting/'
dir_c3s = dir_main+'c3s_hindcasts/region_subsets/'
dir_out = dir_main+'c3s_hindcasts/region_subsets/forecast_probabilities/'

years = np.arange(1993,2016+1,1)
c3s_nmembers = [25,7,25,30]
center_ids = ['ecmwf_sys5','ukmo_sys13','meteofrance_sys6','dwd_sys2']
tercile_names = ['Below normal','Normal','Above normal']
tercile_short_name = ['Dry','Normal','Wet']
# define percentile categories to analyse
percentiles = [100./3.,200./3.]
p_name = 'tercile' # give categories an appropriate name


def do_stuff(month_start,month_end,month_lead,region):
  c3s_file = dir_c3s+'CHIRPS_total_rainfall_'+region+'_s'+mon_string(month_start)+'_e'+mon_string(month_end)+'_lead'+str(month_lead)+'.nc'
  nc_fid_c3s = Dataset(c3s_file,'r')
  all_c3s = np.array(nc_fid_c3s.variables['rfe'])
  c3s_lon = np.array(nc_fid_c3s.variables['longitude'])
  c3s_lat = np.array(nc_fid_c3s.variables['latitude'])
  nc_fid_c3s.close()

  land_sea_mask = np.mean(all_c3s,axis=0).squeeze()
  land_sea_mask[np.isnan(land_sea_mask)==False]=1

  all_p = np.nan*np.zeros((len(years),3,len(all_c3s[0,:,0]),len(all_c3s[0,0,:])))

  for y in np.arange(0,len(years)):
    print 'starting '+str(years[y])+' '+str(datetime.now())
    #get percentiles for CHIRPS
    y_id = np.where(years == years[y])[0]
    tmp_c3s = []
    tmp_c3s = np.delete(all_c3s,y_id,axis=0) # delete year of interest so independant of climatology
    c3s_percentiles = np.nanpercentile(tmp_c3s,percentiles,axis=0)
    curr_c3s = np.copy(all_c3s[y_id,:,:])
    for nn in np.arange(0,len(percentiles)+1):
     if nn == 0:
       all_p[y,nn,:,:] = mask_percentiles(curr_c3s,c3s_percentiles[nn,:,:],'below')
     elif nn == len(percentiles):
       all_p[y,nn,:,:] = mask_percentiles(curr_c3s,c3s_percentiles[nn-1,:,:],'above')
     else:
       all_p[y,nn,:,:] = mask_percentiles(curr_c3s,c3s_percentiles[nn-1:nn+1,:,:],'between')

  nc_outfile = dir_out+'CHIRPS_total_rainfall_'+p_name+'_outcomes_'+region+'_s'+mon_string(month_start)+'_e'+mon_string(month_end)+'_lead'+str(month_lead)+'.nc'
  # save final array as netcdf
  dataset = Dataset(nc_outfile,'w',format='NETCDF4')
  time = dataset.createDimension('time',len(years)) # create time
  perc = dataset.createDimension(p_name,len(percentiles)+1) # create lon
  lat = dataset.createDimension('lat',len(c3s_lat)) # create lat (dims depend on region)
  lon = dataset.createDimension('lon',len(c3s_lon)) # create lon
  # create variables
  p_out = dataset.createVariable('p','d',('time',p_name,'lat','lon'))
  latitudes = dataset.createVariable('latitude','f',('lat',))
  longitudes = dataset.createVariable('longitude','f',('lon',))
  times = dataset.createVariable('years', np.float64, ('time',))

  # Global Attributes (will need modified accordingly)
  dataset.description = 'CHIRPSv2.0 '+p_name+' outcomes of accumulated rainfall from beginning of month '+str(month_start)+' to end of month '+str(month_end)+', 1 degree horizontal resolution subset over '+region
  dataset.history = 'Created ' + tt.ctime(tt.time())
  dataset.source = 'Subset by M. Young'
  # Variable Attributes
  latitudes.units = 'degree_north'
  longitudes.units = 'degree_east'
  p_out.units = ''
  times.units = 'years'
  # times.calendar = 'gregorian'

  # Fill variables with data
  latitudes[:] = c3s_lat
  longitudes[:] = c3s_lon
  p_out[:] = all_p*land_sea_mask
  # Fill in time.
  #date = []
  #date = datetime(year, month, day)
  times[:] = years#date2num(date, units = times.units, calendar = times.calendar)
  dataset.close()
  print 'saved CHIRPS '+str(datetime.now())

if __name__ == "__main__":
   output = do_stuff(int(sys.argv[1]),int(sys.argv[2]),int(sys.argv[3]),sys.argv[4])
