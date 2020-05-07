'''
Combines individual forecast ensemble outcomes (saved in tmp folder)
for certain event (e.g. lower tercile) into forecast probability
for that event.

Saves as netcdf file

Submit this to lotus using the shell script
'submit_save_sforecast_probabilities.sh'

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


dir_obs = '/gws/nopw/j04/ncas_climate_vol2/users/myoung02/datasets/CHIRPS/global_daily/'
dir_main = '/gws/nopw/j04/ncas_climate_vol1/users/vboult/seasonal_forecasting/'
dir_c3s = dir_main+'c3s_hindcasts/region_subsets/'
dir_out = dir_main+'c3s_hindcasts/region_subsets/forecast_probabilities/'
dir_tmp = dir_out+'tmp_probs/'


years = np.arange(1993,2016+1,1)
c3s_nmembers = [25,7,25,30]
center_ids = ['ecmwf_sys5','ukmo_sys13','meteofrance_sys6','dwd_sys2']
tercile_names = ['Below normal','Normal','Above normal']
tercile_short_name = ['Dry','Normal','Wet']
# define percentile categories to analyse
percentiles = [100./3.,200./3.]
p_name = 'tercile' # give categories an appropriate name


c3s_file = dir_c3s+'ecmwf_sys5_total_rainfall_Kenya_s10_e12_lead0.nc'
nc_fid = Dataset(c3s_file, 'r')
c3s_lat = np.array(nc_fid.variables['latitude'][:])  # extract/copy the data
c3s_lon = np.array(nc_fid.variables['longitude'][:])
nc_fid.close()

def do_stuff(model,month_start,month_end,month_lead,region):

  in_files=dir_tmp+'tmp_'+center_ids[model]+'_total_rainfall_'+region+'_s'+mon_string(month_start)+'_e'+mon_string(month_end)+'_lead'+str(month_lead)+'_p_*.npy'

  while len(glob.glob(in_files)) < len(years)*c3s_nmembers[model]:
    print 'dir but not enough files, sleeping...'+str(datetime.now())
    tt.sleep(60*15) # wait

  if len(glob.glob(in_files)) == len(years)*c3s_nmembers[0]:
    print 'All files here '+str(datetime.now())
    ye_list = []
    for a in itt.product(years,np.arange(0,c3s_nmembers[0])):
      ye_list.append(a)

    # load data and save
    all_p = np.nan*np.zeros((len(years),c3s_nmembers[model],len(percentiles)+1,len(c3s_lat),len(c3s_lon)))
    for ye in ye_list:
      y_id = np.where(ye[0]==years)[0]#year index
      e_id = np.where(ye[1]==np.arange(0,c3s_nmembers[0]))# ensemble member
      curr_file = dir_tmp+'tmp_'+center_ids[model]+'_total_rainfall_'+region+'_s'+mon_string(month_start)+'_e'+mon_string(month_end)+'_lead'+str(month_lead)+'_p_'+str(ye[0])+'_'+str(ye[1])+'.npy'
      all_p[y_id,e_id,:,:,:] = np.load(curr_file)
      os.system('rm '+curr_file)
      #print curr_file

    nc_outfile = dir_out+center_ids[model]+'_total_rainfall_'+p_name+'_probabilities_'+region+'_s'+mon_string(month_start)+'_e'+mon_string(month_end)+'_lead'+str(month_lead)+'.nc'
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
    dataset.description = center_ids[model]+' lead '+str(month_lead)+'months '+p_name+' probabilities of accumulated rainfall from beginning of month '+str(month_start)+' to end of month '+str(month_end)+', 1 degree horizontal resolution subset over '+region
    dataset.history = 'Created ' + tt.ctime(tt.time())
    dataset.source = 'Subset by M. Young'
    # Variable Attributes
    latitudes.units = 'degree_north'
    longitudes.units = 'degree_east'
    p_out.units = 'mm'
    times.units = 'years'
    # times.calendar = 'gregorian'

    # Fill variables with data
    latitudes[:] = c3s_lat
    longitudes[:] = c3s_lon
    p_out[:] = np.nanmean(all_p,axis=1)
    # Fill in time.
    #date = []
    #date = datetime(year, month, day)
    times[:] = years#date2num(date, units = times.units, calendar = times.calendar)
    dataset.close()
    print 'saved '+center_ids[model]+' '+str(datetime.now())

if __name__ == "__main__":
   output = do_stuff(int(sys.argv[1]),int(sys.argv[2]),int(sys.argv[3]),int(sys.argv[4]),sys.argv[5])
