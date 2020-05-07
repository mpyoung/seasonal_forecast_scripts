'''
Saves individual seasonal forecast ensemble outcomes into tmp folder
for certain event (e.g. lower tercile).

This needs to be run before running 'save_sforecast_probabilities'
which uses the individual ensemble outcomes to compute the
forecast probabilities for the given event

Individual ensemble member outcomes (1 or 0) are saved into
'.npy' files in specified tmp folder. 

Submit this to lotus using the shell script:
'submit_save_sforecast_outcomes.sh'

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


def ensemble_mask(y,e,years,percentiles,c3s_file):
  '''
  Function generates tercile outcomes (1 or 0) for each ensemble member.
  Tercile outcomes are saved then combined later to generate the
  forecast probabilities.
  '''

  print 'starting '+str(years[y])+', ensemble member '+str(e)+' '+str(datetime.now())
  nc_fid_c3s = Dataset(c3s_file,'r')
  all_c3s = np.array(nc_fid_c3s.variables['rfe'][:,e,:,:])
  nc_fid_c3s.close()
  #get percentiles for ECMWF
  y_id = np.where(years == years[y])[0]
  tmp_c3s = []
  tmp_c3s = np.delete(all_c3s,y_id,axis=0) # delete year of interest so independant of climatology
  c3s_percentiles = np.nanpercentile(tmp_c3s,percentiles,axis=0)
  curr_c3s = np.copy(all_c3s[y_id,:,:])
  curr_p = np.nan*np.zeros((3,len(all_c3s[0,:,0]),len(all_c3s[0,0,:])))
  for nn in np.arange(0,len(percentiles)+1):
   if nn == 0:
     curr_p[nn,:,:] = mask_percentiles(curr_c3s,c3s_percentiles[nn,:,:],'below')
   elif nn == len(percentiles):
     curr_p[nn,:,:] = mask_percentiles(curr_c3s,c3s_percentiles[nn-1,:,:],'above')
   else:
     curr_p[nn,:,:] = mask_percentiles(curr_c3s,c3s_percentiles[nn-1:nn+1,:,:],'between')
  fname_out = 'tmp_'+ os.path.basename(c3s_file).replace('.nc','')+'_p_'+str(years[y])+'_'+str(e)
  np.save(dir_tmp+fname_out,curr_p)

  print 'saved '+str(years[y])+', ensemble member '+str(e)+' '+str(datetime.now())
  return []


def do_stuff(model,month_start,month_end,month_lead,region,y,e):
  c3s_file = dir_c3s+center_ids[model]+'_total_rainfall_'+region+'_s'+mon_string(month_start)+'_e'+mon_string(month_end)+'_lead'+str(month_lead)+'.nc'
  out = ensemble_mask(y,e,years,percentiles,c3s_file)


if __name__ == "__main__":
   output = do_stuff(int(sys.argv[1]),int(sys.argv[2]),int(sys.argv[3]),int(sys.argv[4]),sys.argv[5],int(sys.argv[6]),int(sys.argv[7]))
