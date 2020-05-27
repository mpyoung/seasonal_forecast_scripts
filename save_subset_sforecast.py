'''

M. Young 12/02/2019
'''
from __future__ import division
# from astropy.convolution import convolve
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
from sklearn.metrics import roc_curve, roc_auc_score
import sys
execfile('date_str.py')

dir_main = '/gws/nopw/j04/ncas_climate_vol1/users/vboult/seasonal_forecasting/'
dir_c3s = dir_main+'c3s_hindcasts/global/'
dir_obs = '/gws/nopw/j04/ncas_climate_vol2/users/myoung02/datasets/CHIRPS/global_daily/'
dir_out = dir_main+'c3s_hindcasts/region_subsets/'

# Define region for subsetting
# region = 'Africa'
# latlim = [-35,35]
# lonlim = [-20,52]
region = 'Kenya'
latlim = [-5,5]
lonlim = [33,43]
start_mon = [10,11,12] # OND, ND, D
end_mon = [12,12,12] # OND, ND, D

years = np.arange(1993,2016+1,1)
center_ids = ['ecmwf_sys5','ukmo_sys13','meteofrance_sys6','dwd_sys2'] # modelling center names
c3s_nmembers = [25,7,25,30] # number of ensemble members for each model

# open dummy file to set lon's and lat's
c3s_file = dir_c3s+'c3s_ecmwf_sys5_seas_1993_2016_12.nc'
nc_fid = Dataset(c3s_file, 'r')
c3s_lat1 = np.flipud(np.array(nc_fid.variables['latitude'][:]))  # extract/copy the data
c3s_lon1 = np.array(nc_fid.variables['longitude'][:])
rfe = np.array(nc_fid.variables['tp'][:])
nc_fid.close()
rfe,c3s_lon2 = shiftgrid(180.,rfe,c3s_lon1,start=False)
rfe = []

# Open chirps to get lat and lon
# tmp_file= dir_obs+'chirps-v2.0.1983.days_p05_Africa.nc'
# nc_fid = Dataset(tmp_file, 'r')
# tmp_lat = np.array(nc_fid.variables['latitude'][:])  # extract/copy the data
# tmp_lon = np.array(nc_fid.variables['longitude'][:])
# nc_fid.close()
# lon_id = np.where((tmp_lon >= lonlim[0]) & (tmp_lon <= lonlim[1]))[0]
# lat_id = np.where((tmp_lat >= latlim[0]) & (tmp_lat <= latlim[1]))[0]
# chp_lon = tmp_lon[lon_id]
# chp_lat = tmp_lat[lat_id]
# rg_lon,rg_lat = np.meshgrid(chp_lon,chp_lat)


def do_stuff(model,month_start,month_end,month_lead,region,lonlim,latlim):
  '''

  '''
  # grab rfe at specific latitude and longitude (inlat,inlon)
  c3s_lon_id = np.where((c3s_lon2 >= lonlim[0]) & (c3s_lon2 <= lonlim[1]))[0]
  c3s_lat_id = np.where((c3s_lat1 >= latlim[0]) & (c3s_lat1 <= latlim[1]))[0]
  c3s_lon = c3s_lon2[c3s_lon_id]
  c3s_lat = np.flipud(c3s_lat1[c3s_lat_id])

  print c3s_nmembers[model]
  all_c3s = []
  all_c3s = np.nan*np.zeros((len(years),c3s_nmembers[model],len(c3s_lat),len(c3s_lon)))

  # open c3s data for the given season
  c3s_file = dir_c3s+'c3s_'+center_ids[model]+'_seas_1993_2016_'+mon_string(month_start-month_lead)+'.nc'
  if (os.path.isfile(c3s_file) == True):
    #print 'c3s file exists for '+sea_names[s]
    c3s_nc_fid = []
    c3s_nc_fid = Dataset(c3s_file, 'r')
    t_time = np.array(c3s_nc_fid.variables['time'][:])
    t_units = str(c3s_nc_fid.variables['time'].units)
    # # convert time to date using netCDF4 function
    c3s_dates = nc4.num2date(t_time,t_units)
    lon = np.array(c3s_nc_fid.variables['longitude'][:])
    c3s_rfe = np.array(c3s_nc_fid.variables['tp'][:])
    c3s_nc_fid.close()
    c3s_rfe,lon2 = shiftgrid(180.,c3s_rfe,lon,start=False)
    c3s_tp = c3s_rfe[:,:,c3s_lat_id,:][:,:,:,c3s_lon_id].squeeze()

    for y in np.arange(0,len(years)):
      print years[y]
      # Get forecast data accumulated through time
      r_start = datetime.strptime(str(years[y])+'/'+mon_string(month_start-month_lead)+'/01','%Y/%m/%d')
      f_start = r_start + timedelta(days=30*month_lead)

      if month_start > month_end: # allows to cross year boundary, e.g. DJF
        f_end = r_start + timedelta(days=30*(month_lead+((month_end+month_start)-month_start+1)))
      else:
        f_end = r_start + timedelta(days=30*(month_lead+(month_end-month_start+1)))

      if y == 0:
        print 'lead '+str(month_lead)
        print f_start
        print f_end

      # find start and end dates in forecast
      st_id = np.where(c3s_dates == f_start)[0]
      ed_id = np.where(c3s_dates == f_end)[0]
      if (len(st_id) > 0) & (len(ed_id)>0):
        all_c3s[y,:,:,:] = np.flip((c3s_tp[ed_id,:,:,:]-c3s_rfe[st_id,:,:,:])*1000.,axis=2).squeeze()
      elif (len(st_id) == 0):
        all_c3s[y,:,:,:] = np.flip(c3s_tp[ed_id,:,:,:]*1000.,axis=2).squeeze()
      else:
        print 'error with time subsetting: check forecast start and end dates'
        xxx
        # for nm in np.arange(0,c3s_nmembers[model]):
          # all_c3s[y,nm,:,:] = np.flip((c3s_tp[ed_id,nm,:,:]-c3s_rfe[st_id,nm,:,:])*1000.,axis=1).squeeze()
          # all_c3s[y,nm,:,:] = basemap.interp(tmp,c3s_lon,c3s_lat,rg_lon,rg_lat,order=0)

    nc_outfile = dir_out+center_ids[model]+'_total_rainfall_'+region+'_s'+mon_string(month_start)+'_e'+mon_string(month_end)+'_lead'+str(month_lead)+'.nc'
    # save final array as netcdf
    dataset = Dataset(nc_outfile,'w',format='NETCDF4')
    time = dataset.createDimension('time',len(years)) # create time
    nmem = dataset.createDimension('members',c3s_nmembers[model]) # create lon
    lat = dataset.createDimension('lat',len(c3s_lat)) # create lat (dims depend on region)
    lon = dataset.createDimension('lon',len(c3s_lon)) # create lon
    # create variables
    rfe_out = dataset.createVariable('rfe','d',('time','members','lat','lon'))

    latitudes = dataset.createVariable('latitude','f',('lat',))
    longitudes = dataset.createVariable('longitude','f',('lon',))
    times = dataset.createVariable('years', np.float64, ('time',))

    # Global Attributes (will need modified accordingly)
    dataset.description = center_ids[model]+' lead '+str(month_lead)+' months rainfall accumulation from beginning of month '+str(month_start)+' to end of month '+str(month_end)+', 1 degree horizontal resolution subset over '+region
    dataset.history = 'Created ' + tt.ctime(tt.time())
    dataset.source = 'Subset by M. Young'
    # Variable Attributes
    latitudes.units = 'degree_north'
    longitudes.units = 'degree_east'
    rfe_out.units = 'mm'
    times.units = 'years'
    # times.calendar = 'gregorian'

    # Fill variables with data
    latitudes[:] = c3s_lat#np.flipud(sub_lat)
    longitudes[:] = c3s_lon
    rfe_out[:] = all_c3s#np.flip(all_data,axis=2)# average in mm h then mmd

    # Fill in time.
    #date = []
    #date = datetime(year, month, day)
    times[:] = years#date2num(date, units = times.units, calendar = times.calendar)
    dataset.close()
    print 'saved '+center_ids[model]+' '+str(datetime.now())

if __name__ == "__main__":
   output = do_stuff(int(sys.argv[1]),int(sys.argv[2]),int(sys.argv[3]),int(sys.argv[4]),sys.argv[5],[int(sys.argv[6]),int(sys.argv[7])],[int(sys.argv[8]),int(sys.argv[9])])
