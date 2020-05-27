'''
Saves, subsets and regrids CHIRPSv2.0 precipiation
into format for comparison to C3S seasonal forecast data

Saves as netcdf file

Submit this to lotus using the shell script
'submit_subset_chirps.sh'

M. Young 7/05/2020
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
dir_c3s_region = dir_main+'c3s_hindcasts/region_subsets/'

dir_obs = dir_main+'CHIRPS/'
dir_out = dir_c3s_region

# Define region for subsetting
# region = 'Africa'
# latlim = [-35,35]
# lonlim = [-20,52]
# region = 'Kenya'
# latlim = [-5,5]
# lonlim = [33,43]
# start_mon = [10,11,12] # OND, ND, D
# end_mon = [12,12,12] # OND, ND, D

years = np.arange(1993,2016+1,1)

center_ids = ['ecmwf_sys5','ukmo_sys13','meteofrance_sys6','dwd_sys2'] # modelling center names


# rg_lon,rg_lat = np.meshgrid(chp_lon,chp_lat)


def do_stuff(model,month_start,month_end,month_lead,region,lonlim,latlim):
  '''
  The following inputs are needed to open the seasonal forecast data
  and obtain the appropriate accumulation periods (start date and end date)
  to accumulate the chirpsv2.0 observations over.
  CHIRPSv2.0 is regridded to the same horizontal resolution as the forecast

  Inputs
  model: integer, index of the model (keep fixed as all models will have same format)
  month_start: integer, start month for rainfall accumulation
  month_end: integer, end month for rainfall accumulation
  month_lead: integer, lead time of forecast for accumulation
  region: name of region
  lonlim: tuple [lonmin, lonmax], longitude limits of region
  latlim: tuple [latmin, latmax], latitude limits of region

  Outputs
  Saves accumulated and regridded CHIRPSv2.0 in netcdf file
  '''
  # Open chirps to get lat and lon
  tmp_file= dir_obs+'chirps-v2.0.1983.days_p05_Africa.nc'
  nc_fid = Dataset(tmp_file, 'r')
  tmp_lat = np.array(nc_fid.variables['latitude'][:])  # extract/copy the data
  tmp_lon = np.array(nc_fid.variables['longitude'][:])
  nc_fid.close()
  lon_id = np.where((tmp_lon >= lonlim[0]-5) & (tmp_lon <= lonlim[1]+5))[0]
  lat_id = np.where((tmp_lat >= latlim[0]-5) & (tmp_lat <= latlim[1]+5))[0]
  chp_lon = tmp_lon[lon_id]
  chp_lat = tmp_lat[lat_id]

  # open dummy file to set lon's and lat's
  tmp_file = dir_c3s_region+center_ids[model]+'_total_rainfall_'+region+'_s'+mon_string(month_start)+'_e'+mon_string(month_end)+'_lead'+str(month_lead)+'.nc'
  nc_fid = Dataset(tmp_file, 'r')
  out_lat = np.array(nc_fid.variables['latitude'][:])  # extract/copy the data
  out_lon = np.array(nc_fid.variables['longitude'][:])
  nc_fid.close()
  rg_lon,rg_lat = np.meshgrid(out_lon,np.flip(out_lat))

  all_chp = []
  all_chp = np.nan*np.zeros((len(years),len(out_lat),len(out_lon)))

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
    c3s_nc_fid.close()


    for y in np.arange(0,len(years)):
      print years[y]
      # open chirps data for the equivalent year and regrid to c3s grid
      obs_file = dir_obs+'chirps-v2.0.'+str(years[y])+'.days_p05_Africa.nc'
      if (os.path.isfile(obs_file) == True):
        nc_fid_chp = []
        nc_fid_chp = Dataset(obs_file, 'r')
        t_time = np.array(nc_fid_chp.variables['time'][:])
        t_units = str(nc_fid_chp.variables['time'].units)
        nc_fid_chp.close()
        # convert time to date using netCDF4 function
        chp_dates = nc4.num2date(t_time,t_units)

      # Get chirps data for corresponding forecast dates
      r_start = datetime.strptime(str(years[y])+'/'+mon_string(month_start-month_lead)+'/01','%Y/%m/%d')
      f_start = r_start + timedelta(days=30*month_lead)
      if month_start < month_end:
        f_end = r_start + timedelta(days=30*(month_lead+(month_end-month_start+1)))
      elif month_start > month_end: # allows to cross year boundary, e.g. DJF 
        f_end = r_start + timedelta(days=30*(month_lead+((month_end+month_start)-month_start+1)))

      if y == 0:
        print 'lead '+str(month_lead)
        print f_start
        print f_end
      chp_start_id = np.where(chp_dates == f_start)[0]
      chp_end_id = np.where(chp_dates == f_end)[0]
      chp_time_id = np.arange(chp_start_id,chp_end_id+1)
      nc_fid_chp = Dataset(obs_file, 'r')
      tmp = []
      tmp =np.array(nc_fid_chp.variables['precip'][chp_time_id,:,:][:,lat_id,:][:,:,lon_id],dtype='d').squeeze()
      nc_fid_chp.close()
      tmp[tmp<0] = np.nan
      all_chp[y,:,:] = basemap.interp(np.sum(tmp,axis=0),chp_lon,chp_lat,rg_lon,rg_lat,order=0)

    nc_outfile = dir_out+'CHIRPS_total_rainfall_'+region+'_s'+mon_string(month_start)+'_e'+mon_string(month_end)+'_lead'+str(month_lead)+'.nc'
    # save final array as netcdf
    dataset = Dataset(nc_outfile,'w',format='NETCDF4')
    time = dataset.createDimension('time',len(years)) # create time
    lat = dataset.createDimension('lat',len(out_lat)) # create lat (dims depend on region)
    lon = dataset.createDimension('lon',len(out_lon)) # create lon
    # create variables
    rfe_out = dataset.createVariable('rfe','d',('time','lat','lon'))

    latitudes = dataset.createVariable('latitude','f',('lat',))
    longitudes = dataset.createVariable('longitude','f',('lon',))
    times = dataset.createVariable('years', np.float64, ('time',))

    # Global Attributes (will need modified accordingly)
    dataset.description = 'CHIRPSv2.0 rainfall accumulation from beginning of month '+str(month_start)+' to end of month '+str(month_end)+', regridded to 1 degree horizontal resolution subset over '+region
    dataset.history = 'Created ' + tt.ctime(tt.time())
    dataset.source = 'Subset by M. Young'
    # Variable Attributes
    latitudes.units = 'degree_north'
    longitudes.units = 'degree_east'
    rfe_out.units = 'mm'
    times.units = 'years'
    # times.calendar = 'gregorian'

    # Fill variables with data
    latitudes[:] = out_lat#np.flipud(sub_lat)
    longitudes[:] = out_lon
    rfe_out[:] = all_chp#np.flip(all_data,axis=2)# average in mm h then mmd

    # Fill in time.
    #date = []
    #date = datetime(year, month, day)
    times[:] = years#date2num(date, units = times.units, calendar = times.calendar)
    dataset.close()
    print 'saved CHIRPS '+str(datetime.now())

if __name__ == "__main__":
   output = do_stuff(int(sys.argv[1]),int(sys.argv[2]),int(sys.argv[3]),int(sys.argv[4]),sys.argv[5],[int(sys.argv[6]),int(sys.argv[7])],[int(sys.argv[8]),int(sys.argv[9])])
