'''
Prepare CHIRPS observations for comparison
to ECMWF SEAS5 and other C3S forecasts
c3s data is downloaded using get_c3s_seasonal.py


M. Young 12/02/2019
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
from sklearn.metrics import roc_curve, roc_auc_score

execfile('date_str.py')

dir_c3s = '/gws/nopw/j04/ncas_climate_vol1/users/myoung02/datasets/C3S/hindcasts/'
#dir_obs = '/gws/nopw/j04/ncas_climate_vol2/users/myoung02/datasets/TAMSAT/TAMSAT3/'
# dir_obs = '/gws/nopw/j04/ncas_climate_vol1/users/myoung02/datasets/TAMSAT3_seasonal/'
dir_obs = '/gws/nopw/j04/ncas_climate_vol2/users/myoung02/datasets/CHIRPS/global_daily/'
f_outdir = '/gws/nopw/j04/ncas_climate_vol1/users/myoung02/datasets/spatial_scale_paper/'
dir_out	= '/home/users/myoung02/spatial_scale_paper/'

# Define region for subsetting (extent of TAMSAT3 grid)
region = 'Africa'
latlim = [-35,35]
lonlim = [-20,52]

years = np.arange(1984,2016+1,1)#np.arange(1994,2016+1,1)
c3s_mon_file = ['11','02','05','08']#['12','03','06','09'] # DJF, MAM & JJA
c3s_mon = ['12','03','06','09'] # DJF, MAM & JJA
obs_seas = ['02','05','08','11'] #,'11'] #
sea_names = ['DJF','MAM','JJA','SON']#,'SON']
start_mon = [1,3,6,9]
end_mon = [2,5,8,11]
c3s_nmembers = 25
center_ids = ['ecmwf_sys5','ukmo_sys13','meteofrance_sys6','dwd_sys2']

box_size = 5.0
grid_size= int(np.round(box_size/1.0,0))
kernel = np.ones((grid_size,grid_size))/(grid_size*grid_size)


# Open chirps to get lat and lon
tmp_file= dir_obs+'chirps-v2.0.1983.days_p05_Africa.nc'
nc_fid = Dataset(tmp_file, 'r')
tmp_lat = np.array(nc_fid.variables['latitude'][:])  # extract/copy the data
tmp_lon = np.array(nc_fid.variables['longitude'][:])
nc_fid.close()
lon_id = np.where((tmp_lon >= lonlim[0]) & (tmp_lon <= lonlim[1]))[0]
lat_id = np.where((tmp_lat >= latlim[0]) & (tmp_lat <= latlim[1]))[0]
chp_lon = tmp_lon[lon_id]
chp_lat = tmp_lat[lat_id]

all_obs = np.nan*np.zeros((len(c3s_mon),len(years),len(chp_lat),len(chp_lon)))
# nan_id = np.zeros((len(c3s_mon),len(years)))

# get seasonal totals at leadtime 3 months (90 days)
# and open & regrid TAMSAT seasonal files
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
    for s in [0,1,2,3]:##[1,3]:
      print(sea_names[s])
      rfe_total = []
      if s == 0:
        # c3s file time 3 months ahead
        r_start = datetime.strptime(str(years[y]-1)+'/'+c3s_mon_file[s]+'/01','%Y/%m/%d')
        f_start = r_start + timedelta(days=30)#datetime.strptime(str(years[y]-1)+'/'+c3s_mon[s]+'/01','%Y/%m/%d')
        f_end = r_start + timedelta(days=120)

        # open c3s data for the given season
        # c3s_file = dir_c3s+'c3s_'+center_ids[0]+'_seas_1993_2016_'+c3s_mon_file[s]+'.nc'
        # if (os.path.isfile(c3s_file) == True):
        #   #print 'c3s file exists for '+sea_names[s]
        #   c3s_nc_fid = []
        #   c3s_nc_fid = Dataset(c3s_file, 'r')
        #   t_time = np.array(c3s_nc_fid.variables['time'][:])
        #   t_units = str(c3s_nc_fid.variables['time'].units)
        #   # # convert time to date using netCDF4 function
        #   c3s_dates = nc4.num2date(t_time,t_units)
        #   c3s_id = np.where(c3s_dates == f_end)[0]
        #   c3s_nc_fid.close()

        other_file = dir_obs+'chirps-v2.0.'+str(years[y]-1)+'.days_p05_Africa.nc'
        nc_fid2 = []
        nc_fid2 = Dataset(other_file, 'r')
        other_time = np.array(nc_fid2.variables['time'][:])
        other_units = str(nc_fid2.variables['time'].units)
        # # convert time to date using netCDF4 function
        other_dates = nc4.num2date(other_time,other_units)
        t_id_dec = np.where(other_dates >= f_start)[0]
        rfe2 = []
        rfe2 = np.array(nc_fid2.variables['precip'][t_id_dec,:,:][:,lat_id,:][:,:,lon_id],dtype='d').squeeze()
        nc_fid2.close()
        # loop through daily RFE and regrid to 1 degrees then sum
        rfe2[rfe2 <0] = np.nan
        rfe_sum2 = np.sum(rfe2,axis=0)

        t_id_chp = []
        t_id_chp = np.where((chp_dates >= f_start) & (chp_dates <= f_end))[0]
        nc_fid_chp = []
        nc_fid_chp = Dataset(obs_file, 'r')
        rfe = []
        rfe = np.array(nc_fid_chp.variables['precip'][t_id_chp,:,:][:,lat_id,:][:,:,lon_id],dtype='d').squeeze()
        nc_fid_chp.close()
        rfe[rfe <0] = np.nan
        rfe_sum1 = np.sum(rfe,axis=0)
        rfe_total = rfe_sum1+rfe_sum2
      else:
        # c3s file time 3 months ahead
        r_start = datetime.strptime(str(years[y])+'/'+c3s_mon_file[s]+'/01','%Y/%m/%d')
        f_start = r_start + timedelta(days=30)#datetime.strptime(str(years[y])+'/'+c3s_mon[s]+'/01','%Y/%m/%d')
        f_end = r_start + timedelta(days=120)
        # open c3s data for the given season
        # c3s_file = dir_c3s+'c3s_'+center_ids[0]+'_seas_1993_2016_'+c3s_mon_file[s]+'.nc'
        # if (os.path.isfile(c3s_file) == True):
        #   #print 'c3s file exists for '+sea_names[s]
        #   c3s_nc_fid = []
        #   c3s_nc_fid = Dataset(c3s_file, 'r')
        #   t_time = np.array(c3s_nc_fid.variables['time'][:])
        #   t_units = str(c3s_nc_fid.variables['time'].units)
        #   # # convert time to date using netCDF4 function
        #   c3s_dates = nc4.num2date(t_time,t_units)
        #   c3s_id = np.where(c3s_dates == f_end)[0]
        #   c3s_nc_fid.close()

        t_id_chp = []
        t_id_chp = np.where((chp_dates >= f_start) & (chp_dates <= f_end))[0]
        nc_fid_chp = []
        nc_fid_chp = Dataset(obs_file, 'r')
        rfe = []
        rfe = np.array(nc_fid_chp.variables['precip'][t_id_chp,:,:][:,lat_id,:][:,:,lon_id],dtype='d').squeeze()
        nc_fid_chp.close()
        rfe[rfe <0] = np.nan
        rfe_total = np.sum(rfe,axis=0)
            all_c3s[y,nm,:,:] = basemap.interp(rfe_total,chp_lon,chp_lat,rg_lon,rg_lat,order=0)

      all_obs[s,y,:,:] = np.copy(rfe_total)

for s in np.arange(0,len(sea_names)):
  nc_outfile = f_outdir+'CHIRPS.seasonal.0.05d_'+c3s_mon[s]+'_lead1_1983.nc'
  # save final array as netcdf
  dataset = Dataset(nc_outfile,'w',format='NETCDF4')
  time = dataset.createDimension('time',len(years)) # create time
  lat = dataset.createDimension('lat',len(chp_lat)) # create lat (dims depend on region)
  lon = dataset.createDimension('lon',len(chp_lon)) # create lon
  # create variables
  rfe_out = dataset.createVariable('rfe','d',('time','lat','lon'))
  latitudes = dataset.createVariable('latitude','f',('lat',))
  longitudes = dataset.createVariable('longitude','f',('lon',))
  times = dataset.createVariable('years', np.float64, ('time',))

  # Global Attributes (will need modified accordingly)
  dataset.description = 'CHIRPSv2 seasonal rainfall at raw 0.05 degree resolution'
  dataset.history = 'Created ' + tt.ctime(tt.time())
  dataset.source = 'Subset by M. Young'
  # Variable Attributes
  latitudes.units = 'degree_north'
  longitudes.units = 'degree_east'
  rfe_out.units = 'mm'
  times.units = 'years'
  # times.calendar = 'gregorian'

  # Fill variables with data
  latitudes[:] = chp_lat#np.flipud(sub_lat)
  longitudes[:] = chp_lon
  rfe_out[:] = all_obs[s,:,:,:]#np.flip(all_data,axis=2)# average in mm h then mmd

  # Fill in time.
  #date = []
  #date = datetime(year, month, day)
  times[:] = years#date2num(date, units = times.units, calendar = times.calendar)
  dataset.close()

  '''
  nc_outfile = f_outdir+'ECMWF_SEAS5.seasonal.point1d.area'+str(box_size)+'d_'+c3s_mon[s]+'.nc'
  # save final array as netcdf
  dataset = Dataset(nc_outfile,'w',format='NETCDF4')
  time = dataset.createDimension('time',len(years)) # create time
  nmem = dataset.createDimension('members',c3s_nmembers) # create lon
  lat = dataset.createDimension('lat',len(c3s_lat)) # create lat (dims depend on region)
  lon = dataset.createDimension('lon',len(c3s_lon)) # create lon
  # create variables
  rfe1_out = dataset.createVariable('rfe1','d',('time','members','lat','lon'))
  rfe5_out = dataset.createVariable('rfe5','d',('time','members','lat','lon'))

  latitudes = dataset.createVariable('latitude','f',('lat',))
  longitudes = dataset.createVariable('longitude','f',('lon',))
  times = dataset.createVariable('years', np.float64, ('time',))

  # Global Attributes (will need modified accordingly)
  dataset.description = 'ECMWF seas5 seasonal rainfall (3 month lead) at 1 degree and average over '+str(box_size)+' degree region surrounding each 1 degree pixel '
  dataset.history = 'Created ' + tt.ctime(tt.time())
  dataset.source = 'Subset by M. Young'
  # Variable Attributes
  latitudes.units = 'degree_north'
  longitudes.units = 'degree_east'
  rfe1_out.units = 'mm'
  rfe5_out.units = 'mm'
  times.units = 'years'
  # times.calendar = 'gregorian'

  # Fill variables with data
  latitudes[:] = c3s_lat#np.flipud(sub_lat)
  longitudes[:] = c3s_lon
  rfe1_out[:] = all_c3s[s,:,:,:,:]#np.flip(all_data,axis=2)# average in mm h then mmd
  rfe5_out[:] = c3s_5d[s,:,:,:,:]#np.flip(all_data,axis=2)# average in mm h then mmd

  # Fill in time.
  #date = []
  #date = datetime(year, month, day)
  times[:] = years#date2num(date, units = times.units, calendar = times.calendar)
  dataset.close()
  '''
