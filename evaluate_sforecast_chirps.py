'''
Reads in regional subsets of seasonal forecast data 
and equivalent observed data (from CHIRPS)

Computes and plots some standard evaluation statistics


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
import scipy as sp
import sys
execfile('date_str.py')

def calc_U(y_true, y_score):
  '''
  Computes the U statistic for calculating the
  statistical significance of the ROC Area
  from https://johaupt.github.io/roc-auc/model%20evaluation/Area_under_ROC_curve.html
  '''
  n1 = np.sum(y_true==1)
  n0 = len(y_score)-n1

  ## Calculate the rank for each observation
  # Get the order: The index of the score at each rank from 0 to n
  order = np.argsort(y_score)
  # Get the rank: The rank of each score at the indices from 0 to n
  rank = np.argsort(order)
  # Python starts at 0, but statistical ranks at 1, so add 1 to every rank
  rank += 1
  # If the rank for target observations is higher than expected for a random model,
  # then a possible reason could be that our model ranks target observations higher
  U1 = np.sum(rank[y_true == 1]) - n1*(n1+1)/2
  U0 = np.sum(rank[y_true == 0]) - n0*(n0+1)/2
  # Formula for the relation between AUC and the U statistic
  AUC1 = U1/ (n1*n0)
  AUC0 = U0/ (n1*n0)
  return U1,AUC1,n1, U0,AUC0,n0


dir_main = '/gws/nopw/j04/ncas_climate_vol1/users/vboult/seasonal_forecasting/'
dir_c3s_region = dir_main+'c3s_hindcasts/region_subsets/'
dir_probs = dir_c3s_region+'forecast_probabilities/'
dir_out = dir_main+'c3s_validation_plots/'

years = np.arange(1993,2016+1,1)

center_ids = ['ecmwf_sys5','ukmo_sys13','meteofrance_sys6','dwd_sys2'] # modelling center names
model_names = ['ECMWF-SEAS5','UKMO','Meteofrance','DWD']
percentile_names = ['Below normal','Near-normal','Above normal']
percentile_short_name = ['Dry','Normal','Wet']
# define percentile categories to analyse
percentiles = [100./3.,200./3.]
p_name = 'tercile' # give categories an appropriate name

region='Kenya'
latlim = [-5,6]
lonlim = [33,43]
model=0
month_start=12
month_end=12
month_lead=0

# open CHIRPS file
c3s_file = dir_c3s_region+'CHIRPS_total_rainfall_'+region+'_s'+mon_string(month_start)+'_e'+mon_string(month_end)+'_lead'+str(month_lead)+'.nc'
nc_fid = Dataset(c3s_file, 'r')
chp_lat = np.array(nc_fid.variables['latitude'][:])  # extract/copy the data
chp_lon = np.array(nc_fid.variables['longitude'][:])
chp_rfe = np.array(nc_fid.variables['rfe'][:])
land_sea_mask = np.copy(np.mean(chp_rfe,axis=0))
land_sea_mask[np.isnan(land_sea_mask)==False]=1

# open forecast file
c3s_file = dir_c3s_region+center_ids[model]+'_total_rainfall_'+region+'_s'+mon_string(month_start)+'_e'+mon_string(month_end)+'_lead'+str(month_lead)+'.nc'
nc_fid = Dataset(c3s_file, 'r')
c3s_lat = np.flipud(np.array(nc_fid.variables['latitude'][:]))  # extract/copy the data
c3s_lon = np.array(nc_fid.variables['longitude'][:])
c3s_rfe = np.array(nc_fid.variables['rfe'][:])*land_sea_mask
nc_fid.close()

'''
Calculate some validation statistics
'''

# plot mean and standard devation rainfall
cols = 'PuBu'
cmin = 0
cmax = 1000
cspc = 200
clevs = np.arange(cmin,cmax+cspc,cspc)
clabel1 = 'Mean total rainfall (mm)'
norm = BoundaryNorm(boundaries=clevs,ncolors=256)
lw = 1
size1 =[5,3]
fig = plt.figure(figsize=(size1[0],size1[1]))

plt.subplot(1,2,1)
mymap = Basemap(projection='cyl',resolution='l',\
        llcrnrlat=np.min(c3s_lat),urcrnrlat=np.max(c3s_lat),\
        llcrnrlon=np.min(c3s_lon),urcrnrlon=np.max(c3s_lon))
mymap.drawparallels(np.arange(-90,90,2),labels=[1,0,0,0],labelstyle='+/-')
mymap.drawmeridians(np.arange(0,360,2),labels=[0,0,0,1],labelstyle='+/-')
mymap.drawcoastlines(linewidth=1)
mymap.drawcountries(linewidth=1)
x, y = mymap(*np.meshgrid(c3s_lon,c3s_lat))
uncal=mymap.pcolormesh(x,y,np.nanmean(chp_rfe,axis=0),vmin=cmin,vmax=cmax,norm=norm,cmap=cols)
plt.title('CHIRPS')

plt.subplot(1,2,2)
mymap = Basemap(projection='cyl',resolution='l',\
        llcrnrlat=np.min(c3s_lat),urcrnrlat=np.max(c3s_lat),\
        llcrnrlon=np.min(c3s_lon),urcrnrlon=np.max(c3s_lon))
mymap.drawparallels(np.arange(-90,90,2),labels=[0,0,0,0],labelstyle='+/-')
mymap.drawmeridians(np.arange(0,360,2),labels=[0,0,0,1],labelstyle='+/-')
mymap.drawcoastlines(linewidth=1)
mymap.drawcountries(linewidth=1)
x, y = mymap(*np.meshgrid(c3s_lon,c3s_lat))
mymap.pcolormesh(x,y,np.nanmean(c3s_rfe,axis=(0,1)),vmin=cmin,vmax=cmax,norm=norm,cmap=cols)
plt.title(model_names[model])

plt.suptitle(str(month_start)+'-'+str(month_end))
plt.tight_layout(pad=2.5,w_pad=0.02,h_pad=0.5)
fig.subplots_adjust(right=0.90)
cbar_pos = [0.92, 0.15, 0.025, 0.65] #[left, bottom, width, height]
cbar_ax = fig.add_axes(cbar_pos)
cbar = fig.colorbar(uncal,cax=cbar_ax,label=clabel1,extend='max')

fname_plot = dir_out+center_ids[model]+'_mean_total_rainfall_'+region+'_s'+mon_string(month_start)+'_e'+mon_string(month_end)+'_lead'+str(month_lead)
plt.savefig(fname_plot+'.png',bbox_inches='tight',dpi=300)
#plt.show()
plt.close()

# plot bias
# mean difference between totals from model ensemble mean and the obs
bias = np.nanmean(np.nanmean(c3s_rfe,axis=1)-chp_rfe,axis=0)

cols = 'RdBu'
cmin = -200
cmax = 200
cspc = 50
clevs = np.arange(cmin,cmax+cspc,cspc)
clabel1 = 'Mean total rainfall (mm)'
norm = BoundaryNorm(boundaries=clevs,ncolors=256)

fig = plt.figure(figsize=(4,3))
mymap = Basemap(projection='cyl',resolution='l',\
        llcrnrlat=np.min(c3s_lat),urcrnrlat=np.max(c3s_lat),\
        llcrnrlon=np.min(c3s_lon),urcrnrlon=np.max(c3s_lon))
mymap.drawparallels(np.arange(-90,90,2),labels=[1,0,0,0],labelstyle='+/-')
mymap.drawmeridians(np.arange(0,360,2),labels=[0,0,0,1],labelstyle='+/-')
mymap.drawcoastlines(linewidth=1)
mymap.drawcountries(linewidth=1)
x, y = mymap(*np.meshgrid(c3s_lon,c3s_lat))
# mymap.pcolormesh(x,y,bias,vmin=cmin,vmax=cmax,norm=norm,cmap=cols)
mymap.contourf(x,y,bias,clevs,cmap=cols,extend='both')

plt.title(str(month_start)+'-'+str(month_end)+'\n'+model_names[model]+' - CHIRPS')
plt.colorbar(label='Bias (mm)')
plt.tight_layout()
fname_plot = dir_out+center_ids[model]+'_bias_total_rainfall_'+region+'_s'+mon_string(month_start)+'_e'+mon_string(month_end)+'_lead'+str(month_lead)
plt.savefig(fname_plot+'.png',bbox_inches='tight',dpi=300)
#plt.show()
plt.close()

# plot anomaly correlation coefficient
ens_mean = np.nanmean(c3s_rfe,axis=1)# model ensemble mean
c3s_acc = np.nan*np.zeros((2,len(c3s_lat),len(c3s_lon)))
c3s_anom = np.nan*np.zeros((len(years),len(c3s_lat),len(c3s_lon)))
chp_anom = np.nan*np.zeros((len(years),len(c3s_lat),len(c3s_lon)))

# compute anomalies independant of climatology
for y in np.arange(0,len(years)):
  tmp_clim = []
  tmp_clim = np.nanmean(np.delete(chp_rfe,y,axis=0),axis=0)
  chp_anom[y,:,:] = chp_rfe[y,:,:] - tmp_clim

  tmp_clim = []
  tmp_clim = np.nanmean(np.delete(ens_mean,y,axis=0),axis=0)
  c3s_anom[y,:,:] = ens_mean[y,:,:] - tmp_clim
# compute anomaly correlation coefficient (using pearsonr)
# pearson r doesn't work on 2d array so loop through gridboxes
for i in np.arange(0,len(c3s_lat)):
  for j in np.arange(0,len(c3s_lon)):
    if np.sum(np.isnan(c3s_anom[:,i,j])) < len(years): # check for nan
      c3s_acc[:,i,j] = pearsonr(chp_anom[:,i,j],c3s_anom[:,i,j])


cols = 'RdYlBu'
cmin = 0
cmax = 1
cspc = 0.1
clevs = np.arange(cmin,cmax+cspc,cspc)
norm = BoundaryNorm(boundaries=clevs,ncolors=256)

fig = plt.figure(figsize=(4,3))
mymap = Basemap(projection='cyl',resolution='l',\
        llcrnrlat=np.min(c3s_lat),urcrnrlat=np.max(c3s_lat),\
        llcrnrlon=np.min(c3s_lon),urcrnrlon=np.max(c3s_lon))
mymap.drawparallels(np.arange(-90,90,2),labels=[1,0,0,0],labelstyle='+/-')
mymap.drawmeridians(np.arange(0,360,2),labels=[0,0,0,1],labelstyle='+/-')
mymap.drawcoastlines(linewidth=1)
mymap.drawcountries(linewidth=1)
x, y = mymap(*np.meshgrid(c3s_lon,c3s_lat))
uncal=mymap.pcolormesh(x,y,c3s_acc[0,:,:],vmin=cmin,vmax=cmax,norm=norm,cmap=cols)
# mymap.contourf(x,y,c3s_acc[0,:,:],clevs,cmap=cols,extend='both')
# plot significance hatching
for i in np.arange(0,len(c3s_lat)):
  for j in np.arange(0,len(c3s_lon)):
    if c3s_acc[1,i,j] < 0.05:
      plt.annotate('x',(c3s_lon[j]+0.25,c3s_lat[i]+0.25), xytext=(c3s_lon[j]+0.25,c3s_lat[i]+0.25), fontsize=5, color='k')

plt.title(str(month_start)+'-'+str(month_end)+'\n'+model_names[model]+' - CHIRPS')
plt.colorbar(uncal,label='ACC')
plt.tight_layout()

fname_plot = dir_out+center_ids[model]+'_ACC_total_rainfall_'+region+'_s'+mon_string(month_start)+'_e'+mon_string(month_end)+'_lead'+str(month_lead)
plt.savefig(fname_plot+'.png',bbox_inches='tight',dpi=300)
#plt.show()
plt.close()

# probabilistic skill evaluation
# Load probabilities
c3s_file = dir_probs+'CHIRPS_total_rainfall_'+p_name+'_outcomes_'+region+'_s'+mon_string(month_start)+'_e'+mon_string(month_end)+'_lead'+str(month_lead)+'.nc'
nc_fid = Dataset(c3s_file, 'r')
chp_p = np.array(nc_fid.variables['p'][:])
nc_fid.close()

# open forecast file
c3s_file = dir_probs+center_ids[model]+'_total_rainfall_'+p_name+'_probabilities_'+region+'_s'+mon_string(month_start)+'_e'+mon_string(month_end)+'_lead'+str(month_lead)+'.nc'
nc_fid = Dataset(c3s_file, 'r')
c3s_p = np.array(nc_fid.variables['p'][:])*land_sea_mask
nc_fid.close

# ROC Area
c3s_roc = np.nan*np.zeros((3,len(percentiles)+1,len(c3s_lat),len(c3s_lon)))
for ter in np.arange(0,len(percentiles)+1):
  for i in np.arange(0,len(c3s_lat)):
    for j in np.arange(0,len(c3s_lon)):
      tmp_ref = np.copy(chp_p[:,ter,i,j].flatten())
      tmp_est = np.copy(c3s_p[:,ter,i,j].flatten())
      f_ref = tmp_ref[(np.isnan(tmp_ref)==False) & (np.isnan(tmp_est)==False)]
      f_est = tmp_est[(np.isnan(tmp_ref)==False) & (np.isnan(tmp_est)==False)]
      if (len(np.where(f_ref > 0)[0]) > 0) & (len(np.where(f_ref == 1)[0]) < len(f_ref)):
        u1,auc1,n1,u0,auc0,n0 = calc_U(f_ref,f_est)
        # wikipedia says that the smaller value of u is used to consult signif tables
        # assume U is normally distributed
        # with a mean of
        mean_u = (n0*n1)/2.
        # and a standard dev of:
        std_u = np.sqrt((n0*n1*(n0+n1+1))/12.0)
        z_score = (np.min([u1,u0]) - mean_u)/std_u
        # use z-score to get p-values
        p_values_1side =[]
        p_values_1side = sp.stats.norm.sf(abs(z_score)) #one-sided
        c3s_roc[0,ter,i,j] = np.copy(auc1)
        c3s_roc[1,ter,i,j] = np.copy(p_values_1side)
        c3s_roc[2,ter,i,j] = roc_auc_score(f_ref,f_est) # should be ~same as auc1 (just to double check auc1 method works!)


cols = 'RdYlBu'
cmin = 0
cmax = 1
cspc = 0.1
clevs = np.arange(cmin,cmax+cspc,cspc)
norm = BoundaryNorm(boundaries=clevs,ncolors=256)

size1 =[8,3]
fig = plt.figure(figsize=(size1[0],size1[1]))
for ter in np.arange(0,len(percentiles)+1):
  plt.subplot(1,len(percentiles)+1,ter+1)
  mymap = Basemap(projection='cyl',resolution='l',\
          llcrnrlat=np.min(c3s_lat),urcrnrlat=np.max(c3s_lat),\
          llcrnrlon=np.min(c3s_lon),urcrnrlon=np.max(c3s_lon))
  if ter == 0:
    mymap.drawparallels(np.arange(-90,90,2),labels=[1,0,0,0],labelstyle='+/-')
  mymap.drawparallels(np.arange(-90,90,2),labels=[0,0,0,0],labelstyle='+/-')
  mymap.drawmeridians(np.arange(0,360,2),labels=[0,0,0,1],labelstyle='+/-')
  mymap.drawcoastlines(linewidth=1)
  mymap.drawcountries(linewidth=1)
  x, y = mymap(*np.meshgrid(c3s_lon,c3s_lat))
  uncal=mymap.pcolormesh(x,y,c3s_roc[2,ter,:,:],vmin=cmin,vmax=cmax,norm=norm,cmap=cols)
  for i in np.arange(0,len(c3s_lat)):
    for j in np.arange(0,len(c3s_lon)):
      if c3s_roc[1,ter,i,j] < 0.05:
        plt.annotate('x',(c3s_lon[j]+0.25,c3s_lat[i]+0.25), xytext=(c3s_lon[j]+0.25,c3s_lat[i]+0.25), fontsize=5, color='k')
  plt.title(percentile_names[ter])

plt.tight_layout()
fig.subplots_adjust(right=0.90)
cbar_pos = [0.92, 0.15, 0.025, 0.65] #[left, bottom, width, height]
cbar_ax = fig.add_axes(cbar_pos)
cbar = fig.colorbar(uncal,cax=cbar_ax,label='ROC Area')

fname_plot = dir_out+center_ids[model]+'_ROCArea_'+p_name+'_'+region+'_s'+mon_string(month_start)+'_e'+mon_string(month_end)+'_lead'+str(month_lead)
plt.savefig(fname_plot+'.png',bbox_inches='tight',dpi=300)
#plt.show()
plt.close()
