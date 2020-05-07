import geopandas as gpd
import numpy as np
import shapely
from netCDF4 import Dataset
import sys
# sys.path.insert(0, "/gws/nopw/j04/ncas_climate_vol1/users/vboult/tamsat_alert/ta_code")
# from sm_gridded_utils import *
dir_main = '/gws/nopw/j04/ncas_climate_vol1/users/vboult/seasonal_forecasting/'

shapefile = dir_main+"/GIS/Kenya.shp"

# just loading in tamsat rainfall data as working example
nc = Dataset(dir_main+'/prate_tamsat_all.nc')
rfe = nc.variables["rfe"][:]
lons = nc.variables["longitude"][:]
lats = nc.variables["latitude"][:]

# preprocessing of rfe data
rfe[rfe == np.max(rfe)] = np.nan
# rfe = swapdata(rfe)
# rfe = reshape_hist_data(rfe, 1983)

def mask_4d_data(varin):
  '''
  varin should be a 4d numpy array

  This will:
  1) find the mid point of each 0.25 degree grid cell
  2) check if each mid point falls within the boundary of the shapefile
  3) create a 2d mask based on this
  4) expand dims of mask to fit 4d array
  5) return varin as a masked 4d array
  '''
  # Read in boundary shapefile
  mask_shape = gpd.read_file(shapefile)
  # Find mid-points of grid cells
  lons_mid = lons + 0.125
  lats_mid = lats + 0.125
  # Make meshgrid of mid-points
  [lon2d,lat2d] = np.meshgrid(lons_mid,lats_mid)
  lon2 = lon2d.reshape(-1)
  lat2 = lat2d.reshape(-1)
  # Create empty mask
  mask = []
  # Fill mask with true / false
  for lat, lon in zip(lat2, lon2):
    this_point = shapely.geometry.Point(lon, lat)
    res = mask_shape.contains(this_point)
    mask.append(res.values[0])
  # Reshape mask
  mask = np.array(mask).reshape(lon2d.shape)
  mask = ~mask
  mask = mask.T
  # Fit mask to varin shape
  ### EDIT HERE IF YOU HAVE A NUMPY ARRAY OF DIFFERENT DIMENSIONS
  ### Basically repeat the 2d mask over your required dimensions
  mask3d = mask[:,:,np.newaxis]
  mask3d = np.repeat(mask3d, varin.shape[2], axis = 2)
  mask4d = mask3d[:,:,:,np.newaxis]
  mask4d = np.repeat(mask4d, varin.shape[3], axis = 3)
  # Mask var
  varin_masked = np.ma.masked_array(varin, mask = mask4d)
  # Return masked array
  return varin_masked

# mask data
rfe_masked = mask_4d_data(rfe)

# plot to check it looks right - should look like the shape of kenya
plt.pcolor(rfe_masked[:,:,0,0].T)
plt.show()
