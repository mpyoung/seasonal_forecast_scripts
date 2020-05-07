'''
Script downloads seasonal hindcast data from the
ECMWF/EU Copernicus Climate Change Service Climate Data Store

Note that the python package 'cdsapi' is needed and a cds api
key must be installed locally in order for the forecast downloads
to work. See following link for further details:
https://cds.climate.copernicus.eu/api-how-to

M. Young 7/05/2020
'''
#!/usr/bin/env python
import cdsapi
import datetime
import numpy as np
import sys

c = cdsapi.Client()

execfile('date_str.py')


# outputdir = "/gws/nopw/j04/ncas_climate_vol1/users/myoung02/datasets/C3S/hindcasts/"
outputdir = "/gws/nopw/j04/ncas_climate_vol1/users/vboult/seasonal_forecasting/c3s_hindcasts/global/"
# List the forecast dates # date-month
months = np.arange(1,13,1)#[3]#[3,6,9,12]#np.arange(3,10+1,1)
past_years = np.array(np.arange(1993,2016+1,1),dtype="S")
past_years = []
for y in np.arange(1993,2016+1,1):
  past_years.append(str(y))

# Some notes about downloading from Copernicus:
# Do NOT ask for different start months in the same request when retrieving data from the "daily" datasets.
# It is much more efficient to request one single month (e.g. November) for all the hindcast years (1993 to 2016) than requesting all 12 months from a single year
#for month in months:
'''
Before running check which model to download
then run script individually for each model

**** NOTE (Matt 07/05/2020) ******:
Currently the leadtime_hours in the data I have downloaded
only go out to t= 2880 hrs. I have added extra
hours in this script in case interested in looking at longer lead forecasts.
However, to get these the scripts would have to be run again to
download the longer lead data.

'''
system_ls = [5,13,6,2] #,3]
model_ls = ['ecmwf','ukmo','meteo_france','dwd'] #,"cmcc"]
model_name = ['ecmwf_sys5','ukmo_sys13','meteofrance_sys6','dwd_sys2'] #,"cmcc_sys3"]

def do_stuff(model,month):
  month_str = mon_string(month)
  output_fname = outputdir+"c3s_"+model_name[model]+"_seas_"+past_years[0]+"_"+past_years[-1]+"_"+month_str+".nc"
  print output_fname
  c.retrieve("seasonal-original-single-levels",
            {"format":"netcdf",
             "originating_centre":model_ls[model],
             "system":str(system_ls[model]),
             "variable":"total_precipitation",
             "year":past_years,
             "month":month_str,
             "day":"01",
             "leadtime_hour":["720","1440","2160","2880","3600","4320","5040"]},
             output_fname)

if __name__ == "__main__":
   output = do_stuff(int(sys.argv[1]),int(sys.argv[2]))
