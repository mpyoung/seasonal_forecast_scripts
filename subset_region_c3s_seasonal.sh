#!/bin/bash
#setup nco
export inputdir=/gws/nopw/j04/ncas_climate_vol1/users/vboult/seasonal_forecasting/c3s_hindcasts/global/
export outputdir=/gws/nopw/j04/ncas_climate_vol1/users/vboult/seasonal_forecasting/c3s_hindcasts/region_subsets/
export model_names=('ecmwf_sys5' 'ukmo_sys13' 'meteofrance_sys6' 'dwd_sys2')
export months=($(seq 10 12))
export region_name=Kenya
export lonlim=(33 43)
export latlim=(-5 5)

for t in $(seq 0 0); do
  if [ ${months[$t]} -lt 10 ]; then
    currmonth=0${months[$t]}
  else
    currmonth=${months[$t]}
  fi
  for m in $(seq 0 0); do
    echo ${model_names[$m]}
    infile=${inputdir}c3s_${model_names[$m]}_seas_1993_2016_${currmonth}.nc
    outfile=${outputdir}c3s_${model_names[$m]}_seas_1993_2016_${currmonth}_${region_name}.nc
    ncks -d longitude,${lonlim[0]},${lonlim[1]} -d latitude,${latlim[0]},${latlim[1]} $infile $outfile
    # ncks -d longitude,33,43 -d latitude,-5,5 $infile $outfile

  done
done
