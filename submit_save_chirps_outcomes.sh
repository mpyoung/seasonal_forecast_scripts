#!/bin/bash
export err_dir=/gws/nopw/j04/ncas_climate_vol1/users/vboult/seasonal_forecasting/scripts/bsub_output/

# Specify forecast start, end (for accumulation period)
# and lead time (all in months)
month_start=(10 11 12)
month_end=(12 12 12)
month_lead=0

# Specify region
# min/max longitude and latitude boundaries of region
region='Kenya'

for m in $(seq 0 2); do # loop through month_start/end
  echo ${m}

  echo "
  #!/bin/bash
  #BSUB -o ${err_dir}%J.o
  #BSUB -e ${err_dir}%J.e
  #BSUB -q short-serial
  #BSUB -J matt_py_job
  #BSUB -n 1
  #BSUB -W 2:00
  python save_chirps_outcomes.py ${month_start[$m]} ${month_end[$m]} $month_lead $region 2>> ${err_dir}err_chp_outcomes_${m}.log >> ${err_dir}out_chp_outcomes_${m}.log
  " >> job_save.sh
  bsub < job_save.sh
  rm job_save.sh

done
