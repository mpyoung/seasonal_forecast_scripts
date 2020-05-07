#!/bin/bash
export err_dir=/gws/nopw/j04/ncas_climate_vol1/users/vboult/seasonal_forecasting/scripts/bsub_output/

tmp_dir=/gws/nopw/j04/ncas_climate_vol1/users/vboult/seasonal_forecasting/c3s_hindcasts/region_subsets/forecast_probabilities/tmp_probs/
model_n=('ecmwf_sys5' 'ukmo_sys13' 'meteofrance_sys6' 'dwd_sys2')
years=($(seq 1993 2016))


# Specify forecast start, end (for accumulation period)
# and lead time (all in months)
month_start=(10 11 12)
month_end=(12 12 12)
month_lead=0
# Specify region
# min/max longitude and latitude boundaries of region
region='Kenya'
# model order is 0=ECMWF,1=UKMO,2=Meteofrance,3=DWD
for model in $(seq 0 0); do
  for m in $(seq 0 2); do # loop through month_start/end
    echo ${m}
    for y in $(seq 0 23); do
     for e in $(seq 0 24); do
       file=${tmp_dir}tmp_${model_n[$model]}_total_rainfall_${region}_s${month_start[$m]}_e${month_end[$m]}_lead${month_lead}_p_${years[$y]}_${e}.npy

        if [ -e "$file" ]; then
          echo "File exists ${model_n[$model]}_total_rainfall_${region}_s${month_start[$m]}_e${month_end[$m]}_lead${month_lead}_p_${years[$y]}_${e}"
        else
          echo "running ${model_n[$model]}_total_rainfall_${region}_s${month_start[$m]}_e${month_end[$m]}_lead${month_lead}_p_${years[$y]}_${e}"
          echo "
          #!/bin/bash
          #BSUB -o ${err_dir}%J.o
          #BSUB -e ${err_dir}%J.e
          #BSUB -q short-serial
          #BSUB -J matt_py_job
          #BSUB -n 1
          #BSUB -W 0:30
          python save_sforecast_outcomes.py $model ${month_start[$m]} ${month_end[$m]} $month_lead $region $y $e 2>> ${err_dir}err_outcomes_${model}_${m}_y${y}_e${e}.log >> ${err_dir}out_outcomes_${model}_${m}_y${y}_e${e}.log
          " >> job_save.sh
          bsub < job_save.sh
          rm job_save.sh
        fi
      done
    done
  done
done
