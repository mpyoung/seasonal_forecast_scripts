#!/bin/bash
export err_dir=//gws/nopw/j04/ncas_climate_vol1/users/vboult/seasonal_forecasting/scripts/bsub_output/

for model in $(seq 0 0); do
  for m in $(seq 10 12); do
    echo ${m}

    echo "
    #!/bin/bash
    #BSUB -o ${err_dir}%J.o
    #BSUB -e ${err_dir}%J.e
    #BSUB -q short-serial
    #BSUB -J matt_py_job
    #BSUB -n 1
    #BSUB -W 10:00
    python download_c3s_seasonal.py $model $m 2>> ${err_dir}err_dl${s}.log >> ${err_dir}out_dl${s}.log
    " >> job_dl${s}.sh
    bsub < job_dl${s}.sh
    rm job_dl${s}.sh

  done
done
