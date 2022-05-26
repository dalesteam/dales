#!/bin/bash
#SBATCH -p thin          # Normal queue 
#SBATCH -n 1               # 1 CPU cores 
#SBATCH -t 01:00:00        # 10 hours of computing time 

module load 2021
module load MATLAB/2021a-upd3	
echo "mcc -m DALES_ProcessOutput001.m" | matlab -nodisplay
./DALES_ProcessOutput001