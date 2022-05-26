#!/bin/bash
#SBATCH -p thin          # Normal queue 
#SBATCH -n 1               # 1 CPU cores 
#SBATCH -t 2-00:00:00        # 10 hours of computing time
#SBATCH --mem=24000MB 	# Memory usage per cpu core in MB

module load 2021
module load MATLAB/2021a-upd3
echo "mcc -m DALES_ProcessOutput002.m" | matlab -nodisplay
./DALES_ProcessOutput002