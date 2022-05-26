#!/bin/bash
#SBATCH -p thin          # thin queue 
#SBATCH -n 300              # 300 CPU cores 
#SBATCH -t 1-20:00:00        # 4 days hours of computing time 
#SBATCH -e errorfile.001

module purge
module load 2021
module load netCDF-Fortran/4.5.3-gompi-2021a

srun ../dales4 namoptions.001

