#!/bin/bash
#SBATCH -p thin
#SBATCH -t 04:00:00
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mail-user=daan.kivits@wur.nl
#SBATCH --mail-type=FAIL,END
#SBATCH --mem=60000M

module load 2021
module load Anaconda3/2021.05
module load netCDF-Fortran/4.5.3-gompi-2021a

source /projects/0/ctdas/dkivits/scripts/.bashrc_wrf

python MakeWRFChemiFiles.py >& MakeWRFChemiFiles.log

