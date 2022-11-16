#!/bin/bash
#SBATCH -p thin
#SBATCH -t 120:00:00 
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mail-user=daan.kivits@wur.nl
#SBATCH --mail-type=FAIL,END
#SBATCH --mem=60000M

#module load 2021
#module load Anaconda3/2021.05
#module load netCDF-Fortran/4.5.3-gompi-2021a

source /projects/0/ctdas/dkivits/scripts/.bashrc_wrf

python run_ctdas_withrunner_cluster.py /projects/0/ctdas/dkivits/CTDAS_WRF_runs/testrun/settings_file_test_fwd_nrt.pkl --overwrite-ctdas=True >&log

