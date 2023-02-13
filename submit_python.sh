#!/bin/bash
#SBATCH -p thin
#SBATCH -t 6:00:00 
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mail-user=daan.kivits@wur.nl
#SBATCH --mail-type=FAIL,END
#SBATCH --mem=60000M

module purge
module load 2021
module load Anaconda3/2021.05
module load intel/2021a
module load intel-compilers/2021.2.0
module load impi/2021.2.0-intel-compilers-2021.2.0
module load iimpi/2021a
module load UDUNITS/2.2.28-GCCcore-10.3.0
module load time/1.9-GCCcore-10.3.0
module load netCDF-Fortran/4.5.3-gompi-2021a
#module load netCDF/4.8.0-gompi-2021a

source /projects/0/ctdas/dkivits/scripts/.bashrc_wrf

python run_ctdas_withrunner_cluster.py /projects/0/ctdas/dkivits/CTDAS_WRF_runs/testrun/settings_file_test_fwd_highres.pkl --overwrite-ctdas=True >&log

