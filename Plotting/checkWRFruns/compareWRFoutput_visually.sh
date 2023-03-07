#!/bin/bash
# Daan Kivits, 2023

# The script creates a temporal average mixing ratio over the testrun period by combing the 'wrfout' WRF output files 
# from the different WRF runs. It compares simulation results from four different WRF runs using the original empty flux set; 
# a high-emission alternative of this original flux set, in which the emission flux at the first vertical level has been set to 
# 2 mol m-2 s-1; the near-real time (NRT) that have been created using our 'FillWRFChemiWithZeros.py' script, and the fluxes 
# that have been provided by Friedemann Reum.

# This script takes three directories as user input: The directory containing the WRF output from a WRF run that used the (1) 
# original fluxes, (2) the near-real time (NRT) that have been created using our 'FillWRFChemiWithZeros.py' script, and (3) 
# the fluxes that have been provided by Friedemann Reum (personal communcations, 2022). This WRF output is stored in the 
# following directories: 
    # ORIG fluxes (1): '/projects/0/ctdas/dkivits/CTDAS_WRF_runs/testrun/WRF_n5/run_testcases_fwd/test_fwd_original/'
    # NRT fluxes (2): '/projects/0/ctdas/dkivits/CTDAS_WRF_runs/testrun/WRF_n5/run_testcases_fwd/test_fwd_nrt/'
    # FRI fluxes (3): '/projects/0/ctdas/dkivits/CTDAS_WRF_runs/testrun/WRF_n5/run_testcases_fwd/test_fwd_friedeman n/'

# Before running this script, make sure you ran the 'FillWRFChemiWithHighEmissions.py' script, and the 'E_CO2_001' variable
# is filled with a constant high-emission flux field (equal to 2 mol m-2 s-1 or 7200 mol km-2 hr-1)!

cdo -O mergetime /projects/0/ctdas/dkivits/CTDAS_WRF_runs/testrun/WRF_n5/run_testcases_fwd/$1/wrfout_d01_2* /projects/0/ctdas/dkivits/CTDAS_WRF_runs/testrun/WRF_n5/run_testcases_fwd/$1/wrfout_all.nc
cdo -O mergetime /projects/0/ctdas/dkivits/CTDAS_WRF_runs/testrun/WRF_n5/run_testcases_fwd/$2/wrfout_d01_2* /projects/0/ctdas/dkivits/CTDAS_WRF_runs/testrun/WRF_n5/run_testcases_fwd/$2/wrfout_all.nc
cdo -O mergetime /projects/0/ctdas/dkivits/CTDAS_WRF_runs/testrun/WRF_n5/run_testcases_fwd/$3/wrfout_d01_2* /projects/0/ctdas/dkivits/CTDAS_WRF_runs/testrun/WRF_n5/run_testcases_fwd/$3/wrfout_all.nc

WRFOUT1=/projects/0/ctdas/dkivits/CTDAS_WRF_runs/testrun/WRF_n5/run_testcases_fwd/$1/wrfout_fluxesonly.nc
WRFOUT2=/projects/0/ctdas/dkivits/CTDAS_WRF_runs/testrun/WRF_n5/run_testcases_fwd/$2/wrfout_fluxesonly.nc
WRFOUT3=/projects/0/ctdas/dkivits/CTDAS_WRF_runs/testrun/WRF_n5/run_testcases_fwd/$3/wrfout_fluxesonly.nc

cdo -O select,name=CO2_001,CO2_002,CO2_003,CO2_004 /projects/0/ctdas/dkivits/CTDAS_WRF_runs/testrun/WRF_n5/run_testcases_fwd/$1/wrfout_all.nc $WRFOUT1
cdo -O select,name=CO2_001,CO2_002,CO2_003,CO2_004 /projects/0/ctdas/dkivits/CTDAS_WRF_runs/testrun/WRF_n5/run_testcases_fwd/$2/wrfout_all.nc $WRFOUT2
cdo -O select,name=CO2_001,CO2_002,CO2_003,CO2_004 /projects/0/ctdas/dkivits/CTDAS_WRF_runs/testrun/WRF_n5/run_testcases_fwd/$3/wrfout_all.nc $WRFOUT3

WRFOUT1=/projects/0/ctdas/dkivits/CTDAS_WRF_runs/testrun/WRF_n5/run_testcases_fwd/$1/wrfout_timeavg.nc
WRFOUT2=/projects/0/ctdas/dkivits/CTDAS_WRF_runs/testrun/WRF_n5/run_testcases_fwd/$2/wrfout_timeavg.nc
WRFOUT3=/projects/0/ctdas/dkivits/CTDAS_WRF_runs/testrun/WRF_n5/run_testcases_fwd/$3/wrfout_timeavg.nc

cdo -O timavg /projects/0/ctdas/dkivits/CTDAS_WRF_runs/testrun/WRF_n5/run_testcases_fwd/$1/wrfout_fluxesonly.nc $WRFOUT1
cdo -O timavg /projects/0/ctdas/dkivits/CTDAS_WRF_runs/testrun/WRF_n5/run_testcases_fwd/$2/wrfout_fluxesonly.nc $WRFOUT2
cdo -O timavg /projects/0/ctdas/dkivits/CTDAS_WRF_runs/testrun/WRF_n5/run_testcases_fwd/$3/wrfout_fluxesonly.nc $WRFOUT3

python /projects/0/ctdas/dkivits/scripts/Plotting/functions/check_if_WRFout_isequal.py $1 $2 $3
