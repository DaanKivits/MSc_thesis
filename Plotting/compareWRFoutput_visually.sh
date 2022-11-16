#!/bin/bash

# WRF output for following fluxes:
# ORIG fluxes: '/projects/0/ctdas/dkivits/CTDAS_WRF_runs/testrun/WRF_n5/run_testcases_fwd/test_fwd_original/'
# NRT fluxes: '/projects/0/ctdas/dkivits/CTDAS_WRF_runs/testrun/WRF_n5/run_testcases_fwd/test_fwd_nrt/'
# FRI fluxes: '/projects/0/ctdas/dkivits/CTDAS_WRF_runs/testrun/WRF_n5/run_testcases_fwd/test_fwd_friedemann/'

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

python check_if_WRFout_isequal.py $1 $2 $3
