# Daan Kivits, 2023

# This script fills (already existing) WRFchemi emission files with a constant high-emission flux field of 
# 2 mol m-2 s-1 (or 7200 mol km-2 hr-1). The rest of the fields are set to 0!
# Some points of attention:
    # - A set of WRFchemi emission files already exists and can be copied from another directory (in this case the 'backupdir)

##############################################
########## LOAD NECCESSARY PACKAGES ##########
##############################################
import netCDF4 as nc
import numpy as np
from mpl_toolkits.basemap import Basemap
import glob
import os
import shutil

########################################
########## CREATE DIRECTORIES ##########
########################################
# Check if a backup flux field directory exists, otherwise create it
if not os.path.exists(backupdir):
   shutil.copytree(outdir, backupdir)

# To reset the files to the backup standard, copy backup directory to outdir always
shutil.copytree(backupdir, outdir, dirs_exist_ok=True)

######################################
########## CREATE VARIABLES ##########
######################################
# Give the syntax of the filenames, make sure not to include the template files
name = 'wrfchemi_d01_2015*'

# Create run-specific information
dirname = '/projects/0/ctdas/dkivits/CTDAS_WRF_runs/testrun/input/test_fwd_original/emissions/'
mode = 'fwd'
outdir = dirname + mode
backupdir = dirname + mode + '_backup/'
fluxstring =  sorted(glob.glob(backupdir + name))
variablelist = ['E_CO2','E_CO2_proc_1','E_CO2_proc_2','E_CO2_proc_3','e_CO2_000','E_CO2_001','E_CO2_002','E_CO2_003','E_CO2_004']

# Loop over the files
for files in fluxstring:
    # Defining new filename and progress report
    filename = os.path.basename(files)
    newstring = outdir + '/' + filename
    print('Busy with ... ' + filename)
    
    # Data handling
    wrfchemi_nc = nc.Dataset(newstring,'r+',format='NETCDF4')
    for variable in variablelist:
         if variable == 'E_CO2_001':
            wrfchemi_nc.variables[variable][:,:,:,:] = 0
            wrfchemi_nc.variables[variable][:,0,:,:] = 7200 # calculated as follows: 2e-6 * 3600 (from mol km-2 hr-1 to mol m-2 s-1)
         else:
            wrfchemi_nc.variables[variable][:,:,:,:] = 0
    wrfchemi_nc.close()

print('Done changing the emission fluxes!\n Files have been written to: ' + outdir + ' !')
