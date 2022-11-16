import netCDF4 as nc
import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import time
import glob
import xarray as xr
import os
from datetime import datetime, timedelta, time, date
import shutil
from dateutil.relativedelta import relativedelta

def create_file_from_source(src_file, trg_file):
    src = nc.Dataset(src_file)
    trg = nc.Dataset(trg_file, mode='w')

    # Create the dimensions of the file
    for name, dim in src.dimensions.items():
        trg.createDimension(name, len(dim) if not dim.isunlimited() else None)

    # Copy the global attributes
    trg.setncatts({a:src.getncattr(a) for a in src.ncattrs()})

    # Create the variables in the file
    for name, var in src.variables.items():
        trg.createVariable(name, var.dtype, var.dimensions)

        # Copy the variable attributes
        trg.variables[name].setncatts({a:var.getncattr(a) for a in var.ncattrs()})

        # Copy the variables values (as 'f4' eventually)
        trg.variables[name][:] = src.variables[name][:]

    # Save the file
    trg.close()

# Check if a regridded flux directory exists, otherwise create it
dirname = '/projects/0/ctdas/awoude/NRT/ICOS_OUTPUT/'
fluxdir = '/projects/0/ctdas/dkivits/CTDAS_WRF_runs/testrun/input/test_fwd_nrt/emissions/'
mode = 'fwd'
outdir = fluxdir + mode + '/'
backupdir = fluxdir + mode + '_backup/'

if not os.path.exists(outdir):
    os.mkdir(outdir)

# Define which files to loop over
#fluxstring = []
#for name in fluxfilelist:
#    fluxstring += sorted(glob.glob(dirname + name + '.202206.nc'))
fluxstring = ['/projects/0/ctdas/awoude/NRT/ICOS_OUTPUT/nep.202206.nc']

infiles = sorted(glob.glob(backupdir + '/wrfchemi_d01_2*'))

CTEHR_nc =  nc.Dataset(fluxstring[0],'r')
timevar = CTEHR_nc.variables['time'][:]

# Define names of tracers and flux fields
variablelist = ['E_CO2_001','E_CO2_002','E_CO2_003','E_CO2_004']
fluxfilelist = ['nep','fire','ocean','anthropogenic']
fluxnamelist=['nep','fire','ocean','combustion']
domain = 1
basedate = datetime(2000,1,1,0,0,0,0)

for times in range(0, len(timevar)):
        
        # Loop over CTE_HR time variables and create emission flux files according to these times
        curtime= timevar[times]
        time = basedate + timedelta(seconds=int(curtime)) - relativedelta(years=7)
        print('Busy with ... ' + str(time))
        new_filename = outdir + 'wrfchemi_d0' + str(domain) + '_' + str(time.date()) + '_' + str(time.time())
        old_filename = infiles[times]
        create_file_from_source(old_filename, new_filename)

for file in fluxstring:
    # Progress report
    print('WORKING ON ... ' + file)

    # Open CTE-HR flux file, extract time to loop over at later stage
    CTEHR_nc =  nc.Dataset(file,'r')
    timevar = CTEHR_nc.variables['time'][:]


