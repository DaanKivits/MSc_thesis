import netCDF4 as nc
import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import time
import glob
import xarray as xr
import os
import shutil
from dateutil.relativedelta import relativedelta

# Check if a regridded flux directory exists, otherwise create it
dirname = '/projects/0/ctdas/awoude/NRT/ICOS_OUTPUT/'
outdir = '/projects/0/ctdas/dkivits/CTDAS_WRF_runs/testrun/input/test_fwd_friedemann/emissions/fwd/'
#outdir = '/projects/0/ctdas/dkivits/CTDAS_WRF_runs/testrun/input/test_fwd_nrt/emissions/fwd/'
backupdir = '/projects/0/ctdas/dkivits/CTDAS_WRF_runs/testrun/input/test_fwd_original/emissions/fwd_backup/'

if os.path.exists(outdir):
    shutil.rmtree(outdir)
    shutil.copytree(backupdir, outdir)
else:
    shutil.copytree(backupdir, outdir)
 
# Load WRF grid data
geopath = '/projects/0/ctdas/dkivits/CTDAS_WRF_runs/testrun/WPS/WPS_test_fwd_v01/geo_em.d01.nc'
geo = nc.Dataset(geopath,'r')
wrflat = geo.variables['XLAT_M'][0,:,:]
wrflon = geo.variables['XLONG_M'][0,:,:]
we  = geo.variables['XLAT_M'].shape[2]
sn  = geo.variables['XLAT_M'].shape[1]
lon = (geo.variables['XLONG_U'][0,sn/2,we/2] + geo.variables['XLONG_U'][0,sn/2-1,we/2])/2.
lat = (geo.variables['XLAT_V'][0,sn/2,we/2] + geo.variables['XLAT_V'][0,sn/2,we/2-1])/2.
dx_dom = 100000

# Define basemap to prepare for transformation later
m = Basemap(width=dx_dom*we,height=dx_dom*sn,
            rsphere=(6378137.00,6356752.3142),\
            resolution='l',area_thresh=1000.,projection='lcc',\
            lat_1=45.,lat_2=55,lat_0=lat,lon_0=lon)

# Define names of tracers and flux fields
variablelist = ['E_CO2_001','E_CO2_002','E_CO2_003','E_CO2_004']
fluxnamelist=['nep','fire','ocean','anthropogenic']

# Define which files to loop over
#fluxstring = []
#for name in fluxnamelist:
#    fluxstring += sorted(glob.glob(dirname + name + '.2022*.nc'))

outfiles = sorted(glob.glob(outdir + 'wrfchemi_d01_2*'))

# Define CTE-HR grid variables
lon_bounds = [-15.,35.] # in degrees 
lat_bounds = [33.,72.] # in degrees 
lon_step = 0.2 # in degrees 
lat_step = 0.1 # in degrees 
nx = 15 # gridsize in x dimension (amount of gridcells)
ny = 12 # gridsize in y dimension (amount of gridcells)
interpolation_method = 0 # Default = 1; 0 = nearest neighbour, 1 = bilinear interpolation

fluxstring = ['/projects/0/ctdas/awoude/NRT/ICOS_OUTPUT/nep.202206.nc']
CTEHR_nc =  nc.Dataset(fluxstring[0],'r')
timevar = CTEHR_nc.variables['time'][:]

# Loop over flux files
for file in outfiles[0:4]:
       
    # Status report
    print('Busy with ... ' + file)

    # Read CTE-HR flux file
    flux_nc= nc.Dataset(file,'r+')
   
    # Select correct variable and transform fluxes
    for time in range(0,len(timevar)):
        flux_data = CTEHR_nc.variables['nep'][time,:,:]
        flux_trans = m.transform_scalar(flux_data,lons=np.arange(lon_bounds[0],lon_bounds[1],lon_step),lats=np.arange(lat_bounds[0],lat_bounds[1],lat_step),nx=nx,ny=ny,order=interpolation_method)     
            
        # Make sure all NaN values are converted to zeroes (required by Snellius)     
        flux_trans = np.where(np.isnan(flux_trans), 0, flux_trans)
            
        # Apply a Landmask over the converted flux data:
        flux_trans = np.where(geo.variables['LANDMASK'][0,:,:] != 0, flux_trans, 0)
            
        # Put the transformed CTE-HR fluxes into the NC file
        #  fossilvar[:,:,:,:] = 0
        for variable in variablelist:
            flux_nc.variables[variable][0,:,:,:] = 0.0
            flux_nc.variables[variable][0,1,:,:] = flux_trans * 3600 * 1e6
        # Check if dimensions of transformed
        #print(fossilvar.shape)
        #print(flux_trans.shape)     

    flux_nc.close()
CTEHR_nc.close()
print('Done changing the emission fluxes!\n Files have been written to: ' + outdir + ' !')
