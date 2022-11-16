import netCDF4 as nc
import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import time
import glob
import xarray as xr
import os

# Check if a regridded flux directory exists, otherwise create it
dirname = '/projects/0/ctdas/awoude/NRT/ICOS_OUTPUT/'
outdir = '/projects/0/ctdas/dkivits/DATA/CTE_regriddedfluxes/'
if not os.path.exists(outdir):
    os.mkdir(outdir)

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
variablelist = ['E_CO2_001','E_CO2_002','E_CO2_003, E_CO2_004']
fluxnamelist=['nep','fire','ocean','anthropogenic']

# Define which files to loop over
fluxstring = []
for name in fluxnamelist:
    fluxstring += sorted(glob.glob(dirname + name + '.2022*.nc'))

#fluxstring = ['/projects/0/ctdas/awoude/NRT/ICOS_OUTPUT/fire.202206.nc']

# Show progress
#print('Files to convert:')
#for files in fluxstring:
#    print(str(files))

# Define CTE-HR grid variables
lon_bounds = [-15.,35.] # in degrees 
lat_bounds = [33.,72.] # in degrees 
lon_step = 0.2 # in degrees 
lat_step = 0.1 # in degrees 
nx = 250 # gridsize in x dimension (amount of gridcells)
ny = 390 # gridsize in y dimension (amount of gridcells)
interpolation_method = 0 # Default = 1; 0 = nearest neighbour, 1 = bilinear interpolation

# Loop over flux files
for files in fluxstring:
    # Extract old filename
    filename = os.path.basename(files)
    
    # Create new filenames and status report
    newstring = '/projects/0/ctdas/dkivits/DATA/CTE_regriddedfluxes/' + filename
    print('Busy with ... ' + newstring)
    
    # Prepare NC file
    CTEHR_nc = nc.Dataset(newstring,'w',format='NETCDF4')
    CTEHR_nc.createDimension('Time', None)
    CTEHR_nc.createDimension('emissions_zdim',9) 
    CTEHR_nc.createDimension('south_north', 390)
    CTEHR_nc.createDimension('west_east', 250)
    
    # Read CTE-HR flux file
    flux_nc= nc.Dataset(files,'r')
    
    if 'fire' in filename:
        # Create input field in NC file
        fossilvar = CTEHR_nc.createVariable('E_CO2_003','f4',('Time','emissions_zdim','south_north','west_east'))
        fossilvar.FieldType=int(104)
        fossilvar.MemoryOrder="XYZ"
        fossilvar.description="CO2 emission"
        fossilvar.units="mole/km2/hr"
        fossilvar.stagger="Z"
        fossilvar.coordinates="XTIME XZDIM XLAT XLONG"    
        
        # Select correct variable and transform fluxes
        for time in range(0,len(flux_nc['fire'][:,0,0])):
            flux_data = flux_nc.variables['fire'][time,:,:] 
            flux_trans = m.transform_scalar(flux_data,lons=np.arange(lon_bounds[0],lon_bounds[1],lon_step),lats=np.arange(lat_bounds[0],lat_bounds[1],lat_step),nx=nx,ny=ny,order=interpolation_method)     
         
            # Make sure all NaN values are converted to zeroes (required by Snellius)     
            flux_trans = np.where(np.isnan(flux_trans), 0, flux_trans)
            
            # Apply a Landmask over the converted flux data:
            #flux_trans = np.where(geo.variables['LANDMASK'][1,:,:] != 1, flux_trans, 0)
            
            # Put the transformed CTE-HR fluxes into the NC file
            #  fossilvar[:,:,:,:] = 0
            fossilvar[time,1,:,:] = flux_trans
            
        # Check if dimensions of transformed    
        print('Dimensions of newly created nc file: ' + str(fossilvar.shape))    
    
    if 'nep' in filename:
        # Create input field in NC file
        fossilvar = CTEHR_nc.createVariable('E_CO2_001','f4',('Time','emissions_zdim','south_north','west_east'))
        fossilvar.FieldType=int(104)
        fossilvar.MemoryOrder="XYZ"
        fossilvar.description="CO2 emission"
        fossilvar.units="mole/km2/hr"
        fossilvar.stagger="Z"
        fossilvar.coordinates="XTIME XZDIM XLAT XLONG"    
        
        # Select correct variable and transform fluxes
        for time in range(0,len(flux_nc['nep'][:,0,0])):
            flux_data = flux_nc.variables['nep'][time,:,:]
            flux_trans = m.transform_scalar(flux_data,lons=np.arange(lon_bounds[0],lon_bounds[1],lon_step),lats=np.arange(lat_bounds[0],lat_bounds[1],lat_step),nx=nx,ny=ny,order=interpolation_method)     
            
            # Make sure all NaN values are converted to zeroes (required by Snellius)     
            flux_trans = np.where(np.isnan(flux_trans), 0, flux_trans)
                     
            # Apply a Landmask over the converted flux data:
            #flux_trans = np.where(geo.variables['LANDMASK'][1,:,:] != 1, flux_trans, 0)

            # Put the transformed CTE-HR fluxes into the NC file
            #  fossilvar[:,:,:,:] = 0
            fossilvar[time,1,:,:] = flux_trans
            
        # Check if dimensions of transformed   
        print('Dimensions of newly created nc file: ' + (str(fossilvar.shape)))    
           
    if 'ocean' in filename:
        # Create input field in NC file
        fossilvar = CTEHR_nc.createVariable('E_CO2_002','f4',('Time','emissions_zdim','south_north','west_east'))
        fossilvar.FieldType=int(104)
        fossilvar.MemoryOrder="XYZ"
        fossilvar.description="CO2 emission"
        fossilvar.units="mole/km2/hr"
        fossilvar.stagger="Z"
        fossilvar.coordinates="XTIME XZDIM XLAT XLONG"    
        
        # Select correct variable and transform fluxes
        for time in range(0,len(flux_nc['ocean'][:,0,0])):
            flux_data = flux_nc.variables['ocean'][time,:,:] 
            flux_trans = m.transform_scalar(flux_data,lons=np.arange(lon_bounds[0],lon_bounds[1],lon_step),lats=np.arange(lat_bounds[0],lat_bounds[1],lat_step),nx=nx,ny=ny,order=interpolation_method)     
             
            # Make sure all NaN values are converted to zeroes (required by Snellius)     
            flux_trans = np.where(np.isnan(flux_trans), 0, flux_trans)
            
            # Apply a Landmask over the converted flux data:
            #flux_trans = np.where(geo.variables['LANDMASK'][1,:,:] != 0, flux_trans, 0)
            
            # Put the transformed CTE-HR fluxes into the NC file
            #  fossilvar[:,:,:,:] = 0
            fossilvar[time,1,:,:] = flux_trans
            
        # Check if dimensions of transformed
        print('Dimensions of newly created nc file: ' + str(fossilvar.shape))    
    
    if 'anthropogenic' in filename:
        # Create input field in NC file
        fossilvar = CTEHR_nc.createVariable('E_CO2_004','f4',('Time','emissions_zdim','south_north','west_east'))
        fossilvar.FieldType=int(104)
        fossilvar.MemoryOrder="XYZ"
        fossilvar.description="CO2 emission"
        fossilvar.units="mole/km2/hr"
        fossilvar.stagger="Z"
        fossilvar.coordinates="XTIME XZDIM XLAT XLONG"    
        
        # Select correct variable and transform fluxes
        for time in range(0,len(flux_nc['combustion'][:,0,0])):    
            flux_data = flux_nc.variables['combustion'][time,:,:] 
            flux_trans = m.transform_scalar(flux_data,lons=np.arange(lon_bounds[0],lon_bounds[1],lon_step),lats=np.arange(lat_bounds[0],lat_bounds[1],lat_step),nx=nx,ny=ny,order=interpolation_method)     
             
            # Make sure all NaN values are converted to zeroes (required by Snellius)     
            flux_trans = np.where(np.isnan(flux_trans), 0, flux_trans)
            
            # Apply a Landmask over the converted flux data:
            #flux_trans = np.where(geo.variables['LANDMASK'][1,:,:] != 1, flux_trans, 0)
            
            # Put the transformed CTE-HR fluxes into the NC file
            #  fossilvar[:,:,:,:] = 0
            fossilvar[time,1,:,:] = flux_trans
            
        # Check if dimensions of transformed
        print('Dimensions of newly created nc file: ' + str(fossilvar.shape))    

    CTEHR_nc.close()

print('Done changing the emission fluxes!\n Files have been written to: ' + outdir + ' !')
