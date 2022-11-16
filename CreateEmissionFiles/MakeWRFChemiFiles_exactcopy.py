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

# Check if a regridded flux directory exists, otherwise create it
dirname = '/projects/0/ctdas/awoude/NRT/ICOS_OUTPUT/'
fluxdir = '/projects/0/ctdas/dkivits/CTDAS_WRF_runs/testrun/input/test_fwd_nrt/emissions/'
mode = 'fwd'
outdir = fluxdir + mode + '/'
backupdir = fluxdir + mode + '_backup/'

# Check if a backup flux field directory exists, otherwise create it
#if not os.path.exists(backupdir):
#   shutil.copytree(outdir, backupdir)

# To reset the files to the backup standard, copy backup folder to outdir always
if os.path.exists(outdir):
    shutil.rmtree(outdir)
    os.mkdir(outdir)

#shutil.copytree(backupdir, outdir, dirs_exist_ok=True)

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
domain = 1

# Define basemap to prepare for transformation later
m = Basemap(width=dx_dom*we,height=dx_dom*sn,
            rsphere=(6378137.00,6356752.3142),\
            resolution='l',area_thresh=1000.,projection='lcc',\
            lat_1=45.,lat_2=55,lat_0=lat,lon_0=lon)

# Define names of tracers and flux fields
variablelist = ['E_CO2_001','E_CO2_002','E_CO2_003','E_CO2_004']
fluxfilelist = ['nep','fire','ocean','anthropogenic']
fluxnamelist=['nep','fire','ocean','combustion']
unused_variablelist= ['E_CO2_000', 'E_CO2_proc_1', 'E_CO2_proc_2','E_CO2_proc_3','e_CO2_000']

# Define which files to loop over
fluxstring = []
for name in fluxfilelist:
    fluxstring += sorted(glob.glob(dirname + name + '.202206.nc'))

# Define CTE-HR grid variables
lon_bounds = [-15.,35.] # in degrees 
lat_bounds = [33.,72.] # in degrees 
lon_step = 0.2 # in degrees 
lat_step = 0.1 # in degrees 
nx = 15 # gridsize in x dimension (amount of gridcells)
ny = 12 # gridsize in y dimension (amount of gridcells)
interpolation_method = 0 # Default = 1; 0 = nearest neighbour, 1 = bilinear interpolation
injection_level = 1 # The vertical level at which the emissions fluxes are released
basedate = datetime(2000,1,1,0,0,0,0)

for file in fluxstring:
    
    # Progress report
    print('WORKING ON ... ' + file)

    # Open CTE-HR flux file, extract time to loop over at later stage
    CTEHR_nc =  nc.Dataset(file,'r')
    timevar = CTEHR_nc.variables['time'][:]

    # Debug: print CTEHR_nc contents
    #print(list(CTEHR_nc.variables))

    # Loop over time variable to extract each hour in the CTE-HR files
    for times in range(0, len(timevar)):
        # Loop over CTE_HR time variables and create emission flux files according to these times
        curtime= timevar[times]
        time = basedate + timedelta(seconds=int(curtime)) - relativedelta(years=7)
        print('Busy with ... ' + str(time))
        new_filename = outdir + 'wrfchemi_d0' + str(domain) + '_' + str(time.date()) + '_' + str(time.time())
        
        # Check if file exists, and either append or write depending on that condition
        if not os.path.exists(new_filename):
            wrfchemi_nc = nc.Dataset(new_filename, 'w', format='NETCDF3_CLASSIC')
            fluxname = list(CTEHR_nc.variables.keys())[3]
            variableindex = fluxnamelist.index(fluxname)
            variablename = variablelist[variableindex]
            variable_data = CTEHR_nc.variables[fluxname][times,:,:]
            flux_trans = m.transform_scalar(variable_data,lons=np.arange(lon_bounds[0],lon_bounds[1],lon_step),lats=np.arange(lat_bounds[0],lat_bounds[1],lat_step),nx=nx,ny=ny,order=interpolation_method)
            
           # Fill in newly created emission flux NC files
            wrfchemi_nc.createDimension('Time', None)
            wrfchemi_nc.createDimension('emissions_zdim',9)
            wrfchemi_nc.createDimension('south_north', ny)
            wrfchemi_nc.createDimension('west_east', nx)
            wrfchemi_nc.TITLE = 'OUTPUT FROM PREPROCESSING_FOR V4.1'
            wrfchemi_nc.CEN_LAT = 45.0
            wrfchemi_nc.CEN_LON = 6.0
            wrfchemi_nc.TRUELAT1 = 45.0
            wrfchemi_nc.TRUELAT2 = 45.0
            wrfchemi_nc.MOAD_CEN_LAT = 45.0
            wrfchemi_nc.STAND_LON = 6.0
            wrfchemi_nc.POLE_LAT = 90.0
            wrfchemi_nc.POLE_LON = 0.0
            wrfchemi_nc.MAP_PROJ = 1
            wrfchemi_nc.DX = 100000.0
            wrfchemi_nc.DY = 100000.0
            wrfchemi_nc.MMINLU = 'MODIFIED_IGBP_MODIS_NOAH'
            wrfchemi_nc.NUM_LAND_CAT = 21
            wrfchemi_nc.ISWATER = 17
            wrfchemi_nc.ISLAKE = 21
            wrfchemi_nc.ISICE = 15
            wrfchemi_nc.ISURBAN = 13
            wrfchemi_nc.ISOILWATER = 14
            wrfchemi_nc.GMT = 14
            wrfchemi_nc.JULYR = 2015
            wrfchemi_nc.JULDAY = 152
            wrfchemi_nc.history = 'made manually by converting CTE-HR fluxes to the WRF grid'
            wrfchemi_nc.NCO = 'netCDF Operators version 5.0.0 (Homepage = http://nco.sf.net, Code = http://github.com/nco/nco)'

            # Create the variables that are left empty, to be filled in with CAMS boundary + initial conditions
            for variable in unused_variablelist:
                unusedvariable = wrfchemi_nc.createVariable(variable,'f4',('Time','emissions_zdim','south_north','west_east'))
                unusedvariable.FieldType=int(104)
                unusedvariable.MemoryOrder="XYZ"
                unusedvariable.description="CO2 emission"
                unusedvariable.units="mole/km2/hr"
                unusedvariable.stagger="Z"
                unusedvariable.coordinates="XLAT XLONG"

            # Create the variable that will be filled with CTE-HR fluxes
            tracervariable = wrfchemi_nc.createVariable(variablename,'f4',('Time','emissions_zdim','south_north','west_east'))
            tracervariable.FieldType=int(104)
            tracervariable.MemoryOrder="XYZ"
            tracervariable.description="CO2 emission"
            tracervariable.units="mole/km2/hr"
            tracervariable.stagger="Z"
            tracervariable.coordinates="XLAT XLONG"
            
            # Make sure the transformed fluxes contain only zeroes, and all  NaN values are converted to zeroes (required by Snellius)
            flux_trans = np.where(np.isnan(flux_trans), 0, flux_trans)

            # Apply a Landmask over the converted flux data:
            if variablename == 'ocean':
                flux_trans = np.where(geo.variables['LANDMASK'][0,:,:] != 1, flux_trans, 0)
            else:
                flux_trans = np.where(geo.variables['LANDMASK'][0,:,:] != 0, flux_trans, 0)

            # Fill this variable with zeroes first
            tracervariable[0,:,:,:] = 0.0
            
            # Inject the transformed fluxes at the injection level
            tracervariable[0,injection_level,:,:] = flux_trans * 3600
            
            # Close the file
            wrfchemi_nc.close() 
        
        else: 
            # Append the data
            wrfchemi_nc = nc.Dataset(new_filename, 'a')
            fluxname = list(CTEHR_nc.variables.keys())[3]
            variableindex = fluxnamelist.index(fluxname)
            variablename = variablelist[variableindex]
            variable_data = CTEHR_nc.variables[fluxname][times,:,:]
            flux_trans = m.transform_scalar(variable_data,lons=np.arange(lon_bounds[0],lon_bounds[1],lon_step),lats=np.arange(lat_bounds[0],lat_bounds[1],lat_step),nx=nx,ny=ny,order=interpolation_method)
            
           # Create the variable that will be filled with CTE-HR fluxes
            tracervariable = wrfchemi_nc.createVariable(variablename,'f4',('Time','emissions_zdim','south_north','west_east'))
            tracervariable.FieldType=int(104)
            tracervariable.MemoryOrder="XYZ"
            tracervariable.description="CO2 emission"
            tracervariable.units="mole/km2/hr"
            tracervariable.stagger="Z"
            tracervariable.coordinates="XLAT XLONG"

             # Make sure the transformed fluxes contain only zeroes, and all  NaN values are converted to zeroes (required by Snellius)
            flux_trans = np.where(np.isnan(flux_trans), 0, flux_trans)

            # Apply a Landmask over the converted flux data:
            if variablename == 'ocean':
                flux_trans = np.where(geo.variables['LANDMASK'][0,:,:] != 1, flux_trans, 0)
            else:
                flux_trans = np.where(geo.variables['LANDMASK'][0,:,:] != 0, flux_trans, 0)

            # Fill this variable with zeroes first
            tracervariable[0,:,:,:] = 0.0

            # Inject the transformed fluxes at the injection level
            tracervariable[0,injection_level,:,:] = flux_trans * 3600
 
            # Close the file
            wrfchemi_nc.close()

    CTEHR_nc.close()

