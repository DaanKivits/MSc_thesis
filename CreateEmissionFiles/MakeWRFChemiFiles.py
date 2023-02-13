import netCDF4 as nc
from netCDF4 import stringtochar
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
fluxdir = '/projects/0/ctdas/dkivits/CTDAS_WRF_runs/testrun/input/test_fwd_highres/emissions/'
mode = 'fwd'
outdir = fluxdir + mode + '/'
#backupdir = fluxdir + mode + '_backup/'

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
geopath = '/projects/0/ctdas/dkivits/CTDAS_WRF_runs/testrun/WPS/WPS_test_fwd_highres/geo_em.d01.nc'
geo = nc.Dataset(geopath,'r')
dx_dom = 10000
we = 345
sn = 436
truelat1 = 72
truelat2 = 33
ref_lat = 52.5
ref_lon = 10
domain = 1

# Define basemap to prepare for transformation later
m = Basemap(width=dx_dom*we,height=dx_dom*sn,
            rsphere=(6378137.00,6356752.3142),\
            resolution='l',area_thresh=1000.,projection='lcc',\
            lat_1=truelat1,lat_2=truelat2,lat_0=ref_lat,lon_0=ref_lon)

#m = Basemap(llcrnrlon=-15, llcrnrlat=33, urcrnrlon=35, urcrnrlat=73,
#            rsphere=(6378137.00, 6356752.3142),
#            resolution='l', area_thresh=1000., projection='cyl')
            #lat_0=lat, lon_0=lon)

# Define names of tracers and flux fields
variablelist = ['E_CO2_001','E_CO2_002','E_CO2_003','E_CO2_004']
fluxfilelist = ['nep','fire','ocean','anthropogenic']
fluxnamelist=['nep','fire','ocean','combustion']
unused_variablelist= ['E_CO2', 'E_CO2_proc_1', 'E_CO2_proc_2','E_CO2_proc_3','e_CO2_000']

# Define which files to loop over
fluxstring = []
for name in fluxfilelist:
    fluxstring += sorted(glob.glob(dirname + name + '.201807.nc')) + sorted(glob.glob(dirname + name + '.201806.nc')) + sorted(glob.glob(dirname + name + '.201808.nc'))

# Define CTE-HR grid variables
#lon_bounds = [-15.,35.] # in degrees 
#lat_bounds = [33.,72.] # in degrees 
#lon_step = 0.2 # in degrees 
#lat_step = 0.1 # in degrees 
#nx = 250 # gridsize in x dimension (amount of gridcells)
#ny = 390 # gridsize in y dimension (amount of gridcells)
lon_bounds = [-15.,35.] # in degrees
lat_bounds = [33.,72.] # in degrees
lon_step = 0.2 # in degrees
lat_step = 0.1 # in degrees
nx = 344 # gridsize in x dimension (amount of gridcells)
ny = 435 # gridsize in y dimension (amount of gridcells)

interpolation_method = 0 # Default = 1; 0 = nearest neighbour, 1 = bilinear interpolation
injection_level = 1 # The vertical level at which the emissions fluxes are released
basedate = datetime(2000,1,1,0,0,0,0)

for file in fluxstring:
    # Progress report
    print('WORKING ON ... ' + file)

    # Open CTE-HR flux file, extract time to loop over at later stage
    CTEHR_nc =  nc.Dataset(file,'r')
    timevar = CTEHR_nc.variables['time'][:]
    time_array = np.arange(0, len(timevar),1)

    # Debug: print CTEHR_nc contents
    #print(list(CTEHR_nc.variables))

    # Loop over time variable to extract each hour in the CTE-HR files
    for times in time_array:
        # Loop over CTE_HR time variables and create emission flux files according to these times
        curtime= timevar[times]
        #time = basedate + timedelta(seconds=int(curtime)) - relativedelta(years=7)
        time = basedate + timedelta(seconds=int(curtime))
        print('Busy with ... ' + str(time))
        timestr = str(time.date()) + '_' + str(time.time())
        new_filename = outdir + 'wrfchemi_d0' + str(domain) + '_' + timestr

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
            wrfchemi_nc.createDimension('DateStrLen', 19)
            wrfchemi_nc.TITLE = 'OUTPUT FROM CTE-HR FLUX PREPROCESSING_FOR V4.1'
            wrfchemi_nc.GRIDTYPE = geo.GRIDTYPE
            wrfchemi_nc.DYN_OPT = geo.DYN_OPT
            wrfchemi_nc.CEN_LAT = geo.CEN_LAT
            wrfchemi_nc.CEN_LON = geo.CEN_LON
            wrfchemi_nc.TRUELAT1 = geo.TRUELAT1
            wrfchemi_nc.TRUELAT2 = geo.TRUELAT2
            wrfchemi_nc.MOAD_CEN_LAT = geo.MOAD_CEN_LAT
            wrfchemi_nc.STAND_LON = geo.STAND_LON
            wrfchemi_nc.POLE_LAT = geo.POLE_LAT
            wrfchemi_nc.POLE_LON = geo.POLE_LON
            wrfchemi_nc.corner_lats = geo.corner_lats
            wrfchemi_nc.corner_lons = geo.corner_lons
            wrfchemi_nc.MAP_PROJ = geo.MAP_PROJ
            wrfchemi_nc.DX = geo.DX
            wrfchemi_nc.DY = geo.DY
            wrfchemi_nc.MMINLU = geo.MMINLU
            wrfchemi_nc.NUM_LAND_CAT = geo.NUM_LAND_CAT
            wrfchemi_nc.ISWATER = geo.ISWATER
            wrfchemi_nc.ISLAKE = geo.ISLAKE
            wrfchemi_nc.ISICE = geo.ISICE
            wrfchemi_nc.ISURBAN = geo.ISURBAN
            wrfchemi_nc.ISOILWATER = geo.ISOILWATER
            wrfchemi_nc.grid_id = geo.grid_id
            wrfchemi_nc.parent_id = geo.parent_id
            wrfchemi_nc.i_parent_start = geo.i_parent_start
            wrfchemi_nc.j_parent_start = geo.j_parent_start
            wrfchemi_nc.i_parent_end = geo.i_parent_end
            wrfchemi_nc.j_parent_end = geo.j_parent_end
            wrfchemi_nc.parent_grid_ratio = geo.parent_grid_ratio
            wrfchemi_nc.sr_x = geo.sr_x
            wrfchemi_nc.sr_y = geo.sr_y
            wrfchemi_nc.GMT = 1
            wrfchemi_nc.JULYR = time.year
            wrfchemi_nc.JULDAY = int(time.strftime('%j'))
            wrfchemi_nc.history = 'made manually by converting CTE-HR fluxes to the WRF grid'

            # Create the timestring variable
            timestrvar = wrfchemi_nc.createVariable('Times','S1',('Time', 'DateStrLen'))
            timestr = np.array([timestr],dtype='S19')
            timestrvar[:] = stringtochar(timestr)
            #timestrvar.set_auto_chartostring(True)
            timestrvar._Encoding = 'ascii'
            timestrvar[:] = timestr
            del timestrvar._Encoding

            # Create the variables that are left empty, to be filled in with CAMS boundary + initial conditions
            for variable in unused_variablelist:
                unusedvariable = wrfchemi_nc.createVariable(variable,'f4',('Time','emissions_zdim','south_north','west_east'))
                unusedvariable.FieldType=int(104)
                unusedvariable.MemoryOrder="XYZ"
                unusedvariable.description="CO2 emission"
                unusedvariable.units="mole/km2/hr"
                unusedvariable.stagger="Z"
                unusedvariable.coordinates="XLONG XLAT"
                unusedvariable[0,:,:,:] = 0.0

            # Create the variable that will be filled with CTE-HR fluxes
            tracervariable = wrfchemi_nc.createVariable(variablename,'f4',('Time','emissions_zdim','south_north','west_east'))
            tracervariable.FieldType=int(104)
            tracervariable.MemoryOrder="XYZ"
            tracervariable.description="CO2 emission"
            tracervariable.units="mole/km2/hr"
            tracervariable.stagger="Z"
            tracervariable.coordinates="XLONG XLAT"

            # Snellius does not like zeroes, so fill the unused variables with zeroes too
            tracervariable[0,:,:,:] = 0.0

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
            tracervariable[0,injection_level,:,:] = flux_trans * 3600 * 1e6
            
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
            tracervariable.coordinates="XLONG XLAT"
            
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
            tracervariable[0,injection_level,:,:] = flux_trans * 3600 * 1e6
 
            # Close the file
            wrfchemi_nc.close()

    CTEHR_nc.close()
