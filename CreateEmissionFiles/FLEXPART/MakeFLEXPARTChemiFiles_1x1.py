# Daan Kivits, 2023

# This script creates FLEXPART-compatible flux files with a 1 x 1 degree spatial resolution, and fills these files with CTE-HR
# high-resolution fluxes for the biosphere (NEP and wildfires), ocean, and anthropogenic sectors.
# Some points of attention:
#   - The FLEXPART grid consists of a simulation point in the middle of every gridcell, which is why the lon_bounds and lat_bounds 
#       ranges start at a seemingly odd value.
#   - This script creates monthly FLEXPART flux files and later groups all flux files together into a yearly FLEXPART flux file that
#       that can be used by the FLEXPART processing algorithm.

##############################################
########## LOAD NECCESSARY PACKAGES ##########
##############################################
from dateutil.relativedelta import relativedelta
import shutil
from datetime import datetime, timedelta, time, date
import os
import xarray as xr
import glob
import time
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import netCDF4 as nc
from netCDF4 import stringtochar
import numpy as np

########################################
########## CREATE DIRECTORIES ##########
########################################
# Check if a regridded flux directory exists, otherwise create it
dirname = '/projects/0/ctdas/awoude/NRT/ICOS_OUTPUT/'
outdir = '/projects/0/ctdas/dkivits/DATA/FLEXPART/code/fluxes/CTEHR/upscaled/'
# dirname = '/mnt/c/Users/daank/OneDrive - WageningenUR/WUR/MSc thesis/Data/CTEHR/'
# outdir = '/home/dkivits/Footprints/Flexpart/DICE/fluxes/CTE_HR/'

if not os.path.exists(outdir):
    os.mkdir(outdir)

if os.path.exists(outdir):
    shutil.rmtree(outdir)
    os.mkdir(outdir)

######################################
########## CREATE VARIABLES ##########
######################################
# Define which files to loop over
fluxstring = []
for name in fluxfilelist:
    fluxstring += sorted(glob.glob(dirname + name + '.2018*.nc'))

timelist = np.arange(datetime(int(fluxstring[0][-9:-5]),int(fluxstring[0][-5:-3]),1),datetime(int(fluxstring[-1][-9:-5]),int(fluxstring[-1][-5:-3]),1) + relativedelta(months=+1),timedelta(hours=1)).astype(datetime)

# Define names of tracers and flux fields
variablelist = ['biosphere', 'fire', 'ocean', 'fossil']
fluxfilelist = ['nep', 'fire', 'ocean', 'anthropogenic']
fluxnamelist = ['nep', 'fire', 'ocean', 'combustion']

# Define CTE-HR grid variables
lon_bounds = [-14.9, 35.]  # in degrees
lat_bounds = [33.05, 72.]  # in degrees
flexbounds_lon = [-14.875, 35.]  # in degrees
flexbounds_lat = [33.125, 73.]  # in degrees
lon_step = 0.2  # in degrees
lat_step = 0.1  # in degrees
flex_lon_step = 1  # in degrees
flex_lat_step = 1  # in degrees
flexlat = 39  # gridsize in y dimension (amount of gridcells)
flexlon = 50  # gridsize in x dimension (amount of gridcells)
lat_list = np.arange(lat_bounds[0], lat_bounds[1], lat_step)
lon_list = np.arange(lon_bounds[0], lon_bounds[1], lon_step)
flexlat_list = np.arange(flexbounds_lat[0], flexbounds_lat[1], flex_lat_step)
flexlon_list = np.arange(flexbounds_lon[0], flexbounds_lon[1], flex_lon_step)
interpolation_method = 0 # Default = 1; 0 = nearest neighbour, 1 = bilinear interpolation
basedate = datetime(2000, 1, 1, 0, 0, 0, 0)

# Define basemap to prepare for transformation later
m = Basemap(llcrnrlon=-15, llcrnrlat=33, urcrnrlon=35, urcrnrlat=73,
            rsphere=(6378137.00, 6356752.3142),
            resolution='l', area_thresh=1000., projection='cyl')
            #lat_0=lat, lon_0=lon)

for file in fluxstring:
    # Progress report
    print('WORKING ON ... ' + file)

    # Open CTE-HR flux file, extract time to loop over at later stage
    CTEHR_nc = nc.Dataset(file, 'r')
    timevar = CTEHR_nc.variables['time'][:]
    curtime = basedate + timedelta(seconds=int(timevar[0]))

    flexpart_timestring = 'hours since ' + curtime.strftime('%Y-%m-%d')
    new_filename = outdir + 'CTEHR_FLEXPART_' + \
        str(curtime.month) + '_' + str(curtime.year) + '.nc'

    # Loop over time variable to extract each hour in the CTE-HR files
    for timeindex,time in enumerate(timevar):
        CTE_time = basedate + timedelta(seconds=int(time))
        timedif = int((CTE_time - timelist[0]).total_seconds() / 3600)

        # Check if file exists, and either append or write depending on that condition
        if not os.path.exists(new_filename):
            wrfchemi_nc = nc.Dataset(
                new_filename, 'w', format='NETCDF4')
            fluxname = list(CTEHR_nc.variables.keys())[3]
            variableindex = fluxnamelist.index(fluxname)
            variablename = variablelist[variableindex]
            variable_data = CTEHR_nc.variables[fluxname][timeindex, :, :]
            flux_trans = m.transform_scalar(variable_data, lons=lon_list, lats=lat_list, nx=flexlon, ny=flexlat, order=interpolation_method)

            # Make sure the transformed fluxes contain only zeroes, and all  NaN values are converted to zeroes (required by Snellius)
            flux_trans = np.where(np.isnan(flux_trans), 0, flux_trans)

            # Fill in newly created emission flux NC files
            co2 = wrfchemi_nc.createGroup('co2')
            co2.createDimension('time', len(timevar))
            co2.createDimension('lat', flexlat)
            co2.createDimension('lon', flexlon)
            co2.tracer = "co2"
            co2.units = "micromol/m**2/s"
            co2.timestep = "1h"
            co2.categories = 'biosphere fire ocean fossil'.split()

            # Create the variable that will be filled with CTE-HR fluxes
            timevariable = co2.createVariable(
                'time', 'i4', ('time'))
            timevariable.units = flexpart_timestring
            #timevariable.units = 'hours since 2018-01-01'
            timevariable.calendar = "proleptic_gregorian"
            timevariable[timeindex] = timedif

            latvariable = co2.createVariable(
                'lat', 'f4', ('lat'), fill_value=np.NaN)
            latvariable.units = "degrees"
            latvariable.info = "center of the gridcells"
            for lat in range(0, len(flexlat_list)):
                latvariable[lat] = flexlat_list[lat]

            lonvariable = co2.createVariable(
                'lon', 'f4', ('lon'), fill_value=np.NaN)
            lonvariable.units = "degrees"
            lonvariable.info = "center of the gridcells"
            for lon in range(0, len(flexlon_list)):
                lonvariable[lon] = flexlon_list[lon]

            for variables in variablelist:
                tracervariable = co2.createVariable(
                    variables, 'f4', ('time', 'lat', 'lon'), fill_value=np.NaN)
                tracervariable.units = "micromoles/m2/s"

                # Fill this variable with zeroes first
                tracervariable[:, :, :] = 0.0

            # Inject the transformed fluxes at the injection level
            co2.variables[variablename][timeindex, :, :] = flux_trans * 1e6

            # Close the file
            wrfchemi_nc.close()

        else:
            # Append the data
            wrfchemi_nc = nc.Dataset(new_filename, 'a', format='NETCDF4')
            fluxname = list(CTEHR_nc.variables.keys())[3]
            variableindex = fluxnamelist.index(fluxname)
            variablename = variablelist[variableindex]
            variable_data = CTEHR_nc.variables[fluxname][timeindex, :, :]
            flux_trans = m.transform_scalar(variable_data, lons=lon_list, lats=lat_list, nx=flexlon, ny=flexlat, order=interpolation_method)

            # Make sure the transformed fluxes contain only zeroes, and all  NaN values are converted to zeroes (required by Snellius)
            flux_trans = np.where(np.isnan(flux_trans), 0, flux_trans)

            # Inject the transformed fluxes at the injection level
            co2.variables[variablename][timeindex, :, :] = flux_trans * 1e6
            timevariable[timeindex] = timedif
            
            wrfchemi_nc.close()
    CTEHR_nc.close()

# Combine into a singular yearly file
yearly_filename = outdir + 'CTEHR_FLEXPART_' +  str(curtime.year) + '.nc'
DS = xr.open_mfdataset(os.path.join(outdir, '*.nc'), group = 'co2')
DS.to_netcdf(yearly_filename,group = 'co2')
