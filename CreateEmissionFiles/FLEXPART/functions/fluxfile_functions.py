# Daan Kivits, 2023

# This file contains functions used by the FLEXPART flux file creation scripts.

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
from scipy import stats

def median2d(arr, new_shape):
    """ Function to average any given shape, which is a multiple of the domain size, to the domain size
    Input:
        arr: np.ndarray: original array to be averaged
        new_shape: tuple: shape to be averaged to
    Returns:
        np.ndarray: averaged arr"""
    shape = (new_shape[0], arr.shape[-2] // new_shape[-2],
             new_shape[1], arr.shape[-1] // new_shape[-1])
    if len(arr.shape) == 3: # If time is included:
        shape = (len(arr),) + shape
        
    a = stats.mode(arr.reshape(shape), axis=1)[0]
    b = stats.mode(a, axis=3)[0]
    return b.squeeze()
def get_lu(flux_array):
    """ Function to extract the landuse given any given shape. The shape should be a multiplication of 
    0.05 x 0.05 degrees, so a shape with a 0.1 x 0.2 gridcell size would be possible, but a 0.0825 x 0.125 wouldn't be.   
    Returns:
        returns the landuse array from the landuse dataset  of any given shape overlapping with the 
        extent of this landuse dataset """
    with nc.Dataset('/projects/0/ctdas/NRT/data/SiB/CORINE_PFT_EUROPA_NRT.nc') as ds:
        lu = ds['landuse'][:]
        lu = np.flipud(lu)
    lu = median2d(lu, flux_array.shape[1:])  
    return lu

def exchange_fluxes_based_on_landuse(src_file, trg_file, lu_types):
    """ Function to extract fluxes of two different years, and fill in fluxes of one year 
    to the flux field of the other year for all pixels of a certain landuse type. Plot the difference between
    the input and output flux fields and show the image for each of the variables in variablelist. """
    src = nc.Dataset(src_file, mode = 'r')
    
    ## Create a time variable to loop over
    timevar = src.variables['time'][:]
    fluxname = list(src.variables.keys())[3]

    for timeindex,time in enumerate(timevar):
        trg = nc.Dataset(trg_file, mode='r+')
                
        # Loop over variables of interest and put in lu_specific fluxsets of different year
        if name in variablelist:
            var = trg.variables[fluxname][timeindex,:,:]
            luvar = src.variables[fluxname][timeindex,:,:]
            var[lu == lu_types] = luvar[lu == lu_types]
            trg.variables[fluxname][timeindex,:,:] = var
            dif = trg.variables[fluxname][timeindex,:,:] - src.variables[fluxname][timeindex,:,:]
            
            if timeindex == 2:
                plt.imshow(np.flipud(dif))
                plt.show()

        # Save the file
        trg.close()
    src.close()

def check_mask(src_file, lu_type):
    """ Function to check what area is affected by the landuse mask, by simply putting in a large
    emission flux value of 999 and showing the result. """
    src = nc.Dataset(src_file)
    
    ## Create a time variable to loop over
    timevar = np.arange(0, len(src.variables['time'][:]),1)
    time_array = np.arange(0, len(timevar),1)
    
    for name, var in src.variables.items():
        if name in variablelist:
            for time in time_array:
                data = src.variables[name][time,:,:]
                data[lu == lu_type] = 999
                plt.imshow(np.flipud(data))
                plt.show()

def create_fluxfile_from_source(src_file, trg_file):
    """ Function to copy the variables and (global and variable) attributes of an existing emission flux file
    into a new emission flux file. """
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