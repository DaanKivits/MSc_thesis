# Daan Kivits, 2023

# A simple plotting script used to plot the average CTE-HR flux anomalies for July, August, and September (JAS) of different 
# climate reference periods compared to 2018.

##############################################
########## LOAD NECCESSARY PACKAGES ##########
##############################################
import netCDF4 as nc
import numpy as np
import glob 
import os
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, maskoceans
import argparse
from matplotlib import colors
import numpy.ma as ma
from matplotlib.patches import Path, PathPatch
import matplotlib.ticker
import sys
sys.path.insert(0, '/projects/0/ctdas/dkivits/scripts/Plotting/functions')
import plotting_functions

########################################
##### PARSE COMMAND LINE ARGUMENTS #####
########################################
parser = argparse.ArgumentParser()
parser.add_argument(dest="fluxpath")
parser.add_argument(dest="months")
args = parser.parse_args()
            
months = args.months.split(sep=',')
months = np.arange(int(months[0]),int(months[1])+1,1)

fluxpath = args.fluxpath

# Define which files to loop over
fluxstring = []
fluxstring += sorted(sorted(glob.glob(fluxpath + '*.nc')), key = len)
del(fluxstring[-2])
print(fluxstring)

fluxfile_basename = os.path.basename(fluxstring[0]).split(sep='.')[2]

fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(16,9))
plt.subplots_adjust(left = 0.02, top = 0.95, bottom = 0.06, right=0.97, hspace = 0.12, wspace = -0.5)
axs = ax.ravel()

for i in range(0, len(fluxstring)):
#for i in range(0, len(months)):
    # Define fluxes
    fluxdata = nc.Dataset(fluxstring[i],'r') 
    fluxes_cte = fluxdata.variables['nep'][0,:,:]*1e6
    
    fluxfile_yearname = os.path.basename(fluxstring[i]).split(sep='.')[3].split(sep='_')[1]
    
    axs[i].set_title(fluxfile_yearname, size=12, weight = 'bold')
    
    # Define basemap to prepare for transformation later
    m = Basemap(llcrnrlon=-15, llcrnrlat=33, urcrnrlon=35, urcrnrlat=72,
            rsphere=(6378137.00, 6356752.3142),
            resolution='l', area_thresh=1000., projection='cyl', ax = axs[i])
            #lat_0=lat, lon_0=lon)
    
    cbarticks = np.arange(-3,4,1)
    img = m.imshow(fluxes_cte,cmap='RdBu_r',zorder=1, vmin = cbarticks[0], vmax = cbarticks[-1])
    m.drawcoastlines()
    
    if i == len(fluxstring)-1:
        cbaxes = fig.add_axes([0.19, 0, 0.61, 0.03]) #Add position (left, bottom, width, height)
        #cb = m.colorbar(mappable=img, extend='both', extendrect=True, location="right", ticks = cbarticks)
        cb = fig.colorbar(mappable=img, extend='both', extendrect=True, ticks = cbarticks, cax = cbaxes, orientation = 'horizontal')
        #cb = fig.colorbar(mappable=img, extend='both', extendrect=True, location="bottom", ticks = cbarticks, ax = ax[2,:2])
        cb.set_label(label = r'Biosphere $\bf{CO_{2}}$ flux anomaly ($\bf{\mu}$mol m$^{-2}$ s$^{-1}$)', size = 12, weight = 'bold')
    
    x,y = m(fluxdata.variables['longitude'][:],fluxdata.variables['latitude'][:])

    ##getting the limits of the map:
    x0,x1 = axs[i].get_xlim()
    y0,y1 = axs[i].get_ylim()
    map_edges = np.array([[x0,y0],[x1,y0],[x1,y1],[x0,y1]])

    ##getting all polygons used to draw the coastlines of the map
    polys = [p.boundary for p in m.landpolygons]

    ##combining with map edges
    polys = [map_edges]+polys[:]

    ##creating a PathPatch
    codes = [
    [Path.MOVETO] + [Path.LINETO for p in p[1:]]
    for p in polys
    ]   
    polys_lin = [v for p in polys for v in p]
    codes_lin = [c for cs in codes for c in cs]
    path = Path(polys_lin, codes_lin)
    patch = PathPatch(path,facecolor='white', lw=0)

    ##masking the data:
    axs[i].add_patch(patch)

plt.savefig("/projects/0/ctdas/dkivits/figures/drought/Drought_anomaly_JAS_multiyear_" + fluxfile_yearname + ".png", dpi=300, bbox_inches='tight')
plt.savefig("/projects/0/ctdas/dkivits/figures/drought/Drought_anomaly_JAS_multiyear_" + fluxfile_yearname + ".pdf", dpi=300, bbox_inches='tight')
plt.show()
