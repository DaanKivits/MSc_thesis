# Daan Kivits, 2023

# A simple plotting script used to plots the average CTE-HR biospheric CO2 flux field for 2018.

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
parser.add_argument(dest="fluxfile")
args = parser.parse_args()

#fluxfile_basename = os.path.basename(args.fluxfile).split(sep='.')[3]
fluxfile_basename = os.path.basename(args.fluxfile).split(sep='.')[0]

# Define fluxes
fluxdata = nc.Dataset(args.fluxfile,'r')
fluxes_cte = fluxdata.variables['nep'][0,:,:]* 31536000 * (1/0.083259093974539)

fig, ax = plt.subplots(figsize = set_size(width=483.69687))
plt.subplots_adjust(top = 0.99, bottom = 0.15, left = 0, right = 1)

# Define basemap to prepare for transformation later
m = Basemap(llcrnrlon=-15, llcrnrlat=33, urcrnrlon=35, urcrnrlat=72,
            #rsphere=(6378137.00, 6356752.3142),
            resolution='l', area_thresh=1000., projection='cyl', ax = ax)
            #lat_0=lat, lon_0=lon)

cbarticks = np.arange(-300,400,100)
img = m.imshow(fluxes_cte,cmap='RdBu_r',zorder=1, vmin = cbarticks[0], vmax = cbarticks[-1])
cb = m.colorbar(mappable=img, extend='both', extendrect=True, location="bottom", ticks = cbarticks, label = 'mean biosphere flux (gC m$^{-2}$ yr$^{-1}$)')
m.drawcoastlines()

x,y = m(fluxdata.variables['longitude'][:],fluxdata.variables['latitude'][:])

##getting the limits of the map:
x0,x1 = ax.get_xlim()
y0,y1 = ax.get_ylim()
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
ax.add_patch(patch)

plt.savefig("/projects/0/ctdas/dkivits/figures/drought/2018_yearmean.png", dpi=300)
plt.savefig("/projects/0/ctdas/dkivits/figures/drought/2018_yearmean.pdf", dpi=300)
plt.show()
