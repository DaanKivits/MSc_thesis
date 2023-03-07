# Daan Kivits, 2023

# A simple plotting script used to plot the CORINE-based SiB4 PFT landuse type map over Europe.

##############################################
########## LOAD NECCESSARY PACKAGES ##########
##############################################
from pandas import read_hdf, read_csv
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import datetime
import pandas as pd
import numpy as np
from matplotlib.dates import MO, TU, WE, TH, FR, SA, SU
import seaborn as sns
from scipy.stats import norm
from scipy import stats
import netCDF4 as nc
from mpl_toolkits.basemap import Basemap, maskoceans
from matplotlib import colors
import numpy.ma as ma
from matplotlib.patches import Path, PathPatch
import matplotlib.ticker
import sys
sys.path.insert(0, '/projects/0/ctdas/dkivits/scripts/Plotting/functions')
import plotting_functions

luvar_nc = nc.Dataset('/projects/0/ctdas/awoude/NRT/ICOS_OUTPUT/nep.202208.nc','r', format='NETCDF3_CLASSIC')
fluxname = list(luvar_nc.variables.keys())[3]
array = luvar_nc.variables[fluxname]
lu = get_lu(array)

fig, ax = plt.subplots(figsize = (16,9))
#fig = plt.figure(figsize = (16,9))
plt.subplots_adjust(top = 0.99, bottom = 0.15, left = 0, right = 1)
#plt.tight_layout()
#plt.title("Difference in carbon exchange \n for the months JJA between 2018 and 2017 - 2022")

# Define basemap to prepare for transformation later
m = Basemap(llcrnrlon=-15, llcrnrlat=33, urcrnrlon=35, urcrnrlat=72,
            #rsphere=(6378137.00, 6356752.3142),
            resolution='l', area_thresh=1000., projection='cyl', ax = ax)
            #lat_0=lat, lon_0=lon)

#bounds=[0,1,2,3,4,5,6,7,8]
#bounds = np.arange(-0.5, 19.5, 2)
bounds=[-0.5, 0.5,  1.5,  3.5,  6.5,  9.5, 12.5, 15.5, 17.5, 18.5]
cbarticks = [0,  1,  2.5,  5,  8, 11, 14, 16.5, 18]
#cbarticks = [0.5,  1.5,  3.5,  6.5,  9.5, 12.5, 15.5, 17.5]
#cbarticklabels = [0,  1,  2,  5,  8, 11, 14, 17, 18]
cbarticklabels = ["Undefined", "Desert or \nbare ground",  "Evergreen \nNeedleleaf \nForest",  "Evergreen \nBroadleaf \nForest", "Deciduous \nBroadleaf \nForest",  "Shrublands \n(Non-tundra)", "C3 General", "C3 Crops", "C4 Crops"]

cmap = colors.ListedColormap(["white", "tan", "darkgreen", "limegreen", "gold", "slategrey", "blue", "navy", "red"])
#norm = colors.BoundaryNorm(bounds, 9)
norm = colors.BoundaryNorm(bounds, 9)
img = m.imshow(lu, cmap=cmap, zorder=1, norm = norm)
#img = m.imshow(lu, cmap=cmap, zorder=1)
cb = m.colorbar(mappable=img, location="bottom", ticks = cbarticks, label = 'CORINE / SiB4 PFTs', spacing = 'uniform', norm = norm, boundaries = bounds)
#cb = m.colorbar(mappable=img, location="bottom", ticks = cbarticks, label = 'CORINE / SiB4 PFTs', spacing = 'uniform')
cb.ax.set_xticklabels(cbarticklabels, horizontalalignment = 'center')
#cb.ax.yaxis.set_ticks_position('left')
#cb.ax.yaxis.set_label_position('left')
#cb.ax.yaxis.labelpad = 10
#cb.ax.get_yaxis().labelpad = 10
m.drawcoastlines()

x,y = m(luvar_nc.variables['longitude'][:], luvar_nc.variables['latitude'][:])
#pcol = ax.pcolormesh(x,y,fluxes_cte)

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

#m.drawparallels(np.arange(-80.,81.,10),labels=[True,False,False,False])
#m.drawmeridians(np.arange(-180.,181.,10.),labels=[False,False,False,True])

plt.savefig("/projects/0/ctdas/dkivits/figures/CORINE_SiB4_PFT_landusemap.png", dpi=300)
plt.savefig("/projects/0/ctdas/dkivits/figures/CORINE_SiB4_PFT_landusemap.pdf", dpi=300)
plt.show()
