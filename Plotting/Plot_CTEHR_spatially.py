# Load necessary packages
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

########################################
##### PARSE COMMAND LINE ARGUMENTS #####
########################################

parser = argparse.ArgumentParser()
parser.add_argument(dest="fluxfile")
args = parser.parse_args()

########################################
######## DEFINE FUNCTIONS ##############
########################################
def set_size(width, fraction=1, subplots=(1, 1)):
    """Set figure dimensions to avoid scaling in LaTeX.

    Parameters
    ----------
    width: float
            Document textwidth or columnwidth in pts
    fraction: float, optional
            Fraction of the width which you wish the figure to occupy

    Returns
    -------
    fig_dim: tuple
            Dimensions of figure in inches
    """
    # Width of figure (in pts)
    fig_width_pt = width * fraction

    # Convert from pt to inches
    inches_per_pt = 1 / 72.27

    # Golden ratio to set aesthetic figure height
    # https://disq.us/p/2940ij3
    golden_ratio = (5**.5 - 1) / 2

    # Figure width in inches
    fig_width_in = fig_width_pt * inches_per_pt
    # Figure height in inches
    fig_height_in = fig_width_in * golden_ratio * (subplots[0] / subplots[1])

    fig_dim = (fig_width_in, fig_height_in)

    return fig_dim

class OOMFormatter(matplotlib.ticker.ScalarFormatter):
    def __init__(self, order=0, fformat="%1.1f", offset=True, mathText=True):
        self.oom = order
        self.fformat = fformat
        matplotlib.ticker.ScalarFormatter.__init__(self,useOffset=offset,useMathText=mathText)
    def _set_order_of_magnitude(self):
        self.orderOfMagnitude = self.oom
    def _set_format(self, vmin=None, vmax=None):
        self.format = self.fformat
        if self._useMathText:
             self.format = r'$\mathdefault{%s}$' % self.format
             

fluxfile_basename = os.path.basename(args.fluxfile).split(sep='.')[3]

# Define fluxes
fluxdata = nc.Dataset(args.fluxfile,'r')
fluxes_cte = fluxdata.variables['nep'][0,:,:]*1e6
#fluxes_cte = maskoceans(fluxdata.variables['longitude'][:],fluxdata.variables['latitude'][:], fluxes_cte, resolution='h', grid = 1.25)

fig, ax = plt.subplots(figsize = set_size(width=483.69687))
#fig = plt.figure(figsize = (16,9))
plt.subplots_adjust(top = 0.99, bottom = 0.15, left = 0, right = 1)
#plt.tight_layout()
#plt.title("Difference in carbon exchange \n for the months JJA between 2018 and 2017 - 2022")

# Define basemap to prepare for transformation later
m = Basemap(llcrnrlon=-15, llcrnrlat=33, urcrnrlon=35, urcrnrlat=72,
            #rsphere=(6378137.00, 6356752.3142),
            resolution='l', area_thresh=1000., projection='cyl', ax = ax)
            #lat_0=lat, lon_0=lon)

cbarticks = np.arange(-3,4,)
img = m.imshow(fluxes_cte,cmap='RdBu',zorder=1, vmin = cbarticks[0], vmax = cbarticks[-1])
cb = m.colorbar(mappable=img, extend='both', extendrect=True, location="bottom", ticks = cbarticks, label = 'NEP + fire anomaly ($\mu$mol m$^{-2}$ s$^{-1}$)')
#cb.ax.yaxis.set_ticks_position('left')
#cb.ax.yaxis.set_label_position('left')
#cb.ax.yaxis.labelpad = 10
#cb.ax.get_yaxis().labelpad = 10
m.drawcoastlines()

x,y = m(fluxdata.variables['longitude'][:],fluxdata.variables['latitude'][:])
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

#plt.savefig("/projects/0/ctdas/dkivits/figures/drought/Drought_anomaly_JAS_2018_2017-2022.png", dpi=300)
plt.savefig("/projects/0/ctdas/dkivits/figures/drought/Drought_anomaly_JAS_" + fluxfile_basename + ".png", dpi=300)
plt.show()
