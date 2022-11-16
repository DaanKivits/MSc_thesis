# Load necessary packages
import netCDF4 as nc
import numpy as np
import glob 
import os
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import argparse
from matplotlib.patches import Path, PathPatch
import matplotlib.ticker

########################################
##### PARSE COMMAND LINE ARGUMENTS #####
########################################

parser = argparse.ArgumentParser()
parser.add_argument(dest="wrfout_file_1")
parser.add_argument(dest="wrfout_file_2")
parser.add_argument(dest="wrfout_file_3")
projpath = '/projects/0/ctdas/dkivits/CTDAS_WRF_runs/testrun/WRF_n5/run_testcases_fwd/'

args = parser.parse_args()

def set_size(width, fraction=1, subplots=(1, 1)):
    """Set figure dimensions to avoid scaling in LaTeX.

    Parameters
    ----------
    width: float or string
            Document width in points, or string of predined document type
    fraction: float, optional
            Fraction of the width which you wish the figure to occupy
    subplots: array-like, optional
            The number of rows and columns of subplots.
    Returns
    -------
    fig_dim: tuple
            Dimensions of figure in inches
    """
    if width == 'thesis':
        width_pt = 426.79135
    elif width == 'beamer':
        width_pt = 307.28987
    else:
        width_pt = width

    # Width of figure (in pts)
    fig_width_pt = width_pt * fraction
    # Convert from pt to inches
    inches_per_pt = 1 / 72.27

    # Golden ratio to set aesthetic figure height
    # https://disq.us/p/2940ij3
    golden_ratio = (5**.5 - 1) / 2

    # Figure width in inches
    fig_width_in = fig_width_pt * inches_per_pt
    # Figure height in inches
    fig_height_in = fig_width_in * golden_ratio * (subplots[0] / subplots[1])

    return (fig_width_in, fig_height_in)

def all_equal(iterator):
  try:
     iterator = iter(iterator)
     first = next(iterator)
     return all(np.array_equal(first, rest) for rest in iterator)
  except StopIteration:
     return True

#friedemann_files = sorted(glob.glob('/projects/0/ctdas/dkivits/CTDAS_WRF_runs/testrun/WRF_n5/run_testcases_fwd/test_fwd_friedemann/wrfout_d01_2*'))
#nrt_files = sorted(glob.glob('/projects/0/ctdas/dkivits/CTDAS_WRF_runs/testrun/WRF_n5/run_testcases_fwd/test_fwd_nrt/wrfout_d01_2*'))
#files_lenlist = len(nrt_files)

# Extract the tracer simulation variables
variablelist = ['CO2_001', 'CO2_002', 'CO2_003', 'CO2_004']
fluxfilelist = ['NEP','fire','ocean','anthropogenic']

#Define basemap
geopath = '/projects/0/ctdas/dkivits/CTDAS_WRF_runs/testrun/WPS/WPS_test_fwd_friedemann/geo_em.d01.nc'
geo = nc.Dataset(geopath,'r')
wrflat = geo.variables['XLAT_M'][0,:,:]
wrflon = geo.variables['XLONG_M'][0,:,:]
wrflat_bounds = [wrflat[0,0],wrflat[-1,-1]]
wrflon_bounds = [wrflon[0,0],wrflon[-1,-1]]
we  = geo.variables['XLAT_M'].shape[2]
sn  = geo.variables['XLAT_M'].shape[1]
lon = (geo.variables['XLONG_U'][0,sn/2,we/2] + geo.variables['XLONG_U'][0,sn/2-1,we/2])/2.
lat = (geo.variables['XLAT_V'][0,sn/2,we/2] + geo.variables['XLAT_V'][0,sn/2,we/2-1])/2.
dx_dom = 100000

m = Basemap(llcrnrlon=wrflon_bounds[0], llcrnrlat=wrflat_bounds[0], urcrnrlon=wrflon_bounds[1], urcrnrlat=wrflat_bounds[1],
            #rsphere=(6378137.00, 6356752.3142),
            resolution='l', area_thresh=1000., projection='lcc',
            lat_0=lat, lon_0=lon)
                                                                                                                
# Define fluxes
NRT_nc = nc.Dataset(projpath + args.wrfout_file_1 + "/wrfout_timeavg.nc",'r')
FRI_nc = nc.Dataset(projpath + args.wrfout_file_2 + "/wrfout_timeavg.nc",'r')
ORIG_nc = nc.Dataset(projpath + args.wrfout_file_3 + "/wrfout_timeavg.nc",'r')
fluxes_nrt = NRT_nc.variables[variablelist[0]][0,1,:,:]
fluxes_fri = FRI_nc.variables[variablelist[0]][0,1,:,:]
flux_data_orig = ORIG_nc.variables["CO2_002"][0,1,:,:]
flux_data_7200 = ORIG_nc.variables["CO2_001"][0,1,:,:]
fluxes_list = [fluxes_nrt,fluxes_fri]
fluxes_origlist = [flux_data_orig, flux_data_7200]
titles_list = ["CTE-HR constructed file (flux set 5)","CTE-HR original file (flux set 6)"]
titles_origlist = ["No emissions (flux set 1)", "High emissions (flux set 3)"]

print("Constructed flux field is equal to FRI flux field: " + str(np.ma.allequal(fluxes_nrt,fluxes_fri,fill_value=True)))
print("Original flux field filled with zeroes is equal to original flux set filled with CTE-HR fluxes: " + str(np.ma.allequal(flux_data_orig,fluxes_nrt,fill_value=True)))
# Transform fluxes if necessary (transform_scalar can only downscale, not upscale. So, if a 12*15 grid is provided as input data,
# transform_scalar can't make a 390 * 250 grid out of it, there's not enough info for that in the flux field!
#flux_trans_nrt = m.transform_scalar(fluxes_nrt,lons=np.arange(-15.,35.2,0.2),lats=np.arange(33.,72.1,0.1),nx=250,ny=390,order=0)
#flux_trans_fri = m.transform_scalar(fluxes_fri,lons=np.arange(-15.,35.2,0.2),lats=np.arange(33.,72.1,0.1),nx=250,ny=390,order=0)
        
# Define axis title


#fig = plt.figure(figsize = set_size(width=483.69687, subplots=(2,2)))
fig = plt.figure(figsize = (16,9))
plt.subplots_adjust(top = 0.9, bottom = 0.1, left = 0.01, right = 0.99, wspace = -0.45, hspace = 0.1)
#plt.tight_layout()
#fig.suptitle("Average " + fluxfilelist[0] + " between 01-06-2015 and 04-06-2015")

# Create subfigs
subfigs = fig.subfigures(nrows=2, ncols=1)
for row, subfig in enumerate(subfigs):
    axs = subfig.subplots(nrows=1, ncols=2)

    # Fill subfigs with plots
    for col, ax in enumerate(axs):
        if row == 0:
            ax.set_title(titles_origlist[col],fontsize=8)
        if row == 1:
            ax.set_title(titles_list[col],fontsize=8)

        # Define basemap
        m = Basemap(llcrnrlon=wrflon_bounds[0], llcrnrlat=wrflat_bounds[0], urcrnrlon=wrflon_bounds[1], urcrnrlat=wrflat_bounds[1],
            #rsphere=(6378137.00, 6356752.3142),
            resolution='l', area_thresh=1000., projection='lcc', ax = ax,
            lat_0=lat, lon_0=lon)
        
        if row == 0:
            img = m.imshow(fluxes_origlist[col],cmap='RdBu',zorder=0)
        if row == 1:
            img = m.imshow(fluxes_list[col],cmap='RdBu',zorder=0)
    
        if col == 0:
            cb = m.colorbar(mappable=img,label='NEP CO2 mole fraction (ppm)', extend='both', location='left')
            cb.ax.yaxis.set_ticks_position('left')
            cb.ax.yaxis.set_label_position('left')
            cb.ax.yaxis.labelpad = 10
            cb.ax.tick_params(labelsize=8)
            cb.set_label(label='NEP CO2 mole fraction (ppm)', size='small')
        if col == 1:
            cb = m.colorbar(mappable=img,label='NEP CO2 mole fraction (ppm)', extend='both')
            cb.ax.get_yaxis().labelpad = 10
            cb.ax.tick_params(labelsize=8)
            cb.set_label(label='NEP CO2 mole fraction (ppm)', size='small')

        m.drawcoastlines()

        x,y = m(wrflon,wrflat)

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


        if col == 0:
            m.drawparallels(np.arange(-80.,81.,5.),labels=[False,True,False,False])
            m.drawmeridians(np.arange(-180.,181.,10.),labels=[False,False,False,True])
        else:
            m.drawparallels(np.arange(-80.,81.,5.),labels=[False,False,False,False])
            m.drawmeridians(np.arange(-180.,181.,10.),labels=[False,False,False,True])

        m.drawmapboundary(fill_color='aqua')
        
# Do the same for FRI fluxes
#m = Basemap(width=dx_dom*we,height=dx_dom*sn,
#            rsphere=(6378137.00,6356752.3142),\
#            resolution='l',area_thresh=1000.,projection='lcc',\
#            lat_1=45.,lat_2=55,lat_0=lat,lon_0=lon,
#            ax = axs[1])
#m.imshow(fluxes_fri,cmap='autumn_r',zorder=0)
##m.colorbar(label='NEP CO2 flux [mol m-2 s-1]',extend='both')
#m.drawcoastlines()
##m.fillcontinents(color='coral',lake_color='aqua', alpha=0.2)
#m.drawparallels(np.arange(-80.,81.,20.))
#m.drawmeridians(np.arange(-180.,181.,20.))
#m.drawmapboundary(fill_color='aqua')
        
        # Set basemap back to normal for the other variables
        #m = Basemap(width=dx_dom*we,height=dx_dom*sn,
        #                                rsphere=(6378137.00,6356752.3142),\
        #                                resolution='l',area_thresh=1000.,projection='lcc',\
        #                                lat_1=45.,lat_2=55,lat_0=lat,lon_0=lon)
                                                                                                                        
plt.savefig("/projects/0/ctdas/dkivits/figures/Compare_CTEHR_fluxsets.pdf", dpi=300)
plt.show()
