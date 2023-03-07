# Daan Kivits, 2023

# This is an input script used to run the 'compareWRFoutput_visually.sh' script. It is used to visually compare the WRF output of 
# different WRF runs, and includes the extraction of the WRF output into variables and the plotting of these variables. 

##############################################
########## LOAD NECCESSARY PACKAGES ##########
##############################################
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

########################################
###### DEFINE RELEVANT PARAMETERS ######
########################################
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
fluxes_list = [fluxes_fri,fluxes_nrt]
fluxes_origlist = [flux_data_orig, flux_data_7200]
titles_list = ["CTE-HR (original, 5)", "CTE-HR (constructed, 6)"]
titles_origlist = ["No emissions (2)", "High emissions(3)"]

print("Constructed flux field is equal to FRI flux field: " + str(np.ma.allequal(fluxes_nrt,fluxes_fri,fill_value=True)))
print("Original flux field filled with zeroes is equal to original flux set filled with CTE-HR fluxes: " + str(np.ma.allequal(flux_data_orig,fluxes_nrt,fill_value=True)))

## Transform fluxes if necessary (transform_scalar can only downscale, not upscale. So, if a 12*15 grid is provided as input data,
## transform_scalar can't make a 390 * 250 grid out of it, there's not enough info for that in the flux field!
#flux_trans_nrt = m.transform_scalar(fluxes_nrt,lons=np.arange(-15.,35.2,0.2),lats=np.arange(33.,72.1,0.1),nx=250,ny=390,order=0)
#flux_trans_fri = m.transform_scalar(fluxes_fri,lons=np.arange(-15.,35.2,0.2),lats=np.arange(33.,72.1,0.1),nx=250,ny=390,order=0)
        
fig = plt.figure(figsize = (16,9))
plt.subplots_adjust(top = 0.9, bottom = 0.1, left = 0.01, right = 0.99, wspace = -0.45, hspace = 0.1)

# Create subfigs
subfigs = fig.subfigures(nrows=2, ncols=1)
for row, subfig in enumerate(subfigs):
    axs = subfig.subplots(nrows=1, ncols=2)

    # Fill subfigs with plots
    for col, ax in enumerate(axs):
        if row == 0:
            ax.set_title(titles_origlist[col],fontsize=12, weight = 'bold')
        if row == 1:
            ax.set_title(titles_list[col],fontsize=12, weight = 'bold')

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
            cb.ax.tick_params(labelsize=10)
            cb.set_label(label='NEP CO2 mole fraction (ppm)', size='small')
        if col == 1:
            cb = m.colorbar(mappable=img,label='NEP CO2 mole fraction (ppm)', extend='both')
            cb.ax.get_yaxis().labelpad = 10
            cb.ax.tick_params(labelsize=10)
            cb.set_label(label='NEP CO2 mole fraction (ppm)', size='small')

        # Draw the coastlines in the data
        m.drawcoastlines()

        ##getting the limits of the map:
        x,y = m(wrflon,wrflat)
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
                                                                                                
plt.savefig("/projects/0/ctdas/dkivits/figures/Compare_CTEHR_fluxsets.pdf", dpi=300)
plt.show()
