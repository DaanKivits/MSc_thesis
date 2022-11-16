import netCDF4 as nc
import numpy as np
import glob
import os
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

fluxpath = "/projects/0/ctdas/dkivits/CTDAS_WRF_runs/testrun/WRF_n5/run_testcases_fwd/test_fwd_original/wrfout_timeavg.nc"
Flux = nc.Dataset(fluxpath, "r")
flux_data_orig = Flux.variables["CO2_002"][0,0,:,:]
flux_data_7200 = Flux.variables["CO2_001"][0,0,:,:]

m = Basemap(width=dx_dom*we,height=dx_dom*sn,
                rsphere=(6378137.00,6356752.3142),\
                resolution='l',area_thresh=1000.,projection='lcc',\
                lat_1=45.,lat_2=55,lat_0=lat,lon_0=lon)

m.imshow(flux_data_orig, cmap="autumn_r", zorder=0)
m.colorbar(label="NEP CO2 concentration [ppm]", extend="both")
m.drawcoastlines()
m.drawmapboundary(fill_color="aqua")
plt.savefig("OriginalFluxes_zeroes", dpi=300)

m = Basemap(width=dx_dom*we,height=dx_dom*sn,
                 rsphere=(6378137.00,6356752.3142),\
                 resolution='l',area_thresh=1000.,projection='lcc',\
                 lat_1=45.,lat_2=55,lat_0=lat,lon_0=lon)

m.imshow(flux_data_7200, cmap="autumn_r", zorder=0)
m.colorbar(label="NEP CO2 concentration [ppm]", extend="both")
m.drawcoastlines()
m.drawmapboundary(fill_color="aqua")
plt.savefig("OriginalFluxes_7200", dpi=300)

