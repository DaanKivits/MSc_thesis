# Daan Kivits, 2023

# A simple plotting script used to check what effect the different interpolation methods of the Basemap package have on the 
# re-gridding operation of the CTE-HR fluxes.

##############################################
########## LOAD NECCESSARY PACKAGES ##########
##############################################
import netCDF4 as nc
import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import time

########################################
###### DEFINE RELEVANT PARAMETERS ######
########################################
fluxpath = '/projects/0/ctdas/awoude/NRT/ICOS_OUTPUT/nep.202205.nc'
Flux = nc.Dataset(fluxpath,'r')
flux_data = Flux.variables['nep'][1,:,:] # Extract the first timestep of the dataset

geopath = '/projects/0/ctdas/dkivits/CTDAS_WRF_runs/testrun/WPS/WPS_test_fwd_v01/geo_em.d01.nc'
geo = nc.Dataset(geopath,'r')
wrflat = geo.variables['XLAT_M'][0,:,:]
wrflon = geo.variables['XLONG_M'][0,:,:]
we  = geo.variables['XLAT_M'].shape[2]
sn  = geo.variables['XLAT_M'].shape[1]
lon = (geo.variables['XLONG_U'][0,sn/2,we/2] + geo.variables['XLONG_U'][0,sn/2-1,we/2])/2.
lat = (geo.variables['XLAT_V'][0,sn/2,we/2] + geo.variables['XLAT_V'][0,sn/2,we/2-1])/2.
dx_dom = 100000

#Define basemap
m = Basemap(width=dx_dom*we,height=dx_dom*sn,
            rsphere=(6378137.00,6356752.3142),\
            resolution='l',area_thresh=1000.,projection='lcc',\
            lat_1=45.,lat_2=55,lat_0=lat,lon_0=lon)

flux_trans1 = m.transform_scalar(flux_data,lons=np.arange(-15.,35.,0.2),lats=np.arange(33.,72.,0.1),nx=250,ny=390,order=0)
flux_trans2 = m.transform_scalar(flux_data,lons=np.arange(-15.,35.,0.2),lats=np.arange(33.,72.,0.1),nx=250,ny=390,order=1)
flux_trans_dif = flux_trans2-flux_trans1

plt.figure(figsize=(15,10))
m.imshow(flux_trans_dif,cmap='autumn_r',zorder=0)
m.colorbar(label='NEP CO2 flux [mol m-2 s-1]',extend='both')
m.drawcoastlines()
m.drawparallels(np.arange(-80.,81.,20.))
m.drawmeridians(np.arange(-180.,181.,20.))
m.drawmapboundary(fill_color='aqua')
plt.title('Bilinear versus nearest neighbour interpolation, showing CTE-HR NEP flux anomaly over Europe for ' + str(timestamp))
plt.show()

