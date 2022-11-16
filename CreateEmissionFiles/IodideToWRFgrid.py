import netCDF4 as nc
import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import time

#Load data
iodidepath = '/home/WUR/barte035/WRFChem/AuxillaryDataPreprocessing/Iodide/predicted_iodide_0.125x0.125_Ns_All_Ensemble_members.nc'
Iodide = nc.Dataset(iodidepath,'r')
iodide_data = Iodide.variables['Ensemble_Monthly_mean'][7,:,:] #BE SURE TO SELECT THE RIGHT DATA (Data is monthly, from January (=0) to December (=11))

#Define basemap
geopath = '/home/WUR/barte035/WRFChem/WPS/geo_em.d01.nc'
geo = nc.Dataset(geopath,'r')
wrflat = geo.variables['XLAT_M'][0,:,:]
wrflon = geo.variables['XLONG_M'][0,:,:]
we  = geo.variables['XLAT_M'].shape[2]
sn  = geo.variables['XLAT_M'].shape[1]
lon = (geo.variables['XLONG_U'][0,sn/2,we/2] + geo.variables['XLONG_U'][0,sn/2-1,we/2])/2.
lat = (geo.variables['XLAT_V'][0,sn/2,we/2] + geo.variables['XLAT_V'][0,sn/2,we/2-1])/2.
dx_dom = 30000

plt.figure(figsize=(15,10))
m = Basemap(width=(dx_dom*(we-10))-1,height=(dx_dom*(sn-10))-1,resolution='l',area_thresh=1000.,projection='npstere',boundinglat=57.3,lat_0=90.,lon_0=0.)

#Select correct datapoints
iodide_data = np.delete(iodide_data,np.s_[0:1081],axis=0)

print(iodide_data.shape)
print(np.arange(-180.,180.,0.125).shape)
print(np.arange(60.,90.,0.125).shape)

Iodide_trans = m.transform_scalar(iodide_data,lons=np.arange(-180.,180.,0.125),lats=np.arange(45.,90.,0.125),nx=249,ny=249,order=1)
print(geo.variables['LANDMASK'][0,:,:])
Iodide_trans = np.where(geo.variables['LANDMASK'][0,:,:] != 1, Iodide_trans, np.nan)

print(Iodide_trans[:,0])

#lons=np.arange(-180,180.1,0.1),lats=np.arange(60,90.1,0.1)
m.imshow(Iodide_trans,vmin=30,vmax=130,cmap='autumn_r')
m.colorbar(label='Oceanic [I$_{aq}$$^{-}$] concentration [nM]',extend='both')
m.drawcoastlines(linewidth=1., linestyle='solid', color='k', antialiased=1, zorder=1001)
m.drawparallels(np.arange(-80.,81.,10.),labels=[False,False,False,False],zorder=1001)
m.drawmeridians(np.arange(-180.,181.,20.),labels=[True,False,False,True],zorder=1001)
plt.savefig('IodideMapAugust.png',dpi=300)
time.sleep(5)
plt.show()

Iodide_trans = np.where(np.isnan(Iodide_trans), 0.1, Iodide_trans)

'''
#Add to wrflowinput
wrflowinputpath = '/home/WUR/barte035/WRFChem/AuxillaryDataPreprocessing/wrflowinp_d01_test'
wrflowinput = nc.Dataset(wrflowinputpath,'r+',format='NETCDF4')
print(wrflowinput.dimensions.values())
iodideocean = wrflowinput.createVariable('I_OCEAN','f4',('Time','south_north','west_east'))
iodideocean.units = 'nM'
for i in range(iodideocean.shape[0]):
	iodideocean[i,:,:] = Iodide_trans[:,:]
	print(i)
print(iodideocean.shape)
print(Iodide_trans.shape)
wrflowinput.close()
'''

#Add to wrfinput
wrfinputpath = '/home/WUR/barte035/WRFChem/AuxillaryDataPreprocessing/wrfinput_d01_test'
wrfinput = nc.Dataset(wrfinputpath,'r+',format='NETCDF4')
print(wrfinput.dimensions.values())
iodideocean = wrfinput.createVariable('I_OCEAN','f4',('Time','south_north','west_east'))
iodideocean.FieldType = int(104)
iodideocean.MemoryOrder = "XY "
iodideocean.description = "[I-] oceanic concentration"
iodideocean.units = 'nM'
iodideocean.stagger = ""
iodideocean.coordinates = "XLONG XLAT XTIME"
iodideocean[0,:,:] = Iodide_trans[:,:]
print(iodideocean.shape)
print(Iodide_trans.shape)
wrfinput.close()
