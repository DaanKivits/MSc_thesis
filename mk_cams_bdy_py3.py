#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Routines that read MACC 3D fields output and process it to format
usable to add as initial or boundary conditions in WRF. 

Authors : Ingrid Super (based on wrf_ct_ibcs.py)
Friedemann Reum (additions for tracer ensembles)
Liesbeth Florentie (minor modifications)
"""

import os
import datetime as dtm
import cartopy.crs as ccrs
import cartopy.img_transform as ctrans
import netCDF4 as nc
import numpy as np
import argparse


######### INPUT

# CAMS IN AND OUT PATHS (FIXED SETTINGS)
inpath_cams          = '/projects/0/ctdas/dkivits/DATA/CAMS_concentrations/v18r1/'
infiles_pattern      = 'cams73_v18r1_co2_conc_surface_inst_%s.nc'
infiles_date_pattern = '%Y%m'

outpath_cams         = inpath_cams + '/3d_1x1_concentrations/'

# unit of CO2 in original CAMS concentration files. Options: 'mol mol-1' or 'ppmv'
unit_cams    = 'ppmv'
# WRF-Chem requires background concentration to be in ppmv. Convert if necessary? (for addition to wrfbdy and wrfinput)
unit_convert = False


######### OPTIONS >>>

parser = argparse.ArgumentParser()
parser.add_argument(dest="wrf_dir",       help="Directory to wrfinput and wrfbdy files", type=str)
parser.add_argument(dest="tracer_prefix", help="Prefix of tracers to process. If nmembers==0, full tracer name ", type=str)
parser.add_argument(dest="nmembers",      help="Number of ensemble members to process (all equal)", type=int)
parser.add_argument(dest="start_date",    help="Format: yyyymmdd", type=str)
parser.add_argument(dest="end_date",      help="Format: yyyymmdd", type=str)

args = parser.parse_args()

inpath_wrfinput = args.wrf_dir
tracer_prefix   = args.tracer_prefix
nmembers        = args.nmembers
start_date      = dtm.datetime(int(args.start_date[:4]), int(args.start_date[4:6]), int(args.start_date[6:]))
end_date        = dtm.datetime(int(args.end_date[:4]), int(args.end_date[4:6]), int(args.end_date[6:]))

# Manual settings insteead of argparse:
#inpath_wrfinput = '/home/freum/WRF/4.1.1/WRF_n150/run/'
#nmembers        = 150
#start_date      = dtm.datetime(2015,7,1)
#end_date        = dtm.datetime(2015,8,1)  


# Create output directory if it doesn't exist yet
try:
    os.makedirs(outpath_cams)
except:
    pass


if nmembers == 0:
    # Special case: sample 'tracer_optim', don't add member suffix
    tracers = [tracer_prefix]
    nmembers = 1
else:
    tracers = [tracer_prefix+"_%03d"%n for n in range(nmembers)]

# # For testing:
# tracers = ['CO2_000', 'CO2_001']

# Input variables:
tracers_bg = ['CO2bg'] * len(tracers)

# Trigger the option in add_bg_to_wrfbdy to only process the first
# tracer and copy the results to all others.
# Otherwise, wrfbdy creation would be prohibitely time-consuming.
# Could do the same for wrfinput, but this one only takes a few
# minutes anyway.

#tracers_all_equal = len(set(tracers_bg)) == 1


########## <<< OPTIONS


# # No parallelization - doesn't make sense here because of initial conditions!
# # Get setings for dumb parallel processing
# import sys
# sys.path.append("/home/freum/ctdas/freum_projects/wrfchem/exec/")
# from da.tools.initexit import parse_options
# 
# _, settings = parse_options()
# 
# # Convert strings that represent integers in settings
# for key, val in settings.iteritems():
#     try:
#         settings[key] = int(val)
#     except ValueError:
#         pass
# wrfhelper = WRFChemHelper(settings)
# 
# # Loop over time intervals for files
# # Parallelization process only parts of the original range:
# ndays = int(float((end_date-start_date).total_seconds())/86400.)
# id0, id1 = wrfhelper.get_slicing_ids(ndays, settings['nproc'], settings['nprocs'])
# 
# end_date = start_date + dt.timedelta(days=id1)
# start_date = start_date + dt.timedelta(days=id0)




def calc_pres_fields(input,tstep):
    mf=nc.Dataset(input)
    zlen=len(mf.dimensions['level'])
    acoef=(mf.variables['ap'][1:]+mf.variables['ap'][:-1])/2. # mid level coefs
    bcoef=(mf.variables['bp'][1:]+mf.variables['bp'][:-1])/2. # mid level coefs
    Psurf=mf.variables['Psurf'][:]
    xlen=len(mf.dimensions['longitude'])
    ylen=len(mf.dimensions['latitude'])
    infield=zlen*[Psurf[tstep,:]]
    infield=np.asarray(infield)

    #the sigma pressure levels have a pressure of:
    #P=A*P(0)+B*P(s) where P(0) is a constant (1000 hPa) and P(s) is the surface pressure
    #A and B are coefficients for each layer

    ps=[]
    for k in range(zlen):
       ps.append(acoef[k]+bcoef[k]*infield[k,:,:])
    #these are pressures at half levels, similar to the eta levels in wrfinput

    return ps


def merge_netcdf_fields(input,dsname,tstep1,tstep2,addtracers=False,okdebug=True):
    """ Merge TM5 resolutions """
    # from mpl_toolkits.basemap import Basemap

    readfile=nc.Dataset(input,'r')

    """
    # LF NOT SURE WHY CODE BELOW IS AS IS...

    # get dimensions of dataset in file, also get coordinates of grid
    lons_cams=readfile.variables['longitude'][:]
    im=len(readfile.dimensions['longitude'])
    lats_cams=readfile.variables['latitude'][::-1] # turn south to north
    jm=len(readfile.dimensions['latitude'])
    dx=abs(lons_cams[1]-lons_cams[0])
    dy=abs(lats_cams[1]-lons_cams[0])
    lm=len(readfile.dimensions['level'])

    nlon_target=360
    nlat_target=180
    target=np.zeros((lm,nlat_target,nlon_target),dtype=np.float32)
    # construct its lats and lons arrays
    lons_target=np.arange(im)-180.+0.5
    lats_target=np.arange(jm)-180.+0.5

    # put field to global grid, all elements outside lon and lat are masked, set order=0 to subsample instead of interpolate
    m = Basemap(projection='cyl', resolution='i', llcrnrlat=-90, urcrnrlat = 90, llcrnrlon=-180, urcrnrlon = 180)

    target=[]
    for tstep in range(tstep1,tstep2+1):
       if dsname=='pres':
          infield=calc_pres_fields(input=input,tstep=tstep)
          infield=np.asarray(infield)
          infield = infield[:,::-1,:] # turn south to north 
       else:
          infield=readfile.variables[dsname][:]
          infield=infield[tstep,:,::-1,:] # turn south to north

       outfield=[]
       for l in range(lm):
          outfield.append(m.transform_scalar(infield[l,:,:],lons_cams,lats_cams,nlon_target,nlat_target,order=1)) # bilinear interpolation

       # put into target array
       target.append(outfield)
    """

    # LF ALTERNATIVE TO REGRID TO GLOBAL 1X1

    lons_cams_2d, lats_cams_2d = np.meshgrid(readfile.variables['longitude'][:], readfile.variables['latitude'][:])
    lm = len(readfile.dimensions['level'])
    lons_target = np.linspace(-179.5, 180.5, 360)
    lats_target = np.linspace(-89.5, 90.5, 180)
    lons_target_2d, lats_target_2d = np.meshgrid(lons_target, lats_target)
    source_proj = ccrs.PlateCarree()
    target_proj = ccrs.PlateCarree()

    target = []
    for tstep in range(tstep1,tstep2+1):
        outfield = []
        if dsname=='pres':
            infield = calc_pres_fields(input=input,tstep=tstep)
            infield = np.asarray(infield)
        else:
            infield = readfile.variables[dsname][tstep,:,:,:]

        for l in range(lm):
            outfield.append(ctrans.regrid(infield[l,:,:], lons_cams_2d, lats_cams_2d, source_proj, target_proj, lons_target_2d, lats_target_2d))
        target.append(outfield)

    target = np.asarray(target)
    readfile.close()

    return target


def make_1x1_molefractions(dates=(start_date,end_date), inpath=inpath_cams, outpath=outpath_cams, 
        inpattern=infiles_pattern, inpattern_date=infiles_date_pattern):
    import timedate_py3 as timedate  # this is a module from the pythonlib code from Wouter
    from pylab import date2num
    
    dectime0 = date2num(dtm.datetime(2000,1,1))
    (sd,ed)  = dates
    dtw      = dtm.timedelta(days=1)
    dt       = dtm.timedelta(hours=3)
    days     = timedate.timegen(sd,ed,dtw)
    for day in days:
        saveas = os.path.join(outpath,'3d_molefractions_1x1_%s.nc'%day.strftime('%Y%m%d'))
        if os.path.exists(saveas):
            print('Skipping existing file: %s' % saveas)
            continue
        else:
            print('Creating new file: %s' % saveas)

        tstep1 = (day.timetuple()[2]-1)*8
        tstep2 = tstep1+7
        ncfile = os.path.join(inpath,inpattern%(day.strftime(inpattern_date)))

        # Create background tracer field
        bg_CO2 = merge_netcdf_fields(input=ncfile,dsname='CO2',tstep1=tstep1,tstep2=tstep2)

        # Create 3d pressure field
        press = merge_netcdf_fields(input=ncfile,dsname='pres',tstep1=tstep1,tstep2=tstep2)

        nlon=press.shape[3]
        nlat=press.shape[2]
        nlev=press.shape[1]
            
        saveas = os.path.join(outpath,'3d_molefractions_1x1_%s.nc'%day.strftime('%Y%m%d'))

        # Create NetCDF output file
        fid = nc.Dataset(saveas,'w')

        fid.createDimension('lon', np.int32(nlon))
        fid.createDimension('lat', np.int32(nlat))
        fid.createDimension('level', np.int32(nlev))
        fid.createDimension('ntime', np.int32(8))
        fid.createDimension('ndate', np.int32(6))
              
        var = fid.createVariable('lon', np.float32, ('lon'))
        var[:] = (np.arange(nlon)+0.5)*360./np.float(nlon)-180.
        var = fid.createVariable('lat', np.float32, ('lat'))
        var[:] = (np.arange(nlat)+0.5)*180/np.float(nlat)-90.
        var = fid.createVariable('level', np.float32, ('level'))
        var[:] = np.arange(nlev)+1
            
        var = fid.createVariable('CO2bg', bg_CO2.dtype, ('ntime','level','lat','lon'))
        var[:] = bg_CO2[:]
        fid.variables['CO2bg'].units = 'mol mol-1'
        fid.variables['CO2bg'].long_name='mole_fraction_of_carbon_dioxide_in_air'

        var = fid.createVariable('pressure', press.dtype, ('ntime','level','lat','lon'))
        var[:] = press[:]
        fid.variables['pressure'].units ='Pa'
        fid.variables['pressure'].long_name='pressure_at_center_levels'

        idates=[]
        dates=[]
        for hour in range(0,24,3):
            dd = day + dtm.timedelta(hours=hour)
            idates.append((dd+dt/2).timetuple()[:6])
            dates.append(date2num(dd+dt/2)-dectime0)

        var = fid.createVariable('idate', np.int32, ('ntime','ndate'))
        var[:] = idates[:]
        var = fid.createVariable('date', np.float32, ('ntime'))
        var[:] = dates[:]
        fid.close()

    return None


def zero_tracers_in_wrf_ibcs(inpath,tracer_list=None,okdebug=True):
    """ Routine that zeroes the fields for supplied list of tracers in the WRF initial and boundary conditions - wrfinput & wrfbdy files"""
    if tracer_list == None: 
        print("No tracers supplied, please revise tracer_list in the input parameters")
        return None

    input_paths = [os.path.join(inpath,filename) for filename in os.listdir(inpath) if filename.startswith('wrfinput_d')]
    for each_path in input_paths:
        if os.path.exists(each_path):
            ncfile = nc.Dataset(each_path, mode='r+')
            for each_tracer in tracer_list:
                if each_tracer in list(ncfile.variables.keys()):
                    ncobj = ncfile.variables[each_tracer]
                    zero_tracer = np.zeros(ncobj.shape,float)
                    ncobj[:] = zero_tracer
                else:
                    if okdebug: print('Tracer %s skipped - not in the input file.'%each_tracer)
                    continue
            ncfile.close()
            print("Done zeroing the tracers in %s"%each_path)
        else:
            if okdebug: print("Skipping %s - file path doesn't exist."%each_path)
            continue
    bdy_path = os.path.join(inpath,'wrfbdy_d01')
    end_list = ['BXS', 'BXE', 'BYS', 'BYE', 'BTXS', 'BTXE', 'BTYS', 'BTYE']
    if os.path.exists(bdy_path):
        ncfile = nc.Dataset(bdy_path, mode='r+')
        for each_tracer in tracer_list:
            for each_ending in end_list:
                new_tracer = '_'.join([each_tracer,each_ending])                    
                if new_tracer in list(ncfile.variables.keys()):
                    ncobj = ncfile.variables[new_tracer]
                    zero_tracer = np.zeros(ncobj.shape,float)
                    ncobj[:] = zero_tracer
                else:
                    if okdebug: print('Tracer %s skipped - not in the bdy file.'%new_tracer)
                    continue
        ncfile.close()
        print("Done zeroing the tracers in %s"%bdy_path)
    else:
        if okdebug: print("Skipping %s - file path doesn't exist."%bdy_path)
        pass

    return None


def add_bg_to_wrfinput(inpath_wrf,inpath_cams=outpath_cams,tracer_names=False,bgtracer_names=False,okdebug=True):
    """Routine that interpolates bgCO2 to the wrfinput and wrfbdy files as initial and boundary conditions
    """
    from pylab import logical_and,logical_or,where
    from scipy.interpolate import griddata, interp1d

    # loading the wrf input files one by one
    input_paths = [os.path.join(inpath_wrf,filename) for filename in os.listdir(inpath_wrf) if filename.startswith('wrfinput_d')]

    # calculating the initial conditions first!
    for each_path in input_paths:
        
        ncfile     = nc.Dataset(each_path,mode='r+')
        init_times = ncfile.variables['Times'][:].tostring().decode('UTF-8')
        itime      = dtm.datetime(int(''.join(init_times[:4])),int(''.join(init_times[5:7])),int(''.join(init_times[8:10])),int(''.join(init_times[11:13])),int(''.join(init_times[14:16])),int(''.join(init_times[17:])))
        
        # finding the correct camsfile and index for the time at wrfinput
        cams_filename = os.path.join(inpath_cams,'3d_molefractions_1x1_%s.nc'%itime.date().strftime('%Y%m%d'))
        if not os.path.exists(cams_filename):
            if okdebug: print("Skipping %s - file path to cams data doesn't exist."%cams_filename)
            raise IOError('cams file with input data not available.')
        else:
            camsfile  = nc.Dataset(cams_filename)
            cams_time = camsfile.variables['idate'][:]
            cams_dt   = dtm.datetime(cams_time[1][0],cams_time[1][1],cams_time[1][2],cams_time[1][3],cams_time[1][4],cams_time[1][5]) \
                       - dtm.datetime(cams_time[0][0],cams_time[0][1],cams_time[0][2],cams_time[0][3],cams_time[0][4],cams_time[0][5])
            
            icams = None
            for i_cams in range(cams_time.shape[0]):
                timei_cams = dtm.datetime(cams_time[i_cams][0],cams_time[i_cams][1],cams_time[i_cams][2],cams_time[i_cams][3],cams_time[i_cams][4],cams_time[i_cams][5]) 
                if ~logical_and(itime >= timei_cams - cams_dt/2, itime < timei_cams + cams_dt/2):
                    continue
                else:
                    if okdebug: "WRF itime is found at cams_time with index %s"%i_cams
                    icams = i_cams
                    break
            
            if icams == None:
                if okdebug: "WRF itime was not found in the cams file."
                raise ValueError('Unknown time index')

        # at this point we have the correct cams index for the required wrf input time
        
        print("Loading wrf grid data...") 
        # we need to load wrf grid and pressure levels, interpolate cams data to the wrf grid
        wrflat = ncfile.variables['XLAT'][:]
        wrflon = ncfile.variables['XLONG'][:]
        nx = wrflat.shape[2]
        ny = wrflat.shape[1]
        nt = wrflat.shape[0]

        lonmin = wrflon.min()-3
        lonmax = wrflon.max()+3
        latmin = wrflat.min()-3
        latmax = wrflat.max()+3

        # this differs a bit from the matlab script mostly because python starts to count from 0 and matlab from 1
        Xi = np.arange(0,nx)
        Yi = np.arange(0,ny)
        [Xi,Yi] = np.meshgrid(Xi,Yi)
        Ixy_ic  = where((Xi>=0).flatten())[0].reshape((ny,nx))
      
        # calculating WRF pressure fields at initial time
        eta = ncfile.variables['ZNU'][:]
        nz = eta.shape[1]
        ptop = ncfile.variables['P_TOP'][:]

        MUB     = ncfile.variables['MUB'][:]
        MUP_ic  = ncfile.variables['MU'][:]
        MU_ic   = np.tile((MUB.flatten()[Ixy_ic].reshape(MUP_ic.shape)+MUP_ic )[:,np.newaxis,:,:] ,(1,nz,1,1))        
        outp_ic = np.tile(eta.reshape(1,nz,1,1),(nt,1,ny,nx)) * MU_ic + ptop ### gives 10x higher pressure?!?! check it out

        # converting 2D lat lon to 3D (nz, ny, nz)
        wrflon_3d = np.tile(wrflon.flatten()[Ixy_ic].reshape(Ixy_ic.shape)[np.newaxis,:,:],(nz,1,1))
        wrflat_3d = np.tile(wrflat.flatten()[Ixy_ic].reshape(Ixy_ic.shape)[np.newaxis,:,:],(nz,1,1))

        # interpolation from MACC to WRF for initial conditions
        VAR_ic = np.ones((nt,nz,ny,nx))*np.nan
        print("Making initial conditions for time %s"%itime.isoformat().replace('T','_'))

        # load the MACC grid and find the start index and number of cells needed
        camslat = camsfile.variables['lat'][:]
        camslon = camsfile.variables['lon'][:]
        Ilon = where(logical_and(camslon.flatten()>lonmin,camslon.flatten()<lonmax))[0]
        Ilat = where(logical_and(camslat.flatten()>latmin,camslat.flatten()<latmax))[0]
        Lon1 = Ilon[0]
        Lat1 = Ilat[0]
        LonN = Ilon.size
        LatN = Ilat.size        
        p1 = 0
        pN = 36 # freum: what is this?

        inlat = camsfile.variables['lat'][Lat1:Lat1+LatN]
        inlon = camsfile.variables['lon'][Lon1:Lon1+LonN]
        inp   = camsfile.variables['pressure'][icams,p1:p1+pN,Lat1:Lat1+LatN,Lon1:Lon1+LonN]

        # put all variables in common dimensions (Z, LAT, LON)
        inlon = np.tile(inlon[np.newaxis,np.newaxis,:],(pN,LatN,1))
        inlat = np.tile(inlat[np.newaxis,:,np.newaxis],(pN,1,LonN))

        # interpolating MACC data to WRF grid
        in_nx = inp.shape[2]
        in_ny = inp.shape[1]
        in_nz = inp.shape[0]
        incoord  = np.array([inlon[0].flatten(),inlat[0].flatten()]).T
        outcoord = np.array([wrflon_3d[0].flatten(),wrflat_3d[0].flatten()]).T
        Ptmp     = np.ones((in_nz,ny,nx))*np.nan

        for i,tracer_name in enumerate(bgtracer_names):
            inC     = camsfile.variables[tracer_name][icams,p1:p1+pN,Lat1:Lat1+LatN,Lon1:Lon1+LonN]
            outCtmp = np.ones((in_nz,ny,nx))*np.nan

            # Convert unit to ppmv if required
            if unit_convert and unit_cams == 'mol mol-1':
                inC = inC * 1.e6
                print('%s converted from mol mol-1 to ppmv'%tracer_name)
            elif unit_convert and unit_cams == 'ppmv':
                print('%s already in ppmv, no conversion necessary'%tracer_name)
            elif not unit_convert and unit_cams == 'mol mol-1':
                print('WARNING: %s in mol mol-1!')

            # first a horizontal interpolation at each level of cams
            for iz in range(in_nz):
                Cin = inC[iz].flatten()
                Pin = inp[iz].flatten()

                # ORIGINAL SCRIPT HAS LINEAR INTERPOLATION, doesn't work on maunaloa??
                outCtmp[iz] = griddata(incoord,Cin,outcoord,method='nearest').reshape(ny,nx)
                Ptmp[iz] = griddata(incoord,Pin,outcoord,method='nearest').reshape(ny,nx)

            print("Done with horizontal interpolation %s"%tracer_name)
            
            # secondly, vertical interpolation 
            outC = np.ones((nz,ny,nx))*np.nan
            for iy in range(ny):
                for ix in range(nx): 
                    # interp1d requires the pressure values to be increasing, so we reverse the arrays
                    Pin  = Ptmp[:,iy,ix][::-1]
                    Pout = outp_ic[0,:,iy,ix][::-1].filled(fill_value=np.nan)
                    Cin  = outCtmp[:,iy,ix][::-1]
                    valid_Pout    = where(logical_and(Pout<Pin[-1],Pout>=Pin[0]))[0]
                    invalid_Pout  = where(Pout>=Pin[-1])[0]
                    invalid_Pout2 = where(Pout<Pin[0])[0]
                    f_vert        = interp1d(Pin,Cin,kind='nearest') ## could be linear, works 
                    reverse_Cout  = f_vert(Pout[valid_Pout])    
                    
                    # we fix the values near surface of WRF outside the pressure of MACC to be equal to MACC surface values and reverse back the rest of the interpolated values
                    valid_Pout = np.asarray(valid_Pout)
                    outC[nz-1-invalid_Pout,iy,ix]  = outCtmp[0,iy,ix]
                    outC[nz-1-invalid_Pout2,iy,ix] = 0
                    outC[nz-1-valid_Pout,iy,ix]    = reverse_Cout
            
            print("Done with vertical interpolation %s"%tracer_name)
            
            if where(np.isnan(outC))[0].size>0 :
                print("Nan found in output %s fields, please revise script!"%tracer_name)
                raise ValueError("NaN found in output fields")
            
            VAR_ic[0] = outC
            tracer    = ncfile.variables[tracer_names[i]]
            tracer[:] = VAR_ic
        
        ncfile.close()
        camsfile.close()
        
        print("Initial conditions written in %s"%each_path)
    
    return None


def add_bg_to_wrfbdy(inpath_wrf,inpath_cams=outpath_cams,tracer_names=False,bgtracer_names=False,okdebug=True,copy_from_first=False):
    """
    Routine that interpolates bgCO2 to the wrfinput and wrfbdy files as initial and boundary conditions
    freum: New argument copy_from_first. If True: do the calculations only for the first tracer, and  copy the result to all others.
           This is used for making tracer ensembles for CTDAS.
    """
    from pylab import logical_and,logical_or,where,date2num
    from scipy.interpolate import griddata, interp1d

    if copy_from_first:
        # Reduce tracer list to first entry, save other tracers
        tracer_names_copy_to   = tracer_names[1:]
        bgtracer_names_copy_to = tracer_names[1:]

        tracer_names   = [tracer_names[0]]
        bgtracer_names = [bgtracer_names[0]]
        
    ncfile     = nc.Dataset(os.path.join(inpath_wrf,'wrfinput_d01'))
    bdyfile    = nc.Dataset(os.path.join(inpath_wrf,'wrfbdy_d01'),mode='r+')
    bdy_times  = bdyfile.variables['Times'][:]
    bdy_times2 = [x.tostring() for x in bdy_times]
    btimes     = [dtm.datetime(int(b[:4]),int(b[5:7]),int(b[8:10]),int(b[11:13]),int(b[14:16]),int(b[17:])) for b in bdy_times2]
    dt         = (date2num(btimes[-1])-date2num(btimes[-2]))*86400 # interval between calls in seconds - needed for the tendencies
    btimes.append(btimes[-1]+(btimes[-1]-btimes[-2])) # add one extra time step needed to calculate the tendency of the last time step

    # Load the WRF data
    print("Loading wrf grid data...")
    wrflat = ncfile.variables['XLAT'][:]
    wrflon = ncfile.variables['XLONG'][:]
    eta    = ncfile.variables['ZNU'][:]
    ptop   = ncfile.variables['P_TOP'][:]

    lonmin = wrflon.min()-3
    lonmax = wrflon.max()+3
    latmin = wrflat.min()-3
    latmax = wrflat.max()+3

    nx = wrflat.shape[2]
    
    ny = wrflat.shape[1]
    nt = len(btimes)
    nz = eta.shape[1]
    nb = len(bdyfile.dimensions['bdy_width'])

    # this differs a bit from the matlab script mostly because python starts to count from 0 and matlab from 1
    Xi      = np.arange(0,nx)
    Yi      = np.arange(0,ny)
    [Xi,Yi] = np.meshgrid(Xi,Yi)
    Ixy = {}
#    Ixy['ic'] = where((Xi>=0).flatten())[0].reshape((ny,nx))
    Ixy['BXS'] = where((Xi<nb).flatten())[0].reshape((ny,nb))
    Ixy['BXE'] = where((Xi>=nx-nb).flatten())[0].reshape((ny,nb))
    Ixy['BYS'] = where((Yi<nb).flatten())[0].reshape((nb,nx))
    Ixy['BYE'] = where((Yi>=ny-nb).flatten())[0].reshape((nb,nx))    

    bdys = ['BXS','BXE','BYS','BYE']
    tdys = ['BTXS','BTXE','BTYS','BTYE']

    # calculating WRF pressure fields for the entire period
    print("Calculating WRF pressure fields..")
    MUB  = ncfile.variables['MUB'][0]
    MUP  = {}
    MU   = {}
    outp ={}
    for i in range(len(bdys)):
        bdy = bdys[i]
        MUP[bdy]      = np.zeros((bdyfile.variables['MU_'+bdy].shape[0]+1,bdyfile.variables['MU_'+bdy].shape[1],bdyfile.variables['MU_'+bdy].shape[2]))
        MUP[bdy][:-1] = bdyfile.variables['MU_'+bdy][:]
        MUP[tdys[i]]  = bdyfile.variables['MU_'+tdys[i]][:]
        MUP[bdy][-1]  = MUP[bdy][nt-1]+MUP[tdys[i]][nt-2]*dt 
    
    # pressure state is not available at the last time step, therefore we calculate that one from the tendencies
    for i in range(len(bdys)):
        bdy = bdys[i]
        if bdy in ['BXS','BXE']:
            aMU     = MUB.flatten()[Ixy[bdy]]
            MU[bdy] = np.tile((np.tile(np.transpose(aMU,(1,0))[np.newaxis,:,:],(nt,1,1))+MUP[bdy])[:,:,np.newaxis,:],(1,1,nz,1))
            nh      = ny
            nttmp   = nt
        elif bdy in ['BYS','BYE']:
            aMU     = MUB.flatten()[Ixy[bdy]]
            MU[bdy] = np.tile((np.tile(aMU[np.newaxis,:,:],(nt,1,1))+MUP[bdy])[:,:,np.newaxis,:],(1,1,nz,1))
            nh      = nx
            nttmp   = nt      
        outp[bdy]= np.tile(eta.reshape(1,1,nz,1),(nttmp,nb,1,nh)) * MU[bdy] + ptop
    
    # converting 2D lat lon to 3D (nz, ny, nz)
    wrflon_3d = {}
    wrflat_3d = {}
    for i in range(len(bdys)):
        bdy = bdys[i]
        wrflon_3d[bdy] = np.tile(wrflon.flatten()[Ixy[bdy]][np.newaxis,:,:],(nz,1,1))
        wrflat_3d[bdy] = np.tile(wrflat.flatten()[Ixy[bdy]][np.newaxis,:,:],(nz,1,1))
        
    for j,tracer in enumerate(bgtracer_names):
        # initialize variables to calculate and write
        VAR = {}
        for i in range(len(bdys)):
            bdy = bdys[i]
            tdy = tdys[i]
            if bdy in ['BXS','BXE']:
                VAR[bdy] = np.ones((nt,nb,nz,ny))*np.nan
                VAR[tdy] = np.ones((nt,nb,nz,ny))*np.nan
            elif bdy in ['BYS','BYE']:
                VAR[bdy] = np.ones((nt,nb,nz,nx))*np.nan
                VAR[tdy] = np.ones((nt,nb,nz,nx))*np.nan
        
        # loop over times
        for it in range(nt):
            itime = btimes[it]
            print("Starting calculations for time step %s at %s"%(itime.isoformat(),str(dtm.datetime.utcnow())))
            
            # finding the correct camsfile and index for the time at wrfinput
            cams_filename = os.path.join(inpath_cams,'3d_molefractions_1x1_%s.nc'%itime.date().strftime('%Y%m%d'))
            if not os.path.exists(cams_filename):
                if okdebug: print("Skipping %s - file path to cams data doesn't exist."%cams_filename)
                raise IOError('cams file with input data not available.')
            else:
                camsfile = nc.Dataset(cams_filename)
                cams_time = camsfile.variables['idate'][:]
                cams_dt = dtm.datetime(cams_time[1][0],cams_time[1][1],cams_time[1][2],cams_time[1][3],cams_time[1][4],cams_time[1][5]) \
                          - dtm.datetime(cams_time[0][0],cams_time[0][1],cams_time[0][2],cams_time[0][3],cams_time[0][4],cams_time[0][5])
                icams = None
                for i_cams in range(cams_time.shape[0]):
                     timei_cams = dtm.datetime(cams_time[i_cams][0],cams_time[i_cams][1],cams_time[i_cams][2],cams_time[i_cams][3],cams_time[i_cams][4],cams_time[i_cams][5])
                     if ~logical_and(itime >= timei_cams - cams_dt/2 , itime < timei_cams + cams_dt/2):
                         continue
                     else:
                         if okdebug: "WRF itime is found at cams_time with index %s" %i_cams
                         icams = i_cams
                         break

                if icams == None:
                   if okdebug: "WRF itime was not found in the camsfile."
                   raise ValueError('Unknown time index')

            # at this point we have the correct cams index for the required wrf time step
    
            # load the cams grid and find the start index and number of cells needed
            print("Loading cams grid data...")
            camslat = camsfile.variables['lat'][:]
            camslon = camsfile.variables['lon'][:]
            Ilon = where(logical_and(camslon.flatten()>lonmin,camslon.flatten()<lonmax))[0]
            Ilat = where(logical_and(camslat.flatten()>latmin,camslat.flatten()<latmax))[0]
            Lon1 = Ilon[0]
            Lat1 = Ilat[0]
            LonN = Ilon.size
            LatN = Ilat.size        
            p1 = 0
            pN = 36
        
            inlat = camsfile.variables['lat'][Lat1:Lat1+LatN]
            inlon = camsfile.variables['lon'][Lon1:Lon1+LonN]
            inp = camsfile.variables['pressure'][icams,p1:p1+pN,Lat1:Lat1+LatN,Lon1:Lon1+LonN]
            	
            inlon = np.tile(inlon[np.newaxis,np.newaxis,:],(pN,LatN,1))
            inlat = np.tile(inlat[np.newaxis,:,np.newaxis],(pN,1,LonN))
            
            # interpolating cams data to WRF grid
            print("Interpolating cams to WRF grid..")
            inC = camsfile.variables[tracer][icams,p1:p1+pN,Lat1:Lat1+LatN,Lon1:Lon1+LonN]
            
            # Convert unit to ppmv if required
            if unit_convert and unit_cams == 'mol mol-1':
                inC = inC * 1.e6
                print('%s converted from mol mol-1 to ppmv'%tracer)
            elif unit_convert and unit_cams == 'ppmv':
                print('%s already in ppmv, no conversion necessary'%tracer)
            elif not unit_convert and unit_cams == 'mol mol-1':
                print('WARNING: %s in mol mol-1!')

            for i in range(len(bdys)):
                bdy = bdys[i]
                tdy = tdys[i]
    
                in_nx = inp.shape[2]
                in_ny = inp.shape[1]
                in_nz = inp.shape[0]				
                incoord  = np.array([inlon[0].flatten(),inlat[0].flatten()]).T
                outcoord = np.array([wrflon_3d[bdy][0].flatten(),wrflat_3d[bdy][0].flatten()]).T
                Ptmp     = np.ones((in_nz,wrflon_3d[bdy].shape[1],wrflon_3d[bdy].shape[2]))*np.nan
                outCtmp  = np.ones((in_nz,wrflon_3d[bdy].shape[1],wrflon_3d[bdy].shape[2]))*np.nan

                # first a horizontal interpolation at each level of cams
                for iz in range(in_nz):
                    Zin = inC[iz].flatten()
                    Pin = inp[iz].flatten()
                    # ORIGINAL SCRIPT HAS LINEAR INTERPOLATION, doesn't work on maunaloa??
                    outCtmp[iz] = griddata(incoord,Zin,outcoord,method='nearest').reshape(wrflon_3d[bdy].shape[1],wrflon_3d[bdy].shape[2])
                    Ptmp[iz]    = griddata(incoord,Pin,outcoord,method='nearest').reshape(wrflon_3d[bdy].shape[1],wrflon_3d[bdy].shape[2])
    
                print("Done with horizontal interpolation")

                # secondly, vertical interpolation 
                outC = np.ones((wrflon_3d[bdy].shape[1],nz,wrflon_3d[bdy].shape[2]))*np.nan
                for iy in range(wrflon_3d[bdy].shape[1]):
                    for ix in range(wrflon_3d[bdy].shape[2]): 
                        
                        # interp1d requires the pressure values to be increasing, so we reverse the arrays
                        Pin = Ptmp[:,iy,ix][::-1]
                        if bdy in ['BXS','BXE']:
                           Pout = outp[bdy][it,ix,:,iy][::-1].filled(fill_value=np.nan)
                        elif bdy in ['BYS','BYE']:
                           Pout = outp[bdy][it,iy,:,ix][::-1].filled(fill_value=np.nan)
                        
                        Cin           = outCtmp[:,iy,ix][::-1]
                        valid_Pout    = where(logical_and(Pout<Pin[-1],Pout>=Pin[0]))[0]
                        invalid_Pout  = where(Pout>=Pin[-1])[0]
                        invalid_Pout2 = where(Pout<Pin[0])[0]
                        f_vert        = interp1d(Pin,Cin,kind='nearest') ## could be linear, works 
                        reverse_Cout  = f_vert(Pout[valid_Pout])    

                        # we fix the values near surface of WRF outside the pressure of cams to be equal to cams surface values 
                        # and reverse back the rest of the interpolated values
                        valid_Pout = np.asarray(valid_Pout)
                        outC[iy,nz-1-invalid_Pout,ix]  = outCtmp[0,iy,ix]
                        outC[iy,nz-1-invalid_Pout2,ix] = 0
                        outC[iy,nz-1-valid_Pout,ix]    = reverse_Cout
                
                print("Done with vertical interpolation")
                
                if where(np.isnan(outC))[0].size>0 :
                    print("Nan found in output %s fields, please revise script!"%tracer_name)
                    raise ValueError("NaN found in output fields")
                if bdy in ['BXS','BXE']: 
                    VAR[bdy][it] = np.transpose(outC,(2,1,0))
                elif bdy in ['BYS','BYE']:
                    VAR[bdy][it] = np.transpose(outC,(0,1,2))
                
            camsfile.close()
    
        # calculate tendencies
        print("Calculating tendencies..")
        for i in range(len(bdys)):
            bdy = bdys[i]
            tdy = tdys[i]
            VAR[tdy] = (VAR[bdy][1:nt]-VAR[bdy][:nt-1])/dt
    
        # for 'BXE', 'BYE' we need to reverse the direction of the boundary grid
        VAR['BXE']  = VAR['BXE'][:,::-1,:,:]
        VAR['BTXE'] = VAR['BTXE'][:,::-1,:,:]
        VAR['BYE']  = VAR['BYE'][:,::-1,:,:]
        VAR['BTYE'] = VAR['BTYE'][:,::-1,:,:]
    
        # writing the bc and tendencies in the bdy
        for key in list(VAR.keys()):
            varname = tracer_names[j]+'_'+key
            tracern = bdyfile.variables[varname]
            tracern[:] = VAR[key][:nt-1]
            
    bdyfile.close()
    ncfile.close()

    ### freum >>>
    if copy_from_first:
        # Copy first tracer to all others
        print("Copying results of first tracer in list to all others")
        cmd_template = "ncap2 -A -s '{var_out}={var_in};' {file_name} {file_name}"
    
        # Loop over all variable names
        keys = bdys+tdys
        for key in keys:
            # I open the file in this loop so that processing per key is independent in case something goes wrong.
            bdyfile = nc.Dataset(os.path.join(inpath_wrf,'wrfbdy_d01'),mode='r+')
    
            val_in = bdyfile.variables[tracer_names[0] + "_" + key][:]
            print("Processing " + key + "...")
    
            # Loop over all tracers
            for tracer in tracer_names_copy_to:
                var_name_out = tracer + "_" + key
                print("Making " + var_name_out + "...")
                # Copy variable
                var_out = bdyfile.variables[var_name_out]
                var_out[:] = val_in[:]
    
                # # Version with ncap2: didnt work because a) took forever (1-2min before crash) and
                # # b) failed: returned -11 (no idea what that is) and only wrote the first ~58000
                # # of the ~250000 lines in the ncdump, the rest stayed NaNf.
                # # Maybe the reason why it took so long was that ncap2 makes a temporary copy of the
                # # entire file (9.5 GB in this case) each time it's called.
                # cmd = cmd_template.format(var_out=var_out, var_in=var_in, file_name=file_name)
                # p = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
                # # Wait for process to finish - I dont know what happens if I try to write simultaneously.
                # p.wait()
                # # Check for errors
                # out, err = p.communicate()
                # if p.returncode != 0:
                #     msg = "key %s, tracer %s: returncode: %s" % (key, tracer, str(p.returncode))
                #     print msg
                #     print "Process error output: "
                #     if not err is None:
                #         for line in err:
                #             print line.rstrip()
                #     raise OSError(msg)
                # Variable existed before, so attributes are not affected - done!
    
            bdyfile.close()
    ### <<< freum

    print("Boundary conditions written in %s"%os.path.join(inpath_wrf,'wrfbdy_d01'))

    return None


# ====================
# EXECUTION
# ====================

make_1x1_molefractions(dates=(start_date,end_date), inpath=inpath_cams, outpath=outpath_cams, inpattern=infiles_pattern, inpattern_date=infiles_date_pattern)

#zero_tracers_in_wrf_ibcs(inpath_wrfinput,tracer_list=tracers,okdebug=True)

add_bg_to_wrfinput(inpath_wrf=inpath_wrfinput, inpath_cams=outpath_cams, tracer_names=tracers, bgtracer_names=tracers_bg, okdebug=True)

#add_bg_to_wrfbdy(inpath_wrf=inpath_wrfinput, inpath_cams=outpath_cams, tracer_names=tracers, bgtracer_names=tracers_bg, okdebug=True, copy_from_first=tracers_all_equal)

add_bg_to_wrfbdy(inpath_wrf=inpath_wrfinput, inpath_cams=outpath_cams, tracer_names=tracers, bgtracer_names=tracers_bg, okdebug=True, copy_from_first=False)

