import pickle
import datetime as dt

# -----
# INPUT
# -----

# run identifier
RUNNAME = 'test_fwd_highres'

# path of the projects folder that contains data and scripts
PROJDIR = '/projects/0/ctdas/dkivits/'

# directory in which folders with output will be created for run
RUNDIR  = PROJDIR + 'CTDAS_WRF_runs/testrun/'

# directory in which input specific for run is located
INPUTDIR = RUNDIR + '/input/' + RUNNAME

# directory that contains compiled WRF and WPS code
SRCDIR = PROJDIR + 'WRF_clean/test_4.1.1_compiled_CO2/'

# directory for WRF runs, should be located on scratch
WRFRUNDIR = '/gpfs/scratch1/shared/rdkok/4Daan/WRF_runs/testrun/'

# ensemble members
NMEMBERS_FWD = 5
NMEMBERS_ENS = 150


# -----
# CREATE DICTIONARY WITH SETTINGS
# -----

datadict = {

# GENERAL
'period'            : [dt.datetime(2015, 6, 1, 0, 0), dt.datetime(2015, 6, 4, 0, 0)],
'run_name'          : RUNNAME,

# METEO SETTINGS
# use ERA5 or ERA-Interim meteo
'driving_meteo'     : 'ERA5',
# folder that contains .grb1 and .grb2 meteo data
'grib_dir'          : PROJDIR + '/DATA/ERA5_ml',
# naming convention for date of meteo files (ERA5_<grib_time_format>_<hh>00_ml.grb2 and ERA5_<grib_time_format>_sfc.grb1)
'grib_time_format'  : '%Y%m%d',
# is meteo downloaded for full atmosphere (start=0) or just subset (start=47, other options not implemented)
'meteo_start_level' : 47,

# WPS SETTINGS
# namelist for wps (domain, resolution, wrf settings...)
'fp_nml_wps_src'    : INPUTDIR + '/namelist.' + RUNNAME+'.wps',
# folder with clean but compiled WPS files 
'wps_src_dir'       : SRCDIR + '/WPS/',
# wps_run_dir will be created as copy of wps_src_dir and will contain files for run
'wps_run_dir'       : RUNDIR + '/WPS/WPS_' + RUNNAME,
# required to run WPS (links will be created to fp_ecmwf_coeffs and fp_vtable)
'fp_ecmwf_coeffs'   : SRCDIR + '/ecmwf_coeffs',
'fp_vtable'         : SRCDIR + '/Vtable.ERA5',
'geog_data_path'    : PROJDIR + '/DATA/WPS_GEOG_complete',
# script to add boundary and initial conditions to wrfbdy and wrfinput files
'bc_ic_script'      : PROJDIR + '/scripts/mk_cams_bdy_py3.py',

# WRF SETTINGS
# environment variables to run WRF
'bashrc'            : PROJDIR + '/scripts/.bashrc_wrf',
# namelist for WRF (domain, resolution, wrf settings...)
'fp_nml_input_src'  : INPUTDIR + '/namelist.' + RUNNAME + '.input',
# folders with clean but compiled WRF files, separate compile necessary for different number of tracers!
'wrf_fwd_src_dir'   : SRCDIR + '/WRF_n%s/run' %(NMEMBERS_FWD),
'wrf_ens_src_src_dir':SRCDIR + '/WRF_n%s/run' %(NMEMBERS_ENS),
# run directories, created as copies of above source directories, will contain preprocessed files for run
'wrf_fwd_run_dir'   : RUNDIR + '/WRF_n%s/run_testcases_fwd/%s/' %(NMEMBERS_FWD, RUNNAME),
'wrf_ens_src_dir'   : RUNDIR + '/WRF_n%s/run_testcases_ens/%s/' %(NMEMBERS_ENS, RUNNAME),
# this one is used as run dir for ctdas, it is a copy of wrf_ens_src_dir. This is where WRF will actually be ran! 
'wrf_ens_run_dir'   : WRFRUNDIR + '/WRF_n%s/run_testcases_ens/%s/' %(NMEMBERS_ENS, RUNNAME),
# file that determines WRF output that is written
'fp_iofields'       : INPUTDIR + '/iofields.txt',
# files with fluxes for forward and ensemble run
'wrfchemi_fwd_dir'  : INPUTDIR + '/emissions/fwd',
'wrfchemi_ens_dir'  : INPUTDIR + '/emissions/ens',
'wrfchemi_pattern'  : 'wrfchemi_d??_????-??-??_??:00:00',
'emis_zdim'         : 9,

# WPS AND WRF RUNNER (PREPROCESSING STEP)
'wps_wrf_runner'    : PROJDIR + '/scripts/run_wps_wrf_fwd_withrunner.py',

# CTDAS SETTINGS: GENERAL
'ctdas_src_dir'     :  PROJDIR + '/CTDAS-WRF/',
'ctdas_run_dir'     :  RUNDIR + '/ctdas_runs/' + RUNNAME,
'time.cycle'        :  1,
'time.nlag'         :  2,
'ctdas.ncycles'     :  2,
'nmembers_fwd'      :  NMEMBERS_FWD,
'nmembers_ens'      :  NMEMBERS_ENS,
# wall clock time per ctdas job
'ctdas.ntime'       :  '3:00:00',
# number of threads, depends on observation operator (= WRF) since ctdas is not parallelized
'ctdas.ntasks'      :  1,
# template files for running CTDAS, will be copied to CTDAS run dir and modified based on settings specified here
'ctdas_template_rc'    : INPUTDIR + '/ctdas_templates/ctdas_template.rc',
'ctdas_template_py'    : INPUTDIR + '/ctdas_templates/ctdas_template.py',
'ctdas_template_jb'    : INPUTDIR + '/ctdas_templates/ctdas_template.jb',
'dasystem_template_rc' : INPUTDIR + '/ctdas_templates/dasystem_wrfchem_template.rc',
'wrfchem_template_rc'  : INPUTDIR + '/ctdas_templates/obsoper_wrfchem_template.rc',

# CTDAS SETTINGS: STATE VECTOR
'regionsfile'       :  INPUTDIR + '/regions_'+RUNNAME+'.nc',
'n_emis_proc'       :  0, # 3,
'sigma_scale_list'  :  [], #[0.3, 0.5, 0.5], # should contain as many values as n_emis_proc
# optimize offset in initial/boundary conditions?
'do_offset'         :  False,
'sigma_offset'      :  'no clue yet what this should be...', # only used if do_offset = True

# OBSERVATIONS
'obs_path'          :  INPUTDIR + '/obs/',
'obs_pattern'       :  'xco2_' + RUNNAME + '_pri_truth7.0_<YYYYMMDD>.nc',
# retrieval definition
'level_def'         :  'layer_average',
# xco2_weights currently only holds observation rejection threshold (based on todel-data mismatch)
'xco2_weights_rc'   :  INPUTDIR + '/obs/xco2_weights.rc',
# sampling options: 1 = in center of footprint, value > 1 uses footprint corners to calculate sample as mean of multiple points within footprint
'footprint_samples_dim' :  1,

#'atmos_data'        :  'artificial',
#'dataset'           :  'pseudo',
#'data_version'      :  None,
#'dat_dir'           :  None,
#'dat_spatial_extent':  None,
#'ave'               :  0,
#'data_unc_unit'     :  'absolute',
#'xco2_unc_rnd'      :  0.5,
#'xco2_unc_sys'      :  0.1,
#'unc_version'       :  3.3,
#'unc_sys_corr_len'  :  500,
#'unc_sys_rng_entropy':  308726478643685548151399757865390919267,
#'prior_version'     :  1,
#'dat_cloud_filter'  :  None,
#'sample_coords_file_format'  :  'sample_coords_'+RUNNAME+'_%s_%s.nc',
#'sample_coords_file_tformat' :  '%Y%m%d',
#'sample_coords_file_dir'     :  INPUTDIR + '/sample_coords_'+RUNNAME,
#'obs_pri_pattern'   :  'xco2_'+RUNNAME+'_pri_<YYYYMMDD>.nc',
#'dt_obs_file'       :  '1 day, 0:00:00',
#'stagger_emis'      :  False,
#'wrfchemi_fwd_template_pattern_prefix' :  'wrfchemi_d%02d_template_fwd',
#'wrfchemi_ens_template_pattern_prefix' :  'wrfchemi_d%02d_template_ens',
#'ti_emis'           :  '2015-06-01 00:00:00',
#'te_emis'           :  '2015-06-04 00:00:00',
#'truth_version'     :  7.0,
#'obs_pri_truth_pattern' :  'xco2_'+RUNNAME+'_pri_truth7.0_<YYYYMMDD>.nc',
#'truth_file_format'     :  'sampled_sample_coords_'+RUNNAME+'_%s_%s.nc',
#'members'           :  {'mem_true': 0, 'mem_icbc': 1},
#'ti_dat'            :  '2015-06-01 12:00:00',
#'te_dat'            :  '2015-06-04 12:00:00'

}

# write to file
file = open('%s/settings_file_%s.pkl' %(RUNDIR,RUNNAME),'wb')
pickle.dump(datadict, file)
file.close()

# open file to check settings
#file = open('test_settings_file.pkl','rb')
#data = pickle.load(file)
#file.close()
#print(data)

