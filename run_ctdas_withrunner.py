#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  4 15:17:49 2020

Launch a CTDAS job based on runner file

@author: Friedemann
Modifications by Liesbeth

"""
import os
import sys
import subprocess
import shutil
import re
import argparse
import pickle
import datetime as dtm
import netCDF4 as nc

HOME = os.getenv("HOME")

#####################
##### EXECUTION #####
#####################

#########################
##### Read settings #####
#########################
parser = argparse.ArgumentParser()
parser.add_argument(dest="settings_file", help="Path to settings file (pickled dictionary)")
parser.add_argument('--overwrite-ctdas',dest='overwrite_ctdas', help="Overwrite CTDAS run when directory already exists?", default=False)
#if running_as_cluster_job:
#    parser.add_argument(dest="mem_MB", help="Memory in MB available for job step creation.")

args = parser.parse_args()

with open(args.settings_file, "rb") as f:
    s = pickle.load(f)

# add to path:
sys.path.append(s['ctdas_src_dir'])
import da.tools.rc as rc


######################################
##### Create CTDAS run directory #####
######################################

src_dir      = os.path.join(s["ctdas_src_dir"], "da")
exec_run_dir = os.path.join(s["ctdas_run_dir"], "exec")

if os.path.exists(exec_run_dir):
    if args.overwrite_ctdas:
        print('CTDAS run directory (%s) already exists, files will be overwritten.' %exec_run_dir)
        rsyncopt = '-r'
    else:
        print(exec_run_dir + " exists and overwrite-ctdas is False. Existing files will not be modified.")
        rsyncopt = '-r --ignore-existing'
else:
    os.makedirs(exec_run_dir)
    
print('Preparing CTDAS run in ' + exec_run_dir)

cmd_ = "rsync %s --cvs-exclude --exclude=*nc --exclude=*jb --exclude=*out --exclude=*log {src} {run}/" %rsyncopt
cmd  = cmd_.format(src=src_dir, run=exec_run_dir)
p    = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
p.wait()

cmd_ = "rsync %s --cvs-exclude {src}/analysis/cteco2/*nc {run}/da/analysis/cteco2/" %rsyncopt
cmd  = cmd_.format(src=src_dir, run=exec_run_dir)
p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
p.wait()


##################################
##### Write run and rc files #####
##################################

# - Copy template py, rc and jb files
# - Adjust keys in rc files that are present in settings
# - This is similar to the functionality in clone_ctdas.sh, but allows me to
#   keep the settings files independent of the development version of CTDAS.
# - This could mean some duplicated work if I ever adjust rc files, but that's
#   worth keeping the separation.

ctdas_prefix = os.path.basename(s["ctdas_run_dir"])


########## CTDAS rc file

# Read template
ctdas_rc = rc.RcFile(s["ctdas_template_rc"])

# Modify for this run
ctdas_rc.replace("time.cycle",  s["time.cycle"])
ctdas_rc.replace("time.nlag",   s["time.nlag"])
ctdas_rc.replace("time.start",  s["period"][0])
ctdas_rc.replace("time.finish", s["period"][0] + dtm.timedelta(days = s["time.cycle"]* s["ctdas.ncycles"]))
ctdas_rc.replace("dir.da_run",  s["ctdas_run_dir"])
ctdas_rc.replace("da.resources.ntime",    s["ctdas.ntime"])
ctdas_rc.replace("da.resources.ntasks",   s["ctdas.ntasks"])
ctdas_rc.replace("da.optimizer.nmembers", s["nmembers_ens"])
ctdas_rc.replace("da.obsoperator.rc",     "da/rc/wrfchem/obsoper_wrfchem_%s.rc"%s["run_name"])
ctdas_rc.replace("da.system.rc",          "da/rc/wrfchem/dasystem_wrfchem_%s.rc"%s["run_name"])

# Write modified file
fp_ctdas_rc = os.path.join(exec_run_dir, ctdas_prefix+".rc")
ctdas_rc.WriteFile(fp_ctdas_rc)
del ctdas_rc


########## template.jb

# Read in the file
with open(s["ctdas_template_jb"], 'r') as file :
  jbfile = file.read()

# Replace the templace placeholder by the runname
jbfile = jbfile.replace('template', ctdas_prefix)

# Write the file out again
with open(os.path.join(exec_run_dir, ctdas_prefix+".jb"), 'w') as file:
  file.write(jbfile)


########## dasystem.rc

# Read Template
dasystem_rc = rc.RcFile(s["dasystem_template_rc"])

# path to and naming template of observations
dasystem_rc.replace("obs.column.input.dir", s["obs_path"])
dasystem_rc.replace("obs.column.ncfile",    s["obs_pattern"])

# xco2_weights.rc: file currently only contains rejection threshold for observations
if "xco2_weights_rc" in s:
    dasystem_rc.replace("obs.column.rc", s["xco2_weights_rc"])
else:
    dasystem_rc.replace("obs.column.rc", s["obs_path"] + '/xco2_weights.rc')

# Test: retrieve unperturbed truth (scaling factors should be 1)
#s["obs_pri_truth_pattern"] = s["obs_pri_truth_pattern"].replace("truth5.0", "truth5.3")

# LF: EASIER TO JUST MAKE A DIFFERENT DIRECTORY WITH ARTIFICIAL OBSERVATIONS, AND SELECT WITH OBS_PATH AND OBS_PATTERN
#if s["atmos_data"] == "artificial":
#    dasystem_rc.replace("obs.column.ncfile", s["obs_pri_truth_pattern"])
#elif s["atmos_data"] == "real":
#    dasystem_rc.replace("obs.column.ncfile", s["obs_pattern"])
#else:
#    raise ValueError('Unknown value for s["atmos_data"]')

# Statevector definition (= regionsfile)
dasystem_rc.replace("regionsfile", s["regionsfile"])

# Get number of parameters in statevector from regionsfile
ncf = nc.Dataset(s["regionsfile"])
#nparams_flux = len(ncf.dimensions["parameters"])
nparams_flux = len(ncf.dimensions["parameters"])*s["n_emis_proc"] # overlaying processes
ncf.close()

# If we optimize offset in initial/boundary conditions: add one parameter to state vector
do_offset      = "do_offset" in s.keys() and s["do_offset"]
nparams_offset = int(do_offset)
dasystem_rc.replace("nparameters", nparams_flux + nparams_offset)
if do_offset:
    dasystem_rc.add("sigma_offset", s["sigma_offset"])

# Overlaying parameters
# LF: WHAT IS SIGMA_SCALE?
dasystem_rc.replace("n_emis_proc", s["n_emis_proc"])
for nproc in range(1, s["n_emis_proc"]+1):
    dasystem_rc.add("sigma_scale_%d"%nproc, s["sigma_scale_list"][nproc-1])

# Write modified file
fp_dasystem_rc = os.path.join(exec_run_dir, rc.RcFile(fp_ctdas_rc).get("da.system.rc"))
dasystem_rc.WriteFile(fp_dasystem_rc)
del dasystem_rc


########## wrfchem.rc

# Read template
wrfchem_rc = rc.RcFile(s["wrfchem_template_rc"])

# Modify paths for ens_src run
#wrfchem_rc.replace("source_dir", s["wrf_ens_src_dir"])
#wrfchem_rc.replace("run_dir",    s["wrf_ens_run_dir"])
#wrfchem_rc.replace("obs.column.footprint_samples_dim", s["footprint_samples_dim"])

# Modify paths for fwd_src run
wrfchem_rc.replace("source_dir", s["wrf_fwd_src_dir"])
wrfchem_rc.replace("run_dir",    s["wrf_fwd_run_dir"])
wrfchem_rc.replace("obs.column.footprint_samples_dim", s["footprint_samples_dim"])

# Write modified file
fp_wrfchem_rc = os.path.join(exec_run_dir, rc.RcFile(fp_ctdas_rc).get("da.obsoperator.rc"))
wrfchem_rc.WriteFile(fp_wrfchem_rc)
del wrfchem_rc
print(fp_wrfchem_rc)


########## ctdas.py

# Copy template to destination
fn_ctdas = ctdas_prefix + ".py" 
fp_ctdas = os.path.join(exec_run_dir, fn_ctdas)
shutil.copy2(s["ctdas_template_py"], fp_ctdas)


###########################################
##### Run WPS and WRF in fwd mode #####
###########################################

print('Running WPS and WRF in fwd mode...')
logfile = 'run_wps_wrf_fwd_withrunner.%s.log' %s['run_name']
r = subprocess.run("python %s %s fwd >&%s" %(s['wps_wrf_runner'],args.settings_file,logfile),
                   shell=True,
                   stdout=subprocess.PIPE,
                   stderr=subprocess.STDOUT)
r.check_returncode()
print('Finished execution of run_wps_wrf_fwd_withrunner.py')

ran_real = os.path.exists(os.path.join(s["wrf_fwd_run_dir"], "wrfbdy_d01"))
if ran_real:
    print('WRF boundary fily present, assuming WRF directory (%s) ready for run' %s["wrf_fwd_run_dir"])
else:
    raise IOError("Ensemble source not found. Run run_wps_wrf... in fwd mode.")
    

###########################################
##### Run WPS and WRF in ens_src mode #####
###########################################

#print('Running WPS and WRF in ens_src mode...')
#logfile = 'run_wps_wrf_fwd_withrunner.%s.log' %s['run_name']
#r = subprocess.run("python %s %s ens_src >&%s" %(s['wps_wrf_runner'],args.settings_file,logfile),
#                   shell=True,
#                   stdout=subprocess.PIPE,
#                   stderr=subprocess.STDOUT)
#r.check_returncode()
#print('Finished execution of run_wps_wrf_fwd_withrunner.py')

#ran_real = os.path.exists(os.path.join(s["wrf_ens_src_dir"], "wrfbdy_d01"))
#if ran_real:
#    print('WRF boundary fily present, assuming WRF directory (%s) ready for run' %s["wrf_ens_src_dir"])
#else:
#    raise IOError("Ensemble source not found. Run run_wps_wrf... in ens_src mode.")


############################
##### Submit CTDAS job #####
############################

cwd = os.getcwd()
os.chdir(exec_run_dir)
#cmd = os.path.join(HOME, "submit_ctdas.sh") + " " + fn_ctdas
cmd = 'sbatch %s.jb' %s['run_name']
print('Submit command: ' + cmd)
p = subprocess.Popen(cmd, shell=True,stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
os.chdir(cwd)
o, e = p.communicate()
print(o)
print(e)
