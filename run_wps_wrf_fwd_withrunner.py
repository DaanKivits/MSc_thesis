#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  9 10:13:40 2019
@author: friedemann, modifications by liesbeth

- Runs WPS and WRF by creating fresh run directories for either and external namelists
- based on settings files, e.g. produced by create_settings_file.py

Usage:
"python create_settings_file.py" -> settings_file.pkl
"python run_wps_wrf_fwd_withrunner.py settings_file.pkl"

Note: Re-running without overwriting is not guaranteed to work. For example, WPS programs look
if *any* of their output exists to decide whether or not to run, so if some but
not all files were produced in an attempt and overwrite is False, they have to
be deleted manually first. Also, I'm not sure all namelist adjustments go through
(i.e. mod_levs)
    
Note:
Running as above on login-node doesn't work for 2-month-che runs, because
they take longer than one hour. So I fiddled with submitting this to the cluster
with .exe, and I kept getting a bunch of MPI errors, among them:
    MPID_Init(1276)......: channel initialization failed
Until I added "srun -n1 -N1 --mem_opt=XXXXM" to the subprocess.
This most likely needs to be done for the WPS calls as well. Submitting WRF
worked without, so I implemented this ("submit_cmd") for all subprocess calls
to WPS executables, but not bash scripts (link_grib, cat, submit_exe). Perhaps
the difference is whether an executable was compiled for parallel usage. Then
real.exe would be the only one where this is necessary.
Didn't test it yet.

Slurm usage:

wrf_runner.jb:
#! /bin/bash
#SBATCH -p thin
#SBATCH -t 4:00:00
#SBATCH -n 24
#SBATCH -N 1
#SBATCH --mail-user=<email>
#SBATCH --mail-type=FAIL,END

echo "All output piped to file wrf_runner.log"
module load 2021
module load Python/3.9.5-GCCcore-10.3.0

python ./wrf_runner.py  >& wrf_runner.log   

sbatch --job-name=wrf_runner -o log/slurm-%j.out wrf_runner.jb



"""

# Run CHE WPS for times defined here ("cases") and domains/configurations
# defined in namelists.

# Namelist location: namelist_dir
# Namelist names: "namelist." + cases.keys()[n] + ".wps"

# Before processing, it prepares the WPS directory. Therefore, several paths
# need to be defined here.

import os
import copy
import subprocess
import shutil
import re
import glob
import argparse
import pickle
import datetime as dtm
import pandas as pd
import f90nml


###################
##### OPTIONS #####
###################

nml_time_format = "%Y-%m-%d_%H:%M:%S"

if os.getenv("SLURM_JOB_ID") is not None:
    running_as_cluster_job = True
else:
    running_as_cluster_job = False

if running_as_cluster_job:
    print("Warning: Didn't test if subprocess launching works. See docstring of this file.")


########################################
##### PARSE COMMAND LINE ARGUMENTS #####
########################################

parser = argparse.ArgumentParser()
parser.add_argument(dest="settings_file", help="Path to settings file (pickled dictionary)")
parser.add_argument(dest="mode", choices=['fwd','ens_src'], help="Specify run mode: fwd or ens_src")

#if running_as_cluster_job:
#    parser.add_argument(dest="mem_MB", help="Memory in MB available for job step creation.")

args = parser.parse_args()
with open(args.settings_file, "rb") as f:
    s = pickle.load(f)

# fwd does a forward run, ens_src runs real for ctdas ens source dir (slightly different flux files and few other options)
mode = args.mode
print('Running script in %s mode' %mode.upper())


#######################
##### WPS OPTIONS #####
#######################

do_wps = True

if mode=="ens_src":
    #do_wps = False # Yah never do that
    pass # actually have to for scarbo

if mode=="fwd":
    #do_wps = False # Yah never do that
    pass # actually have to for scarbo


# If do_wps = False, the remaining WPS settings have no effect
overwrite_wps = True

do_ungrib  = True
do_geogrid = True
do_metgrid = True

concat_intermediate = True # For WPS with ERA5 always use True.


####################### 
##### WRF OPTIONS #####
#######################

do_wrf = True

if mode=="ens_src":
    do_wrf = True

if mode=="fwd":
    do_wrf = True

# If do_wrf = False, the remaining WRF settings have no effect
overwrite_real    = True
overwrite_tracers = True
overwrite_wrf     = True
submit_wrf        = True

if mode=="ens_src":
    submit_wrf = False

if mode =="fwd":
    submit_wrf = True

# A command that can be run by subprocess.run on this platform
# for WRF/, 5 days:
wrf_submit_cmd = "/projects/0/ctdas/dkivits/scripts/submit_exe.sh -p thin -t 10:00:00 -n 1 -N 1"
# for WRF/, 1-2 months:
#wrf_submit_cmd = "~/submit_exe.sh -p thin -t 80:00:00 -n 72"
# for WRF/, 5 days months:
#wrf_submit_cmd = "~/submit_exe.sh -t 10:00:00 -n 32"
# for WRF_n150/:
#wrf_submit_cmd = "~/submit_exe.sh -p thin -t 5:00:00 -n 72"
#wrf_submit_cmd = "/projects/0/ctdas/dkivits/scripts/submit_exe.sh -n 1"


######################################
########## HELPER FUNCTIONS ##########
######################################

def clone_dir(src_dir, dst_dir, copy_patterns=[], skip_patterns=[],
        #link_patterns=[],
              overwrite=False):
    """
    Clones a directory. Default: Link all files. copy_patterns and skip_patterns
    can be used to copy or skip specific files instead.
    
    Arguments
    =========
    src_dir: string
    dst_dir: string
    copy_patterns: list of file patterns recognizable by re.match
    skip_patterns: list of file patterns recognizable by re.match
    #link_patterns: None or list of file patterns recognizable by re.match
    overwrite: Checks each file individually.
    """
    
    # Make run directory/check if is empty    
    os.makedirs(dst_dir, exist_ok=True)
                
    # Create one single pattern to check
    copy_full_pattern = "(" + "|".join(copy_patterns) + ")"
    skip_full_pattern = "(" + "|".join(skip_patterns) + ")"
    
    # Get a list of all files to copy/link
    source_files = os.listdir(src_dir)
    
    # Copy/link files to target_dir
    for filename in source_files:
        fp_source = os.path.join(src_dir, filename)
        fp_destination = os.path.join(dst_dir, filename)

        # Do nothing for files to skip
        if re.match(skip_full_pattern, filename):
            continue

        # If exists: remove if overwrite, next file if not
        if os.path.exists(fp_destination):
            if overwrite:
                try:
                    os.remove(fp_destination)
                except IsADirectoryError:
                    # If the error was caused because the destination was a directory
                    shutil.rmtree(fp_destination)
            else:
                continue
            
        # Copy files to copy
        if re.match(copy_full_pattern, filename):
            try:
                shutil.copy2(fp_source, fp_destination)
            except IsADirectoryError:
                # If the error was caused because the source was a directory
                shutil.copytree(fp_source, fp_destination)
            
            # logging.debug("Copied file %s", filename)
            continue

        # Link all other files
        os.symlink(fp_source, fp_destination)
        
    
#####################
##### EXECUTION #####
#####################  
if running_as_cluster_job:
    #submit_cmd = "srun -n1 -N1 --mem=" + str(args.mem_MB) + "M "
    submit_cmd = "srun -n1 -N1 "
else:
    submit_cmd = ""

## Check all paths if they are absolute
#for path in [namelist_dir, wps_dir, wps_out_dir, geog_data_path, grib_dir,
#             fp_vtable, run_src_dir, run_base_dir]:
#    if not os.path.isabs(path):
#        raise OSError("Path " + path + " not absolute.")

# Suffix for log files
log_suffix = "." + s["run_name"] + ".log"
  
if mode=="fwd":
    wrf_src_dir = s["wrf_fwd_src_dir"]
    wrf_run_dir = s["wrf_fwd_run_dir"]
    wrfchemi_dir_key = "wrfchemi_fwd_dir"
    nmembers = s["nmembers_fwd"]

elif mode=="ens_src":
    wrf_src_dir = s["wrf_ens_src_src_dir"] # This is to enable multiple runs at once
    wrf_run_dir = s["wrf_ens_src_dir"]
    wrfchemi_dir_key = "wrfchemi_ens_dir"
    nmembers = s["nmembers_ens"]
    if submit_wrf:
        print("Resetting submit_wrf to False")
        submit_wrf = False

if ("run_wps" in s.keys()) and (s["run_wps"] == do_wps):
    raise ValueError("Settings file run_wps and do_wps don't match.")


#############################
########## RUN WPS ##########
#############################

if not do_wps or (glob.glob(os.path.join(s["wps_run_dir"], "met_em.*")) and not overwrite_wps):
    print("WPS not selected or already ran and overwrite_wps is false, skipping WPS.")

else:    
    # Make directory for WPS output
    print("Making WPS run directory in " + s["wps_run_dir"] + "...")
    
    # Files to skip
    wps_skip_patterns = ["^ungrib&",
                         "^geogrid&",
                         "^metgrib&",
                         "^util&",
                         "^arch&",
                         "^namelist.wps",
                         "^FILE:.*",
                         "^PRES:.*",
                         "^nPRES:.*",
                         "^SFC:.*",
                         "^ML:.*",
                         r"^geo_em\.d??\.nc",
                         r"^met_em\.d??\.nc",
                         r".*\.txt"]
    
    # Files to copy
    wps_copy_patterns = []
    clone_dir(s["wps_src_dir"], s["wps_run_dir"], wps_copy_patterns, wps_skip_patterns, overwrite=overwrite_wps)    
    
    # Link Vtable
    cwd = os.getcwd()
    os.chdir(s["wps_run_dir"])
    
    if not os.path.exists("Vtable") or overwrite_wps:
        try:
            os.remove("Vtable")
        except FileNotFoundError:
            pass
        os.symlink(s["fp_vtable"], "Vtable")
    os.chdir(cwd)
    
    # Link ecmwf_coeffs
    if s["driving_meteo"] in ["ERA-Interim", "ERA5"]:
        cwd = os.getcwd()
        os.chdir(s["wps_run_dir"])
        
        if not os.path.exists("ecmwf_coeffs") or overwrite_wps:
            try:
                os.remove("ecmwf_coeffs")
            except FileNotFoundError:
                pass
            os.symlink(s["fp_ecmwf_coeffs"], "ecmwf_coeffs")
        os.chdir(cwd)
    else:
        raise ValueError("Unknown driving meteo: %s" % s["driving_meteo"])
    
    # Read namelist template
    nml_wps_src = f90nml.read(s["fp_nml_wps_src"])
        
    # Fill namelist template
    print("Making namelist...")
    nml = copy.deepcopy(nml_wps_src)

    # Set dates
    max_dom_wps = int(nml_wps_src["share"]["max_dom"])
    nml["share"]["start_date"] = [s["period"][0].strftime(nml_time_format)]*max_dom_wps
    nml["share"]["end_date"] = [s["period"][1].strftime(nml_time_format)]*max_dom_wps

    # Set directories
    nml["share"]["opt_output_from_geogrid_path"] = s["wps_run_dir"]
    nml["geogrid"]["geog_data_path"] = s["geog_data_path"]
    
    if s["driving_meteo"] == "ERA-Interim":
        fg_name_metgrid = [s["wps_run_dir"] + "/" + x for x in ["FILE", "PRES"]]
    if s["driving_meteo"] == "ERA5":
        #fg_name_metgrid = [s["wps_run_dir"] + "/" + x for x in ["FILE", "nPRES"]]
        fg_name_metgrid = [s["wps_run_dir"] + "/" + x for x in ["SFC", "ML", "nPRES"]]
        
    nml["metgrid"]["fg_name"] = fg_name_metgrid 
    nml["metgrid"]["opt_output_from_metgrid_path"] = s["wps_run_dir"]

    # Write namelist
    fn_nml_dst = "namelist." + s["run_name"] + ".wps"
    fp_nml_dst = os.path.join(s["wps_run_dir"], fn_nml_dst)
    
    if not os.path.exists(fp_nml_dst) or overwrite_wps:
        f90nml.write(nml, fp_nml_dst, force=True)

        # Link namelist to namelist.wps
        fp_nml_run = os.path.join(s["wps_run_dir"], "namelist.wps")
        try:
            os.remove(fp_nml_run)
        except FileNotFoundError:
            pass
        os.symlink(fp_nml_dst, fp_nml_run)

    # Switch to wps_dir
    cwd = os.getcwd()
    os.chdir(s["wps_run_dir"])

    # Run WPS: ungrib, geogrid, metgrid, and any intermediate steps
    if do_ungrib:
        # Link driving meteo (only required time steps) and run ungrib
        # Different for different driving meteos
        if s["driving_meteo"] == "ERA-Interim":
            
            # Speficy pattern for meteo source files to link
            if s["grib_time_format"] == "%Y%m":
                met_times   = pd.date_range(start=s["period"][0], end=s["period"][1], freq="MS")         
                time_stamps = met_times.strftime(s["grib_time_format"]).tolist()
            else:
                raise ValueError("Don't know met_time_format '" + s["grib_time_format"] + "'")
    
            if not glob.glob("FILE*") or overwrite_wps:
                print("Running link_grib...")
                patterns = [os.path.join(s["grib_dir"], "*" + x + "*") for x in time_stamps]
                pattern  = " ".join(patterns)
        
                r = subprocess.run("./link_grib.csh " + pattern,
                                   shell=True,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT)
                r.check_returncode()
                                
                print("Adjusting namelist for ungrib (all data)...")
                nml = f90nml.read(fp_nml_dst)
                nml["ungrib"]["prefix"] = s["wps_run_dir"] + "/FILE"
                f90nml.write(nml, fp_nml_dst, force=True)
                
                print("Running ungrib...")
                log_file = os.path.join(s["wps_run_dir"], "ungrib" + log_suffix)
                r = subprocess.run(submit_cmd + "./ungrib.exe >& " + log_file,
                                   shell=True,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT)
                r.check_returncode()
                # There's a second log file in the run directory... move it to the output directory.
                fp_ungrib_log = os.path.join(s["wps_run_dir"], "ungrib" + log_suffix + ".execdir")
                shutil.move(src="ungrib.log", dst=fp_ungrib_log)
            
        if s["driving_meteo"] == "ERA5":
            
            # Speficy pattern for meteo source files to link
            if s["grib_time_format"] == "%Y%m":
                met_times   = pd.date_range(start=s["period"][0], end=s["period"][1], freq="MS")
                time_stamps = met_times.strftime(s["grib_time_format"]).tolist()
            
            elif s["grib_time_format"] == "%Y%m%d":
                met_times   = pd.date_range(start=s["period"][0], end=s["period"][1], freq="D")
                time_stamps = met_times.strftime(s["grib_time_format"]).tolist()
            
            else:
                raise ValueError("Don't know met_time_format '" + s["grib_time_format"] + "'")
    
            
            if not glob.glob("SFC*") or overwrite_wps:
                print("Running link_grib for surface data...")
                patterns_sfc = [os.path.join(s["grib_dir"], "ERA5_" + x + "*_sfc.grb1") for x in time_stamps]
                pattern_sfc  = " ".join(patterns_sfc)
                r = subprocess.run("./link_grib.csh " + pattern_sfc,
                                   shell=True,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT)
                r.check_returncode()
                
                print("Adjusting namelist for ungrib of surface data...")
                nml = f90nml.read(fp_nml_dst)
                nml["ungrib"]["prefix"] = s["wps_run_dir"] + "/SFC"
                f90nml.write(nml, fp_nml_dst, force=True)
                
                print("Running ungrib for surface data...")
                log_file = os.path.join(s["wps_run_dir"], "ungrib_sfc" + log_suffix)
                r = subprocess.run(submit_cmd + "./ungrib.exe >& " + log_file,
                                   shell=True,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT)
                r.check_returncode()
                # There's a second log file in the run directory... move it to the output directory.
                fp_ungrib_sfc_log = os.path.join(s["wps_run_dir"], "ungrib_sfc" + log_suffix + ".execdir")
                shutil.move(src="ungrib.log", dst=fp_ungrib_sfc_log)
    
            if not glob.glob("ML*")  or overwrite_wps:
                print("Running link_grib for model level data...")
                patterns_ml = [os.path.join(s["grib_dir"], "ERA5_" + x + "*_ml.grb2") for x in time_stamps]
                pattern_ml  = " ".join(patterns_ml)
                r = subprocess.run("./link_grib.csh " + pattern_ml,
                                   shell=True,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT)
                r.check_returncode()
                
                print("Adjusting namelist for ungrib of model level data...")
                nml = f90nml.read(fp_nml_dst)
                nml["ungrib"]["prefix"] = s["wps_run_dir"] + "/ML"
                f90nml.write(nml, fp_nml_dst, force=True)
                
                print("Running ungrib for model level data...")
                log_file = os.path.join(s["wps_run_dir"], "ungrib_ml" + log_suffix)
                r = subprocess.run(submit_cmd + "./ungrib.exe >& " + log_file,
                                   shell=True,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT)
                r.check_returncode()
                # There's a second log file in the run directory... move it to the output directory.
                fp_ungrib_ml_log = os.path.join(s["wps_run_dir"], "ungrib_ml" + log_suffix + ".execdir")
                shutil.move(src="ungrib.log", dst=fp_ungrib_ml_log)
                
            if concat_intermediate and (not glob.glob("FILE*") or overwrite_wps):
                print("Concatenating intermediate model level and surface files...")
                
                # Reference: https://github.com/wrf-model/WPS/issues/113.
                # im = "intermediate"
                time_stamps_im_d = met_times.strftime("%Y-%m-%d").tolist()
                for time_stamp_im_d in time_stamps_im_d:
                    ml_files = glob.glob(os.path.join(s["wps_run_dir"], "ML:" + time_stamp_im_d + "*"))
                    ml_files.sort()
                    time_stamps_im = [os.path.basename(x)[3:] for x in ml_files]
                   
                    for time_stamp_im in time_stamps_im:
                        ml_file  = os.path.join(s["wps_run_dir"], "ML:"+time_stamp_im)
                        sfc_file = os.path.join(s["wps_run_dir"], "SFC:"+time_stamp_im)
                        cat_file = os.path.join(s["wps_run_dir"], "FILE:"+time_stamp_im)
                        r = subprocess.run("cat " + sfc_file + " " + ml_file + " > " + cat_file, shell=True)
                        r.check_returncode()
                        
        if not glob.glob("PRES*") or overwrite_wps:
            print("Adjusting namelist for calc_ecmwf_p...")
            nml = f90nml.read(fp_nml_dst)
            if s["driving_meteo"]=="ERA-Interim" or concat_intermediate:
                nml["metgrid"]["fg_name"] = s["wps_run_dir"] + "/FILE"
            else:
                nml["metgrid"]["fg_name"] = [s["wps_run_dir"] + x for x in ["/SFC", "/ML"]]
            f90nml.write(nml, fp_nml_dst, force=True)
        
            print("Running calc_ecmwf_p...")
            log_file = os.path.join(s["wps_run_dir"], "calc_ecmwf_p" + log_suffix)
            r = subprocess.run(submit_cmd + "./util/calc_ecmwf_p.exe >& " + log_file,
                               shell=True,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT)
            r.check_returncode()
            fp_cep_log = os.path.join(s["wps_run_dir"], "logfile" + log_suffix + ".execdir")
            shutil.move(src="logfile.log", dst=fp_cep_log)
        
        if s["driving_meteo"] == "ERA5":
            
            if (not glob.glob("nPRES*") or overwrite_wps) and (s["meteo_start_level"] == 47):
                # This step is necessary if only levels 47 to 137 of ERA5 are downloaded
                print("Running mod_levs...")
                time_stamps_im_d = met_times.strftime("%Y-%m-%d").tolist()
                for time_stamp_im_d in time_stamps_im_d:
                    ml_files = glob.glob(os.path.join(s["wps_run_dir"], "ML:" + time_stamp_im_d + "*"))
                    ml_files.sort()
                    time_stamps_im = [os.path.basename(x)[3:] for x in ml_files]                    
                    
                    for time_stamp_im in time_stamps_im:
                        in_file  = os.path.join(s["wps_run_dir"], "PRES:"+time_stamp_im)
                        out_file = os.path.join(s["wps_run_dir"], "nPRES:"+time_stamp_im)
                        log_file = "mod_levs_" + time_stamp_im + log_suffix
                        r = subprocess.run(submit_cmd + "./util/mod_levs.exe " + in_file + " " + out_file + " >& " + log_file,
                                           shell=True,
                                           stdout=subprocess.PIPE,
                                           stderr=subprocess.STDOUT)
                        r.check_returncode()
    
    ########## geogrid
    if do_geogrid:
        if not glob.glob(os.path.join(s["wps_run_dir"], "geo_em.*")) or overwrite_wps:
            print("Running geogrid...")
            log_file = os.path.join(s["wps_run_dir"], "geogrid" + log_suffix)
            r = subprocess.run(submit_cmd + "./geogrid.exe >& " + log_file,
                               shell=True,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT)
            r.check_returncode()
            fp_geogrid_log = os.path.join(s["wps_run_dir"], "geogrid" + log_suffix + ".execdir")
            shutil.move(src="geogrid.log", dst=fp_geogrid_log)

    ########## metgrid
    if do_metgrid:
        print("Adjusting namelist for metgrid...")
        nml = f90nml.read(fp_nml_dst)
        nml["metgrid"]["fg_name"] = fg_name_metgrid 
        f90nml.write(nml, fp_nml_dst, force=True)
    
        print("Running metgrid...")
        log_file = "metgrid" + log_suffix
        r = subprocess.run(submit_cmd + "./metgrid.exe >& " + log_file,
                           shell=True,
                           stdout=subprocess.PIPE,
                           stderr=subprocess.STDOUT)
        r.check_returncode()
        fp_metgrid_log = os.path.join(s["wps_run_dir"], "metgrid" + log_suffix + ".execdir")
        shutil.move(src="metgrid.log", dst=fp_metgrid_log)

    # Go back to previous directory.
    os.chdir(cwd)
    print("WPS done.")


################################
########## END OF WPS ##########
################################

##################################
########## START OF WRF ##########
##################################

if not do_wrf or (glob.glob(os.path.join(wrf_run_dir, "wrfout.*")) and not (overwrite_real or overwrite_wrf)):
    print("WRF not selected or already ran and overwrite_real and overwrite_wrf are false, skipping WRF.")

else:
    print("Preparing WRF run in %s" %wrf_run_dir)
    
    # Read namelist template
    nml_input_src = f90nml.read(s["fp_nml_input_src"])

    #################################################
    ########## Check namelist requirements ##########
    #################################################
    
    print("Checking required namelist items...")
    ##iofields
    #if not "iofields_filename" in nml_input_src["time_control"].keys():
    #    raise ValueError("iofields_filename not found.")
    # daily restarts
    if not "restart_interval" in nml_input_src["time_control"].keys():
        raise ValueError("restart_interval not found.")
    
    ########################################
    ########## Make fresh run dir ##########
    ########################################
    print("Making WRF run directory in " + wrf_run_dir + "...")
    
    # The default behaviour is to link everything from the source run directory.
    # However, the files that are modified later need to be copied, and previous
    # output (if any) skipped. These cases are defined here. This should only be
    # changed if the source run directory has files other than those in a fresh WRF/run/.
    
    wrf_copy_patterns = []
    
    # Skip some files
    # wrfchemi-files come from a seperate directory specified in the settings file
    wrf_skip_patterns = ["^wrfchemi_.*",
                         r"^wrfinput_track.*\.txt$",
                         "^wrfout.*",
                         "^wrfrst.*",
                         "wrfinput_d.*",
                         "^wrfbdy_d.*",
                         "^wrffdda_d.*",
                         r"^namelist\.input$",
                         r"^namelist\.output",
                         r"^met_em\..*",
                         ".*~$",
                         "log",
                         "^rsl.error.*",
                         "^rsl.out.*",
                         "sampled_.*",
                         "^restart.*",
                         "^rsl_real$",
                         "^case_selection$"
                         ]
    
    if "fp_iofields" in s.keys():
        fn_iofields = os.path.basename(s["fp_iofields"])
        wrf_skip_patterns.append(fn_iofields)

    clone_dir(wrf_src_dir, wrf_run_dir, wrf_copy_patterns, wrf_skip_patterns, overwrite=overwrite_real or overwrite_wrf)
    
    # Create a log dir
    os.makedirs(os.path.join(wrf_run_dir, "log"), exist_ok=True)
    
    # Change to run directory
    cwd = os.getcwd()
    os.chdir(wrf_run_dir)
    

    ###########################################
    ########## UPDATE namelist.input ##########
    ###########################################
    
    if os.path.exists("namelist.input") and (overwrite_real or overwrite_wrf):
        os.remove("namelist.input")

    print("Making namelist...")
    nml_input = copy.deepcopy(nml_input_src)

    # Set dates in namelist
    max_dom_input = int(nml_input["domains"]["max_dom"])

    time_components = ["year", "month", "day", "hour", "minute", "second"]
    for time_component in time_components:
        nml_input["time_control"]["start_" + time_component] = \
            [getattr(s["period"][0], time_component)]*max_dom_input
        nml_input["time_control"]["end_" + time_component] = \
            [getattr(s["period"][1], time_component)]*max_dom_input

    # Set iofields...
    if "fp_iofields" in s.keys():
        fn_iofields = os.path.basename(s["fp_iofields"])
        nml_input["time_control"]["iofields_filename"] = [fn_iofields]*max_dom_input
    
    # Set chemistry options - already done in template, but for true restart runs it is needed here.
    nml_input["time_control"]["io_form_auxinput5"]    = 2
    nml_input["time_control"]["auxinput5_interval_m"] = [60]*max_dom_input
    nml_input["time_control"]["frames_per_auxinput5"] = [1]*max_dom_input
    nml_input["time_control"]["auxinput5_inname"]     = 'wrfchemi_d<domain>_<date>'
    
    if "vprm_input_dir" in s.keys():
        nml_input["time_control"]["io_form_auxinput15"]       = 2
        nml_input["time_control"]["auxinput15_inname"]        = 'vprm_input_d<domain>_<date>'
        nml_input["time_control"]["auxinput15_interval_m"]    = [1440]*max_dom_input
        nml_input["time_control"]["frames_per_auxinput15"]    = [1]*max_dom_input
        nml_input["chem"]["bioemdt"]         = [30]*max_dom_input # minutes
        nml_input["chem"]["bio_emiss_opt"]   = [17]*max_dom_input
        # not sure if this is necessary:
        nml_input["chem"]["conv_tr_wetscav"] = [0]*max_dom_input 
        
    nml_input["chem"]["chem_opt"]           = [17]*max_dom_input
    nml_input["chem"]["io_style_emissions"] = 2
    nml_input["chem"]["emiss_opt"]          = [17]*max_dom_input

    if "emis_zdim" in s.keys():
        nml_input["chem"]["kemit"] = s["emis_zdim"]
        
    # Check that timestep divides auxinput5_interval (assume no fractional time steps and same interval for all domains)
    intrvl = nml_input["time_control"]["auxinput5_interval_m"]
    if isinstance(intrvl,list):
        intrvl = intrvl[0]
    
    if intrvl*60 % nml_input["domains"]["time_step"] != 0:
        raise ValueError("time_step does not divide auxinput5_interval_m. Reading emissions won't work.")
        
    # Write namelist
    fn_nml_input_dst = "namelist." + s["run_name"] + ".input"
    fp_nml_input_dst = os.path.join(wrf_run_dir, fn_nml_input_dst)
    if not os.path.exists(fp_nml_input_dst) or overwrite_real or overwrite_wrf:
        f90nml.write(nml_input, fp_nml_input_dst, force=True)
        if not os.path.exists("namelist.input") or overwrite_real or overwrite_wrf:
            os.symlink(fn_nml_input_dst, "namelist.input")
        

    #################################
    ########## LINK MET_EM ##########
    #################################
    if not glob.glob(os.path.join(wrf_run_dir, "met_em.*")) or overwrite_real:
        print("Linking met_em.*...")
        source_files = glob.glob(os.path.join(s["wps_run_dir"], "met_em.*"))
        for file in source_files:
            fp_dst = os.path.join(wrf_run_dir, os.path.basename(file))
            if os.path.exists(fp_dst) and overwrite_real:
                os.remove(fp_dst)
            os.symlink(file, fp_dst)
    

    #######################################
    ########## COPY IOFIELDS.TXT ##########
    #######################################
    if "fp_iofields" in s.keys():
        fn_iofields_dst = os.path.basename(s["fp_iofields"])
        fp_iofields_dst = os.path.join(wrf_run_dir, fn_iofields_dst)
        if not os.path.exists(fp_iofields_dst) or overwrite_real or overwrite_wrf:
            print("Copying iofields file " + s["fp_iofields"] + "...")           
            shutil.copy2(s["fp_iofields"], fp_iofields_dst)

    
    ###################################
    ########## LINK WRFCHEMI ##########
    ###################################
    if wrfchemi_dir_key in s.keys():
        fps_source = glob.glob(os.path.join(s[wrfchemi_dir_key], s["wrfchemi_pattern"]))
        if len(fps_source) == 0:
            print("No wrfchemi files found.")
        else:
            print("Linking wrfchemi files from " + s[wrfchemi_dir_key] + "...")
            
        for fp_src in fps_source:
            fn = os.path.basename(fp_src)
            # Check if the file falls in the simulation period
            tstamp = fn[-len("xxxx-xx-xx_xx:xx:xx"):]
            tfile = dtm.datetime.strptime(tstamp, "%Y-%m-%d_%H:%M:%S")
            if s["period"][0] <= tfile and tfile <= s["period"][1]:
                fp_dst = os.path.join(wrf_run_dir, fn)
                if os.path.exists(fp_dst) and overwrite_wrf:
                    os.remove(fp_dst)
                if not os.path.exists(fp_dst):
                    os.symlink(fp_src, fp_dst)


    #####################################
    ########## LINK vprm input ##########
    #####################################
    if "vprm_input_dir" in s.keys():
        fps_source = glob.glob(os.path.join(s["vprm_input_dir"], s["vprm_input_pattern"]))
        if len(fps_source)==0:        
            print("VPRM input not found.")
        else:
            print("Linking vprm input...")
        for fp_src in fps_source:
            fn = os.path.basename(fp_src)
            # Check if the file falls in the simulation period
            tstamp = fn[-len("xxxx-xx-xx_xx:xx:xx"):]
            tfile = dtm.datetime.strptime(tstamp, "%Y-%m-%d_%H:%M:%S")
            if s["period"][0] <= tfile and tfile <= s["period"][1]:
                fp_dst = os.path.join(wrf_run_dir, fn)
                if os.path.exists(fp_dst) and overwrite_wrf:
                    os.remove(fp_dst)
                os.symlink(fp_src, fp_dst)



    ##############################
    ########## RUN REAL ##########
    ##############################
    if not (os.path.exists(os.path.join(wrf_run_dir, "wrfbdy_d01")) or \
            os.path.exists(os.path.join(wrf_run_dir, "wrfinput_d01"))) or \
            overwrite_real:
        print("Running real...")
        log_file =  "real" + log_suffix
        if running_as_cluster_job:
            # DK (2022): Added this part to make real.exe submittable to the cluster
            r = subprocess.run("srun /projects/0/ctdas/dkivits/scripts/submit_real.sh %s >& %s " %(wrf_run_dir, log_file),
                     shell=True,
                     stdout=subprocess.PIPE,
                     stderr=subprocess.STDOUT)
            r.check_returncode()
        else:
            r = subprocess.run("cd %s; source %s; ./real.exe >& %s" %(wrf_run_dir, s['bashrc'], log_file),
                     shell=True,
                     stdout=subprocess.PIPE,
                     stderr=subprocess.STDOUT)
            r.check_returncode()
    
        # Move rsl log files to own directory
        real_log_dir = os.path.join(wrf_run_dir, "rsl_real")
        os.makedirs(real_log_dir, exist_ok=True)
        rsl_files = glob.glob(os.path.join(wrf_run_dir, "rsl.*.0000"))
        for file in rsl_files:
            fp_dst = os.path.join(real_log_dir, os.path.basename(file))
            if os.path.exists(fp_dst):
                os.remove(fp_dst)
            shutil.move(file, fp_dst)
        ran_real = True  
    else:
        print("Output of real exists, overwrite_real=False, skipping running real.")
        ran_real = False
    
    
    ##########################################################
    ########## WRITE TRACERS TO wrfbdy AND wrfinput ##########
    ##########################################################
    if overwrite_tracers or ran_real:
        bdy_script = "python3 %s" %s['bc_ic_script']
                
        bdy_args = wrf_run_dir + " " + \
                   "CO2 " + \
                   str(nmembers) + " " + \
                   s["period"][0].strftime("%Y%m%d") + " " + \
                   s["period"][1].strftime("%Y%m%d")
        
        log_file = 'bc_ic' + log_suffix

        print("Writing tracers to wrfbdy and wrfinput (check %s)..."%log_file)

        r = subprocess.run(bdy_script + " " + bdy_args + " >& " + log_file,
                           shell=True,
                           stdout=subprocess.PIPE,
                           stderr=subprocess.STDOUT)
        r.check_returncode()
        
        # Also write VPRM input tracers --> Not necessary for CTDAS-WRF run!
        if "vprm_input_dir" in s.keys():
            print("Writing vprm tracers to wrfbdy and wrfinput..."%log_file)
            bdy_args_vprm = wrf_run_dir + " " + \
                   "CO2_BIO " + \
                   "0 " + \
                   s["period"][0].strftime("%Y%m%d") + " " + \
                   s["period"][1].strftime("%Y%m%d")

            log_file = 'vprm' + log_suffix

            r = subprocess.run(bdy_script + " " + bdy_args_vprm + " >& " + log_file,
                               shell=True,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT)
            r.check_returncode()
        
        
    else:
        print("real.exe not run and overwrite_tracers = False, skipping writing tracers to wrfbdy and wrfinput.")


    ############################################
    ########## PREPARE AS RESTART RUN ##########
    ############################################
    
    if "restart_run" in s.keys() and s["restart_run"]:
        # Todo:
        # - get wrfrst from wrf src dir
        # - copy tracers in there
        #cp wrfrst_d02_2015-02-09_00\:00\:00 wrfrst_d02.bak
        #ncks -A -C -v CO2_000,CO2_001,CO2_002,CO2_003,CO2_004 wrfinput_d02 wrfrst_d02_2015-02-09_00\:00\:00 
        #ncks -A -v NSTEP_D wrfrst_d02.bak wrfrst_d02_2015-02-09_00\:00\:00 
        #rm  wrfrst_d02.bak
        # - edit namelist
        raise ValueError("to do")

    
    ###########################################
    ########## SUBMIT WRF TO CLUSTER ##########
    ###########################################
    
    log_file =  "wrf" + log_suffix
    wrf_cmd = wrf_submit_cmd + " ./wrf.exe >&" + log_file
    print("Command for submitting WRF: " + wrf_cmd)
    if submit_wrf or (not glob.glob(os.path.join(wrf_run_dir, "wrfout*"))) : # 08/04/2022: changed 'and' into 'or'. 
        # Now logically unsound but at least WRF is submitting!
        print("Submitting WRF...")
        r = subprocess.run(wrf_cmd,
                           shell=True,
                           stdout=subprocess.PIPE,
                           stderr=subprocess.STDOUT)
        r.check_returncode()
    else:
        print("submit_wrf = False or output of wrf exists, skipping.")
    
    # Change back to original directory
    os.chdir(cwd)
    

################################
########## END OF WRF ##########
################################
