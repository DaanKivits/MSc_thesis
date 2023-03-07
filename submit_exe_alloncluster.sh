#!/bin/bash

# Purpose: Submit a generic command and its arguments to a job scheduler
# Author: Friedemann Reum, updated by Liesbeth
# DK, 2023: Changed original submit_python.exe script to include the submission of the real.exe WRF executable

# Usage message
function print_usage_and_exit
{
    echo "Usage: `basename $0` [-p partition_profile] [-n number_of_cores] [-t time_limit] [-N number_of_nodes] command [arguments]"
    exit 0
}


# Parse options
# https://sookocheff.com/post/bash/parsing-bash-script-arguments-with-shopts/
while getopts ":p:n:t:N:" opt; do
  # - getopts is a function. the following "<args>"  is the list of valid
  #   options, the next identifier is the running variable
  # - first ":" disable default error handling, enables specifying the "?" case
  # - second ":" indicates an option with an argument
  case ${opt} in
    p )
      partition_profile=$OPTARG
      # - OPTARG holds the argument of the currently parsed option
      ;;
    n )
      ntasks_arg=$OPTARG
      ;;
    N )
      nnodes_arg=$OPTARG
      ;;
    t )
      tlimit_arg=$OPTARG
      ;;
    \? )
      print_usage_and_exit
      ;;
  esac
done

# Remove options from argument list
shift $((OPTIND -1))

# Parse commands
if [ $# -eq 0 ]
  then
    print_usage_and_exit
fi

# Read name of command from arguments:
full_scriptname=$1

# Remove command from arguments list:
shift 

# For job- and logfiles, Strip of the extension and path if there is
# any from the command name:
# No path:
script_basename=`basename $full_scriptname`
# No path, no extension:
scriptname="${script_basename%.*}"
# Only path:
dir=`dirname $full_scriptname`

# Choose the log directory:
log_dir=$dir/log

# Write a job script

# Choose a default partition profile if it wasn't given in the options:
# https://unix.stackexchange.com/questions/212183/how-do-i-check-if-a-variable-exists-in-an-if-statement
if [[ ! -v partition_profile ]]
then
  partition_profile=thin
  echo "-p not given, using partition_profile '$partition_profile'"
fi

# Set SBATCH options according to partition_profile
case $partition_profile in
    short )
	# 'short' partition:
	partition=thin
	tlimit="01:00:00"
	mem=60000M
	ntasks=32 # less than 1 full node not possible anyway
	;;
    thin )
	# 'normal' partition:
	partition=thin
	tlimit="120:00:00"
	mem=60000M
	ntasks=32
	;;
    staging )
	# 'staging' partition:
	partition=staging
	tlimit="120:00:00" # that's the max limit
	#mem=48000
	ntasks=1
	;;
esac

# If number of cores was specified as argument, overwrite ncore
if [[ -v ntasks_arg ]]
then
  ntasks=$ntasks_arg
fi

# If tlimit was specified as argument, overwrite
if [[ -v tlimit_arg ]]
then
  tlimit=$tlimit_arg
fi


# Command to submit:
# (I don't know whether I should use $*, $@, "$*" or "$@" for the argument list -
# is there even a difference?)
# (Documentation of argument list: http://tldp.org/LDP/abs/html/internalvariables.html#ARGLIST)
full_command="srun $full_scriptname $*  >& $log_dir/$scriptname.\$SLURM_JOB_ID.log"
# Case without parallel computing:
# full_command="$full_scriptname $*  >& $log_dir/$scriptname.\$SLURM_JOB_ID.log"

# Write job submission script
cat > this.job <<EOF
#! /bin/bash
#SBATCH -p $partition
#SBATCH -t $tlimit
#SBATCH -n $ntasks
#SBATCH --mail-user=daan.kivits@wur.nl
#SBATCH --mail-type=FAIL,END

source ~/.bashrc_wrf

source activate daan

echo "Executing: $full_command"

$full_command
EOF

if [[ -v nnodes_arg ]]
then
  mv this.job this.job.tmp
  head -n3 this.job.tmp > this.job
  echo "#SBATCH -N $nnodes_arg" >> this.job
  tail -n7 this.job.tmp >> this.job
  rm this.job.tmp
fi

# Make the job script executable:
chmod u+x this.job

# Make the output directory:
mkdir -p $log_dir

# Choose a jobname
jobname=$scriptname

# Submit the job:
sbatch --job-name=$jobname -o $log_dir/$scriptname.%j.out this.job


