# ctdas_lib.r
# Helper functions for working with CTDAS-data

library(ncdf4)

get_ctdas_input <- function(path, file_pattern, ti, te, vars=NA){
  # Reimplementation of column_obs.add_observations()
  
  get_file_name <- function(pattern, date){
    # Helper function
    date_str <- format(date, "%Y%m%d")
    file_name <- sub(pattern="<YYYYMMDD>",replacement = date_str, x= pattern)
    return(file_name)
  }
  
  dates <- seq(from=ti, to=te, by="day")
  
  
  # Get number of observations
  nobs <- numeric(length(dates)-1)
  for(n in 1:(length(dates)-1)){
    file_name <- get_file_name(file_pattern, dates[n])
    file_path <- file.path(path, file_name)
    if(!file.exists(file_path)){
      next
    }
    ncf <- nc_open(file_path)
    nobs[n] <- ncf$dim$soundings$len
  }
  
  # Select variables (from first existing file)
  if(all(is.na(vars))){
    for(n in 1:(length(dates)-1)){
      file_name <- get_file_name(file_pattern, dates[n])
      file_path <- file.path(path, file_name)
      if(file.exists(file_path)){
        break
      }
    }
    ncf <- nc_open(file_path)
    vars <- names(ncf$var)
    nc_close(ncf)
  }
  vars <- union(vars, c("sounding_id", "date"))
  
  # Allocate output
  dat <- data.frame(matrix(data=NA_real_, nrow=sum(nobs), ncol=length(vars)))
  names(dat) <- vars
  if(is.element("date", vars)){
    dat$date <- as.POSIXct(dat$date, origin=POSIXt_origin, tz="GMT")
  }
  
  # Get dat
  for(n in 1:(length(dates)-1)){
    file_name <- get_file_name(file_pattern, dates[n])
    file_path <- file.path(path, file_name)
    if(!file.exists(file_path)){
      next
    }
    
    ncf <- nc_open(file_path)
    for(var in vars){
      #print(var)
      if(!is.element(var, names(ncf$var))){
        warning(p0("Variable ", var, " not found, skipping."))
        next
      }
      tmp <- ncvar_get(ncf, var, collapse_degen = F)
      if(var=="date"){ # Parse date
        tmp <- as.POSIXct(apply(t(tmp),1,function(x)pct(sprintf("%d-%02d-%02d %02d:%02d:%02d.%d",x[1],x[2],x[3],x[4],x[5],x[6],x[7]))),origin=POSIXt_origin,tz="GMT")
      }
      if(is.null(dim(tmp)) || length(dim(tmp))==1){
        dat[[var]][sum(nobs[1:n]) - (nobs[n]:1) + 1] <- tmp
        # This code was for debugging the first case with 1 obs in a file
        # if(length(tmp) != length(sum(nobs[1:n]) - (nobs[n]:1) + 1)){
        #   print(p0("n = ", n))
        #   print(p0("var = ", var))
        #   stop("look here")
        # }
      } else if(length(dim(tmp))==2){
        if(class(dat[[var]]) != "list"){
          dat[[var]] <- as.list(dat[[var]])
        }
        tmp <- t(tmp)
        for(nn in 1:dim(tmp)[1]){
          dat[[var]][[sum(nobs[1:n]) - (nobs[n]:1)[nn] + 1]] <- tmp[nn,]
        }
      } else {
        warning(p0("Variable '", var, "' has more than two dimensions, don't know how to handle this. Skipping."))
        next
      }
      
    }
    nc_close(ncf)
  }
  return(dat)
}


write_regionsfile <- function(filename, r, geogrid, parameter_names=NULL){
  # Writes regionsfile that can be read by CTDAS
  # Only works with lambert projection right now.
  library(ncdf4)
  # filename = filename
  # r = raster data to write
  # geogrid = list of projection-parameters as in geogrid-section of namelist
  
  if(geogrid$map_proj != "lambert")
    stop("can only deal with Lambert projection")
  
  writeRaster(x = r,
              filename = filename,
              format = "CDF",
              varname="regions",
              varunit="unitless",
              longname="parameter_for_each_gridbox",
              # xname="lon", yname="lat",
              overwrite=F)
  
  # THINGS TO ADD FOR PANOPLY AND OTHERS
  # I.e., cf conventions. see cf-conventions-1.7.pdf
  
  # # add coordinates attribute
  # ncf <- nc_open(filename,write = T)
  # ncatt_put(ncf,"regions","coordinates","easting northing")
  # nc_close(ncf)
  
  # Add standard name to coordinates
  ncf <- nc_open(filename,write = T)
  ncatt_put(ncf,"easting","standard_name","projection_x_coordinate")
  ncatt_put(ncf,"northing","standard_name","projection_y_coordinate")
  nc_close(ncf)
  
  # Add projection information compliant with cf conventions
  ncf <- nc_open(filename,write = T)
  var <- ncvar_def(name="projection",units = "",dim = list(),
                   prec = "integer")
  ncvar_add(ncf,var)
  nc_close(ncf)
  
  ncf <- nc_open(filename,write = T)
  ncatt_put(ncf,"projection","grid_mapping_name","lambert_conformal_conic")
  ncatt_put(ncf,"projection","standard_parallel",c(geogrid$truelat1,geogrid$truelat2))
  ncatt_put(ncf,"projection","latitude_of_projection_origin",geogrid$ref_lat)
  # This problem didn't occur in domv4, so I'm commenting the warnings now.
  #warning("I'm not sure about this part. This is just to have Panoply plot it.")
  #warning("When reading it back in, lat_2 becomes 45?! whatevs, only needs to work roughly right now.")
  ncatt_put(ncf,"projection","longitude_of_central_meridian",geogrid$stand_lon)
  nc_close(ncf)
  
  # # Not necessary for panoply, but for cf compliance:
  # # However, when I do this, then reading it back in with raster (dat_read below) gives:
  # # Warning message:
  # #   In cbind(m[i, ], vals) :
  # #   number of rows of result is not a multiple of vector length (arg 2)
  # ncf <- nc_open(filename,write = T)
  # ncatt_put(ncf,"regions","grid_mapping","projection")
  # nc_close(ncf)
  
  # Region names
  if(!is.null(parameter_names)){
    nparam_dim <- ncdim_def(name = "parameters", units="", vals=1:length(parameter_names), create_dimvar = FALSE)
    char_len_dim <- ncdim_def(name = "char", units="", vals=1:as.integer(max(50, max(nchar(parameter_names)))), create_dimvar = FALSE)
    name_var <- ncvar_def(name = "parameter_names", longname = "Parameter names", units = "", dim = list(char_len_dim, nparam_dim), prec = "char")
    ncf <- nc_open(filename,write = T)
    ncvar_add(ncf,name_var)
    nc_close(ncf)
    ncf <- nc_open(filename,write = T)
    ncvar_put(ncf, "parameter_names", parameter_names)
    nc_close(ncf)
  }
  
  
  # IMPORTANT! FOR CTDAS, FLIP HOW Y COORDINATE IS STORED 
  # I couldn't find an option to do this in writeRaster
  # I previously did this with "cdo invertlat", but it deleted all other variables (projection, parameter_names)
  # Man, screw cdo! It also messed up cat_ncfiles! Gonna use 
  system(p0("ncpdq -O -h -a -northing ", filename, " ", filename))
  
  
  
}

read_ctdas_output <- function(ctdas_path, rc_file, home_rc, obs_vars=NA, get_obs_ens=FALSE){
  # At the moment, reads optimized state vector, observations and regions map
  # ctdas_path = the directory containing the directories exec/ and output/
  # Name of the rc_file used in the ctdas run (located in ctdas_path/exec/)
  # home_rc = the home path as given in the rc_file (need to change this to the local home)
  # obs_vars - only return certain variables in the observation files (reading the profiles is slow)
  # get_obs_ens - get ensemble of simulated observations in addition to means. For tens of
  #               thousands of data points, consider turning this off to save time
  print("Reading settings")
  
  library(reticulate)
  pymod_rc <- reticulate::import_from_path(module="rc",path=file.path(ctdas_path, "exec/da/tools"))
  HOME <- Sys.getenv("HOME")
  
  # Read settings
  output_path <- file.path(ctdas_path,"output")
  exec_path <- file.path(ctdas_path,"exec")
  
  ctdas_rc <- pymod_rc$read(rcfilename = file.path(exec_path, rc_file))
  dasystem_rc <- pymod_rc$read(rcfilename = file.path(exec_path, ctdas_rc$da.system.rc))
  obsop_rc <- pymod_rc$read(rcfilename = file.path(exec_path, ctdas_rc$da.obsoperator.rc))
  
  rc <- list(ctdas=ctdas_rc, dasystem=dasystem_rc, obsop=obsop_rc)
  
  # rc only returns strings, convert some values to their appropriate data types
  rc$ctdas$time.start <- pct(rc$ctdas$time.start)
  rc$ctdas$time.finish <- pct(rc$ctdas$time.finish)
  rc$ctdas$time.cycle <- as.integer(rc$ctdas$time.cycle)
  rc$ctdas$da.optimizer.nmembers <- as.integer(rc$ctdas$da.optimizer.nmembers)
  rc$ctdas$time.nlag <- as.integer(rc$ctdas$time.nlag)
  rc$ctdas$time.cycle <- as.integer(rc$ctdas$time.cycle)
  rc$dasystem$nparameters <- as.integer(rc$dasystem$nparameters)
  
  # Replace home in paths
  paths <- c("dir.da_run", "obs.column.rc", "regionsfile", "obs.column.input.dir")
  rcs <-   c("ctdas",      "dasystem",      "dasystem",    "dasystem")
  for(np in 1:length(paths)){
    rc[[rcs[np]]][[paths[np]]] <- sub(home_rc, HOME, rc[[rcs[np]]][[paths[np]]])
  }
  
  ########## Read statevector and atmospheric data
  # I *THINK* that for posterior emissions and data, I'm looking for the
  # first lag in the respective date, and for the prior data, I'm looking for the
  # last date the data in date - time.cycle*time.nlag because that's when it was
  # initialized.
  
  # First, adjust time.finish according to where we have output
  # (my preprocessing doesn't reach the end of the month + time.nlag, so CTDAS always fails nlag days before its supposed end)
  #rc$ctdas$time.finish <- rc$ctdas$time.finish - as.difftime(rc$ctdas$time.nlag*rc$ctdas$time.cycle - 1, units="days")
  
  # Allocate variables
  n_cycles <- as.numeric(rc$ctdas$time.finish-rc$ctdas$time.start,units="days") / rc$ctdas$time.cycle
  
  param_mean_pri <- array(dim=c(n_cycles, rc$dasystem$nparameters))
  param_mean_opt <- param_mean_pri
  
  param_ens_pri <- array(dim=c(n_cycles, rc$ctdas$da.optimizer.nmembers, rc$dasystem$nparameters))
  param_ens_opt <- param_ens_pri
  
  # Read observations
  # There is more information in the original files, and I can already read them. So I will.
  # This means I'm reimplementing column_obs.get_simulations... maybe I should have learned
  # plotting in python instead...
  if(!is.null(obs_vars)){
    obs_path <- rc$dasystem$obs.column.input.dir
    obs_pattern <- rc$dasystem$obs.column.ncfile
    ti <- rc$ctdas$time.start
    te <- rc$ctdas$time.finish + as.difftime(rc$ctdas$time.nlag*rc$ctdas$time.cycle - 1, units="days")
    print("Getting observations")
    obs <- get_ctdas_input(path = obs_path, file_pattern = obs_pattern, ti=ti, te=te, vars=obs_vars)
    
    # Map for getting id to which a simulation file belongs
    # (iirc looking for names is instantaneous compared to extremely slow sapply)
    id_map <- 1:length(obs$sounding_id)
    names(id_map) <- as.character(obs$sounding_id)
    
    # Initialize simulation variables
    obs$mod_mean_pri <- NA_real_
    obs$mod_mean_opt <- NA_real_
    if(get_obs_ens){
      obs$mod_ens_pri <-  NA_real_
      obs$mod_ens_opt <-  NA_real_
      obs$mod_ens_pri <- as.list(obs$mod_ens_pri)
      obs$mod_ens_opt <- as.list(obs$mod_ens_opt)
    }
  }
  print("Getting parameters")
  
  for(n in 0:(n_cycles-1)){
    date <- rc$ctdas$time.start + as.difftime(rc$ctdas$time.cycle*n,units="days")
    date_str <- strftime(date, "%Y%m%d")
    
    # Read statevector
    file_name <- file.path(output_path,date_str,p0("savestate_", date_str, ".nc"))
    
    # My preprocessing for the first monthly runs doesn't reach the end of
    # the month + time.nlag, so CTDAS always fails nlag days before its supposed end
    # Catch this case here
    if(!file.exists(file_name)){
      break
    }
    
    ncf <- nc_open(file_name)
    # The statevector elements of which optimization is final: always only the first
    nlag_opt <- 1
    if(n+nlag_opt <= n_cycles){
      param_mean_opt[n + nlag_opt, ] <- ncvar_get(ncf, "statevectormean_opt")[, nlag_opt]
      param_ens_opt[n + nlag_opt, , ] <- t(ncvar_get(ncf, "statevectorensemble_opt")[, ,nlag_opt])
    }
    # The statevector elements of which the background is the prior (i.e.,
    # newly initialized ensemble)
    if(n == 0){
      # First cycle: all are initialized
      nlags_pri <- 1:rc$ctdas$time.nlag
    } else {
      # All other cycles: last nlag is initialized
      nlags_pri <- rc$ctdas$time.nlag
    }
    for(nlag_pri in nlags_pri){
      if(n+nlag_pri <= n_cycles){
        param_mean_pri[n + nlag_pri, ] <- ncvar_get(ncf, "statevectormean_prior")[, nlag_pri]
        param_ens_pri[n + nlag_pri, , ] <- t(ncvar_get(ncf, "statevectorensemble_prior")[, , nlag_pri])
      }
    }
    nc_close(ncf)
    
    # Read simulated atmospheric observations
    if(!is.null(obs_vars)){
      file_path <- file.path(output_path,date_str, p0("optimizer.", date_str, ".nc"))
      if(file.exists(file_path)){
        print(p0("Getting simulated observations for lag ", date_str))
        ncf <- nc_open(file_path)
        
        # Get sounding id
        #sid <- ncvar_get(ncf, "sounding_id")
        # Well, it's still called obspack_num... let's just be thankful it's in there at all
        sid <- ncvar_get(ncf, "obspack_num")
        id <- id_map[as.character(sid)]
        
        # # Verify obs is the same in obs and optimizer nc
        # obs_ <- ncvar_get(ncf, "observed")
        # print(p0("max(abs(obs_-obs$obs[id])) = ", max(abs(obs_-obs$obs[id]))))
        
        # Read simulated data
        if(!all(is.na(obs$mod_mean_pri[id]))){
          stop("Zoinks! mod_mean_pri was set before!")
        }
        if(!all(is.na(obs$mod_mean_opt[id]))){
          stop("Zoinks! mod_mean_opt was set before!")
        }
        obs$mod_mean_pri[id] <- ncvar_get(ncf, "modelsamplesmean_prior")
        obs$mod_mean_opt[id] <- ncvar_get(ncf, "modelsamplesmean_optimized")
        
        # Read simulated ensemble
        if(get_obs_ens){
          ens_pri <- t(ncvar_get(ncf, "modelsamplesdeviations_prior"))
          ens_opt <- t(ncvar_get(ncf, "modelsamplesdeviations_optimized"))
          for(nid in 1:length(id)){
            if(!all(is.na(obs$mod_ens_pri[[id[nid]]]))){
              stop("Zoinks! mod_ens_pri was set before!")
            }
            if(!all(is.na(obs$mod_ens_opt[[id[nid]]]))){
              stop("Zoinks! mod_ens_opt was set before!")
            }
            obs$mod_ens_pri[[id[nid]]] <- ens_pri[nid,]
            obs$mod_ens_opt[[id[nid]]] <- ens_opt[nid,]
          }
        }
        
        nc_close(ncf)
      } else {
        print(p0("No optimizer file for lag ", date_str, ". No data for this lag?"))
      }
    
    }
    
    # Save which cycle has an optimized statevector
    n_has_data <- n+nlag_opt
  }
  
  # In case the run failed at some point (see comment on the loop break above), adjust n_cycles
  if(n_has_data<n_cycles){
    n_cycles <- n_has_data
    rc$ctdas$time.finish <- date + as.difftime(rc$ctdas$time.cycle,units="days")
    param_mean_pri <- param_mean_pri[1:n_cycles,,drop=F]
    param_mean_opt <- param_mean_opt[1:n_cycles,,drop=F]
    
    param_ens_pri <-  param_ens_pri[1:n_cycles,,,drop=F]
    param_ens_opt <-  param_ens_opt[1:n_cycles,,,drop=F]
  }
  
  
  ########## Read regions file
  print("Getting regions definition")
  library(ncdf4)
  library(raster)
  data("wrld_simpl",package = "maptools")
  regions_map <- raster(rc$dasystem$regionsfile,varname="regions")
  ncf <- nc_open(rc$dasystem$regionsfile)
  projection(regions_map) <- ncatt_get(ncf, 0, "proj4")$value
  parameter_names <- ncvar_get(ncf, "parameter_names")
  nc_close(ncf)
  
  # Duplicate parameter names for overlaying parametrs
  if(is.element("n_emis_proc", names(rc$dasystem))){
    n_emis_proc <- as.integer(rc$dasystem$n_emis_proc)
    nparams <- length(parameter_names)
    parameter_names <- p0(rep(parameter_names,times=n_emis_proc), "_proc_", rep(1:n_emis_proc, each=nparams))
  }
  
  
  ########## Make dimension names
  ti <- rc$ctdas$time.start
  dt <- rc$ctdas$time.cycle
  dates <- seq(from=ti, length.out=n_cycles,by=as.difftime(dt,units="days"))
  dates_fmt <- strftime(dates,tz="GMT")
  nmembers <- rc$ctdas$da.optimizer.nmembers
  
  dimnames_ens <- list(date=dates_fmt, member=1:nmembers, parameter = parameter_names)
  dimnames(param_ens_pri) <- dimnames_ens
  dimnames(param_ens_opt) <- dimnames_ens
  
  dimnames_mean <- list(date=dates_fmt, parameter = parameter_names)
  dimnames(param_mean_pri) <- dimnames_mean
  dimnames(param_mean_opt) <- dimnames_mean
  
  
  ########## Return all data in one list
  result <- list(rc=rc,
                 n_cycles=n_cycles,
                 param_mean=list(pri=param_mean_pri,opt=param_mean_opt),
                 param_ens=list(pri=param_ens_pri,opt=param_ens_opt),
                 regions=list(map=regions_map, names=parameter_names))
  if(!is.null(obs_vars)){
    result$obs=obs
  }
  print("Done")
  return(result)
}

get_total <- function(ctdas, fp_total_flux){
  # For simulations with overlaying parameters,
  # computes parameters (mean and ensemble) of total fluxes
  # Have to run calc_prior_fluxes.r first, it computes
  # total fluxes per parameter per cycle.
  
  # Get total flux
  tf <- read.csv(fp_total_flux, row.names = 1)
  
  # Set parameters
  do_offset <- is.element("sigma_offset", names(ctdas$rc$dasystem))
  n_emis_proc <- as.integer(ctdas$rc$dasystem$n_emis_proc)
  n_parameters <- ctdas$rc$dasystem$nparameters
  n_flux_parameters <- n_parameters - as.integer(do_offset)
  n_regions <- n_flux_parameters/n_emis_proc
  tmp <- strsplit(ctdas$regions$names[1:n_regions],"_proc_")
  region_names <- sapply(tmp, `[[`, 1)
  
  # Initialize result variables
  mn_arr <- array(dim=c(ctdas$n_cycles, n_regions),
                  dimnames=list(dimnames(ctdas$param_mean$pri)[[1]], region_names))
  ens_arr <- array(dim=c(ctdas$n_cycles, ctdas$rc$ctdas$da.optimizer.nmembers, n_regions),
                   dimnames=list(dimnames(ctdas$param_ens$pri)[[1]], dimnames(ctdas$param_ens$pri)[[2]], region_names))
  ctdas$param_mean_total <- list(pri=mn_arr, opt=mn_arr)
  ctdas$param_ens_total <- list(pri=ens_arr, opt=ens_arr)
  
  ctdas$param_mean_nee <- list(pri=mn_arr, opt=mn_arr)
  ctdas$param_ens_nee <- list(pri=ens_arr, opt=ens_arr)
  
  mn_proc_arr <- array(dim=c(ctdas$n_cycles, n_regions, n_emis_proc),
                  dimnames=list(dimnames(ctdas$param_mean$pri)[[1]], region_names, 1:n_emis_proc))
  ens_proc_arr <- array(dim=c(ctdas$n_cycles, ctdas$rc$ctdas$da.optimizer.nmembers, n_regions, n_emis_proc),
                   dimnames=list(dimnames(ctdas$param_ens$pri)[[1]], dimnames(ctdas$param_ens$pri)[[2]], region_names, 1:n_emis_proc))
  ctdas$param_mean_abs <- list(pri=mn_proc_arr, opt=mn_proc_arr)
  ctdas$param_ens_abs <- list(pri=ens_proc_arr, opt=ens_proc_arr)
  
  # Compute
  for(nr in 1:n_regions){
    # Total flux
    proc_ind <- nr +(0:(n_emis_proc-1))*n_regions
    ctdas$param_mean_total$pri[, nr] <- apply(ctdas$param_mean$pri[,proc_ind]*tf[,proc_ind],1,sum)
    ctdas$param_mean_total$opt[, nr] <- apply(ctdas$param_mean$opt[,proc_ind]*tf[,proc_ind],1,sum)
    
    # This line would work as first step if I add some magic that makes it an array
    # of the right dimension. But it would be badly readable, so I'll use a loop.
    # tmp <- apply(ctdas$param_ens$opt[,,proc_ind],2,`*`, tf[,proc_ind])
    for(nm in 1:dim(ctdas$param_ens$pri)[2]){
      ctdas$param_ens_total$opt[, nm, nr] <- apply(ctdas$param_ens$opt[,nm,proc_ind]*tf[,proc_ind],1,sum)
      ctdas$param_ens_total$pri[, nm, nr] <- apply(ctdas$param_ens$pri[,nm,proc_ind]*tf[,proc_ind],1,sum)
    }
    
    # NEE (copy-paste from above but with proc_ind only covering GPP and RESP)
    proc_ind <- nr +(1:(n_emis_proc-1))*n_regions
    ctdas$param_mean_nee$pri[, nr] <- apply(ctdas$param_mean$pri[,proc_ind]*tf[,proc_ind],1,sum)
    ctdas$param_mean_nee$opt[, nr] <- apply(ctdas$param_mean$opt[,proc_ind]*tf[,proc_ind],1,sum)
    for(nm in 1:dim(ctdas$param_ens$pri)[2]){
      ctdas$param_ens_nee$opt[, nm, nr] <- apply(ctdas$param_ens$opt[,nm,proc_ind]*tf[,proc_ind],1,sum)
      ctdas$param_ens_nee$pri[, nm, nr] <- apply(ctdas$param_ens$pri[,nm,proc_ind]*tf[,proc_ind],1,sum)
    }
    
    # Other processes
    for(ne in 1:n_emis_proc){
      proc_ind <- nr +(ne-1)*n_regions
      ctdas$param_mean_abs$pri[, nr, ne] <- ctdas$param_mean$pri[,proc_ind]*tf[,proc_ind]
      ctdas$param_mean_abs$opt[, nr, ne] <- ctdas$param_mean$opt[,proc_ind]*tf[,proc_ind]
      for(nm in 1:dim(ctdas$param_ens$pri)[2]){
        ctdas$param_ens_abs$opt[, nm, nr, ne] <- ctdas$param_ens$opt[,nm,proc_ind]*tf[,proc_ind]
        ctdas$param_ens_abs$pri[, nm, nr, ne] <- ctdas$param_ens$pri[,nm,proc_ind]*tf[,proc_ind]
      }
    }
    
    
  }
  ctdas$prior_flux <- tf
  return(ctdas)
  
}



get_unc_red <- function(ctdas){
  # Get CTDAS uncertainty reduction
  # ctdas = output from read_ctdas_output
  # the "-1" remove the ensemble member that represents the mean
  ctdas$unc_rnd <- list()
  ctdas$unc_rnd$pri <- apply(ctdas$param_ens$pri[,-1,],c(1,3),sd)
  ctdas$unc_rnd$opt <- apply(ctdas$param_ens$opt[,-1,],c(1,3),sd)
  ctdas$unc_rnd$red <- 1-ctdas$unc_rnd$opt/ctdas$unc_rnd$pri
  return(ctdas)
}

get_unc_red_total_nee_abs <- function(ctdas){
  # Get CTDAS uncertainty reduction of total flux
  # ctdas = output from read_ctdas_output, must have applied get_total
  # the "-1" remove the ensemble member that represents the mean
  ctdas$unc_rnd_total <- list()
  ctdas$unc_rnd_total$pri <- apply(ctdas$param_ens_total$pri[,-1,],c(1,3),sd)
  ctdas$unc_rnd_total$opt <- apply(ctdas$param_ens_total$opt[,-1,],c(1,3),sd)
  ctdas$unc_rnd_total$red <- 1-ctdas$unc_rnd_total$opt/ctdas$unc_rnd_total$pri
  
  ctdas$unc_rnd_nee <- list()
  ctdas$unc_rnd_nee$pri <- apply(ctdas$param_ens_nee$pri[,-1,],c(1,3),sd)
  ctdas$unc_rnd_nee$opt <- apply(ctdas$param_ens_nee$opt[,-1,],c(1,3),sd)
  ctdas$unc_rnd_nee$red <- 1-ctdas$unc_rnd_nee$opt/ctdas$unc_rnd_nee$pri
  
  ctdas$unc_rnd_abs <- list()
  ctdas$unc_rnd_abs$pri <- apply(ctdas$param_ens_abs$pri[,-1,,],c(1,3,4),sd)
  ctdas$unc_rnd_abs$opt <- apply(ctdas$param_ens_abs$opt[,-1,,],c(1,3,4),sd)
  ctdas$unc_rnd_abs$red <- 1-ctdas$unc_rnd_abs$opt/ctdas$unc_rnd_abs$pri
  
  return(ctdas)
}


get_cor <- function(ctdas_run){
  # Compute correlations among prior and posterior state vector elements
  
  
  nparams <- dim(ctdas_run$param_mean$opt)[2]

  ctdas_run$cor_pri <- array(dim=c(ctdas_run$n_cycles, nparams, nparams))
  ctdas_run$cor_opt <- ctdas_run$cor_pri
  ctdas_run$p_pri <- ctdas_run$cor_pri
  ctdas_run$p_opt <- ctdas_run$cor_pri
  
  for(nc in 1:ctdas_run$n_cycles){
    # Get covariance matrices
    cov_pri <- cov(ctdas_run$param_ens$pri[nc,-1,])
    ctdas_run$cor_pri[nc, , ] <- cov2cor(cov_pri)
    cov_opt <- cov(ctdas_run$param_ens$opt[nc,-1,])
    ctdas_run$cor_opt[nc, , ] <- cov2cor(cov_opt)
    # Get p-values
    for(n1 in 1:(nparams-1)){
      for(n2 in (n1+1):nparams){
        ctdas_run$p_pri[nc, n1, n2] <- summary(lm(ctdas_run$param_ens$pri[nc, -1, n1] ~ ctdas_run$param_ens$pri[nc, -1, n2]))$coefficients[2,4]
        ctdas_run$p_opt[nc, n1, n2] <- summary(lm(ctdas_run$param_ens$opt[nc, -1, n1] ~ ctdas_run$param_ens$opt[nc, -1, n2]))$coefficients[2,4]
      }
    }
  }
  
  return(ctdas_run)
}
