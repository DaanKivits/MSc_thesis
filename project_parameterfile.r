# Project parameter file from CTDAS gridded statevector to WRF domain
# only really makes sense for European domains

# Borrowing from learn_cartopy.ipynb


# Get domain parameters

## R library for reading nml files: can't cope with empty lines!
# if(!require("nml",quietly = T)){
#   remotes::install_github("jsta/nml")
# }

# Instead, execute python!
library(reticulate)

f90nml <- import("f90nml")

# fp_namelist <- "/Users/friedemann/data/namelist.wps"
# fp_c <- "~/data/countries_lambert.nc"
# fp_namelist <- "/Users/friedemann/data/WPS/namelist.wps.4.1.1"
# fp_c <- "~/data/ctdas/countries_lambert_4.1.1.nc"
# fp_namelist <- "/Users/friedemann/data/WPS/namelist.wps.4.1.1"




namelist_wps = f90nml$read(fp_namelist)

ggp = namelist_wps[["geogrid"]]
# namelist_input = f90nml.read("/Users/friedemann/data/namelist.input")

if(ggp[["ref_lon"]] != ggp[["stand_lon"]]){
  stop("ref_lon and stand_lon must be the same until I implement otherwise.")
  # They are:
  # ref_lon = center of coarse domain
  # stand_lon = center of projection
  
  # If they are not similar, the domain is not symmetric. I'd have to implement that.
}


# Create crs from namelist parameters:
# (see http://www.pkrc.net/wrf-lambert.html)

library(sp)
library(raster)

# # Get data
# fp = "~/data/ctdas/griddedNHparameters.nc"
# dat = raster(fp, varname="regions")

# I tried what happens when you swap truelat1 and truelat2 - projected points have the same coordinates
# ggp[["truelat1"]] <- 52.5
# ggp[["truelat2"]] <- 52.5

dst_crs = paste0("+proj=lcc ",
                 "+lat_1=", ggp[["truelat1"]], " ",
                 "+lat_2=", ggp[["truelat2"]], " ",
                 "+lat_0=", ggp[["ref_lat"]], " ",
                 "+lon_0=", ggp[["stand_lon"]], " ",
                 "+ellps=sphere ", # WRF globe representation
                 "+R=6370000 ",       # WRF globe representation (apparently has to be in m, see '$ proj -le')
                 "+units=m ", # to avoid confusion - dx and dy in namelist are in m as well
                 "+towgs84=0,0,0 ",
                 "+no_defs ")

# # try what these mean by projecting a ll point
# # to conform with wrf, I think lat_0 and lon_0 is (0,0)
# pts <- data.frame(lat=c(51,52,53, ggp[["ref_lat"]]),
#                     lon=c(7,8,9, ggp[["ref_lon"]]))
# coordinates(pts) <- c("lon","lat")
# crs(pts)="+proj=longlat +ellps=sphere +R=6370000 +units=km"
# pts_dst <- spTransform(pts,dst_crs)
# signif(coordinates(pts_dst),4)
# #       lon       lat
# # [1,] -105.00 -54.51
# # [2,]  -34.22  55.70
# # [3,]   33.46 166.90
# # [4,]    0.00   0.00



# coords1 <- coordinates(pts_dst)
# coords2 <- coordinates(pts_dst)

# So... if I take the WRF manual literally:
# then the translation is as follows:
# PROJ4   WRF
# lat_1   truelat1
# lat_2   truelat2
# lat_0   ref_lat
# lon_0   stand_lon

# wrf-python makes the following correspondence between global attributes in wrf netcdf files and proj4:
# _cartopy = crs.LambertConformal(
#   central_longitude=self.stand_lon,
#   central_latitude=self.moad_cen_lat,
#   standard_parallels=self._std_parallels,
#   globe=self._globe(),
#   cutoff=cutoff)

# In short:
# PROJ4   WRF
# lat_0   MOAD_CEN_LAT
# lon_0   STAND_LON

# This is confusing, because what is CEN_LON then?! Reminder, 
# the namelist and nc files seem to be like this:
# WRF namelist    WRF nc
# ref_lat         MOAD_CEN_LAT (apparently always same as CEN_LAT)
# ref_lon         CEN_LON         
# stand_lon       STAND_LON

# However: CEN_LAT and CEN_LON change with the domains! So it's
# logical that the other two, which stay invariant, refer to the projection!
# Now I think I should know everything, in case I disregard ref_x and ref_y.

if(any(is.element(c("ref_x","ref_y"),names(ggp)))){
  stop("ref_x and ref_y are specified. I don't know how to deal with that.")
}




# The following steps are no longer in learn_cartopy.ipynb
# hgt =  raster("~/data/geo_em.d01.nc",varname="HGT_M")
# Damn, it can't read it: 
# class       : RasterLayer 
# dimensions  : 91, 98, 8918  (nrow, ncol, ncell)
# resolution  : 1, 1  (x, y)
# extent      : 0.5, 98.5, 0.5, 91.5  (xmin, xmax, ymin, ymax)
# coord. ref. : NA 
# data source : /Users/friedemann/data/geo_em.d01.nc 
# names       : HGT_M 
# z-value     : 1 
# zvar        : HGT_M 

# Try eixport and https://earthscience.stackexchange.com/questions/13688/how-to-map-emission-inventory-from-latlon-corrdinate-to-wrf-model-grid

# # From eixport documentation
# library(eixport)
# wrfg <- wrf_grid("~/data/geo_em.d01.nc",type = "geo", epsg=4326)
# wrfg1 <- wrf_grid("~/data/geo_em.d01.nc",type = "geo", epsg=2000)
# # What the fuck! It should recognize the projection, and moreover, how can it be limited to epsg codes?! Hello?!
# 
# # Ok, try the next package
# remotes::install_github("atmoschem/wrftools")
# library(wrftools)


# Create destination raster
# Strategy: One raster for all domains. 
# Method: Raster with the extent of the coarsest domain and the resolution of the finest domain.

# Ratio between coarsest and finest domain:
ratio_fine <- max(ggp[["parent_grid_ratio"]])

# Resolution of finest domain:
res_x = ggp[["dx"]] / ratio_fine
res_y = ggp[["dy"]] / ratio_fine

# Grid cells in coarsest domain (mass grid)
nx_coarse = ggp[["e_we"]][1] - 1
ny_coarse = ggp[["e_sn"]][1] - 1

# Grid cells in fine domain resolution:
nx_fine <- nx_coarse*ratio_fine
ny_fine <- ny_coarse*ratio_fine

# Get extent - only centered domains implemented. This
# would change if ref_x and ref_y are specified, and
# 

xmin <- -res_x*nx_fine/2
xmax <- res_x*nx_fine/2
ymin <- -res_y*ny_fine/2
ymax <- res_y*ny_fine/2

# Create raster
dst_r <- raster(crs=dst_crs,
                xmn=xmin, xmx=xmax, ymn=ymin, ymx=ymax,
                resolution=c(res_x,res_y))

# cc <- coordinates(dst_r)


# Project parameter map
# Problem: only have methods ngb and bilinear
# Solution: same as I did previously: disaggregate - project (in this case, ngb!) - aggregate(modal!)

# How much finer than the destination grid do I want to go?
sampling_factor <- 10
dst_r_fine <- disaggregate(dst_r,fact=sampling_factor)

# Project to fine grid. Here: categorial variable, fine sampled. so use 'ngb'
dat_dst_fine <- projectRaster(from=dat,to=dst_r_fine,method="ngb")

# Aggregate - which is the value that's in there most often
dat_dst <- aggregate(dat_dst_fine,fact=sampling_factor,fun=modal)

# A short check:
# all(unique(dat_dst) %in% unique(dat))

# A better check: plot!

data(wrld_simpl, package="maptools")

clims <- range(values(dat_dst),na.rm=T)

wrld_source <- wrld_to_points(wrld_simpl,dat)
dat_df <- as.data.frame(rasterToPoints(dat))
names(dat_df) <- c("x","y","dat")

# From the world (source) map, plot Europe:
ggplot() +
  geom_raster(data=dat_df,mapping=aes(x=x,y=y,fill=dat)) +
  scale_fill_gradientn(colors=rainbow(5),limits=clims) +
  geom_path(data=wrld_source,mapping=aes(x=long,y=lat,group=group)) +
  theme_blank() + coord_fixed() +
  xlim(-12,32) +
  ylim(35,70)


# Plot the dst map:
wrld_dst <- wrld_to_points(wrld_simpl,dat_dst)
dat_dst_df <- as.data.frame(rasterToPoints(dat_dst))
names(dat_dst_df) <- c("x","y","dat")

ggplot() +
  geom_raster(data=dat_dst_df,mapping=aes(x=x,y=y,fill=dat)) +
  #scale_fill_gradientn(colors=rainbow(5),limits=clims) +
  geom_path(data=wrld_dst,mapping=aes(x=long,y=lat,group=group)) +
  theme_blank() + coord_fixed()

# Woohoo :-)

# I will now simplify this so that I can remove covariances for now.
# It seems to make sense to group the parameter numbers, because they
# appear already to be grouped




# Next step: write the projected parameter file for this domain to file. 
# This must adhere to the same format as the previous one.
# However, dimension names do not seem to be important.
# remap the values to 1:n_values
dat_dst_1n <- dat_dst
values(dat_dst_1n) <- as.integer(as.factor(values(dat_dst)))

# # Get a few settings from the original file
# library(ncdf4)
# ncf <- nc_open(fp)
# attributes <- ncatt_get(ncf,"regions")
# nc_close(ncf)


# Write to file
fp_dst <- "~/data/ctdas/griddedNHparameters_proj.nc"
crs(dat_dst_1n) <- "+proj=lcc +ellps=sphere +lat_1=51.5 +lat_2=52.5 +lat_0=51.5 +lon_0=8.5 +units=m"

write_regionsfile(fp_dst, dat_dst_1n, ggp)


dat_read <- raster(fp_dst, varname="regions")
dat_read





#sapply(id1, function(x) which(x==id2))


# alternatively, I could start right away from a nations map. Like,
# 1 parameter per state.
# See if I can do this quickly.

# Failed attempts
# # rasterize should do the trick
# countries_r <- rasterize(wrld_simpl, dat_dst, "UN")
# countries_r <- rasterize(wrld_simpl, dat_dst)
# 
# 
# xy_min <- extent(dat_dst)[c(1,3)]
# xy_d <- res(dat_dst)
# xy_n <- dim(dat_dst)[2:1] # !
# tg_grid_topo <- GridTopology(cellcentre.offset = c(xy_min[1]+0.5*xy_d[1],xy_min[2]+0.5*xy_d[2]),            
#                              cellsize = xy_d,                                                           
#                              cells.dim = xy_n)                                                          
# grid <- SpatialGrid(grid = tg_grid_topo,proj4string = proj4string(dat_dst))  
# 
# 
# # wrld_dst <- crop_wrld(wrld_simpl, dat_dst)
# wrld_dst <- proj_wrld(wrld_simpl, dat_dst)
# countries_r <- over(grid, wrld_dst)
# # nope
# 
# countries_r <- rasterize(wrld_proj, dat_dst, "UN")
# 
# library(rworldmap)
# countries_g_wrld <- rworldmap::gridCountriesDegreesHalf
# 
# countries_g_proj <- spTransform(countries_g_wrld, crs(dat_dst))
# 
# countries <- over(grid, countries_g_proj)
# # nope
# 
# countries <- rasterize(x=countries_g_proj, y=dat_dst, field="UN" )


# Something that works - start with a grid
# Start with a grid
c_g_wrld <- rworldmap::gridCountriesDegreesHalf
# manually make it a raster - thought that could be done automatically, but can't
c_r_wrld <- raster(c_g_wrld@grid)
crs(c_r_wrld) <- crs(c_g_wrld)
values(c_r_wrld) <- c_g_wrld$ISO_N3


# plot 
dat_wrld <- as.data.frame(rasterToPoints(c_r_wrld))
names(dat_wrld) <- c("x","y","dat")
wrld_wrld <- wrld_to_points(wrld_simpl,c_r_wrld)

ggplot() +
  geom_raster(data=dat_wrld,mapping=aes(x=x,y=y,fill=dat)) +
  #scale_fill_gradientn(colors=rainbow(5),limits=clims) +
  geom_path(data=wrld_wrld,mapping=aes(x=long,y=lat,group=group)) +
  theme_blank() + coord_fixed()

# yes


c_r_proj <- projectRaster(from=c_r_wrld, to=dst_r, method="ngb")

# plot 
dat_smpl_df <- as.data.frame(rasterToPoints(c_r_proj))
names(dat_smpl_df) <- c("x","y","dat")
wrld_dst <- wrld_to_points(wrld_simpl, dst_r)

ggplot() +
  geom_raster(data=dat_smpl_df,mapping=aes(x=x,y=y,fill=dat)) +
  #scale_fill_gradientn(colors=rainbow(5),limits=clims) +
  geom_path(data=wrld_dst,mapping=aes(x=long,y=lat,group=group)) +
  theme_blank() + coord_fixed()

# I fucking have it.

# clean it up, i.e. remove small countries and add them to neighboring countries
# country code definitions
# https://en.wikipedia.org/wiki/ISO_3166-1_numeric
v <- values(c_r_proj)
tv <- table(v)
tv[tv<100]
# -99  12 234 248 352 442 499 504 788 
# 96   2  10  12  95  15  96   1  12 
# names(tv) <- rworldmap::isoToName(names(tv)) # this one doesn't match. whatever.

# -99 = Kosovo - leave it

# 12 = algeria -> ocean
values(c_r_proj)[v==12] <- NA
# 234 = faroe islands -> ocean
values(c_r_proj)[v==234] <- NA
# 248 = Aland islands -> ocean
values(c_r_proj)[v==248] <- NA
# 352 = iceland -> ocean
values(c_r_proj)[v==352] <- NA
# 442 = Luxembourg -> france
values(c_r_proj)[v==442] <- 250
# 504 = morocco -> ocean
values(c_r_proj)[v==504] <- NA
# 788 = tunisia -> ocean
values(c_r_proj)[v==788] <- NA

head(sort(table(values(c_r_proj))))
# New minimum sizes:
# -99 499 705 807   8  56 
# 96  96 152 182 199 228 
# sort of ok

# now divide the ocean manually a bit by coordinates:
# atlantic, baltic, north sea, mediterranean, black sea
# values: use 4-digit numbers, because ISO_3 is only 3-digit numbers
water_code <- c(north_sea=1000, baltic_sea=2000, atlantic=3000, mediterranean=4000, black_sea=5000) + 1000



# get lat-long coordinates of grid cells
coords <- as.data.frame(coordinates(c_r_proj))
coordinates(coords) <- c("x","y")
projection(coords) <- projection(c_r_proj)
ll <- as.data.frame(spTransform(x=coords,CRSobj = "+proj=longlat")@coords)
head(ll) # OK
names(ll) <- c("long","lat")


# save
v <- values(c_r_proj)
# reload
values(c_r_proj) <- v


# north sea
values(c_r_proj)[is.na(values(c_r_proj)) & 
                   is_inside(ll$long,c(-2.5, 9)) & 
                   is_inside(ll$lat, c(48, 60))] <- water_code["north_sea"]


# black sea
values(c_r_proj)[is.na(values(c_r_proj)) & 
                   is_inside(ll$long,c(27, 43)) & 
                   is_inside(ll$lat, c(40, 48))] <- water_code["black_sea"]

# mediterranean (two rectangles, requires black sea to be set)
values(c_r_proj)[is.na(values(c_r_proj)) & 
                   is_inside(ll$long,c(-5.5, 38)) & 
                   is_inside(ll$lat, c(29, 42))] <- water_code["mediterranean"]

values(c_r_proj)[is.na(values(c_r_proj)) & 
                   is_inside(ll$long,c(0, 27)) & 
                   is_inside(ll$lat, c(42, 46))] <- water_code["mediterranean"]

# baltic sea (two rectangles)
values(c_r_proj)[is.na(values(c_r_proj)) & 
                   is_inside(ll$long,c(9, 32)) & 
                   is_inside(ll$lat, c(53, 62))] <- water_code["baltic_sea"]

values(c_r_proj)[is.na(values(c_r_proj)) & 
                   is_inside(ll$long,c(16, 27)) & 
                   is_inside(ll$lat, c(62, 66))] <- water_code["baltic_sea"]


# atlantic (the rest)
values(c_r_proj)[is.na(values(c_r_proj))] <- water_code["atlantic"]



# plot now...
dat_smpl_df <- as.data.frame(rasterToPoints(c_r_proj))
names(dat_smpl_df) <- c("x","y","dat")

ggplot() +
  geom_raster(data=dat_smpl_df,mapping=aes(x=x,y=y,fill=dat)) +
  #scale_fill_gradientn(colors=rainbow(5),limits=clims) +
  geom_path(data=wrld_dst,mapping=aes(x=long,y=lat,group=group)) +
  theme_blank() + coord_fixed()


# Write.
c_r_proj_1n <- c_r_proj
values(c_r_proj_1n) <- as.integer(as.factor(values(c_r_proj)))

dat_smpl_df <- as.data.frame(rasterToPoints(c_r_proj_1n))
names(dat_smpl_df) <- c("x","y","dat")

ggplot() +
  geom_raster(data=dat_smpl_df,mapping=aes(x=x,y=y,fill=dat)) +
  #scale_fill_gradientn(colors=rainbow(5),limits=clims) +
  geom_path(data=wrld_dst,mapping=aes(x=long,y=lat,group=group)) +
  theme_blank() + coord_fixed()



write_regionsfile(filename = fp_c,
                  r = c_r_proj_1n,
                  geogrid = ggp)


# IMPORTANT: CTDAS doesn't understand decreasing latitude axis.
# To make the axis increasing, execute the following commands in the folder where fp_c is stored:

# module load cdo
# cdo invertlat countries_lambert.nc countries_lambert_ctdas.nc

# This doesn't change the data, only how how the axis is stored.

#r <- raster(fp_c)


########## SIMPLIFIED PARAMETERS FOR TESTING ##########
fp_c <- "~/data/ctdas_input/regions_files/quarters1dom_v2.nc" # v1 was: values(dst_r) <- as.integer(coords[,1]>0) + 2*as.integer(coords[,2]>0)

# Just divide the domain in 4 pieces meeting at the center of the domain
# This is for 1 domain
res_x = ggp[["dx"]] 
res_y = ggp[["dy"]]

nx = ggp[["e_we"]][1] - 1
ny = ggp[["e_sn"]][1] - 1

xmin <- -res_x*nx/2
xmax <- res_x*nx/2
ymin <- -res_y*ny/2
ymax <- res_y*ny/2

# Create raster
dst_r <- raster(crs=dst_crs,
                xmn=xmin, xmx=xmax, ymn=ymin, ymx=ymax,
                resolution=c(res_x,res_y))

coords <- coordinates(dst_r)

values(dst_r) <- 1 + as.integer(coords[,1]>0) + 2*as.integer(coords[,2]>0)

write_regionsfile(fp_c, dst_r, ggp)








########## Countries merged, 1 domain ##########
# Started from copying the code from the start of this file, but cutting away all the legacy shit

data(wrld_simpl, package="maptools")

# Instead, execute python!
library(reticulate)
f90nml <- reticulate::import("f90nml")

# Input options
# domain_version <- "3"
domain_version <- "4"


# Output options
path_out <- "~/data/ctdas_input/regions_files/"
fp_out <- p0(path_out, "countries_merged_domv", domain_version,"_2.nc")

             
             
# Get domain definitions
if(domain_version == "3"){
  fp_namelist <- "/Users/friedemann/data/WPS/namelist.wps.4.1.1"
} else if(domain_version == "4"){
  fp_namelist <- "/Users/friedemann/projects/SCARBO/WP17/simulations/namelist.domv4.wps"
} else {
  stop(p0("Unknown domain_version: ", domain_version))
}
             
namelist_wps = f90nml$read(fp_namelist)
max_dom <- namelist_wps$share$max_dom

ggp = namelist_wps[["geogrid"]]

if(ggp[["ref_lon"]] != ggp[["stand_lon"]]){
  stop("ref_lon and stand_lon must be the same until I implement otherwise.")
  # I think I could do this by now, but anyway
}

# Create crs from namelist parameters:
# (see http://www.pkrc.net/wrf-lambert.html)

library(sp)
library(raster)

if(domain_version=="3"){
  dst_crs = paste0("+proj=lcc ",
                   "+lat_1=", ggp[["truelat1"]], " ",
                   "+lat_2=", ggp[["truelat2"]], " ",
                   "+lat_0=", ggp[["ref_lat"]], " ",
                   "+lon_0=", ggp[["stand_lon"]], " ",
                   "+ellps=sphere ", # WRF globe representation
                   "+R=6370000 ",       # WRF globe representation (apparently has to be in m, see '$ proj -le')
                   "+units=m ", # to avoid confusion - dx and dy in namelist are in m as well
                   "+towgs84=0,0,0 ",
                   "+no_defs ")
} else if(domain_version=="4"){
  dst_crs = paste0("+proj=lcc ",
                   "+lat_1=", ggp[["truelat1"]], " ",
                   "+lat_2=", ggp[["truelat2"]], " ",
                   "+lat_0=", ggp[["ref_lat"]], " ",
                   "+lon_0=", ggp[["stand_lon"]], " ",
                   "+ellps=WGS84 ",   # Here I switch to wgs84 because in panoply the domv3 map and landusef from wrfinput doesn't match
                   #"+ellps=sphere ", # WRF globe representation
                   #"+R=6370000 ",       # WRF globe representation (apparently has to be in m, see '$ proj -le')
                   "+units=m ", # to avoid confusion - dx and dy in namelist are in m as well
                   "+towgs84=0,0,0 ",
                   "+no_defs ")
}

if(any(is.element(c("ref_x","ref_y"),names(ggp)))){
  stop("ref_x and ref_y are specified. I don't know how to deal with that.")
}


# Create destination raster
# Strategy: One raster for all domains. 
# Method: Raster with the extent of the coarsest domain and the resolution of the finest domain.

# Ratio between coarsest and finest domain:
ratio_fine <- max(ggp[["parent_grid_ratio"]][1:max_dom])

# Resolution of finest domain:
res_x = ggp[["dx"]] / ratio_fine
res_y = ggp[["dy"]] / ratio_fine

# Grid cells in coarsest domain (mass grid)
nx_coarse = ggp[["e_we"]][1] - 1
ny_coarse = ggp[["e_sn"]][1] - 1

# Grid cells in fine domain resolution:
nx_fine <- nx_coarse*ratio_fine
ny_fine <- ny_coarse*ratio_fine

# Get extent - only centered domains implemented. This
# would change if ref_x and ref_y are specified

xmin <- -res_x*nx_fine/2
xmax <- res_x*nx_fine/2
ymin <- -res_y*ny_fine/2
ymax <- res_y*ny_fine/2

# Create raster
dst_r <- raster(crs=dst_crs,
                xmn=xmin, xmx=xmax, ymn=ymin, ymx=ymax,
                resolution=c(res_x,res_y))

# Project parameter map
# Problem: only have methods ngb and bilinear
# Solution: same as I did previously: disaggregate - project (in this case, ngb!) - aggregate(modal!)

# How much finer than the destination grid do I want to go?
sampling_factor <- 10
dst_r_fine <- disaggregate(dst_r,fact=sampling_factor)

wrld_pts <- wrld_to_points(wrld_simpl,dst_r)

library(raster)
library(sf)

file_shp <- "~/data/ctdas_input/regions_files/ne_10m_admin_0_countries/ne_10m_admin_0_countries.shp"
proj_src <- "~/data/ctdas_input/regions_files/ne_10m_admin_0_countries/ne_10m_admin_0_countries.prj"

shp <- st_read(file_shp)
shp_trans <- st_transform(shp, dst_crs)


# convert the shapefile to a raster based on a standardised background
# raster
countries <- c(
  "Albania"
  ,"Algeria"
  ,"Austria"
  ,"Belarus"
  ,"Belgium"
  ,"Bosnia and Herz."
  ,"Bulgaria"
  ,"Czechia"
  ,"Croatia"
  ,"Denmark"
  ,"Estonia"
  ,"Finland"
  ,"France"
  ,"Germany"
  ,"Greece"
  ,"Hungary"
  ,"Ireland"
  ,"Italy"
  ,"Kosovo"
  ,"Latvia"
  ,"Lithuania"
  ,"Luxembourg"
  ,"Macedonia"
  ,"Moldova"
  ,"Montenegro"
  ,"Morocco"
  ,"Netherlands"
  ,"Norway"
  ,"Poland"
  ,"Portugal"
  ,"Romania"
  ,"Russia"
  ,"Serbia"
  ,"Slovakia"
  ,"Slovenia"
  ,"Spain"
  ,"Sweden"
  ,"Switzerland"
  ,"Turkey"
  ,"Tunisia"
  ,"Ukraine"
  ,"United Kingdom"
  )

id <- sapply(countries, function(c)which(shp$NAME==c))

# Somehow it's possible to export the name directly, but I didn't get that to work
r <- rasterize(shp_trans[id,], dst_r)

# Aggregate countries to groups
if(domain_version=="3"){
  country_groups <- list(Alpine=c("Austria", "Switzerland"),
                         Baltics=c("Estonia","Latvia","Lithuania"),
                         Northeast_domain=c("Finland", "Russia"),
                         Benelux=c("Netherlands", "Belgium", "Luxembourg"),
                         Rom_Mold=c("Romania", "Moldova"),
                         Central_Balkan=c("Slovenia", "Croatia", "Bosnia and Herz.", "Serbia", "Kosovo", "Montenegro", "Macedonia", "Albania"))
} else if(domain_version=="4"){
  country_groups <- list(Alpine=c("Austria", "Switzerland"),
                         Baltics=c("Estonia","Latvia","Lithuania"),
                         Northeast_domain=c("Finland", "Russia"),
                         Southeast_domain=c("Greece", "Turkey"), # Turkey got way too small in this domain
                         Benelux=c("Netherlands", "Belgium", "Luxembourg"),
                         Rom_Mold=c("Romania", "Moldova"),
                         Central_Balkan=c("Slovenia", "Croatia", "Bosnia and Herz.", "Serbia", "Kosovo", "Montenegro", "Macedonia", "Albania"),
                         Northern_Africa = c("Algeria", "Morocco", "Tunisia"))
}
countries_grouped <- countries

for(n_group in 1:length(country_groups)){
  id_group <- which(is.element(countries, country_groups[[n_group]]))
  countries_grouped <- c(setdiff(countries_grouped, country_groups[[n_group]]), names(country_groups)[n_group])
  countries <- c(countries, names(country_groups)[n_group])
  values(r)[is.element(values(r),id_group)] <- length(countries)               
}

# Change value keys to countries_grouped
rg <- r
for (val in unique(values(r))) {
  id <- which(values(r)==val)
  values(rg)[id] <- which(countries_grouped==countries[val])
}

r <- rg

# plot 
df <- as.data.frame(rasterToPoints(r))
df$layer <- factor(countries_grouped[df$layer])

ggplot() +
  geom_raster(data=df,mapping=aes(x=x,y=y,fill=layer)) +
  #scale_fill_gradientn(colors=rainbow(5)) +
  geom_path(data=wrld_pts,mapping=aes(x=long,y=lat,group=group)) +
  theme_blank() + coord_fixed()



# now divide the ocean manually a bit by coordinates: same as original countries-vector
# atlantic, baltic, north sea, mediterranean, black sea
# values: use 4-digit numbers, because ISO_3 is only 3-digit numbers
if(domain_version=="3"){
  water_code <- c(north_sea=1000, baltic_sea=2000, atlantic=3000, mediterranean=4000, black_sea=5000) + 1000
} else if(domain_version=="4"){
  # Black sea goes to Bulgaria
  water_code <- c(north_sea=1000, baltic_sea=2000, atlantic=3000, mediterranean=4000) + 1000
}
# Nope, do it right the first time:
water_code[] <- length(countries_grouped)+1:length(water_code)

regions <- c(countries_grouped, names(water_code))

c_r_proj <- rg
# get lat-long coordinates of grid cells
coords <- as.data.frame(coordinates(c_r_proj))
coordinates(coords) <- c("x","y")
projection(coords) <- projection(c_r_proj)
ll <- as.data.frame(spTransform(x=coords,CRSobj = "+proj=longlat")@coords)
head(ll) # OK
names(ll) <- c("long","lat")


# save
v <- values(c_r_proj)
# reload
values(c_r_proj) <- v


# north sea
values(c_r_proj)[is.na(values(c_r_proj)) & 
                   is_inside(ll$long,c(-2.5, 9)) & 
                   is_inside(ll$lat, c(48, 60))] <- water_code["north_sea"]


# black sea
if(domain_version == "3"){
  black_sea_code <- water_code["black_sea"]
} else if(domain_version == "4"){
  # In domain 4, black sea is tiny/not there so add it to Bulgaria
  black_sea_code <- which(regions=="Bulgaria")
}

values(c_r_proj)[is.na(values(c_r_proj)) & 
                   is_inside(ll$long,c(27, 43)) & 
                   is_inside(ll$lat, c(40, 48))] <- black_sea_code

# mediterranean (two rectangles, requires black sea to be set)
values(c_r_proj)[is.na(values(c_r_proj)) & 
                   is_inside(ll$long,c(-5.5, 38)) & 
                   is_inside(ll$lat, c(29, 42))] <- water_code["mediterranean"]

values(c_r_proj)[is.na(values(c_r_proj)) & 
                   is_inside(ll$long,c(0, 27)) & 
                   is_inside(ll$lat, c(42, 46))] <- water_code["mediterranean"]

# baltic sea (two rectangles)
values(c_r_proj)[is.na(values(c_r_proj)) & 
                   is_inside(ll$long,c(9, 32)) & 
                   is_inside(ll$lat, c(53, 62))] <- water_code["baltic_sea"]

values(c_r_proj)[is.na(values(c_r_proj)) & 
                   is_inside(ll$long,c(16, 27)) & 
                   is_inside(ll$lat, c(62, 66))] <- water_code["baltic_sea"]


# atlantic (the rest)
values(c_r_proj)[is.na(values(c_r_proj))] <- water_code["atlantic"]



# plot now...
df <- as.data.frame(rasterToPoints(c_r_proj))

df$layer[!is_inside(df$layer,c(27,31),incl = T)] <- NA # Water only
df$layer[df$layer!=which(regions=="Germany")] <- NA # Single regions only
df$layer[df$layer!=which(regions=="Bulgaria")] <- NA # Single regions only
df$layer[df$layer!=which(regions=="Southeast_domain")] <- NA # Single regions only

df$layer <- factor(regions[df$layer])

# df$layer[df$layer != "Turkey"] <- NA

ggplot() +
  geom_raster(data=df,mapping=aes(x=x,y=y,fill=layer),alpha=0.5) +
  geom_path(data=wrld_pts,mapping=aes(x=long,y=lat,group=group)) +
  guides(fill=guide_legend(title="Region")) +
  theme_blank() + coord_fixed()


fp_plot <- sub(".nc$",".png",fp_out)
ggsave(filename=fp_plot, width=8,height=5,dpi=300)


# How big are certain countries or regions?

length(which(df$layer=="Germany"))
# [1] 275
length(which(df$layer=="Central_Balkan"))
# [1] 222
length(which(df$layer=="France"))
# [1] 427
length(which(df$layer=="Italy"))
# [1] 238

# Well, do it right!
tmp <- sort(table(df$layer))
data.frame(tmp)
# Var1 Freq
# 1           Denmark   36
# 2          Slovakia   37
# 3           Benelux   55
# 4           Ireland   55
# 5           Czechia   59
# 6          Portugal   68
# 7           Hungary   72
# 8            Turkey   85
# 9          Bulgaria   86
# 10           Alpine   94
# 11        black_sea   95
# 12           Greece  100
# 13          Baltics  133
# 14          Belarus  159
# 15           Norway  179
# 16   United Kingdom  188
# 17         Rom_Mold  212
# 18   Central_Balkan  222
# 19            Italy  238
# 20           Poland  245
# 21          Germany  275
# 22           Sweden  285
# 23          Ukraine  294
# 24       baltic_sea  352
# 25            Spain  404
# 26        north_sea  422
# 27           France  427
# 28 Northeast_domain  771
# 29    mediterranean 1039
# 30         atlantic 2231


# Write.


write_regionsfile(filename = fp_out,
                  r = c_r_proj,
                  geogrid = ggp,
                  parameter_names = regions)



