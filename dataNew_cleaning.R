
#library
suppressMessages(library(lubridate))
suppressMessages(library(fields))
suppressMessages(library(ncdf4)) #for .nc file type
suppressMessages(library(raster))
suppressMessages(library(viridis))
suppressMessages(library( scales))

#color
suppressMessages(library(grDevices))

#timing library
suppressMessages(library(bench))
suppressMessages(library(tictoc))

#functions:
setwd("~/COS_LK")
source("DF_LK/cleaning_functions.R")

# bounding regions (wider)

Lat_min <- -50
Lat_max <- 10
Lon_min <- 90
Lon_max <- 180
bounds <- data.frame(Lon_min, Lon_max, Lat_min, Lat_max )

full_domain <- rbind(c(90,-50), 
                     c(180,10))

# bounding regions (narrow)
aus.lat.range <- c(-48, -10) 
aus.lon.range <- c(140, 155)
reduced_bounds <- data.frame(Lon_min = aus.lon.range[1], 
                             Lon_max = aus.lon.range[2],
                             Lat_min = aus.lat.range[1],
                             Lat_max = aus.lat.range[2])

reduced_domain <- rbind(c(140,-48), 
                        c(155,-10))


# tropomi

##first
setwd("~/COS_LK/DF_LK")
file.exists('DataNEW/S5P_RPRO_L2__CO_____20191230T013047_20191230T031217_11459_03_020400_20220903T073807.nc')

#first overpass (pass 1)
trop1.nc <- nc_open('DataNEW/S5P_RPRO_L2__CO_____20191230T013047_20191230T031217_11459_03_020400_20220903T073807.nc')

#second overpass (pass 2) #E Aus
trop2.nc <- nc_open('DataNEW/S5P_RPRO_L2__CO_____20191230T031217_20191230T045347_11460_03_020400_20220903T073808.nc')

#third overpass (pass 3) #E Aus
trop3.nc <- nc_open('DataNEW/S5P_RPRO_L2__CO_____20191230T045347_20191230T063517_11461_03_020400_20220903T073808.nc')

#fourth overpass (pass 4)
trop4.nc <- nc_open('DataNEW/S5P_RPRO_L2__CO_____20191230T063517_20191230T081647_11462_03_020400_20220903T073809.nc')

#data prep/checks
#nc_time <- ncvar_get(trop1.nc, "PRODUCT/time_utc")

##data 1 (pass 1)
trop1_co <- ncvar_get(trop1.nc, "PRODUCT/carbonmonoxide_total_column") #lots of NA
trop1_lonbond <- ncvar_get(trop1.nc, "PRODUCT/SUPPORT_DATA/GEOLOCATIONS/longitude_bounds")
trop1_latbond <- ncvar_get(trop1.nc, "PRODUCT/SUPPORT_DATA/GEOLOCATIONS/latitude_bounds")

##data 2 (pass 2)
trop2_co <- ncvar_get(trop2.nc, "PRODUCT/carbonmonoxide_total_column") #lots of NA
trop2_lonbond <- ncvar_get(trop2.nc, "PRODUCT/SUPPORT_DATA/GEOLOCATIONS/longitude_bounds")
trop2_latbond <- ncvar_get(trop2.nc, "PRODUCT/SUPPORT_DATA/GEOLOCATIONS/latitude_bounds")

trop2_co_correct <- ncvar_get(trop2.nc, "PRODUCT/carbonmonoxide_total_column_corrected") #lots of NA


##data 3 (pass 3)
trop3_co <- ncvar_get(trop3.nc, "PRODUCT/carbonmonoxide_total_column") #lots of NA
trop3_lonbond <- ncvar_get(trop3.nc, "PRODUCT/SUPPORT_DATA/GEOLOCATIONS/longitude_bounds")
trop3_latbond <- ncvar_get(trop3.nc, "PRODUCT/SUPPORT_DATA/GEOLOCATIONS/latitude_bounds")

trop3_co_correct <- ncvar_get(trop3.nc, "PRODUCT/carbonmonoxide_total_column_corrected") #lots of NA
trop3_prec <- ncvar_get(trop3.nc, "PRODUCT/carbonmonoxide_total_column_precision") #lots of NA

##data 4 (pass 4)
trop4_co <- ncvar_get(trop4.nc, "PRODUCT/carbonmonoxide_total_column") #lots of NA
trop4_lonbond <- ncvar_get(trop4.nc, "PRODUCT/SUPPORT_DATA/GEOLOCATIONS/longitude_bounds")
trop4_latbond <- ncvar_get(trop4.nc, "PRODUCT/SUPPORT_DATA/GEOLOCATIONS/latitude_bounds")


#TODO fix the bounding on this, the reduced_bounds doesn't seem to be working
tic()
trop1_poly <- tropomi_polygons(trop1_co, trop1_lonbond, trop1_latbond, reduced_bounds)
toc()

tic()
trop2_poly <- tropomi_polygons(trop2_co, trop2_lonbond, trop2_latbond, reduced_bounds)
toc()

tic()
trop2_poly_correct <- tropomi_polygons(trop2_co_correct , trop2_lonbond, trop2_latbond, reduced_bounds)
toc()

tic()
trop3_poly <- tropomi_polygons(trop3_co, trop3_lonbond, trop3_latbond, reduced_bounds)
toc()

tic()
trop3_poly_correct <- tropomi_polygons(trop3_co_correct , trop3_lonbond, trop3_latbond, reduced_bounds)
toc()

tic()
trop3_poly_prec <- tropomi_polygons(trop3_prec , trop3_lonbond, trop3_latbond, reduced_bounds)
toc()

#note: below values are 
test_data <- trop3_poly_correct$data
test_prec <- trop3_poly_prec$data

#trop1_poly currently giving NAs
#combing data

polyGroup_trop <- c(trop2_poly$polygons, trop3_poly$polygons)
col_co <- c(trop2_poly$data, trop3_poly$data)
convert_co_trop <- col_co * 6.02214e+19
scale_lim <- range(convert_co_trop)
zero_vals <- which(convert_co_trop <= 0)
length(zero_vals)

if (length(zero_vals) > 0) {
  convert_co_trop <- convert_co_trop[-zero_vals]
  polyGroup_trop <- polyGroup_trop[-zero_vals]
}


save(convert_co_trop, polyGroup_trop, file = "trop_data_aus.rda")


#corrected for banding:

polyGroup_trop_correct <- c(trop2_poly_correct$polygons, trop3_poly_correct$polygons)
col_co_corr <- c(trop2_poly_correct$data, trop3_poly_correct$data)
convert_co_trop_corr <- col_co_corr * 6.02214e+19
scale_lim <- range(convert_co_trop_corr)
zero_vals <- which(convert_co_trop_corr <= 0)
length(zero_vals)

if (length(zero_vals) > 0) {
  convert_co_trop <- convert_co_trop[-zero_vals]
  polyGroup_trop <- polyGroup_trop[-zero_vals]
}

save(convert_co_trop_corr, col_co_corr, 
     polyGroup_trop_correct, file = "trop_data_corrected.rda")

#repeat for full area

#TODO fix the bounding on this, the reduced_bounds doesn't seem to be working
tic()
trop1_polyFull <- tropomi_polygons(trop1_co, trop1_lonbond, trop1_latbond, bounds)
toc()

tic()
trop2_polyFull <- tropomi_polygons(trop2_co, trop2_lonbond, trop2_latbond, bounds)
toc()

tic()
trop3_polyFull <- tropomi_polygons(trop3_co, trop3_lonbond, trop3_latbond, bounds)
toc()

tic()
trop4_polyFull <- tropomi_polygons(trop4_co, trop4_lonbond, trop4_latbond, bounds)
toc()

#trop1_poly currently giving NAs
#combing data

polyGroup_tropFull <- c(trop1_polyFull$polygons, trop2_polyFull$polygons,
                    trop3_polyFull$polygons, trop4_polyFull$polygons)
col_co <- c(trop1_polyFull$data, trop2_polyFull$data,
            trop3_polyFull$data, trop4_polyFull$data)
convert_tropFull <- col_co *6.022140857e+19
scale_lim <- range(convert_tropFull)
zero_vals <- which(convert_tropFull <= 0)
length(zero_vals)

if (length(zero_vals) > 0) {
  convert_tropFull <- convert_tropFull[-zero_vals]
  polyGroup_tropFull <- polyGroup_tropFull[-zero_vals]
}

save(polyGroup_tropFull, convert_tropFull, file = "trop_data_full.rda")

##VIS test
#jet color pallet
colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                             "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))(length(convert_co_trop))

colors <- colors[rank(convert_co_trop)]
#polyNew <- polyGroups[[rank(col_co)]]

plot( full_domain, type="n",xlab = "Longitude", ylab = "Latitude",
      main = "TROPOMI - Total Column CO")
for ( k in seq_along(polyGroup_trop)) {
  polygon(polyGroup_trop[[k]], col=colors[k], border=NA)
}
world(add = TRUE, lwd = 2)

#For MOPITT data
##first mop
setwd("~/COS_LK/DF_LK")
file.exists('DataNEW/MOP02J-20191230-L2V19.9.3.he5')
#first day (day 1)
mop1.nc <- nc_open('DataNEW/MOP02J-20191230-L2V19.9.3.he5')

#second day (day 2)
mop2.nc <- nc_open("DataNEW/MOP02J-20191231-L2V19.9.3.he5")

mop_co <- ncvar_get(mop1.nc, "HDFEOS/SWATHS/MOP02/Data Fields/RetrievedCOTotalColumn")
mop_fix <- ncvar_get(mop1.nc, "HDFEOS/SWATHS/MOP02/Data Fields/APrioriCOTotalColumn")


mop_lat <-  ncvar_get(mop1.nc, "HDFEOS/SWATHS/MOP02/Geolocation Fields/Latitude")
mop_lon <-  ncvar_get(mop1.nc, "HDFEOS/SWATHS/MOP02/Geolocation Fields/Longitude")

#not sure if needed but we'll add it
mop_aprior <- ncvar_get(mop1.nc, "HDFEOS/SWATHS/MOP02/Data Fields/APrioriCOTotalColumn")

mop1_co <- mop_co[1, ] #selects for data

#mop1_co[1]
#mop_aprior[1]

mop1_polyFull <- mopitt_polygons(mop1_co, mop_lon, mop_lat, bounds)

mop1_coFull <- mop1_polyFull$data
polyGroups_mopFull <- c(mop1_polyFull$polygons)
zero_vals <- which(mop1_coFull <= 0)
length(zero_vals)

if (length(zero_vals) > 0) {
  mop1_coFull <- mop1_coFull[-zero_vals]
  polyGroups_mopFull <- polyGroups_mopFull[-zero_vals]
}

save(mop1_coFull, polyGroups_mopFull, file = "mop_data_full.rda")

save(mop1_col, polyGroups_mop, file = "mop_data_aus.rda")

##VIS test
#jet color pallet
colors_mop <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                             "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))(length(mop1_col))

colors_mop <- colors_mop[rank(mop1_col)]
#polyNew <- polyGroups[[rank(col_co)]]

plot( reduced_domain, type="n",xlab = "Longitude", ylab = "Latitude",
      main = "MOPITT - Total Column CO")
for ( k in seq_along(polyGroups_mop)) {
  polygon(polyGroups_mop[[k]], col=colors_mop[k], border=NA)
}
world(add = TRUE, lwd = 2)


#mop1_co[1]
#mop_aprior[1]

#TODO: fix this for everything later.
save(mop1_col, polyGroups_mop, file = "mop_data_test.rda")
save(col_co, polyGroups, file = "trop_data_test.rda")


#even more reduced domain:
# bounding regions (narrow)
aus.lat.range <- c(-38, -33) 
aus.lon.range <- c(142, 152)
reduced_bounds <- data.frame(Lon_min = aus.lon.range[1], 
                             Lon_max = aus.lon.range[2],
                             Lat_min = aus.lat.range[1],
                             Lat_max = aus.lat.range[2])



tic()
trop2_polyRed <- tropomi_polygons(trop2_co, trop2_lonbond, trop2_latbond, reduced_bounds)
toc()

tic()
trop3_polyRed <- tropomi_polygons(trop3_co, trop3_lonbond, trop3_latbond, reduced_bounds)
toc()

polyGroup_tropRed <- c(trop2_polyRed$polygons)
col_co <- c(trop2_polyRed$data)
co_tropRed <- col_co *6.022140857e+19

save(co_tropRed, polyGroup_tropRed, file = "reduced_trop.rda")

mop1_polyRed <- mopitt_polygons(mop1_co, mop_lon, mop_lat, reduced_bounds)

co_mopRed <- mop1_polyRed$data
polyGroups_mopRed <- c(mop1_polyRed$polygons)

save(co_mopRed, polyGroups_mopRed, file = "reduced_mop.rda")

