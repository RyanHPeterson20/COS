#updated visualizations for JSM 2024

#libraries
suppressMessages(library(LatticeKrig))

#color
suppressMessages(library(grDevices))

#data
setwd("~/COS_LK")
load("DF_LK/trop_data_aus.rda")
load("DF_LK/mop_data_aus.rda")

load("DF_LK/trop_data_full.rda")
load("DF_LK/mop_data_full.rda")


# bounding regions 
##full
Lat_min <- -50
Lat_max <- 10
Lon_min <- 90
Lon_max <- 180
bounds <- data.frame(Lon_min, Lon_max, Lat_min, Lat_max )

full_domain <- rbind(c(90,-50), 
                     c(180,10))

##reduced
aus.lat.range <- c(-48, -10) 
aus.lon.range <- c(140, 155)
reduced_bounds <- data.frame(Lon_min = aus.lon.range[1], 
                             Lon_max = aus.lon.range[2],
                             Lat_min = aus.lat.range[1],
                             Lat_max = aus.lat.range[2])

reduced_domain <- rbind(c(140,-48), 
                        c(155,-10))


#data mess

temp_index <- which.max(convert_tropFull)
temp_trop <- convert_tropFull[-which(convert_tropFull >= 1e+19)]
temp_polyTrop <- polyGroup_tropFull[-which(convert_tropFull >= 1e+19)]


sorted_vec <- sort(convert_tropFull, decreasing = TRUE)
sorted_vec[1:10]

sorted_mop <- sort(mop1_coFull, decreasing = TRUE)
sorted_mop[1:10]

hist(convert_tropFull)
hist(mop1_coFull)

log_trop <- log(convert_tropFull)

#range
#zlim <- range(z,na.rm = TRUE) #base zlim setup



#TODO: adapt this for all data, and 
zlim_trop <- range(c(temp_trop), na.rm = TRUE)

#color setup
test_trop <- as.matrix(temp_trop)


nlevel <- 64
midpoints<- seq( zlim_trop[1], zlim_trop[2],,nlevel)
delta<- (midpoints[2]- midpoints[1])/2
# nlevel +1 breaks with the min and max as midpoints 
# of the first and last bins.
breaks <- c( midpoints[1]- delta, midpoints + delta)

zcol_trop <- drape.color(test_trop, col = tim.colors(nlevel), midpoint = FALSE, zlim = zlim_trop, 
                         transparent.color = "white", breaks = breaks)$color.index

#zcol_trop[which.max(convert_tropFull), 1]

#corrected viz
plot( full_domain, type="n",xlab = "Longitude", ylab = "Latitude",
      main = "TROPOMI - Total Column CO")
for (i in 1:length(temp_polyTrop)) {
  pcol <- c(zcol_trop[i, 1])
  polygon(temp_polyTrop[[i]], col = pcol, border=NA)
}
world(add = TRUE, lwd = 2)

# draw each poly with different color including the border
# if the border color has not been specified.
# this will avoid missing some space on some output devices.
# one can also crank down width of border lines to avoid rounded corners

#TODO: repeat above but based on the limits defined by the fused and predicted data
#repeat for the reduced range

##base visualizations

##tropomi
#jet color pallet
colors_trop <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                             "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))(length(convert_tropFull))

colors_trop <- colors_trop[rank(convert_tropFull)]
#polyNew <- polyGroups[[rank(col_co)]]

plot( full_domain, type="n",xlab = "Longitude", ylab = "Latitude",
      main = "TROPOMI - Total Column CO")
for ( k in seq_along(polyGroup_tropFull)) {
  polygon(polyGroup_tropFull[[k]], col=colors_trop[k], border=NA)
}
world(add = TRUE, lwd = 2)

##mopitt
#jet color pallet
colors_mop <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                                  "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))(length(mop1_coFull))

colors_mop <- colors_mop[rank(mop1_coFull)]
#polyNew <- polyGroups[[rank(col_co)]]

plot( full_domain, type="n",xlab = "Longitude", ylab = "Latitude",
      main = "MOPITT - Total Column CO")
for ( k in seq_along(polyGroups_mopFull)) {
  polygon(polyGroups_mopFull[[k]], col=colors_mop[k], border=NA)
}
world(add = TRUE, lwd = 2)



