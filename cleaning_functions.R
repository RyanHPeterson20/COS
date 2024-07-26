

tropomi_polygons <- function(total_col_co, test_lonbond, test_latbond, bounds){
  
  #lots of manual movements (fix later)
  trop_lonbond_loc1 <- test_lonbond[1,,]
  trop_lonbond_loc2 <- test_lonbond[2,,]
  trop_lonbond_loc3 <- test_lonbond[3,,]
  trop_lonbond_loc4 <- test_lonbond[4,,]
  
  trop_latbond_loc1 <- test_latbond[1,,]
  trop_latbond_loc2 <- test_latbond[2,,]
  trop_latbond_loc3 <- test_latbond[3,,]
  trop_latbond_loc4 <- test_latbond[4,,]
  
  #get locations in matrix of Non_NA
  trop_col_co <- total_col_co[which(!is.na(total_col_co))]
  
  #select for pixel locations not associated with NA
  trop_lonbond_loc1 <- trop_lonbond_loc1[which(!is.na(total_col_co))]
  trop_lonbond_loc2 <- trop_lonbond_loc2[which(!is.na(total_col_co))]
  trop_lonbond_loc3 <- trop_lonbond_loc3[which(!is.na(total_col_co))]
  trop_lonbond_loc4 <- trop_lonbond_loc4[which(!is.na(total_col_co))]
  
  trop_latbond_loc1 <- trop_latbond_loc1[which(!is.na(total_col_co))]
  trop_latbond_loc2 <- trop_latbond_loc2[which(!is.na(total_col_co))]
  trop_latbond_loc3 <- trop_latbond_loc3[which(!is.na(total_col_co))]
  trop_latbond_loc4 <- trop_latbond_loc4[which(!is.na(total_col_co))]
  
  #delete blocks outsides of bounds
  N <- length(trop_col_co)
  delete_list <- NA
  for (i in 1:N) {
    if(trop_lonbond_loc2[i] < bounds$Lon_min | trop_latbond_loc1[i] < bounds$Lat_min | 
       trop_latbond_loc1[i] >  bounds$Lat_max | trop_lonbond_loc4[i] > bounds$Lon_max ){
      delete_list<- c( delete_list,i)
    }
  }
  delete_list <- delete_list[-1]
  
  #delete step
  trop_col_co <- trop_col_co[-delete_list]
  
  trop_lonbond_loc1 <- trop_lonbond_loc1[-delete_list]
  trop_lonbond_loc2 <- trop_lonbond_loc2[-delete_list]
  trop_lonbond_loc3 <- trop_lonbond_loc3[-delete_list]
  trop_lonbond_loc4 <- trop_lonbond_loc4[-delete_list]
  
  trop_latbond_loc1 <- trop_latbond_loc1[-delete_list]
  trop_latbond_loc2 <- trop_latbond_loc2[-delete_list]
  trop_latbond_loc3 <- trop_latbond_loc3[-delete_list]
  trop_latbond_loc4 <- trop_latbond_loc4[-delete_list]
  
  #get polygroups
  polyGroup <- NA
  N <- length(trop_lonbond_loc1)
  #N <- 50000
  for(i in 1:N){
    newPoly <- rbind(c(trop_lonbond_loc1[i], trop_latbond_loc1[i]),
                     c(trop_lonbond_loc2[i], trop_latbond_loc2[i]),
                     c(trop_lonbond_loc3[i], trop_latbond_loc3[i]),
                     c(trop_lonbond_loc4[i], trop_latbond_loc4[i]),
                     c(trop_lonbond_loc1[i], trop_latbond_loc1[i]))
    polyGroup<- c( polyGroup,list(newPoly))
  }
  
  polyGroup <- polyGroup[-1]
  test_col <- trop_col_co[1:N]
  return(list(polygons = polyGroup, data = test_col))
}


mopitt_polygons <- function(mopitt_co, mopitt_lon, mopitt_lat, bounds) {
  
  mopitt_co <- mopitt_co[which(!is.na(mopitt_co))]
  
  mopitt_lon <- mopitt_lon[which(!is.na(mopitt_co))]
  mopitt_lat <- mopitt_lat[which(!is.na(mopitt_co))]
  
  N <- length(mopitt_co)
  delete_list_2 <- NA
  for (j in 1:N) {
    if(mopitt_lon[j] < bounds$Lon_min | mopitt_lat[j] < bounds$Lat_min | 
       mopitt_lat[j] > bounds$Lat_max | mopitt_lon[j] > bounds$Lon_max) {
      delete_list_2<- c( delete_list_2,j)
    }
  }
  delete_list_2 <- delete_list_2[-1]
  
  mopitt_co <- mopitt_co[-delete_list_2]
  
  mopitt_lon <- mopitt_lon[-delete_list_2]
  mopitt_lat <- mopitt_lat[-delete_list_2]
  
  lat_kmdeg <- 1 / 110.574 #lat km in degrees
  lat_adj <- 11 * lat_kmdeg
  
  N <- length(mopitt_lon)
  mop_lon_loc1 <- NA
  mop_lon_loc2 <- NA
  mop_lon_loc3 <- NA
  mop_lon_loc4 <- NA
  
  mop_lat_loc1 <- NA
  mop_lat_loc2 <- NA
  mop_lat_loc3 <- NA
  mop_lat_loc4 <- NA
  
  for(k in 1:N){
    temp_lon_adj <- 11 * (1 / (111.320 * cos(mopitt_lat[k] * pi / 180)))
    mop_lon_loc1 <- c(mop_lon_loc1, mopitt_lon[k]-temp_lon_adj)
    mop_lon_loc2 <- c(mop_lon_loc2, mopitt_lon[k]+temp_lon_adj)
    mop_lon_loc3 <- c(mop_lon_loc3, mopitt_lon[k]+temp_lon_adj)
    mop_lon_loc4 <- c(mop_lon_loc4, mopitt_lon[k]-temp_lon_adj)
    
    mop_lat_loc1 <- c(mop_lat_loc1, mopitt_lat[k]-lat_adj)
    mop_lat_loc2 <- c(mop_lat_loc2, mopitt_lat[k]-lat_adj)
    mop_lat_loc3 <- c(mop_lat_loc3, mopitt_lat[k]+lat_adj)
    mop_lat_loc4 <- c(mop_lat_loc4, mopitt_lat[k]+lat_adj)
  }
  
  mop_lon_loc1 <- mop_lon_loc1[-1]
  mop_lon_loc2 <- mop_lon_loc2[-1]
  mop_lon_loc3 <- mop_lon_loc3[-1]
  mop_lon_loc4 <- mop_lon_loc4[-1]
  
  mop_lat_loc1 <- mop_lat_loc1[-1]
  mop_lat_loc2 <- mop_lat_loc2[-1]
  mop_lat_loc3 <- mop_lat_loc3[-1]
  mop_lat_loc4 <- mop_lat_loc4[-1]
  
  
  polyGroup_mop <- NA
  N <- length(mop_lon_loc1)
  #N <- 50000
  for(i in 1:N){
    newPoly <- rbind(c(mop_lon_loc1[i], mop_lat_loc1[i]),
                     c(mop_lon_loc2[i], mop_lat_loc2[i]),
                     c(mop_lon_loc3[i], mop_lat_loc3[i]),
                     c(mop_lon_loc4[i], mop_lat_loc4[i]),
                     c(mop_lon_loc1[i], mop_lat_loc1[i]))
    polyGroup_mop<- c( polyGroup_mop,list(newPoly))
  }
  
  
  polyGroup_mop <- polyGroup_mop[-1]
  test_col_mop <- mopitt_co[1:N]
  
  test_col_mop <- test_col_mop[order(mopitt_co)]
  polyGroup_mop <- polyGroup_mop[order(mopitt_co)]
  
  return(list(polygons = polyGroup_mop, data = test_col_mop))
}
