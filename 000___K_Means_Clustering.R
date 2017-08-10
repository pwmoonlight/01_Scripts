
###############################################################################################################
###############################################################################################################
###############################################################################################################
######################### Script by Peter Moonlight, Tiina Sarkinen et al. 2017 ###############################
###############################################################################################################
###############################################################################################################
###############################################################################################################

bg <- lapply(list.files(path="Y:/South America GIS/Brasil/Brazil Masked BIOCLIM", pattern="*.tif$", full.names = T)[[1]], raster)

prec <- lapply(list.files(path="Y:/World/Future Climate Scenarios/Current/prec_30s_bil", pattern="*.bil$", full.names = T), raster)

prec <- for(x in 12:length(prec)){
  prec[[x]] <- crop(prec[[x]], bg)
  prec[[x]] <- mask(prec[[x]], bg)
  writeRaster(prec[[x]], file=paste("Y:/South America GIS/Brasil/Brazil Masked BIOCLIM/Monthly_Precipitation/", prec[[x]]@data@names, ".tif", sep=""))
}

tmax <- lapply(list.files(path="Y:/World/Future Climate Scenarios/Current/tmax_30s_bil", pattern="*.bil$", full.names = T), raster)

tmax <- for(x in 1:12){
  tmax[[x]] <- crop(tmax[[x]], bg)
  tmax[[x]] <- mask(tmax[[x]], bg)
  writeRaster(tmax[[x]], file=paste("Y:/South America GIS/Brasil/Brazil Masked BIOCLIM/Monthly_tMax/", tmax[[x]]@data@names, ".tif", sep=""))
}

tmin <- lapply(list.files(path="Y:/World/Future Climate Scenarios/Current/tmin_30s_bil", pattern="*.bil$", full.names = T), raster)

tmin <- for(x in 8:12){
  tmin[[x]] <- crop(tmin[[x]], bg)
  tmin[[x]] <- mask(tmin[[x]], bg)
  writeRaster(tmin[[x]], file=paste("Y:/South America GIS/Brasil/Brazil Masked BIOCLIM/Monthly_tmin/", tmin[[x]]@data@names, ".tif", sep=""))
}




### Prepare the working space
### -------------------------

setwd("C:/000_Modeling_Working_Directory_000")
source(paste(getwd(), "/01_Scripts/Minor_Modules/1___Prepare_Working_Space.R", sep=""))



### Load Current Data ###
### ----------------- ###

require(raster)

bg <- lapply(list.files(path="Y:/South America GIS/Brasil/Brazil Masked BIOCLIM", pattern="*.tif$", full.names = T), raster)
bg <- stack(bg[1:2])
sample <- raster(nrows=bg@layers[[1]]@nrows/20, ncols=bg@layers[[1]]@ncols/20, ext=extent(bg))

# Read in the Current Precipition Data 
current_prec <- lapply(list.files(path="Y:/South America GIS/Brasil/Brazil Masked BIOCLIM/Monthly_Precipitation", pattern="*.tif$", full.names = T), raster)
current_prec <- lapply(current_prec, resample, sample)


Future_Climate_Scenarios <- list.dirs("Y:/South America GIS/Brasil/Future Climate Scenarios", full.names = F, recursive=F)

for(x in 1:length(Future_Climate_Scenarios)){
  writeLines(paste("Reading Precipitation Data for", Future_Climate_Scenarios[[x]], "..."))
  if(file.exists(paste("Y:/South America GIS/Brasil/Future Climate Scenarios/", Future_Climate_Scenarios[[x]], "/45pr50", sep=""))){
    writeLines(paste("...45 Scenario"))
    files <- list.files(paste("Y:/South America GIS/Brasil/Future Climate Scenarios/", Future_Climate_Scenarios[[x]], "/45pr50", sep=""), pattern="*.tif", full.names=T)
    files <- lapply(files, raster)
    files <- lapply(files, resample, sample)
    for(y in 1:length(files)){
      writeLines(paste("   ...", y, sep=""))
      files[[y]] <- overlay(files[[y]], current_prec[[y]], fun=function(x,y){return(x/y)})
      files[[y]] <- mean(files[[y]]@data@values, na.rm=T)[1]
    }
    results[x,1] <- mean(c(files[[1]][1], files[[2]][1], files[[3]][1], files[[4]][1], files[[5]][1], files[[6]][1], files[[7]][1], files[[8]][1], files[[9]][1], files[[10]][1], files[[11]][1], files[[12]][1]))
  }
  if(file.exists(paste("Y:/South America GIS/Brasil/Future Climate Scenarios/", Future_Climate_Scenarios[[x]], "/85pr50", sep=""))){
    writeLines(paste("...85 Scenario"))
    files <- list.files(paste("Y:/South America GIS/Brasil/Future Climate Scenarios/", Future_Climate_Scenarios[[x]], "/85pr50", sep=""), pattern="*.tif", full.names=T)
    files <- lapply(files, raster)
    files <- lapply(files, resample, sample)
    for(y in 1:length(files)){
      writeLines(paste("   ...", y, sep=""))
      files[[y]] <- overlay(files[[y]], current_prec[[y]], fun=function(x,y){return(x/y)})
      files[[y]] <- mean(files[[y]]@data@values, na.rm=T)[1]
    }
    results[x,4] <- mean(c(files[[1]][1], files[[2]][1], files[[3]][1], files[[4]][1], files[[5]][1], files[[6]][1], files[[7]][1], files[[8]][1], files[[9]][1], files[[10]][1], files[[11]][1], files[[12]][1]))
  }
}


write.csv(results, file="k_means_clustering.csv")





current_tm <- lapply(list.files(path="Y:/South America GIS/Brasil/Brazil Masked BIOCLIM/Monthly_tMin", pattern="*.tif$", full.names = T), raster)
current_tm <- lapply(current_tm, resample, sample)

for(x in 1:length(Future_Climate_Scenarios)){
  writeLines(paste("Reading Temperature Minimum for", Future_Climate_Scenarios[[x]], "..."))
  if(file.exists(paste("Y:/South America GIS/Brasil/Future Climate Scenarios/", Future_Climate_Scenarios[[x]], "/45tn50", sep=""))){
    writeLines(paste("...45 Scenario"))
    files <- list.files(paste("Y:/South America GIS/Brasil/Future Climate Scenarios/", Future_Climate_Scenarios[[x]], "/45tn50", sep=""), pattern="*.tif", full.names=T)
    files <- lapply(files, raster)
    files <- lapply(files, resample, sample)
  for(y in 1:length(files)){
      writeLines(paste("   ...", y, sep=""))
      files[[y]] <- overlay(files[[y]], current_prec[[y]], fun=function(x,y){return(x-y)})
      files[[y]] <- mean(files[[y]]@data@values, na.rm=T)[1]
    }
    results[x,3] <- mean(c(files[[1]][1], files[[2]][1], files[[3]][1], files[[4]][1], files[[5]][1], files[[6]][1], files[[7]][1], files[[8]][1], files[[9]][1], files[[10]][1], files[[11]][1], files[[12]][1]))
  }
  if(file.exists(paste("Y:/South America GIS/Brasil/Future Climate Scenarios/", Future_Climate_Scenarios[[x]], "/85tn50", sep=""))){
    writeLines(paste("...85 Scenario"))
    files <- list.files(paste("Y:/South America GIS/Brasil/Future Climate Scenarios/", Future_Climate_Scenarios[[x]], "/85tn50", sep=""), pattern="*.tif", full.names=T)
    files <- lapply(files, raster)
    files <- lapply(files, resample, sample)
    for(y in 1:length(files)){
      writeLines(paste("   ...", y, sep=""))
      files[[y]] <- overlay(files[[y]], current_prec[[y]], fun=function(x,y){return(x-y)})
      files[[y]] <- mean(files[[y]]@data@values, na.rm=T)[1]
    }
    results[x,6] <- mean(c(files[[1]][1], files[[2]][1], files[[3]][1], files[[4]][1], files[[5]][1], files[[6]][1], files[[7]][1], files[[8]][1], files[[9]][1], files[[10]][1], files[[11]][1], files[[12]][1]))
  }
}


write.csv(results, file="k_means_clustering.csv")



current_tx <- lapply(list.files(path="Y:/South America GIS/Brasil/Brazil Masked BIOCLIM/Monthly_tMax", pattern="*.tif$", full.names = T), raster)
current_tx <- lapply(current_tx, resample, sample)

for(x in 1:length(Future_Climate_Scenarios)){
  writeLines(paste("Reading Temperature Minimum for", Future_Climate_Scenarios[[x]], "..."))
  if(file.exists(paste("Y:/South America GIS/Brasil/Future Climate Scenarios/", Future_Climate_Scenarios[[x]], "/45tx50", sep=""))){
    writeLines(paste("...45 Scenario"))
    files <- list.files(paste("Y:/South America GIS/Brasil/Future Climate Scenarios/", Future_Climate_Scenarios[[x]], "/45tx50", sep=""), pattern="*.tif", full.names=T)
    files <- lapply(files, raster)
    files <- lapply(files, resample, sample)
    for(y in 1:length(files)){
      writeLines(paste("   ...", y, sep=""))
      files[[y]] <- overlay(files[[y]], current_prec[[y]], fun=function(x,y){return(x-y)})
      files[[y]] <- mean(files[[y]]@data@values, na.rm=T)[1]
    }
    results[x,2] <- mean(c(files[[1]][1], files[[2]][1], files[[3]][1], files[[4]][1], files[[5]][1], files[[6]][1], files[[7]][1], files[[8]][1], files[[9]][1], files[[10]][1], files[[11]][1], files[[12]][1]))
  }
  if(file.exists(paste("Y:/South America GIS/Brasil/Future Climate Scenarios/", Future_Climate_Scenarios[[x]], "/85tx50", sep=""))){
    writeLines(paste("...85 Scenario"))
    files <- list.files(paste("Y:/South America GIS/Brasil/Future Climate Scenarios/", Future_Climate_Scenarios[[x]], "/85tx50", sep=""), pattern="*.tif", full.names=T)
    files <- lapply(files, raster)
    files <- lapply(files, resample, sample)
    for(y in 1:length(files)){
      writeLines(paste("   ...", y, sep=""))
      files[[y]] <- overlay(files[[y]], current_prec[[y]], fun=function(x,y){return(x-y)})
      files[[y]] <- mean(files[[y]]@data@values, na.rm=T)[1]
    }
    results[x,5] <- mean(c(files[[1]][1], files[[2]][1], files[[3]][1], files[[4]][1], files[[5]][1], files[[6]][1], files[[7]][1], files[[8]][1], files[[9]][1], files[[10]][1], files[[11]][1], files[[12]][1]))
  }
}


write.csv(results, file="k_means_clustering.csv")















results <- as.data.frame(matrix(nrow=19, ncol=6))
colnames(results) <- c("45_prec_change", "45_tmax_change", "45_tmin_change", "85_prec_change", "85_tmax_change", "85_tmin_change")
Future_Climate_Scenarios <- list.dirs("Y:/South America GIS/Brasil/Future Climate Scenarios", full.names = F, recursive=F)
rownames(results) <- Future_Climate_Scenarios


