###############################################################################################################
  ###############################################################################################################
    ###############################################################################################################
      ######################### Script by Peter Moonlight, Tiina Sarkinen et al 2017 ################################
    ###############################################################################################################
  ###############################################################################################################
###############################################################################################################

writeLines(paste("\nLoading Background Data:\n"))

Future_Climate_Scenarios <- list.dirs("Y:/South America GIS/Brasil/Future Climate Scenarios", full.names = F, recursive=F)
dir.create("12_Future_Projections", showWarnings = F)
lapply(Future_Climate_Scenarios, function(x){dir.create(paste("12_Future_Projections/", x, sep=""), showWarnings = F)})

Future_Climate_Data <- lapply(1:length(Future_Climate_Scenarios), function(x){lapply(1:2, function(x){})})

for(x in 1:length(Future_Climate_Scenarios)){
  writeLines(paste("Reading Data for", Future_Climate_Scenarios[[x]], "..."))
  if(file.exists(paste("Y:/South America GIS/Brasil/Future Climate Scenarios/", Future_Climate_Scenarios[[x]], "/45bi50", sep=""))){
    writeLines(paste("...45 Scenario"))
    files <- list.files(paste("Y:/South America GIS/Brasil/Future Climate Scenarios/", Future_Climate_Scenarios[[x]], "/45bi50", sep=""), pattern="*3_degrees.tif", full.names=T)
    files <- files[PCA_bioclim]
    files <- lapply(files, raster)
    for(z in 1:length(files)){
      files[[z]]@data@names <- bg_bioclim[[z]]@data@names
    }
    Future_Climate_Data[[x]][[1]] <- stack(files)
  }
  if(file.exists(paste("Y:/South America GIS/Brasil/Future Climate Scenarios/", Future_Climate_Scenarios[[x]], "/85bi50", sep=""))){
    writeLines(paste("...85 Scenario"))
    files <- list.files(paste("Y:/South America GIS/Brasil/Future Climate Scenarios/", Future_Climate_Scenarios[[x]], "/85bi50", sep=""), pattern="*3_degrees.tif", full.names=T)
    files <- files[PCA_bioclim]
    files <- lapply(files, raster)
    for(z in 1:length(files)){
      files[[z]]@data@names <- bg_bioclim[[z]]@data@names
    }
    Future_Climate_Data[[x]][[2]] <- stack(files)
  }
}

rm(files, x)
