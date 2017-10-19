###############################################################################################################
  ###############################################################################################################
    ###############################################################################################################
      ######################### Script by Peter Moonlight, Tiina Sarkinen et al 2017 ################################
    ###############################################################################################################
  ###############################################################################################################
###############################################################################################################

dir.create("03_Modelling/04_Species_To_Model_Distribution_Data/Cleaned_2", showWarnings = F)


### -------------------------------
### Read in Polygons for each state

require(rgdal)
BRA_ADM <- list.files("000_GIS_LAYERS/BRA_Adm_2/", pattern="[.]shp$", full.names=T)
BRA_ADM <- readOGR(BRA_ADM[1])

BRA_ADM_States <- toupper(iconv(BRA_ADM@data$NM_ESTADO,, to='ASCII//TRANSLIT'))
BRA_ADM_States[26] <- "RO"

# Change the state names into standard two letter codes
States <- read.table("000_Species/States.txt", header=F, sep="\t", fill=T, quote="")
for(x in 1:length(BRA_ADM_States)){BRA_ADM_States[x] <- as.character(States[which(States[,1] == BRA_ADM_States[x]),2])}



### ----------------------------------------------------------------
### Load in a function for turning the state data into a raster file

FunR<-function(r){
  ext<-raster(extent(r), nrow=(extent(r)[4]-extent(r)[3])/0.01, ncol=(extent(r)[2]-extent(r)[1])/0.01)
  crs(ext)<-crs(r)
  D<-rasterize(r,ext,field=1,update=T)
}

# Turn each polygon into a raster
StateRasters <- lapply(1:length(BRA_ADM_States), function(x) FunR(BRA_ADM[x,]))

# Merge multiple Rasters for individual States
for(x in length(BRA_ADM_States):1){
  index <- which(BRA_ADM_States[-x] == BRA_ADM_States[x])
  if(length(index)==1){
    StateRasters[[index]] <- merge(StateRasters[[x]], StateRasters[[index]], tolerance=0.5)
    BRA_ADM_States <- BRA_ADM_States[-x]
    StateRasters <- StateRasters[-x]
  }
  if(length(index)>1){
    for(y in 1:length(StateRasters[index])){
      StateRasters[[index[y]]] <- merge(StateRasters[[x]], StateRasters[[index[y]]], tolerance=0.5)
    }
    BRA_ADM_States <- BRA_ADM_States[-x]
    StateRasters <- StateRasters[-x]
  }
}


### -----------------------------------
### Check Species Against These Rasters

for(x in 1:length(species)){
  distribution_data <- read.csv(file=paste("03_Modelling/04_Species_To_Model_Distribution_Data/Cleaned_1/", species[[x]], ".csv", sep=""), header=T, stringsAsFactors = F)[,-1]
  dims <- dim(distribution_data)[1]
  if(length(distribution_data[,1])>0){
    for(y in length(distribution_data[,1]):1){
      if(!is.na(distribution_data[y,10])>0){
        if(is.na(distribution_data[y,1])){distribution_data <- distribution_data[-y,]}
        if(!is.na(distribution_data[y,1])){
          if(distribution_data[y,10] %in% BRA_ADM_States){
            if(nchar(distribution_data[y,10])>0){
              index <- extract(StateRasters[[which(BRA_ADM_States==distribution_data[y,10])]],SpatialPoints(distribution_data[y,2:1]))
              if(is.na(index)){distribution_data <- distribution_data[-y,]}
              if(!is.na(index)){
                if(index == 0){distribution_data <- distribution_data[-y,]}
              }
            }
          }
        }
      }
    }
    if(dim(distribution_data)[1]<dims){writeLines(paste(species[x], "now has", dims-dim(distribution_data[1]),"fewer points..."))}
    write.csv(distribution_data, file=paste("03_Modelling/04_Species_To_Model_Distribution_Data/Cleaned_2/", species[[x]], ".csv", sep=""))
  }
}

rm(dims, distribution_data, States, BRA_ADM_States, StateRasters, BRA_ADM)