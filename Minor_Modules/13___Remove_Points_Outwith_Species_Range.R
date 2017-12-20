###############################################################################################################
  ###############################################################################################################
    ###############################################################################################################
      ######################### Script by Peter Moonlight, Tiina Sarkinen et al 2017 ################################
    ###############################################################################################################
  ###############################################################################################################
###############################################################################################################

dir.create("03_Modelling/04_Species_To_Model_Distribution_Data/Cleaned_5", showWarnings = F)




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

mask <- raster("000_GIS_LAYERS/brazil_mask_minus_centroids.tif")
FunR<-function(r){
  ext<-raster(extent(mask), nrow=nrow(mask), ncol=ncol(mask))
  crs(ext)<-crs(r)
  D<-rasterize(r,ext,field=1,update=T)
}

# Turn each polygon into a raster
StateRasters <- lapply(1:length(BRA_ADM_States), function(x) FunR(BRA_ADM[x,]))

# Read in the Reflora Data
Reflora_Species <- read.csv("000_Species/Angiospermas_Flora do Brasil 2020_25 set 2017_Working.csv")



for(x in 1:length(species)){
  if(!dir.exists(paste("03_Modelling/04_Species_To_Model_Distribution_Data/Cleaned_5/", species[[x]], sep=""))){
    dir.create(paste("03_Modelling/04_Species_To_Model_Distribution_Data/Cleaned_5/", species[[x]], sep=""), showWarnings = F)
    
    distribution_data <- read.csv(file=paste("03_Modelling/04_Species_To_Model_Distribution_Data/Cleaned_4/", species[[x]], ".csv", sep=""), header=T, stringsAsFactors = F)[,-1]
    dims <- dim(distribution_data)[1]
    if(length(distribution_data[,1])>0){
      
      Reflora_Index  <- which(Reflora_Species$Species == species[[x]])
      Reflora_States <- which(Reflora_Species[Reflora_Index,20:46]>0)
      Reflora_States <- (colnames(Reflora_Species[20:46])[Reflora_States])
      
      Reflora_State_List <- StateRasters[which(BRA_ADM_States %in% Reflora_States)]
      
      # If there are no recorded states for the species, remove that species
      if(length(Reflora_State_List) == 0){
        distribution_data <- distribution_data[-c(1:length(distribution_data[,1])),]
      }
    
      # If there is 1 recorded state for the species, remove records outwith that state
      if(length(Reflora_State_List) == 1){
        index <- which(is.na(extract(Reflora_State_List[[1]], SpatialPoints(distribution_data[,2:1]))))
        if(length(index)>0){
          distribution_data <- distribution_data[-index,]
        }
      }
    
      # If there are many recorded states for the species, remove records outwith those states
      if(length(Reflora_State_List) > 1){
        Reflora_State_List$fun <- max
        Reflora_State_List <- do.call(mosaic, Reflora_State_List)
        index <- which(is.na(extract(Reflora_State_List, SpatialPoints(distribution_data[,2:1]))))
        if(length(index)>0){
          distribution_data <- distribution_data[-index,]
        }
      }

    }
    if(dim(distribution_data)[1]<dims){writeLines(paste(species[x], "now has", dims-length(distribution_data[,1]),"fewer points..."))}
    if(length(distribution_data[,1])>0){
      write.csv(distribution_data, file=paste("03_Modelling/04_Species_To_Model_Distribution_Data/Cleaned_5/", species[[x]], ".csv", sep=""))
    }
  }
}
