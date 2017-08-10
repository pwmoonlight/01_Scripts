###############################################################################################################
  ###############################################################################################################
    ###############################################################################################################
      ######################### Script by Peter Moonlight, Tiina Sarkinen et al 2017 ################################
    ###############################################################################################################
  ###############################################################################################################
###############################################################################################################


dir.create("12_Future_Projections", showWarnings=F)
lapply(Future_Climate_Scenarios, function(x){dir.create(paste("12_Future_Projections/", x, sep=""), showWarnings = F)})

writeLines(paste("...Predicting Model Using Future Climate Variables"))

chosen_models_45 <- c(1, 6, 8, 9, 12, 13)
chosen_models_85 <- c(6, 8, 11, 12, 13)


for(z in 1:length(Future_Climate_Data)){
  if(z %in% chosen_models_45){
    writeLines(paste("   ...Working With ", Future_Climate_Scenarios[[z]], "--- 45"))
    if(class(Future_Climate_Data[[z]][[1]])=="RasterStack"){
      predict(model, Future_Climate_Data[[z]][[1]], filename=paste("12_Future_Projections/", Future_Climate_Scenarios[[z]], "/", species[[x]], "_45.tif", sep="")) 
    }
  }
  if(z %in% chosen_models_85){
    writeLines(paste("   ...Working With ", Future_Climate_Scenarios[[z]], "--- 85"))
    if(class(Future_Climate_Data[[z]][[2]])=="RasterStack"){
      predict(model, Future_Climate_Data[[z]][[2]], filename=paste("12_Future_Projections/", Future_Climate_Scenarios[[z]], "/", species[[x]], "_85.tif", sep="")) 
    }
  }
}
