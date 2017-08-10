
###############################################################################################################
  ###############################################################################################################
    ###############################################################################################################
      ######################### Script by Peter Moonlight, Tiina Sarkinen et al 2017 ################################
    ###############################################################################################################
###############################################################################################################

dir.create("98_PCA", showWarnings=F)

writeLines(paste("...Reading in distribution data"))

species_distribution_data <- list.files("09_Species_To_Model_Scale_Corrected_Distribution_Data", pattern="*.csv", full.names = T)
species_distribution_data <- lapply(species_distribution_data, read.csv)
species_distribution_data <- do.call(rbind, species_distribution_data)

writeLines(paste("...Reading in background data"))

bg <- lapply(list.files(path="Y:/South America GIS/Brasil/Brazil_Masked_GIS_Layers", pattern="*.tif$", full.names = T), raster)
bg <- stack(bg)

writeLines(paste("...Extracting background data values"))

values <- extract(bg, species_distribution_data[,3:2])
values <- values[complete.cases(values),]
write.csv(values, "98_PCA/Background_variables_for_collectionpoints.csv", quote=F, row.names=F)

writeLines(paste("...Standardising background data values"))
values_std <- scale(values)
write.csv(values_std, "98_PCA/Background_variables_for_collectionpoints_standardised.csv", quote=F, row.names=F)

writeLines(paste("...Finding Correlations"))
correlation <- cor(values)
write.csv(correlation, "98_PCA/correlations.csv", quote=F)

writeLines(paste("...Running PCA"))
Bioclim_pca <- princomp(values, cor=TRUE)
biplot(Bioclim_pca)
plot(Bioclim_pca)
pca_loadings <- loadings(Bioclim_pca)
write.csv(pca_loadings, "98_PCA/pca_loadings.csv")
pca_summary <- summary(Bioclim_pca)
vars <- pca_summary$sdev^2
vars <- vars/sum(vars)
results <- rbind("Standard deviation" = pca_summary$sdev, "Proportion of variance" = vars,"Cumulative Proportion" = cumsum(vars))
write.csv(results, "98_PCA/pca_summary.csv", quote=F, row.names=T)


rm(species_distribution_data, values, values_std, correlation, Bioclim_pca, pca_loadings, pca_summary, vars, results, bg)

bg <- lapply(list.files(path="Y:/South America GIS/Brasil/Brazil_Masked_GIS_Layers", pattern="*.tif$", full.names = T), raster)

writeLines(paste("\n\n\n####################################################", sep=""))
writeLines(paste("###              THE PCA HAS NOW RUN             ###", sep=""))
writeLines(paste("### The Data is available in the '04_PCA' Folder ###", sep=""))
writeLines(paste("###     Manually Remove Extraneous Variables     ###", sep=""))
writeLines(paste("####################################################", sep=""))
writeLines(paste("\n"))




