
###############################################################################################################
###############################################################################################################
###############################################################################################################
######################### Script by Peter Moonlight, Tiina Sarkinen et al. 2017 ###############################
###############################################################################################################
###############################################################################################################
###############################################################################################################




### Prepare the working space
### -------------------------

setwd("E:/000_Modeling_Working_Directory_000")
source(paste(getwd(), "/01_Scripts/Minor_Modules/1___Prepare_Working_Space.R", sep=""))

require(raster)

### Read in the nordeste mask
### -------------------------

nordeste <- raster("000_GIS_LAYERS/nordeste.tif")
nordeste <- raster("000_GIS_LAYERS/nordeste.grd")


### Read in the Models
### ------------------

species <- gsub(".tif$", "", list.files("03_Modelling/12a_Thresholded_Models_Masked_Nordeste/", pattern="*.tif$", full.names=F))
species.level.dir <- lapply(1:length(species), function(x){raster(paste("03_Modelling/12a_Thresholded_Models_Masked_Nordeste/", as.character(species)[[x]], ".tif", sep=""))})
SDM.b <- stack(species.level.dir)


### Determine the site localities
### -----------------------------

sites <- which(!is.na(SDM.b[[1]][1:length(SDM.b[[1]])]))
sites <- as.data.frame(xyFromCell(SDM.b, sites))
rownames(sites) <- which(!is.na(SDM.b[[1]][1:length(SDM.b[[1]])]))

plot(nordeste)
plot(SpatialPoints(sites[,1:2]), add=T)

### Expand the sites matrix to house the results
### --------------------------------------------
sites[,3:(length(species)+2)] <- NA
colnames(sites)[3:(length(species)+2)] <- species


### Populate the results
### --------------------


writeLines(paste("...0%", sep=""))
for(x in 1:nlayers(SDM.b)){
  percents <- round(nlayers(SDM.b)*seq(1, 0, -0.01))
  if(x %in% percents){
    writeLines(paste("...", round(x/nlayers(SDM.b)*100, 0), "% at ", Sys.time(), sep=""))
  }
  sites[,(x+2)] <- extract(SDM.b[[x]], SpatialPoints(sites[,1:2]))
}

site_locs <- sites[,1:2]
sites <- sites[,-c(1:2)]
sites <- t(sites)


write.csv(sites, "07_Virtual_Plot_Checklists/cell_results.csv")





# a) Load all of the required spreadsheets into R
# Obs: I used the read.csv function in here, but there is a limit to the number of rows and columns that csv files can have. If data is being lost because of this, please load tables as text files

spp <- read.csv("07_Virtual_Plot_Checklists/cell_results.csv", sep=",", head=TRUE, row.names=1)
#colnames(spp) <- gsub("\\.", "", colnames(spp))
colnames(spp) <- gsub("X", "", colnames(spp))

sums <- colSums(spp)
#spp <- spp[,-which(sums==0)]


sppxsites <- as.data.frame(matrix(ncol=2))[-1,]
colnames(sppxsites) <- c("CellID", "SppID")

for(x in 1:length(colnames(spp))){
  index <- which(spp[,x] == 1)
  
  temp <- as.data.frame(matrix(ncol=2), nrow=length(index))
  colnames(temp) <- c("CellID", "SppID")
  temp[1:length(index),1] <- colnames(spp)[[x]]
  temp[,2] <- rownames(spp)[index]
  
  sppxsites <- rbind(sppxsites, temp)
}


write.csv(sppxsites, "07_Virtual_Plot_Checklists/sppxcells.csv")

#Spp by sites/grid-cells matrix (a correspondance matrix with two columns: sites and species)
sppxsites <- read.csv("07_Virtual_Plot_Checklists/sppxsites.csv", sep=",", head=TRUE)[,-1]
dim(sppxsites)
head(sppxsites)


#Checking the sites/grid-cells
sites_miss <- sites[which(!sites$AreaID %in% sppxsites$AreaID),]
dim(sites_miss)
sites_miss
#Should be none


#the opposite
sites_miss2 <- sppxsites[which(!sppxsites$AreaID %in% sites$AreaID),]
dim(sites_miss2)
#should also be none


#Now checking species
spp_miss <- spp[which(!rownames(spp) %in% sppxsites$SppID),]
dim(spp_miss)
#Should be none

#the opposite
spp_miss2 <- sppxsites[which(!sppxsites$SppID %in% rownames(spp)),]
dim(spp_miss2)
#should also be none

# b) Build a species by site/grid-cell matrix (rows with species names and columns with sites/grid-cells names)
#In here we'll create a species occurrence matrix (presence X abscence) through a loop function. Depending on the size, this might take a while to run.

#In case you have a difference in number of species included in each one of your matrices, follow the object created bellow instead of using "sppxsites"
#sppxareas_trim <- sppxareas[which(sppxsites$AreaID %in% sites$AreaID),] #Making sure that only species in both matrices are being included 
#length(sppxareas_trim$AreaID)

sites_sub <- unique(sppxsites$AreaID)
spp_sub <- unique(sppxsites$SppID)
spp_sub2 <- spp[which(rownames(spp) %in% sppxsites$SppID),]

spp_commat <- matrix(0,length(sites_sub),length(spp_sub))
for (i in 1:nrow(spp_commat)){
  temp_sites <- sppxsites[which(sppxsites$AreaID==sites_sub[i]),]
  spp_commat[i,which(spp_sub%in%temp_sites$SppID)] <- 1
  print(i)
}

rownames(spp_commat) <- as.character(sites$AreaCode[match(sites_sub,sites$AreaID)])
colnames(spp_commat) <- as.character(rownames(spp)[match(spp_sub,rownames(spp))])
dim(spp_commat)
spp_commat[1:10,1:10]

# c) Building the first clusters and identifying the rogue sites/grid-cells by suing the roguenarok approach (external to R)
#The clusters are built based on a pairwise distance matrix. There is a plethora of metrics that can be used to construct such matrix and the method needs to be chosen according to
#the analyses' objectives. In here, the simpson distance metric was used since it is the best distance for clustering according to species turnover (one of the components of beta diversity).

#Installing required packages
install.packages("recluster")
library(recluster)
#Building the clusters
cluster1 <- recluster.cons(spp_commat, tr=100, p=0.5) #100 clusters should be enough, but feel free to make more clusters if necessary
#Extracting the consensus tree
constree1 <- cluster1$cons
write.tree(constree1,"07_Virtual_Plot_Checklists/100trees_cons.tre")
plot(constree1) # This is gonna look ugly in R, especially if working with a lot of data. However, it's a good way of checking for the cluster's overall topology
#We'll need all of the clusters used to create the consensus tree. The following lines are a loop function to save all of the clusters created into a list that can be read by other softwares
write.tree(cluster1$trees[[1]],"07_Virtual_Plot_Checklists/100trees.tre")
for (i in 2:length(cluster1$trees)){
  write.tree(cluster1$trees[[i]],"07_Virtual_Plot_Checklists/100trees.tre",append=T)
}

#Due to the number of sites/grid-cells being included, it is very commom for the cluster made above to have no resolution at all,
#especially if you change the threshold for the consensus tree. The next section uses the Roguenarok approach to determine
#which sites are constantly changing position from one group to the other among the clusters produced above.
#These will have to be removed from the dataset in order for the clustering approach to work.

#The site is the following "http://rnr.h-its.org/". All that needs to be done is upload the clusters under the "Bootstrap Tree Set" in the website and then give your e-mail address and press upload.
#The site will return a table of sites and values that can be saved/copied elsewhere. We recomend the removal of all sites that have a value more than 0.
#Create a table with the sites/grid-cells that are to be removed ("rogue sites" from now on) and load that into R. 

# d) Removing sites/grid-cells that are interfering with cluster resolution ("politomies")

#Reading the table listing the rogue sites that will be removed from the dataset.
outliers = read.table("07_Virtual_Plot_Checklists/rogue_taxa_ntt10.txt", header = F, sep = "")
spp_commat_norogue= spp_commat[-which(rownames(spp_commat)%in%as.character(outliers[,3])),]
dim(spp_commat_norogue)

#Removing spp with no occurrences
spp_commat_norogue_trim <- spp_commat_norogue[,which(!colSums(spp_commat_norogue) == 0)]
dim(spp_commat_norogue_trim)

#Always important to see if there is a pattern behind the removed sites/grid-cells. Firstly, it is important to look at
#species richness numbers in order to see if these sites are species poor.
spp_commat_rogue <- spp_commat[which(rownames(spp_commat)%in%outliers[,3]),]
dim(spp_commat_rogue)
spp_commat_rogue <- spp_commat_rogue[,which(colSums(spp_commat_rogue) > 0)]
dim(spp_commat_rogue)
spp_rich_rogue <- rowSums(spp_commat_rogue)
sort(spp_rich_rogue)

#It is also a good idea to see where these sites are falling in the map.
rogue_sites <- subset(sites, sites$AreaCode %in% rownames(spp_commat_rogue))
dim(rogue_sites)
plot(SDM.b[[1]])
points(rogue_sites$Long10,rogue_sites$Lat10)

#By removing sites, some species will automaticly no longer have any occurrence on the dataset and therefore needs to be removed.

# e) Building the second cluster and identifying the main groups in the cluster

cluster_norogue <- recluster.cons(spp_commat_norogue_trim, tr=100, p=0.5)
constree_norogue <- cluster_norogue$cons
write.tree(constree_norogue, "07_Virtual_Plot_Checklists/constree_norogue.tre")
#This approach will usually dramatically increases the resolution of the clusters, but it might not. 
#This will be dealt with later on.

#Obs: the following lines can be used to change sites/grid-cells names if necessary. In the case bellow, columns with extra
#information on the sites were used in order to name the sites in a way that would make group identification much easier.

#sites_norogue <- sites[which(sites$AreaCode %in% rownames(spp_commat_norogue_trim)),]
#dim(sites_norogue)
#sites_norogue$AreaCode <- as.character(sites_norogue$AreaCode)
#sites_metada_temp <- sites_norogue[match(constree_norogue$tip.label,sites_norogue$AreaCode),]
#dim(sites_metada_temp)
#colnames(sites_metada_temp)
#sites_metada_temp[1:10,1:10]
#sites_metada_temp$AreaCode <- as.character(sites_metada_temp$AreaCode)
#sites_metada_temp$Domain <- as.character(sites_metada_temp$Domain)
#let's also include EcoUnit this time
#sites_metada_temp$EcoUnit <- as.character(sites_metada_temp$EcoUnit)
#sites_metada_temp$Vegetation.type <- as.character(sites_metada_temp$Vegetation.type)

#new_tips <- vector("character",length=nrow(sites_metada_temp))
#for (i in 1:length(new_tips)){
#  new_tips[i] <- paste(sites_metada_temp[i,c(3,14,16,25)],collapse="_")
#}
#length(new_tips)
#new_tips[1:6]
#constree_norogue_fulltips <- constree_norogue
#constree_norogue_fulltips$tip.label <- new_tips
#write.tree(constree_norogue_fulltips, "constree_norogue_fulltips.tre")

###Creating vector assigning each site/grid-cell to a group/biome. This vector will be called cluster_membership
install.packages("phytools")
library(phytools)

#The first step is to make R acknowledge the presence of politomies in the cluster. As said earlier, the changes of getting
#a cluster without "politomies" are very slim, especially if you are working with Savannas (in my experience).
constree_norogue_nopoly <- di2multi(constree_norogue)

sort(constree_norogue_nopoly$tip.label) #Just to make sure everything is ok with the labels.

#Atlantic forest 1
atlantic_node <- findMRCA(constree_norogue_nopoly, c("AtlBA080","AtlBA058")) #Find the MRCA of a given group
atlantic_tips <- getDescendants(constree_norogue_nopoly,atlantic_node) # Get the identity of all nodes and tips stemming from the above MRCA
atlantic_tips <- na.omit(constree_norogue_nopoly$tip.label[atlantic_tips]) # Remove the nodes and keeps the tips
length(atlantic_tips) # Important to check the lengh of each one of the groups.

#Atlantic forest 2
atlantic_node_2 <- findMRCA(constree_norogue_nopoly, c("AtlSE003","AtlAL020")) #Find the MRCA of a given group
atlantic_tips_2 <- getDescendants(constree_norogue_nopoly,atlantic_node_2) # Get the identity of all nodes and tips stemming from the above MRCA
atlantic_tips_2 <- na.omit(constree_norogue_nopoly$tip.label[atlantic_tips_2]) # Remove the nodes and keeps the tips
length(atlantic_tips_2) # Important to check the lengh of each one of the groups.

#Atlantic forest 3
atlantic_node_3 <- findMRCA(constree_norogue_nopoly, c("CerMG015","AtlMG049")) #Find the MRCA of a given group
atlantic_tips_3 <- getDescendants(constree_norogue_nopoly,atlantic_node_3) # Get the identity of all nodes and tips stemming from the above MRCA
atlantic_tips_3 <- na.omit(constree_norogue_nopoly$tip.label[atlantic_tips_3]) # Remove the nodes and keeps the tips
length(atlantic_tips_3) # Important to check the lengh of each one of the groups.

#amazon 1
amazon_tips <- "AmzMA008"
length(amazon_tips)

#mix
mix_node <- findMRCA(constree_norogue_nopoly, c("AtlMA004","CerPI026"))
mix_tips <- getDescendants(constree_norogue_nopoly, mix_node)
mix_tips <- na.omit(constree_norogue_nopoly$tip.label[mix_tips])
length(mix_tips)

#mix 2
mix_node_2 <- findMRCA(constree_norogue_nopoly, c("AtlBA094","CerMG039"))
mix_tips_2 <- getDescendants(constree_norogue_nopoly, mix_node_2)
mix_tips_2 <- na.omit(constree_norogue_nopoly$tip.label[mix_tips_2])
length(mix_tips_2)

#caatinga
caatinga_node <- findMRCA(constree_norogue_nopoly, c("AtlCE030","CaaCE034"))
caatinga_tips <- getDescendants(constree_norogue_nopoly, caatinga_node)
caatinga_tips <- na.omit(constree_norogue_nopoly$tip.label[caatinga_tips])
length(caatinga_tips)

#cerrado
cerrado_node <- findMRCA(constree_norogue_nopoly, c("CerMG152","CerPI049"))
cerrado_tips <- getDescendants(constree_norogue_nopoly, cerrado_node)
cerrado_tips <- na.omit(constree_norogue_nopoly$tip.label[cerrado_tips])
length(cerrado_tips)

#Remember to sum all of the lenghts obtained here. The sum must be equal to the total number of tips. Pay special attential to groups with big politomies. 
length(constree_norogue_nopoly$tip.label)
length(amazon_tips) + length(caatinga_tips) + length(mix_tips_2) + length(mix_tips) + length(atlantic_tips_3) + length(atlantic_tips_2) + length(atlantic_tips) + length(cerrado_tips)

#Preparing the cluster_membership vector
#all_tips
all_tips <- c(amazon_tips, caatinga_tips, mix_tips_2, mix_tips, atlantic_tips_3, atlantic_tips_2, atlantic_tips, cerrado_tips)
cluster_membership <- vector("character",length(all_tips))
cluster_membership[which(all_tips%in% amazon_tips)] <- "AMZ"
cluster_membership[which(all_tips%in% caatinga_tips)] <- "CAA"
cluster_membership[which(all_tips%in% mix_tips_2)] <- "MIX_2"
cluster_membership[which(all_tips%in% mix_tips)] <- "MIX_1"
cluster_membership[which(all_tips%in% atlantic_tips_3)] <- "ATL_3"
cluster_membership[which(all_tips%in% atlantic_tips_2)] <- "ATL_2"
cluster_membership[which(all_tips%in% atlantic_tips)] <- "ATL_1"
cluster_membership[which(all_tips%in% cerrado_tips)] <- "CER"
length(cluster_membership)
class(cluster_membership)

#The cluster_membership vector needs to be added to the sites matrix after removing the rogue sites.
sites_norogue <- sites[-which(sites$AreaCode%in%outliers[,3]),]

#Matching the cluster_membership vector with the sites
cluster_membership <- cluster_membership[match(sites_norogue$AreaCode,all_tips)] # The two objects need to be matched up first so they'll be bound correctly to each other

#Binding cluster membership to areas_lowtrop_1000m_notemperate_norogue
sites_norogue_membership <- cbind(sites_norogue, cluster_membership)
unique(sites_norogue_membership$cluster_membership) # Making sure that there are no groups/biomes missing
dim(sites_norogue_membership)
head(sites_norogue_membership)
unique(sites_norogue_membership$cluster_membership)

# g) Using the vector created above to map the results

#The codes in here were written without using ggmap and ggplot2, but feel free to use those or others. #OldSchool 
install.packages("maps")
library(maps)
install.packages("maptools")
library(maptools)


nordeste <- raster("000_GIS_LAYERS/nordeste.tif")
plot(nordeste)
map.axes()
#Amazon group
points(sites_norogue_membership$Long10[which(sites_norogue_membership$cluster_membership == "AMZ")]
       ,sites_norogue_membership$Lat10[which(sites_norogue_membership$cluster_membership == "AMZ")],pch=24,col=rgb
       (t(col2rgb("chartreuse4"))/255,alpha=1), bg=rgb(t(col2rgb("chartreuse4"))/255))
#Cerrado group
points(sites_norogue_membership$Long10[which(sites_norogue_membership$cluster_membership == "CER")]
       ,sites_norogue_membership$Lat10[which(sites_norogue_membership$cluster_membership == "CER")],pch="O",col=rgb
       (t(col2rgb("gray53"))/255,alpha=1), bg=rgb(t(col2rgb("gray53"))/255, alpha=1))
#Mix Group 1
points(sites_norogue_membership$Long10[which(sites_norogue_membership$cluster_membership == "MIX_1")]
       ,sites_norogue_membership$Lat10[which(sites_norogue_membership$cluster_membership == "MIX_1")],pch=15,col=rgb
       (t(col2rgb("blue"))/255,alpha=1), bg=rgb(t(col2rgb("blue"))/255, alpha=1))
#Mix Group 1
points(sites_norogue_membership$Long10[which(sites_norogue_membership$cluster_membership == "MIX_2")]
       ,sites_norogue_membership$Lat10[which(sites_norogue_membership$cluster_membership == "MIX_2")],pch=15,col=rgb
       (t(col2rgb("turquoise"))/255,alpha=1), bg=rgb(t(col2rgb("turquoise"))/255, alpha=1))
#Atlantic Group 3
points(sites_norogue_membership$Long10[which(sites_norogue_membership$cluster_membership == "ATL_3")]
       ,sites_norogue_membership$Lat10[which(sites_norogue_membership$cluster_membership == "ATL_3")],pch=19,col=rgb
       (t(col2rgb("saddlebrown"))/255,alpha=1), bg=rgb(t(col2rgb("saddlebrown"))/255))
#Atlantic Group 2
points(sites_norogue_membership$Long10[which(sites_norogue_membership$cluster_membership == "ATL_2")]
       ,sites_norogue_membership$Lat10[which(sites_norogue_membership$cluster_membership == "ATL_2")],pch=19,col=rgb
       (t(col2rgb("red"))/255,alpha=1), bg=rgb(t(col2rgb("red"))/255))
#Atlantic Group 1
points(sites_norogue_membership$Long10[which(sites_norogue_membership$cluster_membership == "ATL_1")]
       ,sites_norogue_membership$Lat10[which(sites_norogue_membership$cluster_membership == "ATL_1")],pch=19,col=rgb
       (t(col2rgb("tomato"))/255,alpha=1), bg=rgb(t(col2rgb("tomato"))/255))
#Caatinga group
points(sites_norogue_membership$Long10[which(sites_norogue_membership$cluster_membership == "CAA")]
       ,sites_norogue_membership$Lat10[which(sites_norogue_membership$cluster_membership == "CAA")],pch=25,col=rgb
       (t(col2rgb("purple"))/255,alpha=1), bg=rgb(t(col2rgb("purple"))/255, alpha=1))

# h) Performing a silhouette analysis in the cluster in order to identify sites/grid-cells that are floristically closer to another group/biome

#Cluster, like many other multivariate methods, by coercing multi-dimensional data into bi-dimensional data. the main consequence
#that needs to be addressed is that some sites/grid-cells might actually be floristically closer to a different group/biome
#to which the site/grid-cell has been assigned to.
#The silhouette analysis will identify sites/grid-cells that have been assigned to a group, but are actually more similar to a different group.
#This is done by using the pairwise distance matrix used to create the cluster. The function will then estimate the means
#of each one of the groups that have been identified and will compare each distance measure to these means. According to this
#comparison, the funcion will then determine to which group a given site-grid cell would have better been placed.
#These results will be converted into a second cluster_membership vector and then will be added to the sites_norogue_membership dataframe. 
#In overrall terms, sites that have been correctly assigned to a group will have a negative silhouette value and sites that have been wrongly assigned
#to a group will have negative silhouette values.

cluster_groups <- as.integer(sites_norogue_membership$cluster_membership) # the function requires a interger vector with the groupings
plot(match(rownames(spp_commat_norogue_trim),sites_norogue_membership$AreaCode)) # checking if they have been correctly matched
cluster_groups <- cluster_groups[match(rownames(spp_commat_norogue_trim),sites_norogue_membership$AreaCode)] # These two need to be matched up
cluster_distmat <- recluster.dist(spp_commat_norogue_trim) # Obtaining the pairwise distance matrix used to produce the cluster/
install.packages("cluster") # this packages is a part of basic R, so it's usually not necessary to install it. Only do it if you're having problems loading its library
library(cluster) # the silhouette function belongs to this package
#Getting the values
silh_output <- silhouette(cluster_groups,cluster_distmat)
plot(silh_output,col=c("chartreuse4","blue","gray53","black","saddlebrown")) #This is a very interesting graph to look at. It gives an idea of how well or poorly defined the groups are. In here, I used the same colors I used in the maps
#Making a new vector and binding that into the main sites/grid-cells dataframe
silh_output_values <- as.data.frame(cbind(silh_output[,1],silh_output[,2],silh_output[,3]))
silh_output_values$original_group <- cluster_groups
silh_output_values$AreaCode <- rownames(spp_commat_norogue_trim)
misclass_sites <- subset(silh_output_values, silh_output_values$V3 < 0)
dim(misclass_sites)
misclass_sites_full <- subset(sites_norogue_membership, sites_norogue_membership$AreaCode %in% misclass_sites$AreaCode) # You can do the opposite and make a dataframe containing only that have been correctly assigned to a group.
dim(misclass_sites_full)
colnames(misclass_sites_full)
misclass_sites_full$AreaCode
write.csv(misclass_sites_full, "07_Virtual_Plot_Checklists/misclass_sites_full.csv")

# j) Mapping the sites/grid-cells identified by the silhouette analysis

plot(nordeste)
map.axes()
#Atlantic group
points(misclass_sites_full$Long10[which(misclass_sites_full$cluster_membership == "AMZ")]
       ,misclass_sites_full$Lat10[which(misclass_sites_full$cluster_membership == "AMZ")],pch=24,col=rgb
       (t(col2rgb("chartreuse4"))/255,alpha=1), bg=rgb(t(col2rgb("chartreuse4"))/255))
#Cerrado group
points(misclass_sites_full$Long10[which(misclass_sites_full$cluster_membership == "CER")]
       ,misclass_sites_full$Lat10[which(misclass_sites_full$cluster_membership == "CER")],pch="O",col=rgb
       (t(col2rgb("gray53"))/255,alpha=1), bg=rgb(t(col2rgb("gray53"))/255, alpha=1))
#Mix 1 group
points(misclass_sites_full$Long10[which(misclass_sites_full$cluster_membership == "MIX_1")]
       ,misclass_sites_full$Lat10[which(misclass_sites_full$cluster_membership == "MIX_1")],pch=15,col=rgb
       (t(col2rgb("blue"))/255,alpha=1), bg=rgb(t(col2rgb("blue"))/255, alpha=1))
#Mix 2 group
points(misclass_sites_full$Long10[which(misclass_sites_full$cluster_membership == "MIX_2")]
       ,misclass_sites_full$Lat10[which(misclass_sites_full$cluster_membership == "MIX_2")],pch=15,col=rgb
       (t(col2rgb("turquoise"))/255,alpha=1), bg=rgb(t(col2rgb("turquoise"))/255, alpha=1))
#Atl 3 group
points(misclass_sites_full$Long10[which(misclass_sites_full$cluster_membership == "ATL_3")]
       ,misclass_sites_full$Lat10[which(misclass_sites_full$cluster_membership == "ATL_3")],pch=19,col=rgb
       (t(col2rgb("red"))/255,alpha=1), bg=rgb(t(col2rgb("red"))/255))
#Atl 2 group
points(misclass_sites_full$Long10[which(misclass_sites_full$cluster_membership == "ATL_2")]
       ,misclass_sites_full$Lat10[which(misclass_sites_full$cluster_membership == "ATL_2")],pch=19,col=rgb
       (t(col2rgb("purple"))/255,alpha=1), bg=rgb(t(col2rgb("purple"))/255))
#Atl 1 group
points(misclass_sites_full$Long10[which(misclass_sites_full$cluster_membership == "ATL_1")]
       ,misclass_sites_full$Lat10[which(misclass_sites_full$cluster_membership == "ATL_1")],pch=19,col=rgb
       (t(col2rgb("saddlebrown"))/255,alpha=1), bg=rgb(t(col2rgb("saddlebrown"))/255))
#Caatinga group
points(misclass_sites_full$Long10[which(misclass_sites_full$cluster_membership == "CAA")]
       ,misclass_sites_full$Lat10[which(misclass_sites_full$cluster_membership == "CAA")],pch=25,col=rgb
       (t(col2rgb("purple"))/255,alpha=1), bg=rgb(t(col2rgb("purple"))/255, alpha=1))

######
#More details on the information presented in here can be found at: https://ourcodingclub.github.io/2017/03/21/data-clustering.html


