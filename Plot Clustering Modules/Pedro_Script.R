### Routine for Cluster analysis
#Authors: Pedro Miranda, Kyle Dexter and Danilo Neves

### Main objectives

# a) Load all of the required spreadsheets into R
# b) Build a species by site/grid-cell matrix (rows with species names and columns with sites/grid-cells names)
# c) Building the first cluster and identifying the rogue sites/grid-cells by suing the roguenarok approach (external to R)
# d) Removing sites/grid-cells that are interfering with cluster resolution ("politomies")
# e) Building the second cluster and identifying the main groups in the cluster
# f) Creating a vector with the groups that the sites/grid-cells have been assigned to (as factor)
# g) Using the vector created above to map the results
# h) Performing a silhouette analysis in the cluster in order to identify sites/grid-cells that are floristically closer to another group/biome
# i) Creating a second vector with the groups that the sites/grid-cells have been assigned to according to the silhuette analysis
# j) Mapping the sites/grid-cells identified by the silhouette analysis

# a) Load all of the required spreadsheets into R
# Obs: I used the read.csv function in here, but there is a limit to the number of rows and columns that csv files can have. If data is being lost because of this, please load tables as text files

#Spp matrix (matrix with all species in the dataset)
spp <- read.csv("Species_ntt10.csv", sep=",", head=TRUE)
#spp <- read.table("Species_ntt10.csv", sep=",", head=TRUE)
head(spp)
dim(spp)

#Sites/grid-cells matrix (matrix with all sites/grid-cells in the dataset)
sites <- read.csv("Areas_ntt10.csv", sep=",", head=TRUE) 
#sites <- read.table("Areas_ntt10.csv", sep=",", head=TRUE) 
dim(areas)
head(areas)

#Spp by sites/grid-cells matrix (a correspondance matrix with two columns: sites and species)
sppxsites <- read.csv("Species-Area_ntt10.csv", sep=",", head=TRUE)
dim(sppxareas)
head(sppxareas)

#Checking if there are any mismatches (species and sites that are present in one matrix, but not in the other)
#Obs: in NeoTropTree, all sites and species have been assigned a number. This number is under the AreaID and SppID columns

#Checking the sites/grid-cells
sites_miss <- sitess[which(!sitess$AreaID %in% sppxsites$AreaID),]
dim(sites_miss)
sites_miss
#Should be none

#the opposite
sites_miss2 <- sppxsites[which(!sppxsites$AreaID %in% sites$AreaID),]
dim(areas_miss2)
#should also be none

#Now checking species
spp_miss <- spp[which(!spp$SppID %in% sppxareas$SppID),]
dim(spp_miss)
#Should be none

#the opposite
spp_miss2 <- sppxareas[which(!sppxareas$SppID %in% spp$SppID),]
dim(spp_miss2)
#should also be none

# b) Build a species by site/grid-cell matrix (rows with species names and columns with sites/grid-cells names)
#In here we'll create a species occurrence matrix (presence X abscence) through a loop function. Depending on the size, this might take a while to run.

#In case you have a difference in number of species included in each one of your matrices, follow the object created bellow instead of using "sppxsites"
#sppxareas_trim <- sppxareas[which(sppxsites$AreaID %in% sites$AreaID),] #Making sure that only species in both matrices are being included 
#length(sppxareas_trim$AreaID)

sites_sub <- unique(sppxsites$AreaID)
spp_sub <- unique(sppxsites$SppID)
spp_sub2 <- spp[which(spp$SppID%in%sppxsites$SppID),]

spp_commat <- matrix(0,length(sites_sub),length(spp_sub))
for (i in 1:nrow(spp_commat)){
  temp_sites <- sppxsites[which(sppxsites$AreaID==sites_sub[i]),]
  spp_commat[i,which(spp_sub%in%temp_sites$SppID)] <- 1
  print(i)
}

rownames(spp_commat) <- as.character(sites$AreaCode[match(sites_sub,areas$AreaID)])
colnames(spp_commat) <- as.character(spp$Species.code[match(spp_sub,spp$SppID)])
dim(spp_commat)
spp_commat[1:6,1:6]

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
plot(constree1) # This is gonna look ugly in R, especially if working with a lot of data. However, it's a good way of checking for the cluster's overall topology
#We'll need all of the clusters used to create the consensus tree. The following lines are a loop function to save all of the clusters created into a list that can be read by other softwares
write.tree(cluster1$trees[[1]],"100trees.tre")
for (i in 2:length(cluster1$trees)){
  write.tree(cluster1$trees[[i]],"100trees.tre",append=T)
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
outliers = read.table("rogue_taxa_ntt10.txt", header = F, sep = "")
spp_commat_norogue= spp_commat[-which(rownames(spp_commat)%in%outliers[,1]),]
dim(spp_commat_norogue)

#Removing spp with no occurrences
spp_commat_norogue_trim <- spp_commat_norogue[,which(!colSums(spp_commat_norogue) == 0)]
dim(spp_commat_norogue_trim)

#Always important to see if there is a pattern behind the removed sites/grid-cells. Firstly, it is important to look at
#species richness numbers in order to see if these sites are species poor.
spp_commat_rogue <- spp_commat[which(rownames(spp_commat)%in%outliers[,1]),]
dim(spp_commat_rogue)
spp_commat_rogue <- spp_commat_rogue[,which(colSums(spp_commat_rogue) > 0)]
dim(spp_commat_rogue)
spp_rich_rogue <- rowSums(spp_commat_rogue)
sort(spp_rich_rogue)

#It is also a good idea to see where these sites are falling in the map.
rogue_sites <- subset(sites, sites$AreaCode %in% rownames(spp_commat_rogue))
dim(rogue_sites)
library(maps)
map(xlim=c(-83,-32),ylim=c(-35,6))
map.axes()
points(rogue_sites$Long10,rogue_sites$Lat10)

#By removing sites, some species will automaticly no longer have any occurrence on the dataset and therefore needs to be removed.

# e) Building the second cluster and identifying the main groups in the cluster

cluster_norogue <- recluster.cons(spp_commat_norogue_trim, tr=100, p=0.5)
constree_norogue <- cluster_norogue$cons
write.tree(constree_norogue, "constree_norogue.tre")
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

#Atlantic forest
atlantic_node <- findMRCA(constree_lowtrop_norogue_trim, c("AtlSP182","AtlPE001")) #Find the MRCA of a given group
atlantic_tips <- getDescendants(constree_lowtrop_norogue_trim,atlantic_node) # Get the identity of all nodes and tips stemming from the above MRCA
atlantic_tips <- na.omit(constree_lowtrop_norogue_trim$tip.label[atlantic_tips]) # Remove the nodes and keeps the tips
length(atlantic_tips) # Important to check the lengh of each one of the groups.

#sdtf
sdtf_node <- findMRCA(constree_lowtrop_norogue_trim, c("CaaPB029","CerMS061"))
sdtf_tips <- getDescendants(constree_lowtrop_norogue_trim, sdtf_node)
sdtf_tips <- na.omit(constree_lowtrop_norogue_trim$tip.label[sdtf_tips])
length(sdtf_tips)

#amazon
amazon_node <- findMRCA(constree_lowtrop_norogue_trim, c("AmzPA158","AmzAM304"))
amazon_tips <- getDescendants(constree_lowtrop_norogue_trim, amazon_node)
amazon_tips <- na.omit(constree_lowtrop_norogue_trim$tip.label[amazon_tips])
length(amazon_tips)

#chaco
chaco_node <- findMRCA(constree_lowtrop_norogue_trim, c("ChaAR026","ChaAR054"))
chaco_tips <- getDescendants(constree_lowtrop_norogue_trim, chaco_node)
chaco_tips <- na.omit(constree_lowtrop_norogue_trim$tip.label[chaco_tips])
length(chaco_tips)

#cerrado
cerrado_node <- findMRCA(constree_lowtrop_norogue_trim, c("CerMG239","AmzRO056"))
cerrado_tips <- getDescendants(constree_lowtrop_norogue_trim, cerrado_node)
cerrado_tips <- na.omit(constree_lowtrop_norogue_trim$tip.label[cerrado_tips])
length(cerrado_tips)

#Remember to sum all of the lenghts obtained here. The sum must be equal to the total number of tips. Pay special attential to groups with big politomies. 

#Preparing the cluster_membership vector
#all_tips
all_tips <- c(atlantic_tips, chaco_tips, amazon_tips, sdtf_tips, cerrado_tips)
cluster_membership <- vector("character",length(all_tips))
cluster_membership[which(all_tips%in% atlantic_tips)] <- "Atlantic Forest"
cluster_membership[which(all_tips%in% chaco_tips)] <- "Chaco"
cluster_membership[which(all_tips%in% amazon_tips)] <- "Amazon Forest"
cluster_membership[which(all_tips%in% sdtf_tips)] <- "SDTF"
cluster_membership[which(all_tips%in% cerrado_tips)] <- "Cerrado"
length(cluster_membership)
class(cluster_membership)

#The cluster_membership vector needs to be added to the sites matrix after removing the rogue sites.
sites_norogue <- sites[-which(sites$AreaCode%in%outliers[,1]),]

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


map(xlim=c(-83,-25),ylim=c(-35,6))
map.axes()
#Atlantic group
points(sites_norogue_membership$Long10[which(sites_norogue_membership$cluster_membership == "Atlantic Forest")]
       ,sites_norogue_membership$Lat10[which(sites_norogue_membership$cluster_membership == "Atlantic Forest")],pch=24,col=rgb
       (t(col2rgb("chartreuse4"))/255,alpha=1), bg=rgb(t(col2rgb("chartreuse4"))/255))
#Cerrado group
points(sites_norogue_membership$Long10[which(sites_norogue_membership$cluster_membership == "Cerrado")]
       ,sites_norogue_membership$Lat10[which(sites_norogue_membership$cluster_membership == "Cerrado")],pch="O",col=rgb
       (t(col2rgb("gray53"))/255,alpha=1), bg=rgb(t(col2rgb("gray53"))/255, alpha=1))
#Amazon group
points(sites_norogue_membership$Long10[which(sites_norogue_membership$cluster_membership == "Amazon Forest")]
       ,sites_norogue_membership$Lat10[which(sites_norogue_membership$cluster_membership == "Amazon Forest")],pch=15,col=rgb
       (t(col2rgb("blue"))/255,alpha=1), bg=rgb(t(col2rgb("blue"))/255, alpha=1))
#SDTF group
points(sites_norogue_membership$Long10[which(sites_norogue_membership$cluster_membership == "SDTF")]
       ,sites_norogue_membership$Lat10[which(sites_norogue_membership$cluster_membership == "SDTF")],pch=19,col=rgb
       (t(col2rgb("saddlebrown"))/255,alpha=1), bg=rgb(t(col2rgb("saddlebrown"))/255))
#Chaco group
points(sites_norogue_membership$Long10[which(sites_norogue_membership$cluster_membership == "Chaco")]
       ,sites_norogue_membership$Lat10[which(sites_norogue_membership$cluster_membership == "Chaco")],pch=25,col=rgb
       (t(col2rgb("black"))/255,alpha=1), bg=rgb(t(col2rgb("black"))/255, alpha=1))

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
write.csv(misclass_sites_full, "misclass_sites_full.csv")

# j) Mapping the sites/grid-cells identified by the silhouette analysis

map(xlim=c(-83,-25),ylim=c(-35,6))
map.axes()
#Atlantic group
points(misclass_sites_full$Long10[which(misclass_sites_full$cluster_membership == "Atlantic Forest")]
       ,misclass_sites_full$Lat10[which(misclass_sites_full$cluster_membership == "Atlantic Forest")],pch=24,col=rgb
       (t(col2rgb("chartreuse4"))/255,alpha=1), bg=rgb(t(col2rgb("chartreuse4"))/255))
#Cerrado group
points(misclass_sites_full$Long10[which(misclass_sites_full$cluster_membership == "Cerrado")]
       ,misclass_sites_full$Lat10[which(misclass_sites_full$cluster_membership == "Cerrado")],pch="O",col=rgb
       (t(col2rgb("gray53"))/255,alpha=1), bg=rgb(t(col2rgb("gray53"))/255, alpha=1))
#Amazon group
points(misclass_sites_full$Long10[which(misclass_sites_full$cluster_membership == "Amazon Forest")]
       ,misclass_sites_full$Lat10[which(misclass_sites_full$cluster_membership == "Amazon Forest")],pch=15,col=rgb
       (t(col2rgb("blue"))/255,alpha=1), bg=rgb(t(col2rgb("blue"))/255, alpha=1))
#SDTF group
points(misclass_sites_full$Long10[which(misclass_sites_full$cluster_membership == "SDTF")]
       ,misclass_sites_full$Lat10[which(misclass_sites_full$cluster_membership == "SDTF")],pch=19,col=rgb
       (t(col2rgb("saddlebrown"))/255,alpha=1), bg=rgb(t(col2rgb("saddlebrown"))/255))
#Chaco group
points(misclass_sites_full$Long10[which(misclass_sites_full$cluster_membership == "Chaco")]
       ,misclass_sites_full$Lat10[which(misclass_sites_full$cluster_membership == "Chaco")],pch=25,col=rgb
       (t(col2rgb("black"))/255,alpha=1), bg=rgb(t(col2rgb("black"))/255, alpha=1))

######
#More details on the information presented in here can be found at: https://ourcodingclub.github.io/2017/03/21/data-clustering.html






