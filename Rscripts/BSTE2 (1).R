#07/12/2018 - Project on Environemntal Spaces within Comunities of Dry Forests in Nicaragua
#Project with Jesus Gomes Zurita, based on Gisela's Thesis 
##
install.packages('dismo')
install.packages('maptools')
install.packages('hypervolume')
install.packages('nicheROVER')
install.packages('raster')
install.packages('alphahull')
install.packages('magick')
library(dismo)
library(raster)
library(maptools)
library(hypervolume)
library(nicheROVER)
library(alphahull)
library(magick)

#Preparing variables -------------------------------------------------------
#Testing Correlations
install.packages('Hmisc')
library(Hmisc)
library(raster)
var_list<- stack(list.files('D:\\Jesus_data\\BiosWC_ENVIROM\\Nicaragua_30s\\Nicaragua', pattern='.asc', full.names = T))
bio1 <-raster('D:\\Jesus_data\\BiosWC_ENVIROM\\Nicaragua_30s\\Nicaragua\\PCAs\\Comp1.asc')
NicaraguaXY <- rasterToPoints(bio1, fun=NULL, spatial=FALSE) #head(NicaraguaXY)                    # transform to points to create a table of values of points 
NicaraguaXY <- cbind(NicaraguaXY[,1], NicaraguaXY[,2])
colnames(NicaraguaXY) <- c('long', 'lat')
NicaraguaValuesBios <- extract(NicaraguaXY,var_list)
sp <- SpatialPoints(NicaraguaXY) 
Nicaragua_variables <- extract(var_list, sp)
library(Hmisc)
cor<-rcorr(as.matrix(Nicaragua_variables))
correlation<-data.frame(cor$r)
write.csv(correlation,'D:\\Jesus_data\\BiosWC_ENVIROM\\Thin_data\\Nicaragua_correlation.csv')
# graphics - GOTTA MAKE THIS WORK PRINT CORRELATION OF VARIABLES TO INCLUDE IN POWERPOINT
install.packages('corrgram')
library(corrgram)
cor_fig<-corrgram(Nicaragua_variables, order=TRUE, lower.panel=panel.shade,
                  upper.panel=panel.pie, text.panel=panel.txt,
                  main="All variables")
dev.copy(cor_fig,"D:\\Jesus_data\\PlantCoverage\\VariableCorrelation\\all_cor.png",width=8,height=6,units="in",res=300)
dev.off()

install.packages('corrplot')
library(corrplot)
corrplot(cor(Nicaragua_variables))


#Stack of PCs from Nicaragua-------------------------------------------
install.packages('rgdal')
install.packages('raster')
install.packages('maptools')
library(maptools)
library(rgdal)
library(raster)
Comp1<-raster("D:\\Jesus_data\\BiosWC_ENVIROM\\ENVIROM_only\\Comp1.asc")
Comp2<-raster("D:\\Jesus_data\\BiosWC_ENVIROM\\ENVIROM_only\\Comp2.asc")
Comp3<-raster("D:\\Jesus_data\\BiosWC_ENVIROM\\ENVIROM_only\\Comp3.asc")
Comp4<-raster("D:\\Jesus_data\\BiosWC_ENVIROM\\ENVIROM_only\\Comp4.asc")
Comp5<-raster("D:\\Jesus_data\\BiosWC_ENVIROM\\ENVIROM_only\\Comp5.asc")
Comp6<-raster("D:\\Jesus_data\\BiosWC_ENVIROM\\ENVIROM_only\\Comp6.asc")
#Comp6<-raster("D:\\Jesus_data\\BiosWC_ENVIROM\\Worldclim2_only\\Comp6.asc")

st <- stack(Comp1,Comp2,Comp3,Comp4,Comp5, Comp6) #4 first components explain 95%
plot(st)
head(st)
#Getting XYs from each area
library(sp)
library(rgeos)
library(raster)
library(rworldmap)
install.packages('sf')
library(sf)
Esteli_5 <- st_read("D:\\Jesus_data\\NicheA\\Buffers_shp\\Esteli_5.shp") #head)
#r=mask(Comp1,Esteli_50)
#writeRaster(r, filename="D:\\Jesus_data\\NicheA\\Buffers_shp\\XYs\\Granada_20.asc", format="ascii",overwrite=TRUE)
#rivas=readShapePoly("C:\\Users\\mari1\\Desktop\\Jesus\\data\\Rivas.shp")
writeRaster(r, 'D:\\Jesus_data\\NicheA\\Buffers_shp\\XYs\\Rivas_50.asc')
x <- raster("D:\\Jesus_data\\NicheA\\Buffers_shp\\XYs\\Rivas_50.asc")
data_matrix <- rasterToPoints(x) #head(data_matrix)
write.csv(data_matrix, file = "D:\\Jesus_data\\NicheA\\Buffers_shp\\XYs\\Rivas_50.csv") 

#Creating tables with PCs values for each pixel of each area
xy <- read.csv("D:\\Jesus_data\\NicheA\\Buffers_shp\\XYs\\Managua_20.csv") 
xy_values <- extract(st, xy[,2:3]) #str(xy_values)
write.csv(xy_values, file = "D:\\Jesus_data\\BiosWC_ENVIROM\\Worldclim2_only\\XY_values\\Managua_20.csv") #just Pcs values
#xy_values1 <- cbind(xy[2:3], xy_values);   #Pc values and coordinates head(xy_values1)
#write.csv(xy_values1, file = "Rivas_values.csv")


##Analysis: Niche Overlap and Distance##-----------------------------
#1. Hypervolume ---------------------------------------------------------
install.packages('alphahull')
install.packages('hypervolume')
install.packages('rgl')
install.packages('later')
library(later)
library(rgl)
library(alphahull)
library(hypervolume)

data <- read.csv('D:\\Jesus_data\\BiosWC_ENVIROM\\Worldclim2_only\\XY_values\\All_20.csv')  #na.omit(data) #Pcs of morphological data
species_list = unique(data$Species)  
num_species = length(species_list)
hv_data_list = new("HypervolumeList")
hv_data_list@HVList = vector(mode="list",length=num_species)  #creating the list where all the values will be stored

# compute hypervolumes for each species #I DONT USE THE LOG OF SPECIES VALUES - DIFFERS FROM MORPHO SPACE FOR THAT !
for (i in 1:num_species)
{
  this_species = subset(data, Species==species_list[i])
  # keep the trait data
  this_species_log <- this_species[,2:ncol(this_species)]
  na.omit(this_species_log)
  print(this_species_log)
  # make a hypervolume using auto-bandwidth
  hv_data_list@HVList[[i]] <- hypervolume_gaussian(this_species_log, name=as.character(species_list[i]),verbose=FALSE)
}
# compute all pairwise overlaps
overlap = matrix(NA, nrow=num_species, ncol=num_species)
dimnames(overlap)=list(species_list, species_list)
for (i in 1:num_species)
{
  for (j in i:num_species)
  {
    if (i!=j)
    {
      # compute set operations on each pair
      this_set = hypervolume_set(hv_data_list@HVList[[i]], hv_data_list@HVList[[j]], check.memory=FALSE)
      # calculate a Sorensen overlap index (2 x shared volume / sum of |hv1| + |hv2|)
      overlap[i,j] = hypervolume_overlap_statistics(this_set)["sorensen"]
    }
  }   
}

#save a table with the amount of overlap. This will give you numbers and allow you to interpret the level of overlpat of hypervolumes
overlap
# kde.bandwidth=estimate_bandwidth(this_species_log),reps=10000, quantile=0
#hv_annulus <- hypervolume_gaussian(data_annulus, kde.bandwidth=0.05,name='annulus',samples.per.point=1)

write.csv(overlap, file = "D:\\Jesus_data\\BiosWC_ENVIROM\\Worldclim2_only\\Hypervolume\\OverlapAreas_20.csv")

# show all hypervolumes
# Plot hypervolume matrix plot
plot(hv_data_list,npmax.random=500,darkfactor=0.5,cex.legend=0.8,cex.names=1, xjust= 2, yjust= 2) #plot 3d pairwise=F

#3D hypervolume image
plot(hv_data_list,npmax.random=5000, show.3d=TRUE)
hypervolume_save_animated_gif('D:\\Jesus_data\\Hypervolume\\Results\\Hypervlume_20.gif')
library(rgl)
rgl.snapshot('D:\\Jesus_data\\Hypervolume\\Results\\Hypervlume_20.png')


#movie3d(spin3d(),duration=5,movie='gif_5',dir='./').
#animated version
hypervolume_save_animated_gif(image.size = 400,
                              axis = c(0, 0, 1), rpm = 4, duration = 15, fps = 10,
                              file.name = "movie", directory.output = "D:\\Jesus_data\\Hypervolume\\Results\\OverlapAreas_10km.gif")

#2. Niche Analyst-----------------------------------------
#Transposing results
x1 <- read.csv("D:\\Jesus_data\\NicheA\\N_NewVar\\NicheOverlap\\overlap_result5km.csv")
as.data.frame(x1)
dist_NicheA<- dist.quant(t(x1), 1, diag = TRUE, upper = TRUE) 
write.csv(t(x1), "D:\\Jesus_data\\NicheA\\N_NewVar\\NicheOverlap\\overlap_NicheA_NewVar.csv")
##Distance from every oint to every point - creating distributions of dictances-----------
#Distance analysis
##XYs -------------------------------------------------------------------------------------------------------------
library(maptools)
library(rgdal)
library(raster)
library(ade4)
area <- readOGR("C:\\Users\\mari1\\Desktop\\Jesus\\data\\Nicaragua.shp") #plot(area)
#XYs from Esteli, Granada, Managua
x1 <- read.csv("D:\\Jesus_data\\NicheA\\Esteli_10.csv") #head(x1)
x1 <- cbind(x1[,2], x1[,3])
colnames(x1) <- c('long', 'lat')
x2 <- read.csv("D:\\Jesus_data\\NicheA\\Granada_10.csv") #head(x2)
x2 <- cbind(x2[,2], x2[,3])
colnames(x2) <- c('long', 'lat')
x3 <- read.csv("D:\\Jesus_data\\NicheA\\Managua_10.csv") #head(x3)
x3 <- cbind(x3[,2], x3[,3])  
colnames(x3) <- c('long', 'lat')
x4 <- read.csv("D:\\Jesus_data\\NicheA\\Rivas_10.csv") #head(x4)
x4 <- cbind(x4[,2], x4[,3])
colnames(x4) <- c('long', 'lat')
head(x4)
##transformando esri > raster e extraindo valores para pointos nessa layer ##head of the coordernates
#get comps rom files - didi in the first part of teh code
#extract every component value to each area
#AREAS #Esteli: 1 #Managua: 2  #Granada: 3 #Rivas: 4
StackValue_1 <- data.frame(extract(st,x1)) #esteli #Extract XY for every component in the stack str(StackValue)      ##extract the values for every point of m4 from the layer of altitude
StackValue_2 <- data.frame(extract(st,x2)) #granada
StackValue_3 <- data.frame(extract(st,x3)) #managua
StackValue_4 <- data.frame(extract(st,x4)) #rivas

#3. Euclidean Distance Matrix----------------
Comp1<-raster("D:\\Jesus_data\\BiosWC_ENVIROM\\Nicaragua_30s\\Nicaragua\\PCAs\\Comp1.asc")
Comp2<-raster("D:\\Jesus_data\\BiosWC_ENVIROM\\Nicaragua_30s\\Nicaragua\\PCAs\\Comp2.asc")
Comp3<-raster("D:\\Jesus_data\\BiosWC_ENVIROM\\Nicaragua_30s\\Nicaragua\\PCAs\\Comp3.asc")
Comp4<-raster("D:\\Jesus_data\\BiosWC_ENVIROM\\Nicaragua_30s\\Nicaragua\\PCAs\\Comp4.asc")
Comp5<-raster("D:\\Jesus_data\\BiosWC_ENVIROM\\Nicaragua_30s\\Nicaragua\\PCAs\\Comp5.asc")
Comp6<-raster("D:\\Jesus_data\\BiosWC_ENVIROM\\Nicaragua_30s\\Nicaragua\\PCAs\\Comp6.asc")
Comp7<-raster("D:\\Jesus_data\\BiosWC_ENVIROM\\Nicaragua_30s\\Nicaragua\\PCAs\\Comp7.asc")
Comp8<-raster("D:\\Jesus_data\\BiosWC_ENVIROM\\Nicaragua_30s\\Nicaragua\\PCAs\\Comp8.asc")
#Comp6<-raster("D:\\Jesus_data\\BiosWC_ENVIROM\\Worldclim2_only\\Comp6.asc")
st <- stack(Comp1,Comp2, Comp3, Comp4, Comp5, Comp6, Comp7, Comp8)#8 first components explain 95%
plot(st)
head(st)
#Get column means 
install.packages('pdist')
library(pdist)
Esteli_10 <- colMeans(StackValue_1, na.rm=T, dim=1)
Granada_10 <- colMeans(StackValue_2, na.rm=T, dim=1)
Managua_10 <- colMeans(StackValue_3, na.rm=T, dim=1)
Rivas_10 <- colMeans(StackValue_4, na.rm=T, dim=1)
a<-as.data.frame(Esteli_10)
b<-as.data.frame(Granada_10)
c<-as.data.frame(Managua_10)
d<-as.data.frame(Rivas_10)
matrix_10 <- cbind(a, b, c, d) #str(dframe) 
t_matrix_10 <- t(matrix_10)
euc_dist_10 <- pdist(t(matrix_10))
euc_dist_10_mtx <-as.matrix(euc_dist_10)
write.csv(euc_dist_10_mtx, "D:\\Jesus_data\\BiosWC_ENVIROM\\Env_Dist_II\\Euclid_10.csv")                
#Gower Distance--------------------------
install.packages('cluster')
library(cluster)
## Example 1 in ref:
##  Dissimilarities using Euclidean metric and without standardization
gower_dist_10 <- daisy(t_matrix_10, metric = "gower", stand = FALSE)
gower_dist_10_mtx <-as.matrix(gower_dist_10)
gower_dist_10
write.csv(gower_dist_10_mtx, "D:\\Jesus_data\\BiosWC_ENVIROM\\Env_Dist_II\\Gower_10.csv")   
#Gower's coefficient (1971), expressed as a dissimilarity, 
# implies that a particular standardisation will be applied to each 
#, and the "distance" between two units is the sum of all the 
#variable-specific distances

#Standardized distances- Raw variables -------------------------
setwd('D:\\Jesus_data\\BiosWC_ENVIROM\\Nicaragua_30s\\Nicaragua')
var1<-raster("current_30arcsec_annualPET.asc")
var2<-raster("current_30arcsec_aridityIndexThornthwaite.asc")
var3<-raster("current_30arcsec_PETColdestQuarter.asc")
#var4<-raster("current_30arcsec_PETDriestQuarter.asc")
var5<-raster("current_30arcsec_PETWarmestQuarter.asc")
var6<-raster("current_30arcsec_PETWettestQuarter.asc")
#var7<-raster("wc2.0_bio_30s_02.asc")
var8<-raster("wc2.0_bio_30s_15.asc")
var_st <-stack(var1,var2, var3, var5, var6, var8)
plot(var_st)
install.packages('robustHD')
library(robustHD)

stand <-scale(var_st, center = TRUE, scale = TRUE) #scale stack of variables
StackValue_1_st <- data.frame(extract(stand,x1)) #extract values to each area from the scaled stack
na.omit(StackValue_1_st)
StackValue_2_st <- data.frame(extract(stand,x2)) #granada
na.omit(StackValue_2_st)
StackValue_3_st <- data.frame(extract(stand,x3)) #managua
na.omit(StackValue_3_st)
StackValue_4_st <- data.frame(extract(stand,x4)) #rivas
na.omit(StackValue_4_st)
library(pdist)
Esteli_10_st <- colMeans(StackValue_1_st, na.rm=TRUE, dim=1)
Granada_10_st <- colMeans(StackValue_2_st, na.rm=TRUE, dim=1)
Managua_10_st <- colMeans(StackValue_3_st, na.rm=TRUE, dim=1)
Rivas_10_st <- colMeans(StackValue_4_st, na.rm=TRUE, dim=1)
a_st<-as.data.frame(Esteli_10_st)
b_st<-as.data.frame(Granada_10_st)
c_st<-as.data.frame(Managua_10_st)
d_st<-as.data.frame(Rivas_10_st)
matrix_10_st <- cbind(a_st, b_st, c_st, d_st) #str(dframe) 
matrix_10_st <- t(matrix_10_st)
euc_dist_10_st <- pdist(matrix_10_st)
euc_dist_10_st <-as.matrix(euc_dist_10_st)
write.csv(euc_dist_10_st, "D:\\Jesus_data\\BiosWC_ENVIROM\\Env_Dist_II\\Euclid_10_st_6.csv")

#Gower Distance--------------------------
install.packages('cluster')
library(cluster)
## Example 1 in ref:
##  Dissimilarities using Euclidean metric and without standardization
gower_dist_10_st <- daisy(matrix_10_st, metric = "gower", stand = FALSE)
gower_dist_10_st

#Gower's coefficient (1971), expressed as a dissimilarity, 
# implies that a particular standardisation will be applied to each 
#, and the "distance" between two units is the sum of all the 
#variable-specific distances



#Preparing csv files to compare component Values of each area
dframe <- cbind(StackValue_1[,1], StackValue_2[,1], StackValue_3[,1], StackValue_4[,1]) #str(dframe)  
colnames(dframe) <- c('Esteli', 'Managua','Granada','Rivas')
write.csv(dframe, "D:\\Jesus_data\\BiosWC_ENVIROM\\Euclidean_distance\\5Km\\Euclid_Comp1.csv")                

#Dist.Quant function: # centroid distance? Not really----------------------------------------------------------------
#1 - Canonical distance, 2- Joreskog; 3-mahalanobis distance
install.packages('stats')
install.packages('reshape')
library(ade4)
library(stats)
library(reshape)
dframe_dist <- read.csv('D:\\Jesus_data\\BiosWC_ENVIROM\\Euclidean_distance\\5Km\\Euclid_Comp1.csv')
dist_Comp1 <- dist.quant(t(dframe_dist), 2, diag = TRUE, upper = TRUE) #calculating distances eucl.
m <- as.matrix(dist_Comp1) #transform distance object to matrix
m2 <- melt(m)[melt(upper.tri(m))$value,] # melt so we have 3 columns
names(m2) <- c("From", "To", "distance") #add heading 
write.csv(m2, "D:\\Jesus_data\\BiosWC_ENVIROM\\Euclidean_distance\\EuclidDist_Canonical_Comp7.csv")

#Euclidean: Network distance --------------------------------------
install.packages('qgraph')
install.packages('stringi')
install.packages('network')
library(stringi)
library(network)
library(qgraph)
dist_mi <- 1/dist_Comp1 # one over, as qgraph takes similarity matrices as input
jpeg('D:\\Jesus_data\\BiosWC_ENVIROM\\Euclidean_distance\\Euc_5km_.jpg', width=1000, height=1000, unit='px')
qgraph(dist_Comp1, title= 'Distance PC1 Buffer 5Km',edge.labels= T, edge.label.bg=T, edge.label.cex=1.5, labels= colnames(dframe), label.scale= T, weighted=T, cut=0.00001, label.color= 'blue', layout='spring', vsize=10, edge.color = "dark grey",legend.cex = 0.5, vsize = 1 )

# Euclidean: dist function ------------------------------------------------------
install.packages('stats')
install.packages('reshape')
library(stats)
library(reshape)
library(rgdal)
#DataFrame of Comp1 of each area_EUCLIDEAN OF ONLY ONE COMP
library(ade4)
area <- readOGR("C:\\Users\\mari1\\Desktop\\Jesus\\data\\Nicaragua.shp") #plot(area)
x1 <- read.csv("D:\\Jesus_data\\NicheA\\Esteli_20.csv") #head(x1)
x1 <- cbind(x1[,2], x1[,3])
colnames(x1) <- c('long', 'lat')
x2 <- read.csv("D:\\Jesus_data\\NicheA\\Granada_20.csv") #head(x2)
x2 <- cbind(x2[,2], x2[,3])
colnames(x2) <- c('long', 'lat')
x3 <- read.csv("D:\\Jesus_data\\NicheA\\Managua_20.csv") #head(x3)
x3 <- cbind(x3[,2], x3[,3])  
colnames(x3) <- c('long', 'lat')
x4 <- read.csv("D:\\Jesus_data\\NicheA\\Rivas_20.csv") #head(x4)
x4 <- cbind(x4[,2], x4[,3])
colnames(x4) <- c('long', 'lat')
#extract values for areas XYs
StackValue_1 <- extract(st,x1) #Esteli  #  str(StackValue_1)      ##extract the values for every point of m4 from the layer of altitude
StackValue_2 <- extract(st,x2) #Managua
StackValue_3 <- extract(st,x3) #Granada
StackValue_4 <- extract(st,x4) #Rivas
dframe <- cbind(StackValue_1[,1], StackValue_2[,1], StackValue_3[,1], StackValue_4[,1]) #str(dframe)  
colnames(dframe) <- c('Esteli', 'Managua','Granada','Rivas')
#write.csv(dframe, "D:\\Jesus_data\\Euclidean_distance\\Comp1_Euclidean.csv")                
dis_euc<-dist(t(dframe), method = "euclidean", diag = FALSE, upper = FALSE)
library(qgraph)
qgraph(dis_euc, title= 'Distance 20Km buffer',edge.labels= T, edge.label.bg=T, edge.label.cex=1.5, labels= colnames(dframe), label.scale= T, cut=0.00001, label.color= 'blue', layout='spring', vsize=10, edge.color = "dark grey",legend.cex = 0.5, vsize = 1 )
m <- as.matrix(dis_euc)
m2 <- melt(m)[melt(upper.tri(m))$value,]
names(m2) <- c("From", "To", "distance")
write.csv(m2, "D:\\Jesus_data\\BiosWC_ENVIROM\\Thin_data\\Euclidean_distance\\Network_Dist\\Euclidean_20km_Comp1.csv")


# Euclidean: all cells and comps + histograms------------------------------------
library(raster)
x1 <- read.csv("D:\\Jesus_data\\NicheA\\Esteli_20.csv") #head(x1)
x1 <- cbind(x1[,2], x1[,3])
x2 <- read.csv("D:\\Jesus_data\\NicheA\\Granada_20.csv") #head(x2)
x2 <- cbind(x2[,2], x2[,3])
x3 <- read.csv("D:\\Jesus_data\\NicheA\\Managua_20.csv") #head(x3)
x3 <- cbind(x3[,2], x3[,3])  
x4 <- read.csv("D:\\Jesus_data\\NicheA\\Rivas_20.csv") #head(x4)
x4 <- cbind(x4[,2], x4[,3])  
StackValue_1 <- data.frame(extract(st,x1)) #esteli #Extract XY for every component in the stack str(StackValue)      ##extract the values for every point of m4 from the layer of altitude
StackValue_2 <- data.frame(extract(st,x2)) #granada
StackValue_3 <- data.frame(extract(st,x3)) #managua
StackValue_4 <- data.frame(extract(st,x4)) #rivas
Est_10 <- as.matrix(StackValue_1) #all comps for area Esteli head(Est_5)
colnames(Est_20) <- c('Esteli_comp1', 'Esteli_comp2','Esteli_comp3','Esteli_comp4', 'Esteli_comp5', 'Esteli_comp6', 'Esteli_comp7', 'Esteli_comp8')
Gra_20 <- as.matrix(StackValue_2) #all comps for area Granada head(Gra_5)
colnames(Gra_20) <- c('Granada_comp1', 'Granada_comp2','Granada_comp3','Granada_comp4', 'Granada_comp5', 'Granada_comp6', 'Granada_comp7', 'Granada_comp8')
Man_20 <- as.matrix(StackValue_3) #all comps for area Managua
colnames(Man_20) <- c('Managua_comp1', 'Managua_comp2','Managua_comp3','Managua_comp4', 'Managua_comp5', 'Managua_comp6', 'Managua_comp7', 'Managua_comp8')
Riv_20 <- as.matrix(StackValue_4) #all comps for area Rivas
colnames(Riv_20) <- c('Rivas_comp1', 'Rivas_comp2','Rivas_comp3','Rivas_comp4', 'Rivas_comp5', 'Rivas_comp6', 'Rivas_comp7', 'Rivas_comp8')
#Euclidean paiwise distance - X to Y - ALL COMPS AT ONCE
install.packages('pdist')
library(pdist) #Euclidean paiwise distance - X to Y
d_R_M_20<-pdist(Riv_20, Man_10) #Euclidean paiwise distance between rows X to rows Y
hist(as.numeric(d_R_M_20@dist)) #range(d_M_E_5@dist, na.rm = FALSE)

#Creating overlapping densities 20, 10, 5
#Used this to create the graphs
library(ggplot2)
a_5<- as.data.frame(d_G_E_5@dist)  #range(d_G_E_5@dist, na.rm = FALSE)
b_5 <- as.data.frame(d_M_E_5@dist) #range(d_M_E_5@dist, na.rm = FALSE)
c_5 <- as.data.frame(d_M_G_5@dist)
d_5 <- as.data.frame(d_R_E_5@dist)
e_5 <- as.data.frame(d_R_G_5@dist)
f_5 <- as.data.frame(d_R_M_5@dist)

a_10<- as.data.frame(d_G_E_10@dist)  #range(f_10, na.rm = FALSE)
b_10 <- as.data.frame(d_M_E_10@dist) #range(d_M_E_10@dist, na.rm = FALSE)
c_10 <- as.data.frame(d_M_G_10@dist)
d_10 <- as.data.frame(d_R_E_10@dist)
e_10 <- as.data.frame(d_R_G_10@dist)
f_10 <- as.data.frame(d_R_M_10@dist)

a_20<- as.data.frame(d_G_E_20@dist)  #range(f_20, na.rm = FALSE)
b_20 <- as.data.frame(d_M_E_20@dist) #range(d_M_E_20@dist, na.rm = FALSE)
c_20 <- as.data.frame(d_M_G_20@dist)
d_20 <- as.data.frame(d_R_E_20@dist)
e_20 <- as.data.frame(d_R_G_20@dist)
f_20 <- as.data.frame(d_R_M_20@dist)

#Histograms of euclidean distances
library(ggplot2)
barfill <- "lightblue"
barlines <- "darkblue"

ggplot(f_20, aes(d_R_M_20@dist)) +
  geom_histogram(colour=barlines,fill=barfill, binwidth=0.3)+
  theme_classic() +
  coord_cartesian(xlim = c(0, 25)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))


























##More variables 
#Percentage of tree coverage-DEM-Anual Evapotranspiration - Anual Aridity - Primary productivity
library(raster)
r <- raster("D:\\Jesus_data\\PlantCoverage\\gm_ve_v1.tif")
library(maptools)
library(rgdal)
area <- readOGR("C:\\Users\\mari1\\Desktop\\Jesus\\data\\Nicaragua.shp") #plot(area)
c <-mask(r, area)
writeRaster(c, "D:\\Jesus_data\\PlantCoverage\\PlantCov_Nic.asc", format="ascii")
c<-raster("D:\\Jesus_data\\PlantCoverage\\PlantCov_Nic.asc")
plot(c)


#extract values to all xy 
library(raster)
bio1<- raster("D:\\Jesus_data\\Bios\\Nicaragua_30s\\Nicaragua\\bio_1.asc")
bio2<- raster("D:\\Jesus_data\\Bios\\Nicaragua_30s\\Nicaragua\\bio_2.asc")
bio3<- raster("D:\\Jesus_data\\Bios\\Nicaragua_30s\\Nicaragua\\bio_3.asc")
bio4<- raster("D:\\Jesus_data\\Bios\\Nicaragua_30s\\Nicaragua\\bio_4.asc")
bio5<- raster("D:\\Jesus_data\\Bios\\Nicaragua_30s\\Nicaragua\\bio_5.asc")
bio6<- raster("D:\\Jesus_data\\Bios\\Nicaragua_30s\\Nicaragua\\bio_6.asc")
bio7<- raster("D:\\Jesus_data\\Bios\\Nicaragua_30s\\Nicaragua\\bio_7.asc")
#bio8<- raster("D:\\Jesus_data\\Bios\\Nicaragua_30s\\Nicaragua\\bio_8.asc")
#bio9<- raster("D:\\Jesus_data\\Bios\\Nicaragua_30s\\Nicaragua\\bio_9.asc")
bio10<- raster("D:\\Jesus_data\\Bios\\Nicaragua_30s\\Nicaragua\\bio_10.asc")
bio11<- raster("D:\\Jesus_data\\Bios\\Nicaragua_30s\\Nicaragua\\bio_11.asc")
bio12<- raster("D:\\Jesus_data\\Bios\\Nicaragua_30s\\Nicaragua\\bio_12.asc")
bio13<- raster("D:\\Jesus_data\\Bios\\Nicaragua_30s\\Nicaragua\\bio_13.asc")
bio14<- raster("D:\\Jesus_data\\Bios\\Nicaragua_30s\\Nicaragua\\bio_14.asc")
bio15<- raster("D:\\Jesus_data\\Bios\\Nicaragua_30s\\Nicaragua\\bio_15.asc")
bio16<- raster("D:\\Jesus_data\\Bios\\Nicaragua_30s\\Nicaragua\\bio_16.asc")
bio17<- raster("D:\\Jesus_data\\Bios\\Nicaragua_30s\\Nicaragua\\bio_17.asc")
#bio18<- raster("D:\\Jesus_data\\Bios\\Nicaragua_30s\\Nicaragua\\bio_18.asc")
#bio19<- raster("D:\\Jesus_data\\Bios\\Nicaragua_30s\\Nicaragua\\bio_19.asc")
alt<- raster("D:\\Jesus_data\\Bios\\Nicaragua_30s\\Nicaragua\\alt.asc")
StackRawBios <- stack(bio1, bio2, bio3, bio4, bio5, bio6, bio7, bio10, bio11, bio12, bio13, bio14, bio15, bio16, bio17,alt)
NicaraguaXY <- rasterToPoints(bio1, fun=NULL, spatial=FALSE) #head(NicaraguaXY)                    # transform to points to create a table of values of points 
NicaraguaXY <- cbind(NicaraguaXY[,1], NicaraguaXY[,2])
colnames(NicaraguaXY) <- c('long', 'lat')
NicaraguaValuesBios <- rasterToPoints(StackRawBios, fun=NULL, spatial=FALSE) #head(NicaraguaValuesBios)





