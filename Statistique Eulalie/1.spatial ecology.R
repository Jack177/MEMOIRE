
######################## Introductory lines

# set the working directory, which is the folder containing the files
# in my case:
setwd("formation stat")

# install the different packages that you will need
install.packages("rgdal")
install.packages("ape")
install.packages("vegan")
install.packages("SpatialTools")
install.packages("betapart")
install.packages("pvclust")
install.packages("mapr")
install.packages("rgbif")
install.packages("dismo")
install.packages("mapplots")
install.packages("BiodiversityR")
install.packages("corrplot")
install.packages("MuMIn")
install.packages("caret")
install.packages("car")
install.packages("pbkrtest")
install.packages("ggmap")
install.packages("gstat")
install.packages("tidyverse")


####################### Rank-abundance curves

library("BiodiversityR")

bees <-read.table("Composition_pollinators.txt", header=T)
RankAbun.1 <- rankabundance(bees)
RankAbun.1
rankabunplot(RankAbun.1, scale='abundance', addit=FALSE, specnames=c(1,2,3))


# ######################## Rarefaction
# 
 install.packages("iNEXT")
 library(iNEXT) # This is the package developped by Chao et al. for rarefaction. 
# # See : https://cran.r-project.org/web/packages/iNEXT/vignettes/Introduction.html
# # Reference papers :
# # Why they use completeness rather than abundance for rarefaction :
# # https://www.researchgate.net/publication/235713090_Coverage-based_rarefaction_and_extrapolation_Standardizing_samples_by_completeness_rather_than_size
# # Rarefaction of all Hill numbers :
# # https://www.researchgate.net/publication/273219459_Rarefaction_and_extrapolation_with_Hill_numbers_A_framework_for_sampling_and_estimation_in_species_diversity_studies
# 
 bees <-read.table("Composition_pollinators.txt", header=T)
 rownames(bees) <- c(1:40)
# # The code :
# # "communaute" is a dataframe filled with numbers only. It is written with columns as species and 
 #rows as stations. Be sure to save the name of the stations as rownames.
# # communaute <- BDD_sec[,-1]# remove the columns with data that are not numbers
# # rownames(communaute)=BDD_sec$TOPO # save the stations in rownames
# 
# #### iNEXT : computing Hill numbers and Hill curves.
# # The iNEXT function needs the transposed table. q is the list of Hill numbers we need 
 #(0 = species richness, 1 is linked to shannon and 2 to simpson)
# # We have to precise that the table contains abundance data (could be incidence).
 test <- iNEXT(t(bees), q = c(0,1,2), datatype ="abundance") 
 # test is a list of 3 elements
 test$DataInfo # summary of abundance, species richness, coverage and number of singletons, doubletons, tripletons, etc
 test$iNextEst # estimation of coordinates for the rarefaction/accumulation curve. 
 #Completeness in X and hill number in Y (there is a list of coordinates for each hill order)
 test$AsyEst # Asymptotic diversity estimates for each order. Sites are labelled with letters.
# 
# # estimateD is the rarefaction function. It works with Hill numbers 0, 1, 2.
# # Transposed table is used. By default datatype ="abundance" (could be incidence)
# # We have to precise which metric is used for completeness (coverage is better than abundance). Their definition of coverage 
 #is in reference papers.
# # level = NULL by default (not written here). This way, the rarefaction use the smallest value of completeness among sites. 
 #Could be determined by user.
 hill <- estimateD(t(bees), base="coverage") # m is the sample size for the reference level of completeness 
 #(we used coverage as completeness metric, thus sample size is not constant)
# # SC = sample coverage (should be roughly equal to the minimum coverage among sites)
# # the table says if the data is interpolated (rarefaction), observed (for sites with the lowest coverage) 
 #or extrapolated (should not happen as we chose to lower each station coverage to the lowest one)
# #produce an ugly table, need some transformations (one station = 3 rows)
# #values for hill numbers are given with upper and lower estimate (95% CI)
 hillmatrix0= hill[hill[,"order"]==0,]
 hillmatrix1= hill[hill[,"order"]==1,]
 hillmatrix2= hill[hill[,"order"]==2,]
 hillmatrix = hillmatrix0[,c(1,2,4)]
 hillmatrix$SC0 = test$DataInfo$SC# retrieve the initial coverage from the table computed by iNEXT()
 hillmatrix[,c("H0r","H0rU","H0rL")]=hillmatrix0[,c("qD","qD.UCL","qD.LCL")]
 hillmatrix[,c("H1r","H1rU","H1rL")]=hillmatrix1[,c("qD","qD.UCL","qD.LCL")]
 hillmatrix[,c("H2r","H2rU","H2rL")]=hillmatrix2[,c("qD","qD.UCL","qD.LCL")]
rownames(hillmatrix)=rownames(communaute) # put the labels on the rows


#################################### SPATIAL ANALYSIS

######################## Moran's I
# for simple variables, string of values: each site has one value, such as species richness
# need to give the spatial coordinates, X and Y
# Tells you if there is spatial autocorrelation in the variable of interest.
library(ape)

database <- read.table("lineal.txt",header=T)

zone.dists <- as.matrix(dist(cbind(database$X, database$Y)))
zone.dists.inv <- 1/zone.dists
diag(zone.dists.inv) <- 0

Moran.I(database$Flower_abundance, zone.dists.inv)
Moran.I(database$Honeybee_rate, zone.dists.inv)



######################## Spatial representation variables

# bubble plot

library(rgdal)
library(gstat)

database<-read.table("lineal.txt", header=T)

coords <- SpatialPoints(database[, c("X", "Y")], proj4string = CRS("+proj=longlat"))
plots <- SpatialPointsDataFrame(coords, database)
ddll <- spTransform(plots, CRS("+proj=longlat"))
pts <- as.data.frame(coordinates(ddll))
names(pts) <- c("lon", "lat")

print(bubble(plots, "Honeybee_rate", maxsize = 5,key.entries = 4*(1:5),col="blue"))#autecorrelation  
print(bubble(plots, "Flower_abundance", maxsize = 5,key.entries = 4*(1:5),col="blue"))#pas-autecorrelation 

# spatial grid colours interpolation

library(sp)
library(tidyverse)
library(ggmap)
library(gstat)

database<-read.table("lineal.txt", header=T)

x.range <- range(database$X)
y.range <- range(database$Y)
x<-seq(x.range[1], x.range[2], length.out=20)
y<-seq(y.range[1], y.range[2], length.out=20)
grd<-expand.grid(x,y)
coordinates(database) = ~X+Y
coordinates(grd) <- ~ Var1+Var2
gridded(grd) <- TRUE
proj4string(database) <- CRS("+proj=longlat +datum=WGS84")
proj4string(grd) <- CRS("+proj=longlat +datum=WGS84")
plot(grd, cex=1.5)

dat.idw <- idw(formula=Honeybee_rate ~ 1, locations = database, newdata = grd, idp = 2.0)
plot(dat.idw)

# bubble chart categories

library(sp)
library(ape)
library(raster)
library(ggmap)
library(mapr)
library(rgbif)
library(dismo)
library(mapplots)
library(tidyr)
library(vegan)
library(SpatialTools)
library(ggplot2)
library(dplyr)

bees <-read.table("Composition_pollinators.txt", header=T)
selectedbees <- bees[,2:5]
database<-read.table("lineal.txt", header=T)
selectedbees$X <- database$X
selectedbees$Y <- database$Y
selectedbees$PLOT <- database$PLOT

selectedbees <- selectedbees %>%
  gather(species, abundance, c("Rodanthidium_sticticum","Anthophora_acervorum","Pseudophilotes_panoptes","Sarcophagidae_4")) %>%
  arrange(PLOT, X, Y, species)

xlim <- c(1.82,1.95)
ylim <- c(41.25,41.32)
xyz <- make.xyz(selectedbees$X,selectedbees$Y,selectedbees$abundance,selectedbees$species)
col <- c('red','green','blue','yellow')
basemap(xlim, ylim,bg='white')
draw.pie(xyz$x, xyz$y, xyz$z, radius = 0.005, col=col)
legend.pie(1.93,41.26,labels=unique(selectedbees$species), radius=0.005, bty="n", col=col,cex=0.8, label.dist=1.3)

legend.z <- round(max(rowSums(xyz$z,na.rm=TRUE)),0)
legend.bubble(1.83,41.30,z=legend.z,round=1,maxradius=0.008,bty="n",txt.cex=0.7)
text(1.83,41.32,"Abundance",cex=0.8)

