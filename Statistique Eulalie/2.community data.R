
######################## Introductory lines

# set the working directory, which is the folder containing the files
# in my case:
setwd("formation stat")


######################## Beta-diversity

# QUALITATIVE DATA: PRESENCE-ABSENCE

library(betapart)

bees.qual <-read.table("Composition_pollinators.txt", header=T)
bees.qual[bees.qual > 0] <- 1

bees.beta<-beta.pair(bees.qual, index.family="sor")
bees.beta

# beta.sim: dissimilarity matrix accounting for spatial turnover (replacement)
# beta.sne: dissimilarity matrix accounting for nestedness-resultant dissimilarity
# beta.sor: dissimilarity matrix accounting for total dissimilarity

# QUANTITATIVE DATA: ABUNDANCES

bees <-read.table("Composition_pollinators.txt", header=T)
bray.part(bees)


######################## Mantel correlogram
# for complex variables, such as beta-diversity
# need to give the spatial coordinates, X and Y
# Tells you if there is spatial autocorrelation in the variable of interest.

library(vegan)
library(SpatialTools)
library(betapart)

database <-read.table("Lineal.txt", header=T)
bees <-read.table("Composition_pollinators.txt", header=T)

bees.bray<-bray.part(bees)

d <- as.matrix(database[,2:3])
f<-dist1(d)
d.dist<-as.dist(f*100)

mite.correlog <- mantel.correlog(bees.bray$bray, D.geo=d.dist, nperm=999)
mite.correlog
plot(mite.correlog)


######################### PCoA (Metric Multi-Dimensional Scaling)

library(vegan)
library(SpatialTools)
library(betapart)

database <-read.table("Lineal.txt", header=T)
bees <-read.table("Composition_pollinators.txt", header=T)

bees.bray<-bray.part(bees)

res <- pcoa(bees.bray)
#Error in D^2 : argument non num?rique pour un op?rateur binaire
res$values
biplot(res)


########################## NMDS (Non-metric Multi-Dimensional Scaling)alo moins precis, juste pour montrer

library(vegan)

# We want the function to find a solution, if not we can increase the 
# number of dimensions. number of dimensions has to be 2, max 3. 
# In any case stress has to be less than 0.2
# autotransform = F, you don't want the function to transform the data

set.seed(2)

bees <-read.table("Composition_pollinators.txt", header=T)
bees_NMDS <- metaMDS(bees, k=2,trymax=100,autotransform = F) # k are dimensions
bees_NMDS$stress

flowers<-read.table("Composition_flowers.txt",header=T)
flowers_NMDS <- metaMDS(flowers, k=2,trymax=100,autotransform = F)
flowers_NMDS$stress

flowersscores <- as.data.frame(scores(flowers_NMDS))
names(flowersscores) <- c("Flowers1","Flowers2")
beesscores <- as.data.frame(scores(bees_NMDS))
flowersscores$Bees1 <- beesscores$NMDS1
flowersscores$Bees2 <- beesscores$NMDS2

ordiplot(bees_NMDS,type="n")
orditorp(flowers_NMDS,display="sites",col="red",air=0.01)
orditorp(bees_NMDS,display="sites",cex=1.25,air=0.01)

 treat=c(rep("Treatment1",20),rep("Treatment2",20))
 ordiplot(bees_NMDS,type="n")
 ordihull(bees_NMDS,groups=treat,draw="polygon",col="grey90",label=F)
 orditorp(bees_NMDS,display="species",col="red",air=0.01)
 orditorp(bees_NMDS,display="sites",col=c(rep("green",5),rep("blue",5)),
          air=0.01,cex=1.25)

########################## Cluster dendrogram
# Ward Hierarchical Clustering with Bootstrapped p values

library(pvclust)

bees <-read.table("Composition_pollinators.txt", header=T)
tbees<-t(bees)

fit <- pvclust(tbees, method.hclust="ward",method.dist="euclidean")
# it takes some minutes, depending on the database
plot(fit) 
pvrect(fit, alpha=.95)


########################## Mantel test

library(vegan)
library(betapart)
library(SpatialTools)

database<-read.table("lineal.txt", header=T)
flowers<-read.table("Composition_flowers.txt",header=T)
bees <-read.table("Composition_pollinators.txt", header=T)

# matrix of geographic distances in km between sites
d <- as.matrix(database[,2:3])
f<-dist1(d)
d.dist<-as.dist(f*100)

quantitative.bees<-bray.part(bees)
quantitative.flowers <-bray.part(flowers)
apisrate <- dist(database$Honeybee_rate)

mantel(quantitative.bees$bray, quantitative.flowers$bray, method = "pearson", permutations = 999, na.rm = FALSE)

mantel.partial(quantitative.bees$bray, apisrate, quantitative.flowers$bray,method = "pearson", permutations = 999, na.rm = FALSE)
#qd veut elever l'effet d'un 3eme elemnt sur les 2 autres

# Check beta-diversity qualitattive nestedness and turnover, effect of environmental variables
library("betapart")
bees.qual <-read.table("Composition_pollinators.txt", header=T)
bees.qual[bees.qual > 0] <- 1
flowers.qual <-read.table("Composition_flowers.txt", header=T)
flowers.qual[flowers.qual > 0] <- 1

bees.beta<-beta.pair(bees.qual, index.family="sor")
flowers.beta<-beta.pair(flowers.qual, index.family="sor")

mantel(bees.beta$beta.sne, flowers.beta$beta.sor, method = "pearson", permutations = 999, na.rm = FALSE) # nestedness
mantel(bees.beta$beta.sim, apisrate, method = "pearson", permutations = 999, na.rm = FALSE) # turnover


################################ Fourth corner analysis


# # 4th corner analysis (Legendre's version)
# 
# # Reference papers :
# # First version of the analysis in this paper (useful to understand the second paper) :
# # http://biol09.biol.umontreal.ca/numecol/Reprints/4th-corner_paper.pdf
# # current version of the analysis used here (explanation of the different models):
# # https://www.researchgate.net/publication/250076886_Testing_the_species_traits_environment_relationships_The_fourth-corner_problem_revisited
# 
#  install.packages("ade4")
# # The code :
 #  library("ade4")
# 
# bees <-read.table("Composition_pollinators.txt", header=T)
#  species_traits <-read.table("species_traits2.txt", header=T)
#  sites_environment <-read.table("sites_environment.txt", header=T)
# 
#  fourth <- fourthcorner(
#    tabR = sites_environment, # Data on landscape (rows= stations, columns = environmental variables)
#    tabL = bees, # table species (columns) x stations (rows)
#    tabQ = species_traits, # traits (species = rows, columns = traits)
#    modeltype = 6, # kind of permutations of columns/rows applied to tabL. This is the best version (explained in reference paper 2)
#    p.adjust.method.G = "none",
#    p.adjust.method.D = "none",
#    nrepet = 999) # number of permutations. Should be really high, the higher the better
# 
# # Correction for multiple testing, here using FDR
# fourth.adj <- p.adjust.4thcorner(
#   fourth,
#   p.adjust.method.G = "fdr",
#   p.adjust.method.D = "fdr",
#   p.adjust.D = "global")
# 
# # Plot
# plot(fourth.adj, alpha = 0.05, stat = "D2")
# 
# # Three stats can be computed :
# # D2 = correlation 
# # D = homogeneity of each category (for qualitative variables) 
# # G is an  anova like stat for qualitative variables
