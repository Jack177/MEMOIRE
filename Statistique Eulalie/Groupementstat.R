library('ggpubr')
library("iNEXT")
library("BiodiversityR")
library("rgdal")
library("ape")
library("vegan")
library("SpatialTools")
library("betapart")
library("pvclust")
library("mapr")
library("rgbif")
library("dismo")
library("mapplots")
library("corrplot")
library("MuMIn")
library("caret")
library("pbkrtest")
library("car")
library("ggmap")
library("gstat")
library("tidyverse")
library("ggplot2")
library(pheatmap)
library("FactoMineR")
library('factoextra')

cbp2 <- c( "#56B4E9", "#009E73",
           "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#E69F00","#8B6914","#999999",
           '#8B1C62','#32CD32',"#0000EE","#FF407A")



##################################### Data combiné buffer

bufferreese <- read.csv ( "bufferalexr.csv",header = T, sep = ",")
bufferruelle <- read.csv ( "ttbuffer.csv",header = T, sep = ",")
buffermartin <- read.csv ( "buffermartin.csv",header = T, sep = ",")
bufferalexl <- read.csv ( "buffeuralexl.csv",header = T, sep = ",")


buffertt <- bind_rows(bufferreese, bufferruelle, buffermartin, bufferalexl)


buffertt <- buffertt %>% filter(grepl('250',path))
buffertt<- select (buffertt,-LCCS,-path )
buffertt <- distinct(buffertt)
rownames(buffertt) <- buffertt[,7]
buffertt <- buffertt[,-7]


##########################################ACP
PCA(buffertt, scale.unit = FALSE, ncp = 3, graph = TRUE)


library(FactoMineR)
# Compute PCA with ncp = 3

PCA(buffertt, scale.unit = FALSE, ncp = 6, graph = TRUE)
res.pca <-PCA(buffertt, scale.unit = FALSE, ncp = 6, graph = FALSE)

# Compute hierarchical clustering on principal components

HCPC(res.pca, graph = TRUE)
res.hcpc <- HCPC(res.pca, graph = TRUE)

fviz_dend(res.hcpc, 
          cex = 0.7,                     # Label size
          palette = "jco",               # Color palette see ?ggpubr::ggpar
          rect = TRUE, rect_fill = TRUE, # Add rectangle around groups
          rect_border = "jco",           # Rectangle color
          labels_track_height = 2      # Augment the room for labels
)

fviz_cluster(res.hcpc,
             repel = TRUE,            # Avoid label overlapping
             show.clust.cent = TRUE, # Show cluster centers
             palette = "jco",         # Color palette see ?ggpubr::ggpar
             ggtheme = theme_minimal(),
             main = "Factor map"
)

