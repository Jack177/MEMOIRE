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
library("reshape2")
library("reshape") 
library("SoDA")
library("dplyr")

cbp2 <- c( "#56B4E9", "#009E73",
           "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#E69F00","#8B6914","#999999",
           '#8B1C62','#32CD32',"#0000EE","#FF407A")



##################################### Data combiné abeille, station, plante, trait
reese <- read.csv ( "reesett.csv",header = T, sep = ";")
ruelle <- read.csv ( "beeplanttraitstat.csv",header = T, sep = ";")
ruelle$trait_phénott_ID <-  as.character(ruelle$trait_phénott_ID)
martin <- read.csv ( "beestatplanttrait_mart.csv",header = T, sep = ";")
martin$Mellifère <-  as.logical(martin$Mellifère)
alexl <- read.csv ( "beestatplanttrait_alexl.csv",header = T, sep = ";")
alexl$Indigène <-  as.logical(alexl$Indigène)
alexl$ITD <-  as.character(alexl$ITD)

station <- bind_rows(reese, ruelle, martin, alexl)

type <- read.csv ( "listetype.csv",header = T, sep = ";")
stationtt <- left_join(station,type, by = "STATCODE")

stationtt<-stationtt %>% filter(SPECTAXPRIO != "Andrena  sp.")
stationtt <-stationtt %>% filter(SPECTAXPRIO != "Lasioglossum  sp.")

stationtt$SPECTAXPRIO[stationtt$SPECTAXPRIO%in% c("Bombus pascuorum floralis")]<-"Bombus pascuorum"
stationtt$SPECTAXPRIO[stationtt$SPECTAXPRIO%in% c("Bombus terrestris terrestris")]<-"Terrestribombus"
stationtt$SPECTAXPRIO[stationtt$SPECTAXPRIO%in% c("Bombus terrestris")]<-"Terrestribombus"
stationtt$SPECTAXPRIO[stationtt$SPECTAXPRIO%in% c("Bombus lucorum")]<-"Terrestribombus"
stationtt$SPECTAXPRIO[stationtt$SPECTAXPRIO  %in% c("Bombus pascuorum floralis")]<-"Bombus pascuorum"
stationtt$SPECTAXPRIO[stationtt$SPECTAXPRIO  %in% c("Bombus terrestris terrestris")]<-"Terrestribombus"
stationtt$SPECTAXPRIO[stationtt$SPECTAXPRIO  %in% c("Bombus (Bombus)  sp.")]<-"Terrestribombus"

stationtt$SPECTAXPRIO[stationtt$SPECTAXPRIO  %in% c("Bombus campestris campestris")]<-"Bombus campestris"
stationtt$SPECTAXPRIO[stationtt$SPECTAXPRIO  %in% c("Halictus confusus")]<-"Seladonia confusa"
stationtt$SPECTAXPRIO[stationtt$SPECTAXPRIO  %in% c("Hoplosmia spinulosa")]<-"Osmia spinulosa"
stationtt$SPECTAXPRIO[stationtt$SPECTAXPRIO  %in% c("osmia spinulosa")]<-"Osmia spinulosa"
stationtt$SPECTAXPRIO[stationtt$SPECTAXPRIO  %in% c("Halictus tumulorum")]<-"Seladonia tumulorum"
stationtt$SPECTAXPRIO[stationtt$SPECTAXPRIO  %in% c("halictus tumulorum")]<-"Seladonia tumulorum"
stationtt$SPECTAXPRIO[stationtt$SPECTAXPRIO  %in% c("Chalicodoma ericetorum")]<-"Megachile ericetorum "
stationtt$SPECTAXPRIO[stationtt$SPECTAXPRIO  %in% c("Bombus hortorum hortorum")]<-"Bombus hortorum"
listebee <- select( stationtt,SPECTAXPRIO)
listebee <- distinct(listebee)

###matrice

beestattab <- select (stationtt,SPECTAXPRIO,N,TOPO)

beestattab <- aggregate(N~SPECTAXPRIO+TOPO, data = beestattab, sum)
beestattab <- xtabs(N~TOPO+SPECTAXPRIO,beestattab)
beestattab <-type.convert(beestattab)

################################################stationnarité du second ordre pour moran


coord <- read.csv ( "centroidettmémoire1.csv",header = T, sep = ",")
row.names(coord)= coord$AAStation

xygeo <-geoXY(coord$x, coord$y, unit = 1000)
row.names(xy)= coord$AAStation


plot(xy[,1], xy[,2], asp = 1)


modèle linéaire (régression de la variable réponse sur ses 
                 coordonnées géographiques)

graphhill <- select(hillmatrix,H0r,H1r,H2r)
bee.abond <-DataInfotab%>% select("S.obs", "SC")
rownames(bee.abond)=coord$AAStation 
names(bee.abond) <- c("spbeeobs","beecoverage")

moran <- bind_cols(graphhill,bee.abond)

coord <- select (coord, -path,-layer,-AAStation)


bufferreese <- read.csv ( "bufferalexr.csv",header = T, sep = ",")
bufferruelle <- read.csv ( "ttbuffer.csv",header = T, sep = ",")
bufferoublié <- read.csv ( "cockeri.csv",header = T, sep = ",")
buffermartin <- read.csv ( "buffermartin.csv",header = T, sep = ",")
bufferalexl <- read.csv ( "buffeuralexl.csv",header = T, sep = ",")

buffertt1 <- bind_rows(bufferreese, bufferruelle, buffermartin, bufferalexl)
buffertt1 <- buffertt1 %>% filter(grepl('250',path))
buffertt1<- select (buffertt1,-path )

buffertt <- bind_rows(buffertt1, bufferoublié)
buffertt<- select (buffertt,-LCCS )

buffertt <- distinct(buffertt)
write.csv2(buffertt, here::here("data_output","buffertt.csv"),row.names = FALSE ) 

buffertt$layer[buffertt$layer %in% c("station_ancienne_gare250")]<-"Ancienne Gare"
buffertt$layer[buffertt$layer  %in% c("station_bois_havré250")]<-"Bois d'Havré"
buffertt$layer[buffertt$layer  %in% c("station_camp-a-cayaux250")]<-"Camp-à-cayaux"
buffertt$layer[buffertt$layer  %in% c("station_cascade250")]<-"Cascade d'Hyon"
buffertt$layer[buffertt$layer  %in% c("station_chateau_havré250")]<-"Chateau d'Havré"
buffertt$layer[buffertt$layer  %in% c("station_cimetierre_spienne250")]<-"Spiennes cimetière"
buffertt$layer[buffertt$layer  %in% c("station_epargne250")]<-"Epargne - UMons"
buffertt$layer[buffertt$layer  %in% c("station_gare250")]<-"Gare"
buffertt$layer[buffertt$layer  %in% c("station_gd_large250")]<-"Grand Large"
buffertt$layer[buffertt$layer  %in% c("station_géothermia250")]<-"Géothermia - IDEA"
buffertt$layer[buffertt$layer  %in% c("station_haine250")]<-"La Haine"
buffertt$layer[buffertt$layer  %in% c("station_jardin_suspendu250")]<-"Jardin suspendu"
buffertt$layer[buffertt$layer  %in% c("station_mont_panisel250")]<-"Mont-Panisel"
buffertt$layer[buffertt$layer  %in% c("station_moulin250")]<-"Ancien moulin"
buffertt$layer[buffertt$layer  %in% c("station_notre_dame_petit_nimy250")]<-"Notre dame du petit Nimy"
buffertt$layer[buffertt$layer  %in% c("station_omya250")]<-"carrière omya/ le Caufour"
buffertt$layer[buffertt$layer  %in% c("station_parc_obourg250")]<-"Parc d'Obourg"
buffertt$layer[buffertt$layer  %in% c("station_pemh_obourg250")]<-"PEMH Obourg - IDEA"
buffertt$layer[buffertt$layer  %in% c("station_pont_prince250")]<-"Pont du prince"
buffertt$layer[buffertt$layer  %in% c("station_prés_du_village250")]<-"Prés du village"
buffertt$layer[buffertt$layer  %in% c("station_ronveaux250")]<-"Ronveaux"
buffertt$layer[buffertt$layer  %in% c("station_silex250")]<-"Silex"
buffertt$layer[buffertt$layer  %in% c("station_social250")]<-"Siège social - UMons"
buffertt$layer[buffertt$layer  %in% c("station_st_waudru250")]<-"Sainte-Waudru"
buffertt$layer[buffertt$layer  %in% c("station_tilou250")]<-"Tilou"
buffertt$layer[buffertt$layer  %in% c("station_trouille250")]<-"La Trouille"
buffertt$layer[buffertt$layer  %in% c("station_village_abeille250")]<-"Village des abeilles - UMons"
buffertt$layer[buffertt$layer  %in% c("Station_Terril Bayemont-Saint Charles")]<-"Terril Bayemont-Saint-Charles"
buffertt$layer[buffertt$layer  %in% c("Station_Terril Belle vue du huite")]<-"Terril Belle vue du huit"
buffertt$layer[buffertt$layer  %in% c("Station_Terril Ciply st1")]<-"Terril de Ciply site 1"
buffertt$layer[buffertt$layer  %in% c("Station_Terril Dix-Huit")]<-"Terril du dix-huit"
buffertt$layer[buffertt$layer  %in% c("Station_Terril Ferrand")]<-"Terril du Ferrand"
buffertt$layer[buffertt$layer  %in% c("Station_Terril Grand Buisson")]<-"Terril du Grand Buisson"
buffertt$layer[buffertt$layer  %in% c("Station_Terril Hensie")]<-"Terril d'Hensies"
buffertt$layer[buffertt$layer  %in% c("Station_Terril HÃ©ribus st1")]<-"Terril de l'Héribus"
buffertt$layer[buffertt$layer  %in% c("Station_Terril Naye-Ã -Bois")]<-"Terril Naye-à-bois"
buffertt$layer[buffertt$layer  %in% c("Station_Terril Pays-Bas")]<-"Terril n°8 Pays-bas"
buffertt$layer[buffertt$layer  %in% c("Station_Terril Quesnoy st1")]<-"Terril du Quesnoy site 1"
buffertt$layer[buffertt$layer  %in% c("Station_Terril Rieux du Coeur")]<-"Terril Rieu-du-Coeur"
buffertt$layer[buffertt$layer  %in% c("Station_Terril SacrÃ©-FranÃ§ais st1")]<-"Terril Sacré-Français site 1"
buffertt$layer[buffertt$layer  %in% c("Station_Terril Saint Antoine st2")]<-"Terril Saint-Antoine site 2"
buffertt$layer[buffertt$layer  %in% c("Station_Terril Saint Antoine st3")]<-"Terril Saint-Antoine site 3"
buffertt$layer[buffertt$layer  %in% c("Station_Terril Sept st1")]<-"Terril du sept huit stat 1"
buffertt$layer[buffertt$layer  %in% c("Station_Terril Sept st3")]<-"Terril du sept-huit stat 3"
buffertt$layer[buffertt$layer  %in% c("Station_Terril Trazegnies st1")]<-"Terril n°5 de Trazegnie site 1"
buffertt$layer[buffertt$layer  %in% c("Nom_MLM10")]<-"Planoi Site 2"
buffertt$layer[buffertt$layer  %in% c("Nom_MLM14")]<-"Nocarcentre Site 2"
buffertt$layer[buffertt$layer  %in% c("Nom_MLM17")]<-"Nocarcentre Site 5"
buffertt$layer[buffertt$layer  %in% c("Nom_MLM2")]<-"Bruyère Site 2"
buffertt$layer[buffertt$layer  %in% c("Nom_MLM7")]<-"Vertbois Site 2"
buffertt$layer[buffertt$layer  %in% c("Vieill_Haine")]<-"Vieille Haine"
buffertt$layer[buffertt$layer  %in% c("Bois_de_Bon-Secours")]<-"Bois de Bon-Secours"
buffertt$layer[buffertt$layer  %in% c("Bois_de_Wadelincourt")]<-"Bois de Wadelincourt"
buffertt$layer[buffertt$layer  %in% c("Chemin_de_Roucourt")]<-"Chemin de Roucourt"
buffertt$layer[buffertt$layer  %in% c("Chemin_de_Trainage")]<-"Chemin du Trainage"
buffertt$layer[buffertt$layer  %in% c("Friche_des_Vignobles")]<-"Friche des Vignobles"
buffertt$layer[buffertt$layer  %in% c("Marais_de_Douvrais")]<-"Marais de Douvrain Ouest"
buffertt$layer[buffertt$layer  %in% c("Mer_de_Sable")]<-"Mer de Sable"
buffertt$layer[buffertt$layer  %in% c("Mont_Ostènes")]<-"Mont Ostènes"
buffertt$layer[buffertt$layer  %in% c("Parc_de_Jemappes")]<-"Parc de Jemappes"
buffertt$layer[buffertt$layer  %in% c("Parc_des_5_Rocs")]<-"Parc des 5 rocs"
buffertt$layer[buffertt$layer  %in% c("Pré_à_Parchon")]<-"Pré à Parchon"
buffertt$layer[buffertt$layer  %in% c("Rue_du_Bois")]<-"Rue du Bois"
buffertt$layer[buffertt$layer  %in% c("Rue_de_Carne")]<-"Rue du Carme"
buffertt$layer[buffertt$layer  %in% c("Rue_de_Castillon")]<-"Rue du Castillon"
buffertt$layer[buffertt$layer  %in% c("Rue_Jean_Winance")]<-"Rue Jean Winance"
buffertt$layer[buffertt$layer  %in% c("Rue_La-dessous")]<-"Rue là-dessous"

row.names(buffertt)= buffertt$layer


abond <- iNEXT(t(beestattab), q = c(0,1,2), datatype ="abundance") 

DataInfotab <-type.convert(abond$DataInfo)
iNextEsttab <-type.convert(abond$iNextEst)
AsyEsttab <-type.convert(abond$AsyEst)

hill <- estimateD(t(beestattab), base="coverage")
hillmatrix0= hill[hill[,"order"]==0,]
hillmatrix1= hill[hill[,"order"]==1,]
hillmatrix2= hill[hill[,"order"]==2,]
hillmatrix = hillmatrix0[,c(1,2,4)]
hillmatrix$SC0 = abond$DataInfo$SC
hillmatrix[,c("H0r","H0rU","H0rL")]=hillmatrix0[,c("qD","qD.UCL","qD.LCL")]
hillmatrix[,c("H1r","H1rU","H1rL")]=hillmatrix1[,c("qD","qD.UCL","qD.LCL")]
hillmatrix[,c("H2r","H2rU","H2rL")]=hillmatrix2[,c("qD","qD.UCL","qD.LCL")]
rownames(hillmatrix)=rownames(beestattab) 

graphhill <- select(hillmatrix,H0r,H1r,H2r)
bee.abond <-DataInfotab%>% select("S.obs", "SC")
rownames(bee.abond)=coord$AAStation 
names(bee.abond) <- c("spbeeobs","beecoverage")
layer<-rownames_to_column(graphhill)
layer <- rename(layer,c('rowname'='layer'))

rownames(type) <- type[,1]
type <- buffertt[,-1]

abondance <- aggregate(N~TOPO, data = stationtt, sum)
rownames(abondance) <- type[,1]
rownames(type) <- type[,1]

 
envi <- full_join(buffertt,layer)
rownames(envi)=envi$layer

envi <- full_join(envi,abondance, by = c("layer"="TOPO"))
envi <- full_join(envi,type, by = c("layer"  ="nom"))

coord <- select (coord, -path,-layer)
envi <- full_join(envi,coord, by = c( "layer"="AAStation"))
envi <- select(envi, -layer)

################################################Corrélogrammes de Moran
# Load the required packages
library('ape')
library('spdep')
library('ade4')
library('adegraphics')
library('adespatial')
# Source additional functions(files must be in the working directory)
source ('http://www.davidzeleny.net/anadat-r/doku.php/en:numecolr:plot.links?do=export_code&codeblock=1')
source("plot.links.R")
source ('http://www.davidzeleny.net/anadat-r/doku.php/en:numecolr:sr.value?do=export_code&codeblock=1')
source("sr.value.R")
source ('https://raw.githubusercontent.com/zdealveindy/anadat-r/master/scripts/NumEcolR2/quickMEM.R')
source("quickMEM.R")
source ('https://raw.githubusercontent.com/zdealveindy/anadat-r/master/scripts/NumEcolR2/scalog.R')
source("scalog.R")

xy<-select(envi, x, y)

zone.dists <- as.matrix(dist(cbind(envi$x, envi$y)))
zone.dists.inv <- 1/zone.dists
diag(zone.dists.inv) <- 0
# Transform the data
mite.h <- decostand (dis, "hellinger")
mite.xy.c <- scale(xygeo, center = TRUE, scale = FALSE)
## Univariate spatial correlogram (based on Moran's I)Search for neighbours of all points within a radius of 0.7 m and multiples (i.e., 0 to 0.7 m, 0.7 to 1.4 m and so on).
plot.links(mite.xy.c, thresh = 500)
nb1 <- dnearneigh(as.matrix(mite.xy.c), 0,500)
summary(nb1)
# Correlogram of substrate density
subs.dens <- envi[ ,7]
subs.correlog <-
  sp.correlogram(nb1, 
                 subs.dens, 
                 order = 2, 
                 method = "I", 
                 zero.policy = TRUE)
print(subs.correlog, p.adj.method = "holm")
plot(subs.correlog)
     
     subs.dens <- envi[ ,8]
     subs.correlog <-
       sp.correlogram(nb1, 
                      subs.dens, 
                      order = 2, 
                      method = "I", 
                      zero.policy = TRUE)
     print(subs.correlog, p.adj.method = "holm")
     plot(subs.correlog)
################################################stationnarité du second ordre pour mantel

analyse des redondances (régressions multiples de la table de contingence sur les coordonnées géographiques). 
Le résultat obtenu étant non significatif (p = 0.647), nous pouvons considérer que la condition est respectée.


################################################Corrélogramme de Mantel  + correctton holm

mite.h.det <- resid(lm(as.matrix(mite.h) ~ .,data =as.data.frame(xy)))
mite.h.D1 <- dist(mite.h.det)
mite.correlog <- mantel.correlog(mite.h.D1,XY = xy, nperm = 999)
summary(mite.correlog)
# Plot the Mantel correlogram
plot(mite.correlog)

mantel.correlog(dis,zone.dists.inv , n.class=0, break.pts=NULL, 
                cutoff=TRUE, r.type="spearman", nperm=999, mult="holm", progressive=TRUE)
"plot"(x, alpha=0.05, ...)

#####mantel
  
mite.hel <- decostand(dis, "hellinger")

# Detrend the species data by regression on the site coordinates
mite.hel.resid <- resid(lm(as.matrix(mite.hel) ~ ., data=xy))

# Compute the detrended species distance matrix
mite.hel.D <- dist(mite.hel.resid)
mite.correlog2 <- mantel.correlog(mite.hel.D,zone.dists.inv, cutoff=FALSE, 
                                  r.type="spearman", nperm=999, mult="holm")
summary(mite.correlog2)
mite.correlog2
plot(mite.correlog2)


######################## Moran's I
# for simple variables, string of values: each site has one value, such as species richness
# need to give the spatial coordinates, X and Y
# Tells you if there is spatial autocorrelation in the variable of interest.

library(ape)

zone.dists <- as.matrix(dist(cbind(envi$x, envi$y)))
zone.dists.inv <- 1/zone.dists
diag(zone.dists.inv) <- 0

Moran.I(envi$H0r, zone.dists.inv,na.rm = TRUE)#pasok
Moran.I(envi$H1r, zone.dists.inv,na.rm = TRUE)#pasok
Moran.I(envi$H2r, zone.dists.inv,na.rm = TRUE)#pasok
Moran.I(envi$N, zone.dists.inv,na.rm = TRUE)#ok
plot(MOran.I)

subs.dens <- envi[ ,7]
subs.correlog <-
  sp.correlogram(nb1, 
                 subs.dens, 
                 order = 3, 
                 method = "I", 
                 zero.policy = TRUE)
print(subs.correlog, p.adj.method = "holm")
plot(subs.correlog)

######################## Spatial representation variables

# bubble plot

library(rgdal)
library(gstat)

database<-read.table("lineal.txt", header=T)

coords <- SpatialPoints(envi[, c("x", "y")], proj4string = CRS("+proj=longlat"))
plots <- SpatialPointsDataFrame(coords, envi)
ddll <- spTransform(plots, CRS("+proj=longlat"))
pts <- as.data.frame(coordinates(ddll))
names(pts) <- c("lon", "lat")

print(bubble(plots, "H0r", maxsize = 72,key.entries = 4*(1:72),col="blue"))#autecorrelation  
print(bubble(plots, "H1r", maxsize = 72,key.entries = 4*(1:725),col="blue"))#pas-autecorrelation 

# spatial grid colours interpolation

library(sp)
library(tidyverse)
library(ggmap)
library(gstat)


x.range <- range(envi$x)
y.range <- range(envi$y)
x<-seq(x.range[1], x.range[2], length.out=20)
y<-seq(y.range[1], y.range[2], length.out=20)
grd<-expand.grid(x,y)
coordinates(envi) = ~x+y
coordinates(grd) <- ~ Var1+Var2
gridded(grd) <- TRUE
proj4string(envi) <- CRS("+proj=longlat +datum=Lambert")
proj4string(grd) <- CRS("+proj=longlat +datum=Lambert72")
plot(grd, cex=1.5)

dat.idw <- idw(formula=H0r ~ 1, locations = envi, newdata = grd, idp = 2.0)
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

selectedbees <- beesstattab[,2:5]
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


