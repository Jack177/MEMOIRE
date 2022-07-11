 

library("bipartite")
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
library("reshape2")

cbp1 <- c( "#56B4E9", "#009E73",
           "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#E69F00","#999999")
cbp2 <- c( "#56B4E9", "#009E73",
           "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#E69F00","#8B6914","#999999",
           '#8B1C62','#32CD32',"#0000EE","#FF407A")


#################relations-plante-abeille
####data
prairie <- filter(stationtt, type == "Prairie")
bee2019 <- read.csv ( "beeréserve2019.csv",header = T, sep = ";")
bee2020 <- read.csv ( "beeréserve2020.csv",header = T, sep = ";")
bee2021 <- read.csv ( "réservetraitplante.csv",header = T, sep = ";")
bee1920 <- bind_rows(bee2019, bee2020)

library(magrittr)
mel <- bee2021 %>% select(Mellifère)
mel%<>% mutate_each(funs(if(is.character(.)) as.logical(.) else .))
bee2021 <- bee2021 %>% select(-Mellifère)
bee2021 <- bind_cols( bee2021, mel)

beetotal <- bind_rows(bee1920, bee2021)
beetotal$SPECTAXPRIO[beetotal$SPECTAXPRIO  %in% c("Bombus pascuorum floralis")]<-"Bombus pascuorum"
beetotal$SPECTAXPRIO[beetotal$SPECTAXPRIO  %in% c("Bombus terrestris terrestris")]<-"Bombus terrestris"
beetotal$CONDGR2[beetotal$CONDGR2  %in% c("")]<-"Vol/Sol"
beetotal$CONDTAXPRIO[beetotal$CONDTAXPRIO  %in% c("")]<-"Vol/Sol"

### création d'une matrice
beetotaly <-beetotal %>%filter(CONDTAXPRIO != 'Vol/Sol' )

beepla <- select(beetotaly,CONDTAXPRIO,SPECTAXPRIO,N)
beemat <- xtabs(N~SPECTAXPRIO+CONDTAXPRIO,beepla)

## réseau plant-polli
plotweb(beemat,text.rot=90, arrow="down", col.interaction="grey", 
        y.lim=c(-1,2.5),col.high="light green",col.low="light yellow")

#intensité desinteractions
low.abun = round(runif(dim(beemat)[1],1,20)) 
names(low.abun) <- rownames(beemat)
plotweb(beemat,text.rot=90,low.abun=low.abun, arrow="up", bor.col.interaction="red3",
        col.interaction="red", y.width.low=0.05, y.width.high=0.05, col.high="blue", col.low="green3")
#sous forme de grille
visweb(beemat,type="diagonal",box.col = "")

library("hrbrthemes")
library("GGally")
library("viridis")

ggparcoord(beemat,
           columns = 33:34, groupColumn = 2
) 
ggparcoord(beetotal,
           columns = 33:45, groupColumn = 13, order = "anyClass",
           showPoints = TRUE, 
           title = "Parallel Coordinate Plot for the Iris Data",
           alphaLines = 0.3) +
  scale_color_viridis(discrete=TRUE) +
  theme_ipsum()+
  theme(
    plot.title = element_text(size=10)
  )

#abondance

RankAbun.1 <- rankabundance(beestattab)
RankAbun.1
rankabunplot(RankAbun.1, scale='abundance', addit=FALSE, specnames=c(1,2,3))
RankAbun.1 <-type.convert(RankAbun.1)

###rarefaction , hill

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

#######comparaison occupations du sol 

occupatio1000nom <- occupatio %>% filter(grepl('1000',layer)) 
occupatio1000nom$layer[occupatio1000nom$layer %in% c("station_ancienne_gare1000")]<-"Ancienne Gare"
occupatio1000nom$layer[occupatio1000nom$layer  %in% c("station_bois_havré1000")]<-"Bois d'Havré"
occupatio1000nom$layer[occupatio1000nom$layer  %in% c("station_camp-a-cayaux1000")]<-"Camp-à-cayaux"
occupatio1000nom$layer[occupatio1000nom$layer  %in% c("station_cascade1000")]<-"Cascade d'Hyon"
occupatio1000nom$layer[occupatio1000nom$layer  %in% c("station_chateau_havré1000")]<-"Chateau d'Havré"
occupatio1000nom$layer[occupatio1000nom$layer  %in% c("station_cimetierre_spienne1000")]<-"Spiennes cimetière"
occupatio1000nom$layer[occupatio1000nom$layer  %in% c("station_epargne1000")]<-"Epargne - Umons"
occupatio1000nom$layer[occupatio1000nom$layer  %in% c("station_gare1000")]<-"Gare"
occupatio1000nom$layer[occupatio1000nom$layer  %in% c("station_gd_large1000")]<-"Grand Large"
occupatio1000nom$layer[occupatio1000nom$layer  %in% c("station_géothermia1000")]<-"Géothermia - IDEA"
occupatio1000nom$layer[occupatio1000nom$layer  %in% c("station_haine1000")]<-"La Haine"
occupatio1000nom$layer[occupatio1000nom$layer  %in% c("station_jardin_suspendu1000")]<-"Jardin suspendu"
occupatio1000nom$layer[occupatio1000nom$layer  %in% c("station_mont_panisel1000")]<-"Mont-Panisel"
occupatio1000nom$layer[occupatio1000nom$layer  %in% c("station_moulin1000")]<-"Ancien moulin"
occupatio1000nom$layer[occupatio1000nom$layer  %in% c("station_notre_dame_petit_nimy1000")]<-"Notre dame du petit Nimy"
occupatio1000nom$layer[occupatio1000nom$layer  %in% c("station_omya1000")]<-"carrière omya"
occupatio1000nom$layer[occupatio1000nom$layer  %in% c("station_parc_obourg1000")]<-"Parc d'Obourg"
occupatio1000nom$layer[occupatio1000nom$layer  %in% c("station_pemh_obourg1000")]<-"PEMH Obourg - IDEA"
occupatio1000nom$layer[occupatio1000nom$layer  %in% c("station_pont_prince1000")]<-"Pont du prince"
occupatio1000nom$layer[occupatio1000nom$layer  %in% c("station_prés_du_village1000")]<-"Prés du village"
occupatio1000nom$layer[occupatio1000nom$layer  %in% c("station_ronveaux1000")]<-"Ronveaux"
occupatio1000nom$layer[occupatio1000nom$layer  %in% c("station_silex1000")]<-"Silex"
occupatio1000nom$layer[occupatio1000nom$layer  %in% c("station_social1000")]<-"Siège social"
occupatio1000nom$layer[occupatio1000nom$layer  %in% c("station_st_waudru1000")]<-"Sainte-Waudru"
occupatio1000nom$layer[occupatio1000nom$layer  %in% c("station_tilou1000")]<-"Tilou"
occupatio1000nom$layer[occupatio1000nom$layer  %in% c("station_trouille1000")]<-"La Trouille"
occupatio1000nom$layer[occupatio1000nom$layer  %in% c("station_village_abeille1000")]<-"Village des abeilles"
names(occupatio1000nom) <- c("Sol_artificialisé","Sol nu","Bois","Sol agricol","Espace ouvert","Eau","Nom")
occupatio1000nom <-occupatio1000nom %>% mutate_if(is.numeric, ~round(., 1))
occupatio1000nom <- occupatio1000nom%>% arrange(Sol_artificialisé)

occupatio1000nom$Sol_artificialisé <- as.vector(occupatio1000nom$Sol_artificialisé)
ordre<-factor(occupatio1000nom$Nom,levels = rev(unique(occupatio1000nom$Nom)),ordered=TRUE)
ordre <- rep(ordre, times = 6)

data_long1 <- melt(occupatio1000nom, id.vars = c("Nom"))
data_long1 <- data_long1%>% mutate(reorder(Nom,ordre))
data_longordre <- as.data.frame(data_long1)


ggplot(occupatio1000nom,aes(reorder(Nom,Sol_artificialisé), y =Sol_artificialisé, fill= 'Sol_artificialisé')) +
  geom_bar(position = "stack",stat="identity")+
  xlab("Station") +
  ylab("Pourcentage")  +
  coord_flip()+
  scale_y_continuous(breaks = round(seq(min(0), max(100), by = 5),1))+
  scale_fill_manual('Sol_artificialisé',values = c("#009E73"))+
  labs(
    title = "Pourcentage de sol artificialisé dans un rayon de 1km de chaque station ",
    subtitle = "Stations échantillonnées en 2020",
    caption = "Occupation du sol donnée par la carte Lifewatch-WB geodatabase (v3.14) et le programme QGIS (version 3.1.3).
    Donnée légèrement differente pour la réserve, du fait du changement d'emplacement du centroïde.
    Réalisé grâce à RStudio,Version 1.3.1093") +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(face = ),
    plot.caption = element_text(face = "italic"),
    legend.position = "none")

#################################### SPATIAL ANALYSIS

######################## Moran's I regarder comment faire 
# for simple variables, string of values: each site has one value, such as species richness
# need to give the spatial coordinates, X and Y
# Tells you if there is spatial autocorrelation in the variable of interest.
moran <- select (stationtt,SPEC.TAXPRIO,N,TOPO)
stationtt <- read.csv ( "toutmémoire.csv",header = T, sep = ";")


zone.dists <- as.matrix(dist(cbind(stationtt $LATI, stationtt $LONG)))
zone.dists.inv <- 1/zone.dists
diag(zone.dists.inv) <- 0

Moran.I(database$Flower_abundance, zone.dists.inv)
Moran.I(database$Honeybee_rate, zone.dists.inv)

######################## Spatial representation variables

# bubble plot

database<-read.table("lineal.txt", header=T)

coords <- SpatialPoints(database[, c("X", "Y")], proj4string = CRS("+proj=longlat"))
plots <- SpatialPointsDataFrame(coords, database)
ddll <- spTransform(plots, CRS("+proj=longlat"))
pts <- as.data.frame(coordinates(ddll))
names(pts) <- c("lon", "lat")

print(bubble(plots, "Honeybee_rate", maxsize = 5,key.entries = 4*(1:5),col="blue"))#autecorrelation  
print(bubble(plots, "Flower_abundance", maxsize = 5,key.entries = 4*(1:5),col="blue"))#pas-autecorrelation 