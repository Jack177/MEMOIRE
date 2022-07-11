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
library("bipartite")

cbp2 <- c( "#56B4E9", "#009E73",
           "#0072B2", "#D55E00", "#CC79A7", "#E69F00","#8B6914","#999999",
           '#8B1C62','#32CD32',"#0000EE","#FF407A","#F0E442")
setwd("C:/Users/lalie/Desktop/mémoire/1mémoire-stat")

##################################### Data combiné abeille, station, plante, trait

reese <- read.csv ( "reesett.csv",header = T, sep = ";")
ruelle <- read.csv ( "beeplanttraitstat.csv",header = T, sep = ";")
ruelle$trait_phénott_ID <-  as.character(ruelle$trait_phénott_ID)
martin <- read.csv ( "beestatplanttrait_mart.csv",header = T, sep = ";")
martin$Mellifère <-  as.logical(martin$Mellifère)
alexl <- read.csv ( "beestatplanttrait_alexl.csv",header = T, sep = ";")
alexl$Indigène <-  as.logical(alexl$Indigène)
alexl$ITD <-  as.character(alexl$ITD)

stationtt <- bind_rows(reese, ruelle, martin, alexl)

type <- read.csv ( "listetype.csv",header = T, sep = ";")
stationtt <- left_join(station,type, by = "STATCODE")

stationtt<-stationtt %>% filter(SPECTAXPRIO != "Andrena  sp.")
stationtt <-stationtt %>% filter(SPECTAXPRIO != "Lasioglossum  sp.")
stationtt <-stationtt %>% filter(SPECTAXPRIO != "Osmia latreillei")

stationtt$SPECTAXPRIO[stationtt$SPECTAXPRIO%in% c("Bombus pascuorum floralis")]<-"Bombus pascuorum"
stationtt$SPECTAXPRIO[stationtt$SPECTAXPRIO%in% c("Bombus terrestris terrestris")]<-"Terrestribombus sp."
stationtt$SPECTAXPRIO[stationtt$SPECTAXPRIO%in% c("Bombus terrestris")]<-"Terrestribombus sp."
stationtt$SPECTAXPRIO[stationtt$SPECTAXPRIO%in% c("Bombus lucorum")]<-"Terrestribombus sp."
stationtt$SPECTAXPRIO[stationtt$SPECTAXPRIO  %in% c("Bombus pascuorum floralis")]<-"Bombus pascuorum"
stationtt$SPECTAXPRIO[stationtt$SPECTAXPRIO  %in% c("Bombus terrestris terrestris")]<-"Terrestribombus sp."
stationtt$SPECTAXPRIO[stationtt$SPECTAXPRIO  %in% c("Bombus (Bombus)  sp.")]<-"Terrestribombus sp."

stationtt$SPECTAXPRIO[stationtt$SPECTAXPRIO  %in% c("Bombus campestris campestris")]<-"Bombus campestris"
stationtt$SPECTAXPRIO[stationtt$SPECTAXPRIO  %in% c("Halictus confusus")]<-"Seladonia confusa"
stationtt$SPECTAXPRIO[stationtt$SPECTAXPRIO  %in% c("Hoplosmia spinulosa")]<-"Osmia spinulosa"
stationtt$SPECTAXPRIO[stationtt$SPECTAXPRIO  %in% c("osmia spinulosa")]<-"Osmia spinulosa"
stationtt$SPECTAXPRIO[stationtt$SPECTAXPRIO  %in% c("Halictus tumulorum")]<-"Seladonia tumulorum"
stationtt$SPECTAXPRIO[stationtt$SPECTAXPRIO  %in% c("halictus tumulorum")]<-"Seladonia tumulorum"
stationtt$SPECTAXPRIO[stationtt$SPECTAXPRIO  %in% c("Chalicodoma ericetorum")]<-"Megachile ericetorum"
stationtt$SPECTAXPRIO[stationtt$SPECTAXPRIO  %in% c("Bombus hortorum hortorum")]<-"Bombus hortorum"

#########matrice ABEILLE

bee <-select (stationtt,SPECTAXPRIO,N,TOPO)
beestattab <- select (stationtt,SPECTAXPRIO,N,TOPO)

beestattab <- aggregate(N~SPECTAXPRIO+TOPO, data = beestattab, sum)
beestattab <- xtabs(N~TOPO+SPECTAXPRIO,beestattab)
beestattab <-type.convert(beestattab)
beestattab <-as.data.frame(beestattab)

############################matrice environnement 

type <- read.csv ( "listetype.csv",header = T, sep = ";")
row.names(type)<- type$nom
type <- select(type, - STATCODE)
type$type <-as.factor(type$type)


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
rownames(buffertt) <- buffertt[,7]

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

abondance <- aggregate(N~TOPO, data = stationtt, sum)
rownames(abondance) <- abondance$TOPO

envibee <- full_join(buffertt,layer)
rownames(envibee)=envibee$layer

envi <- full_join(envibee,abondance, by = c("layer"="TOPO"))
envi <- full_join(envi,type, by = c("layer"  ="nom"))
rownames(envi)=envi$layer
envi <- filter(envi,layer != "Aérodrome de Maubray")
envi <- select(envi, -layer)

#########################################################################poisson
library('MuMIn')
library('caret')
library('statmod') 
library('lmtest')
library('AER')
library("MASS")


#####I choose the full coverage for my analyses

################################################################abondance, binomial négatif

bin <- glm.nb(N~ artificial+sol_nu_tot+agricol_to+espace_ouv+sp.florale.Totale, data = envi)
summary(bin) #regarder ce qui a des étoile = significatif
car::vif(bin)###vif trop grand avce tout a enlever Bois et eau car fonctionne mal

########################Analyse de déviance

lrtesta<- lrtest(bin)

################################graph diagnostique

qresiduals<- qresiduals(bin) 
qqnorm(qresiduals)
shapiro.test(qresiduals)

### Model selection
dd <- dredge(bin,extra="adjR^2")
ddd <- subset(dd, delta < 2) # select the ones with value of AICc less than 2 points in difference
subset(dd, delta < 2)

### Model averaging
avgmod.95delta2 <- model.avg(ddd) 
summary(avgmod.95delta2) 
lrtest<- lrtest(avgmod.95delta2)

##########################################################################H0r, poisson
#####transformer H0r pour qu'il soit entier
envi$H0r<- round(envi$H0r, digit=0)

poissonN <- glm(H0r~ artificial+sol_nu_tot+agricol_to+espace_ouv +sp.florale.Totale, data = envi, family = poisson())

summary(poissonN) 
car::vif(poissonN)

###########################analyse dispersion
dispersiontest(poissonN) #p-value doit etre  plus grande que 0,05, si 1=OK
  
########################Analyse de déviance
lrtest0<- lrtest(poissonN)
  
################################graph diagnostique

qresiduals<- qresiduals(poissonN) 
qqnorm(qresiduals)
shapiro.test(qresiduals)
ggqqplot(qresiduals)
  
###########selcetion modele
options(na.action = "na.fail")
dd <- dredge(poissonN)
    
### Model selection
dd <- dredge(poissonN,extra="adjR^2")
ddd <- subset(dd, delta < 2) # select the ones with value of AICc less than 2 points in difference
subset(dd, delta < 2)
glm(ddd)

### Model averaging
avgmod.95delta2 <- model.avg(ddd) 
summary(avgmod.95delta2)
lrtest<- lrtest(avgmod.95delta2)
    
##########################################################################H1r
#####transformer H1r pour qu'il soit entier

envi$H1r<- round(envi$H1r, digit=0)
poissonN1 <- glm(H1r~artificial+sol_nu_tot+agricol_to+espace_ouv +sp.florale.Totale, data = envi,family = poisson())
summary(poissonN1) 
car::vif(poissonN1)
  
###########################analyse dispersion
dispersiontest(poissonN1) 
  
########################Analyse de déviance
  lrtest1<- lrtest(poissonN1)
  
################################graph diagnostique
  
qresiduals<- qresiduals(poissonN1) 
qqnorm(qresiduals)
shapiro.test(qresiduals)
ggqqplot(qresiduals)
  
###########selcetion modele
options(na.action = "na.fail")
  
### Model selection
dd <- dredge(poissonN1,extra="adjR^2")
ddd <- subset(dd, delta < 2) # select the ones with value of AICc less than 2 points in difference
subset(dd, delta < 2)
  
### Model averaging
avgmod.95delta2 <- model.avg(ddd) 
summary(avgmod.95delta2)   
  
##########################################################################H2r
#####transformer H0r pour qu'il soit entier
envi$H2r<- round(envi$H2r, digit=0)
poissonN2 <- glm(H2r~sol_nu_tot+(artificial)+agricol_to+espace_ouv +sp.florale.Totale+eau_perman, data = envi,family = poisson())
summary(poissonN2) #regarder ce qui a des étoile = significatif
car::vif(poissonN2)###vif trop grand avce tout a enlever Bois
  
###########################analyse dispersion
dispersiontest(poissonN2) #p-value doit etre  plus grande que 0,05, si 1=OK
  
########################Analyse de déviance
lrtest2<- lrtest(poissonN2)
  
################################graph diagnostique
  
qresiduals<- qresiduals(poissonN2) 
qqnorm(qresiduals)
shapiro.test(qresiduals)# p au dessus de 0.05
ggqqplot(qresiduals)

###########selcetion modele
options(na.action = "na.fail")

### Model selection
dd <- dredge(poissonN2,extra="adjR^2")
ddd <- subset(dd, delta < 2) # select the ones with value of AICc less than 2 points in difference
subset(dd, delta < 2)
  
### Model averaging
avgmod.95delta2 <- model.avg(ddd) 
summary(avgmod.95delta2) # select the conditional average    
  
 
 

