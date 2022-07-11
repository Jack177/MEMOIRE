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
stationtt <- left_join(stationtt,type, by = "STATCODE")

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

beestattab <- select (stationtt,SPECTAXPRIO,N,TOPO)
beestattab <- filter(beestattab, SPECTAXPRIO != "Andrena propinqua")

beestattab <- aggregate(N~SPECTAXPRIO+TOPO, data = beestattab, sum)
beestattab <- xtabs(N~TOPO+SPECTAXPRIO,beestattab)
beestattab <-type.convert(beestattab)
beestattab <-as.data.frame(beestattab)

############################matrice environneent 

type <- read.csv ( "listetype.csv",header = T, sep = ";")
row.names(type)<- type$nom
type <- select(type, - STATCODE)
type <- select(type, - nom, -couleur,-sp.florale.Totale )
type$type <-as.factor(type$type)
type$sp.florale.Totale <-as.numeric(type$sp.florale.Total)

############################matrice trait 

trait <-  read.csv ( "trait.csv",header = T, sep = ";")
listebee <- select(stationtt, SPECTAXPRIO,N)
listebee <- filter(listebee, SPECTAXPRIO != "Andrena propinqua")
listebee <- aggregate(N ~SPECTAXPRIO, listebee, sum)
listebee <- select (listebee, - N)
names(trait) <- c("espece", "ITD", "STI","SCI","nesting","lectisme","sociality","statut_IUCN")

trait <- left_join(listebee,trait, by = c("SPECTAXPRIO"= "espece"))
row.names(trait)<- trait$SPECTAXPRIO
trait <- select(trait,-SPECTAXPRIO,-statut_IUCN)

trait$STI <- as.numeric(trait$STI)
trait$SCI <- as.numeric(trait$SCI)
trait$lectisme<- as.factor(trait$lectisme)
trait$nesting<- as.factor(trait$nesting)
trait$sociality<- as.factor(trait$sociality)

#####################data sans Aérodrom de maubray
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
stationtt <- left_join(stationtt,type, by = "STATCODE")

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

statmaub <- filter(stationtt, TOPO != "Aérodrome de Maubray")

##################### matrice abeille sans Aérodrome 

beesmaub <- select (statmaub,SPECTAXPRIO,N,TOPO)
beesmaub <- filter( beesmaub,SPECTAXPRIO!="Andrena propinqua")


beesmaub <- aggregate(N~SPECTAXPRIO+TOPO, data = beesmaub, sum)
beesmaub <- xtabs(N~TOPO+SPECTAXPRIO,beesmaub)
beesmaub <-type.convert(beesmaub)

############################## matrice enviroment sans Aérodrome, occupation du sol

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

rownames(buffertt) <- buffertt[,7]
buffertt <- buffertt[,-7]
buffertt <- select(buffertt,-"bois_total", -"eau_perman")

###################matrice trait ans Aérodrome

traitmaub <-  read.csv ( "trait.csv",header = T, sep = ";")
listemaub <- select(statmaub, SPECTAXPRIO,N)
listemaub <- filter( listemaub,SPECTAXPRIO!= "Andrena propinqua")
listemaub <- aggregate(N ~SPECTAXPRIO, listemaub, sum)
listemaub <- select (listemaub, - N)
names(traitmaub) <- c("espece", "ITD", "STI","SCI","nesting","lectisme","sociality","statut_IUCN")

traitmaub <- left_join(listemaub,traitmaub, by = c("SPECTAXPRIO"= "espece"))
row.names(traitmaub)<- traitmaub$SPECTAXPRIO

traitmaub$STI <- as.numeric(traitmaub$STI)
traitmaub$SCI <- as.numeric(traitmaub$SCI)

traitmaub$nesting <- as.factor(traitmaub$nesting)
traitmaub$lectisme <- as.factor(traitmaub$lectisme)
traitmaub$sociality <- as.factor(traitmaub$sociality)

traitmaub <- select(traitmaub,-SPECTAXPRIO,-statut_IUCN)

################################ Fourth corner analysis
#################################lien trait-occupations du sol, pas utiliser Aérodrome
# # The code :
library("ade4")

beesmaub <- as.data.frame(beesmaub)
buffertt <- as.data.frame(buffertt)
traitmaub<- as.data.frame(traitmaub)

#####Division de la matrice trait pour les valeur numéric et les valeurs facteurs
traitmaubnum<- select(traitmaub,ITD,SCI,STI)
traitmaubcha<- select(traitmaub,nesting,lectisme,sociality)
# 

fourth <- fourthcorner(tabR = buffertt, tabL = beesmaub,tabQ = traitmaub,modeltype = 6, 
                       p.adjust.method.G = "none",p.adjust.method.D = "none",nrepet = 4999) 

# # Correction for multiple testing, here using FDR

fourth.adj <- p.adjust.4thcorner(fourth,p.adjust.method.G = "holm",p.adjust.method.D = "holm", p.adjust.D = "global")

  # # Plot

plot(fourth.adj, alpha = 0.05, stat = "D2")

  # # Three stats can be computed :
  # D2 = correlation 
  # D = homogeneity of each category (for qualitative variables) 
  # G is an  anova like stat for qualitative variables
  
 ################################################################lien entre trait et type de milieu

library("ade4")
# 
beestattab
type
trait
#####Division de la matrice trait pour les valeur numéric et les valeurs facteurs
traitnum <- select(trait, ITD,SCI,STI)
traitfac <- select(trait, -ITD,-SCI,-STI)
  # 
  
fourth <- fourthcorner(tabR = type, tabL = beestattab,tabQ = trait,modeltype = 6,p.adjust.method.G = "none",p.adjust.method.D = "none",nrepet = 4999) 
  
  # # Correction for multiple testing, here using FDR
  
fourth.adj <- p.adjust.4thcorner(fourth,p.adjust.method.G = "holm",p.adjust.method.D = "holm", p.adjust.D = "global")
  
  # # Plot
plot(fourth.adj, alpha = 0.05, stat = "D") 
plot(fourth.adj, alpha = 0.05, stat = "D2")
plot(fourth.adj, alpha = 0.05, stat = "G")

  # # Three stats can be computed :
  # D2 = correlation 
  # D = homogeneity of each category (for qualitative variables) 
  # G is an  anova like stat for qualitative variables
  