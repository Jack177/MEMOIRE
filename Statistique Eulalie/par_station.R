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
library("reshape") 


cbp2 <- c( "#56B4E9", "#009E73",
           "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#E69F00","#8B6914","#999999",
           '#8B1C62','#32CD32',"#0000EE","#FF407A")

#######Data
réserve <- read.csv ("réservetraitplante.csv",header = T, sep = ";")

####difference de famille par station par individu
station <- select (réserve, TOPO, SPECGR2,N)

station <- aggregate(N ~ TOPO +SPECGR2 , data =station, sum)

ggplot(station, aes(SPECGR2, y = N, fill = TOPO)) +
  geom_col(position = "dodge")+
  xlab("Famille") +
  ylab("Nombre")  +
  #☻coord_flip()+
  scale_fill_manual('Station',values = cbp2)+
  scale_y_continuous(breaks = round(seq(min(0), max(100), by = 5),1))+
  labs(
    title = "Nombre d'individus par familles capturés sur chaque station",
    subtitle = "Division des familles par station",
    caption = "Calcul réalisé grâce aux données de 2021, sur un total de 279 individus. Réalisé grâce à RStudio,Version 1.3.1093") +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(face = ),
    plot.caption = element_text(face = "italic"))

####difference de famille par station par sp

stationin <- select (réserve, TOPO, SPECGR2,SPECTAXPRIO)

stationin <- distinct(stationin)
stationin <- stationin %>% group_by(TOPO)%>% count(SPECGR2)

ggplot(stationin, aes(SPECGR2, y = n, fill = TOPO)) +
  geom_col(position = "dodge")+
  xlab("Famille") +
  ylab("Nombre")  +
  #☻coord_flip()+
  scale_fill_manual('Station',values = cbp2)+
  scale_y_continuous(breaks = round(seq(min(0), max(100), by = 5),1))+
  labs(
    title = "Nombre d'espèces par familles capturées sur chaque station",
    subtitle = "Division des familles par station",
    caption = "Calcul réalisé grâce aux données de 2021, sur un total de 34 espèces 
    + andrena sp., nomada sp. et le groupe des Terrestribombus.Réalisé grâce à RStudio,Version 1.3.1093") +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(face = ),
    plot.caption = element_text(face = "italic"))

####difference de genre par station par individu
stationgen <- select (réserve, TOPO, SPECGEN,N)

stationgen <- aggregate(N ~ TOPO +SPECGEN , data =stationgen, sum)

ggplot(stationgen, aes(SPECGEN, y = N, fill = TOPO)) +
  geom_col(position = "dodge")+
  xlab("Genre") +
  ylab("Nombre")  +
  #☻coord_flip()+
  scale_fill_manual('Station',values = cbp2)+
  scale_y_continuous(breaks = round(seq(min(0), max(100), by = 5),1))+
  labs(
    title = "Nombre d'individus par genre capturé sur chaque station",
    subtitle = "Division des genre par station",
    caption = "Calcul réalisé grâce aux données de 2021, sur un total de 279 individus.
    Réalisé grâce à RStudio,Version 1.3.1093") +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(face = ),
    plot.caption = element_text(face = "italic"))

####difference de genre par station par sp
stationspgen <- select (réserve, TOPO, SPECGEN,SPECTAXPRIO)

stationspgen <- distinct(stationspgen)
stationspgen <- stationspgen %>% group_by(TOPO)%>% count(SPECGEN)

ggplot(stationspgen, aes(SPECGEN, y = n, fill = TOPO)) +
  geom_col(position = "dodge")+
  xlab("Genre") +
  ylab("Nombre")  +
  #☻coord_flip()+
  scale_fill_manual('Station',values = cbp2)+
  scale_y_continuous(breaks = round(seq(min(0), max(100), by = 5),1))+
  labs(
    title = "Nombre d'espèces par genres capturées sur chaque station",
    subtitle = "Division des genres par station",
    caption = "Calcul réalisé grâce aux données de 2021, sur un total de 34 espèces 
    + andrena sp., nomada sp. et le groupe des Terrestribombus.Réalisé grâce à RStudio,Version 1.3.1093") +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(face = ),
    plot.caption = element_text(face = "italic"))

####difference d'sp par station par individu
stationsp <- select (réserve, TOPO, SPECTAXPRIO,N)

stationsp  <- aggregate(N ~ TOPO +SPECTAXPRIO , data =stationsp , sum)

ggplot(stationsp, aes(SPECTAXPRIO, y = N, fill = TOPO)) +
  geom_col(position = "dodge")+
  xlab("Espèce") +
  ylab("Nombre")  +
  coord_flip()+
  scale_fill_manual('Station',values = cbp2)+
  scale_y_continuous(breaks = round(seq(min(0), max(100), by = 5),1))+
  labs(
    title = "Nombre d'individus par espèces capturés sur chaque station",
    subtitle = "Division des espèces par station",
    caption = "Calcul réalisé grâce aux données de 2021, sur un total de 279 individus réparti en 34 espèces 
    + andrena sp., nomada sp. et le groupe des Terrestribombus. Réalisé grâce à RStudio,Version 1.3.1093") +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(face = ),
    plot.caption = element_text(face = "italic"))

####difference de nombre d'individu capturé par station 
stationnbr <- select (réserve, TOPO, SPECTAXPRIO,N)

stationnbr  <- aggregate(N ~ TOPO  , data =stationnbr , sum)

ggplot(stationnbr, aes(TOPO, y = N, fill = TOPO)) +
  geom_col(position = "dodge")+
  xlab("Station") +
  ylab("Nombre")  +
  coord_flip()+
  scale_y_continuous(breaks = round(seq(min(0), max(300), by = 8),1))+
  theme(legend.position = "none")+
  scale_fill_manual('Occupation du sol',values = cbp2)+
  labs(
    title = "Nombre d'individus capturés par station",
    subtitle = ,
    caption = "Calcul réalisé grâce aux données de 2021, sur un total de 279 individus.
    Réalisé grâce à RStudio,Version 1.3.1093") +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(face = ),
    plot.caption = element_text(face = "italic"))

####difference de nombre d'sp capturé par station 

stationsptot <- select (réserve, TOPO,N)
stationsptot <- distinct(stationsptot)
stationsptot  <- aggregate(N ~ TOPO , data =stationsptot , sum)

ggplot(stationsptot, aes(TOPO, y = N, fill = TOPO)) +
  geom_col(position = "dodge")+
  xlab("Station") +
  ylab("Nombre")  +
  scale_y_continuous(breaks = round(seq(min(0), max(100), by = 4),1))+
  theme(legend.position = "none")+
  scale_fill_manual('Occupation du sol',values = cbp2)+
  labs(
    title = "Nombre d'espèces capturées par station",
    subtitle = ,
    caption = "Calcul réalisé grâce aux données de 2021, sur un total de 34 espèces 
    + andrena sp., nomada sp. et le groupe des Terrestribombus. Réalisé grâce à RStudio,Version 1.3.1093") +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(face = ),
    plot.caption = element_text(face = "italic"))

######matrice
réserve <- read.csv ("réservetraitplante.csv",header = T, sep = ";")
réservetab <- select (réserve,SPECTAXPRIO,N,TOPO)

réservetab <- aggregate(N~SPECTAXPRIO+TOPO, data = réservetab, sum)
réservetab <- xtabs(N~TOPO+SPECTAXPRIO,réservetab)
réservetab <-type.convert(réservetab)


####diversité shannon par station 
Shannon <- diversity(réservetab, index = "shannon", MARGIN = 1, base = exp(1))
view(Shannon)
Shannon <-type.convert(Shannon)
fisher.alpha(réservetab, MARGIN = 1)
specnumber(réservetab, groups, MARGIN = 1)

#### boxplot shannon et t-test
shannon <- diversity(réservetab, index = "shannon", MARGIN = 1, base = exp(1))
view(Shannon)
shannon <-as.data.frame(shannon)
shannon <- bind_cols(shannon,TOPO)

ggplot(shannon,aes(TOPO,shannon))+
geom_boxplot()

t.test(data = crabs, rear ~ sex,
       alternative = "two.sided", conf.level = 0.95, var.equal = TRUE)

##### boxlot abondance et t-test
stationnbr <- select (réserve, TOPO, SPECTAXPRIO,N)


stationnbr  <- aggregate(N ~ TOPO  , data =stationnbr , sum)
TOPO <- select( stationnbr,TOPO)

ggplot(stationnbr,aes(TOPO,N))+
  geom_boxplot()

####boxplot richesse specifique et t-test
stationsp <- select (réserve, TOPO, SPECTAXPRIO,N)
stationsp  <- aggregate(N ~ TOPO +SPECTAXPRIO , data =stationsp , sum)

