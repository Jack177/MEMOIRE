
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
setwd("C:/Users/lalie/Desktop/mÈmoire/1mÈmoire-stat")
##################################### Data combinÈ abeille, station, plante, trait
####abeille
reese <- read.csv ( "reesett.csv",header = T, sep = ";")
ruelle <- read.csv ( "beeplanttraitstat.csv",header = T, sep = ";")
ruelle$trait_phÈnott_ID <-  as.character(ruelle$trait_phÈnott_ID)
martin <- read.csv ( "beestatplanttrait_mart.csv",header = T, sep = ";")
martin$MellifËre <-  as.logical(martin$MellifËre)
alexl <- read.csv ( "beestatplanttrait_alexl.csv",header = T, sep = ";")
alexl$IndigËne <-  as.logical(alexl$IndigËne)
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

stationtt$CONDGR2[stationtt$CONDGR2  %in% c("")]<-"Vol/Sol"
stationtt$CONDTAXPRIO[stationtt$CONDTAXPRIO  %in% c("")]<-"Vol/Sol"
stationtt$Lectisme[stationtt$Lectisme  %in% c("Polylectic b")]<-"Polylectic"
stationtt$Lectisme[stationtt$Lectisme  %in% c("Polylectic b")]<-"Polylectic"
stationtt$Lectisme[stationtt$Lectisme  %in% c("")]<-"Non connu"
stationtt$Sociality[stationtt$Sociality  %in% c("Primitivelyeusocial")]<-"Eusocial primitif"
stationtt$Sociality[stationtt$Sociality  %in% c("Social parasite")]<-"Parasite social"
stationtt$Sociality[stationtt$Sociality  %in% c("Solitary")]<-"Solitaire"
stationtt$Sociality[stationtt$Sociality  %in% c("Solitary+Primitivelyeusocial")]<-"Solitaire"
stationtt$Sociality[stationtt$Sociality  %in% c("Socialparasite")]<-"Parasite social"
stationtt$Nesting[stationtt$Nesting  %in% c("Carder","Renter:Existing_cavities_above_ground","Renter:Existingcavitiesaboveground", "Above")]<-"Au-dessus du sol"
stationtt$Nesting[stationtt$Nesting  %in% c("Excavator:Ground","Underground")]<-"Dans le sol"
stationtt$Nesting[stationtt$Nesting  %in% c("Wood/Existing cavities","Plant_stems/Existing_cavities","Plantstems/Existingcavities","Existingcavities")]<-"Mix"
stationtt$Nesting[stationtt$Nesting  %in% c("parasite","Social_parasite","Socialparasite","Social parasite")]<-"Parasite"
stationtt$Nesting[stationtt$Nesting  %in% c("")]<-"Non connu"

nombre <-  aggregate(N ~ SPECTAXPRIO, data =stationtt, sum)
trait <-read.table("trait_phÈno.csv", header=T, sep = ";")
beetrait <- left_join(nombre, trait, by = c("SPECTAXPRIO"="Ô..espece"))
###########plante
plantruelle <- read.csv ( "plantetoutstation.csv",header = T, sep = ";")
plantruelle$ID <-  as.character(plantruelle$ID)
plantruelle$StationID <-  as.character(plantruelle$StationID)
plantruelle$Temperature <-  as.character(plantruelle$Temperature)
plantereeese <- read.csv ( "plantestatalextt.csv",header = T, sep = ";")

plante <- bind_rows(plantruelle,plantereeese)

###matrice bee
prairie <- filter(stationtt, type == "Prairie")
beestattab <- select (prairie,SPECTAXPRIO,N,TOPO)

beestattab <- aggregate(N~SPECTAXPRIO+TOPO, data = beestattab, sum)
beestattab <- xtabs(N~TOPO+SPECTAXPRIO,beestattab)
beestattab <-type.convert(beestattab)


milieustattab <- select (stationtt,SPECTAXPRIO,N,type)

milieustattab <- aggregate(N~SPECTAXPRIO+type, data = milieustattab, sum)
milieustattab <- xtabs(N~type+SPECTAXPRIO,milieustattab)
milieustattab <-type.convert(milieustattab)

#####################################################################################lectisme 
#########par espece
#####pour tout

lecty <- select(prairie,Lectisme,SPECTAXPRIO)
lecty <- distinct(lecty)
lecty <- count(lecty,Lectisme)

lecty <- lecty %>% mutate(prop = n/ 85 *100) 
lecty <-lecty %>% mutate_if(is.numeric, ~round(., 1))


######pour le nombre d'individu
#####pour tout

lectytt <- select(prairie,Lectisme,N)
lectytt <- aggregate(N ~ Lectisme, data =lectytt, sum)

lectytt <- lectytt %>% mutate(prop = N/ sum(N) *100) 
lectytt <-lectytt %>% mutate_if(is.numeric, ~round(., 1))


###################################################################################### socialit√© 
############par sp
######pour tout
prairie$Sociality[prairie$Sociality  %in% c("Socialparasite")]<-"Social parasite"
prairie$Sociality[prairie$Sociality  %in% c("Social parasite")]<-"Social parasite"
prairie$Sociality[prairie$Sociality  %in% c("Parasite social")]<-"Social parasite"

social <- select(prairie,Sociality,SPECTAXPRIO)
social <- distinct(social)
social <- count(social,Sociality)
social<- arrange(social,n)

social <- social %>% mutate(prop = n/ 85 *100) 
social <-social %>% mutate_if(is.numeric, ~round(., 1))


#########  par individu
####pour tout

socialind <- select(prairie,Sociality,N)
socialind <- aggregate(N ~ Sociality, data =socialind, sum)
socialind<- arrange(socialind,N)

socialind <- socialind %>% mutate(prop = N/ sum(N) *100) 
socialind <-socialind %>% mutate_if(is.numeric, ~round(., 1))


################################################################################### nidifications
###########par sp
####pour tout

nid <- select(prairie,Nesting, SPECTAXPRIO)

nid <- distinct(nid)
nid <- count(nid,Nesting)
nid <- arrange(nid,n)

nid <- nid %>% mutate(prop = n/ 85 *100) 
nid <-nid %>% mutate_if(is.numeric, ~round(., 1))


#########par individu 
##### pour tout

nidi <- select(prairie,Nesting, N)
nidi <- aggregate(N ~ Nesting, data =nidi, sum)
nidi <- arrange(nidi,N)

nidi <- nidi %>% mutate(prop = N/ sum(N) *100) 
nidi <-nidi %>% mutate_if(is.numeric, ~round(., 1))

###################################################################################statut
#######################pour tout les station
####################par sp
statut<- select(prairie,SPECTAXPRIO,Statut_IUCN)
statut <- distinct(statut)
statut$Statut_IUCN[statut$Statut_IUCN  %in% c("lc")]<-"LC"

statut <- count(statut,Statut_IUCN)
statut <- arrange(statut,n)

statut <- statut %>% mutate(prop = n/ sum(n) *100) 
statut <-statut %>% mutate_if(is.numeric, ~round(., 1))

statut$Statut_IUCN[statut$Statut_IUCN  %in% c("LC")]<-"Least Concern"
statut$Statut_IUCN[statut$Statut_IUCN  %in% c("VU")]<-"Vulnerable"
statut$Statut_IUCN[statut$Statut_IUCN  %in% c("NT")]<-"Near Threatened"
statut$Statut_IUCN[statut$Statut_IUCN  %in% c("DD")]<-"Data Deficient"
statut$Statut_IUCN[statut$Statut_IUCN  %in% c("CR")]<-"Critically Endangered"
statut$Statut_IUCN[statut$Statut_IUCN  %in% c("RE")]<-"Regionally extinct"
statut$Statut_IUCN[statut$Statut_IUCN  %in% c("EN")]<-" Endangered"
statut$Statut_IUCN[statut$Statut_IUCN  %in% c("")]<-"Non connu"

####Par individu

statuti<- select(prairie,SPECTAXPRIO,Statut_IUCN,N)
statuti$Statut_IUCN[statuti$Statut_IUCN  %in% c("lc")]<-"LC"
statuti <- aggregate(N ~ Statut_IUCN, data =statuti, sum)
statuti <- arrange(statuti,N)

statuti <- statuti %>% mutate(prop = N/ sum(N) *100) 
statuti <-statuti %>% mutate_if(is.numeric, ~round(., 1))


# statut en danger par famille par espece
install.packages( "here")
library("here")
dir.create("data_output") # pour crÈer le dossier data_output
# Exportation
write.csv2(statfam, here::here("data_output","liste_sp.csv"),row.names = FALSE ) 

statfam <- select (prairie,Statut_IUCN,SPECGR2,SPECTAXPRIO)
statfam$Statut_IUCN[statfam$Statut_IUCN  %in% c("lc")]<-"LC"


statfam$Statut_IUCN[statfam$Statut_IUCN  %in% c("LC")]<-"Least Concern"
statfam$Statut_IUCN[statfam$Statut_IUCN %in% c("VU")]<-"Vulnerable"
statfam$Statut_IUCN[statfam$Statut_IUCN  %in% c("NT")]<-"Near Threatened"
statfam$Statut_IUCN[statfam$Statut_IUCN  %in% c("DD")]<-"Data Deficient"
statfam$Statut_IUCN[statfam$Statut_IUCN  %in% c("CR")]<-"Critically Endangered"
statfam$Statut_IUCN[statfam$Statut_IUCN  %in% c("RE")]<-"Regionally extinct"
statfam$Statut_IUCN[statfam$Statut_IUCN  %in% c("EN")]<-"Endangered"
statfam$Statut_IUCN[statfam$Statut_IUCN  %in% c("")]<-"Non connu"
statfam$Statut_IUCN[statfam$Statut_IUCN  %in% c("lc")]<-"LC"

statfam <-statfam %>% filter(SPECTAXPRIO != 'Andrena  sp.')
statfam <- distinct (statfam)
statfam <- count( statfam,Statut_IUCN,SPECGR2 )
statfam <- aggregate(n ~ Statut_IUCN +SPECGR2, data =statfam, sum)
statfam <- statfam %>% group_by(SPECGR2) %>% mutate(sum(n))
statfam <- statfam  %>% rename (c( "nombre" = "sum(n)"))
statfam <- xtabs(n~SPECGR2+Statut_IUCN ,statfam)
statfam <-type.convert(statfam)


statfam <- select (prairie,Statut_IUCN,type,SPECTAXPRIO)
statfam$Statut_IUCN[statfam$Statut_IUCN  %in% c("lc")]<-"LC"
statfam <- filter (statfam, Statut_IUCN !="LC")
statfam <- filter (statfam, Statut_IUCN !="DD")
statfam <- filter (statfam, Statut_IUCN !="")

statfam$Statut_IUCN[statfam$Statut_IUCN  %in% c("LC")]<-"Least Concern"
statfam$Statut_IUCN[statfam$Statut_IUCN %in% c("VU")]<-"Vulnerable"
statfam$Statut_IUCN[statfam$Statut_IUCN  %in% c("NT")]<-"Near Threatened"
statfam$Statut_IUCN[statfam$Statut_IUCN  %in% c("DD")]<-"Data Deficient"
statfam$Statut_IUCN[statfam$Statut_IUCN  %in% c("CR")]<-"Critically Endangered"
statfam$Statut_IUCN[statfam$Statut_IUCN  %in% c("RE")]<-"Regionally extinct"
statfam$Statut_IUCN[statfam$Statut_IUCN  %in% c("EN")]<-"Endangered"
statfam$Statut_IUCN[statfam$Statut_IUCN  %in% c("")]<-"Non connu"
statfam$Statut_IUCN[statfam$Statut_IUCN  %in% c("lc")]<-"LC"

statfam <-statfam %>% filter(SPECTAXPRIO != 'Andrena  sp.')
statfam <- distinct (statfam)
statfam <- count( statfam,Statut_IUCN,type )
statfam <- aggregate(n ~ Statut_IUCN +type, data =statfam, sum)
statfam <- statfam %>% group_by(type) %>% mutate(sum(n))
statfam <- statfam  %>% rename( c("sum(n)" = "nombre"))

statfam <- xtabs(n~type+Statut_IUCN ,statfam)
statfam <-type.convert(statfam)

ggplot(statfam,aes(reorder(SPECGR2,nombre),n,fill = Statut_IUCN))   +
  geom_col(color="white" )+ 
  xlab("Famille") +
  ylab("Nombre d'espËces")+  
  coord_flip()+
  scale_fill_manual("Statut IUCN",values = cbp2)+
  scale_y_continuous(breaks = round(seq(min(0), max(240), by = 1),1))+
  labs(
    title = "Nombre d'espËces en danger par famille",
    subtitle = "",
    caption = "") +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(face = "bold"),
    plot.caption = element_text(face = "italic"))

##################################### voir ou se trouve les en danger

stattopo <- select (stationtt,Statut_IUCN,SPECGR2,SPECTAXPRIO,TOPO,type)
stattopo <- filter(stattopo, type == "Prairie")

stattopo$Statut_IUCN[stattopo$Statut_IUCN  %in% c("lc")]<-"LC"
stattopo <- filter (stattopo, Statut_IUCN !="LC")
stattopo <- filter (stattopo, Statut_IUCN !="DD")
stattopo <- filter (stattopo, Statut_IUCN !="")
stattopo <- filter (stattopo, Statut_IUCN !="NT")

stattopo$Statut_IUCN[stattopo$Statut_IUCN  %in% c("LC")]<-"Least Concern"
stattopo$Statut_IUCN[stattopo$Statut_IUCN %in% c("VU")]<-"Vulnerable"
stattopo$Statut_IUCN[stattopo$Statut_IUCN  %in% c("NT")]<-"Near Threatened"
stattopo$Statut_IUCN[stattopo$Statut_IUCN  %in% c("DD")]<-"Data Deficient"
stattopo$Statut_IUCN[stattopo$Statut_IUCN  %in% c("CR")]<-"Critically Endangered"
stattopo$Statut_IUCN[stattopo$Statut_IUCN  %in% c("RE")]<-"Regionally extinct"
stattopo$Statut_IUCN[stattopo$Statut_IUCN  %in% c("EN")]<-"Endangered"
stattopo$Statut_IUCN[stattopo$Statut_IUCN  %in% c("")]<-"Non connu"
stattopo$Statut_IUCN[stattopo$Statut_IUCN  %in% c("lc")]<-"LC"

stattopo <-stattopo %>% filter(SPECTAXPRIO != 'Andrena  sp.')
stattopo <- distinct (stattopo)
stattopo <- count( stattopo,Statut_IUCN,SPECGR2,TOPO )
stattopo <- aggregate(n ~ Statut_IUCN +TOPO, data =stattopo, sum)
stattopo <- stattopo %>% group_by(TOPO) %>% mutate(sum(n))
stattopo <- rename(stattopo,c('sum(n)'='nombre'))

couleur <-stattopo$couleur
ggplot(stattopo,aes(reorder(TOPO,nombre),n,fill = Statut_IUCN))   +
  geom_col(color="white", position = "dodge" )+ 
  xlab("Station") +
  ylab("Nombre d'espËces")+  
  coord_flip()+
  scale_fill_manual("Statut IUCN",values = c("#FF3030","#FFD700","purple","#76EE00"))+
  geom_text(aes(label = n,y = n+0.1), size =3.5,position=position_dodge(width=1))+
  scale_y_continuous(breaks = round(seq(min(0), max(240), by = 1),1))+
  labs(
    title = "Nombre d'espËces en danger par station",
    subtitle = "",
    caption = "") +
  theme(
    plot.title = element_text(hjust = 0.5, size=18),
    plot.subtitle = element_text(face = ),
    plot.caption = element_text(face = "italic", size = 9),
    axis.title.x =element_text(size = 15 ), 
    axis.title.y =element_text(size = 15),
    axis.text.x=element_text(size = 12),
    axis.text.y=element_text(size = 13),
    legend.text =element_text(size = 12),
    legend.title =element_text(size = 14))

stattopo <- select (stationtt,Statut_IUCN,N,SPECTAXPRIO,type)

stattopo$Statut_IUCN[stattopo$Statut_IUCN  %in% c("lc")]<-"LC"
stattopo <- filter (stattopo, Statut_IUCN !="LC")
stattopo <- filter (stattopo, Statut_IUCN !="DD")
stattopo <- filter (stattopo, Statut_IUCN !="")
stattopo <- filter (stattopo, Statut_IUCN !="NT")


stattopo$Statut_IUCN[stattopo$Statut_IUCN  %in% c("LC")]<-"Least Concern"
stattopo$Statut_IUCN[stattopo$Statut_IUCN %in% c("VU")]<-"Vulnerable"
stattopo$Statut_IUCN[stattopo$Statut_IUCN  %in% c("NT")]<-"Near Threatened"
stattopo$Statut_IUCN[stattopo$Statut_IUCN  %in% c("DD")]<-"Data Deficient"
stattopo$Statut_IUCN[stattopo$Statut_IUCN  %in% c("CR")]<-"Critically Endangered"
stattopo$Statut_IUCN[stattopo$Statut_IUCN  %in% c("RE")]<-"Regionally extinct"
stattopo$Statut_IUCN[stattopo$Statut_IUCN  %in% c("EN")]<-"Endangered"
stattopo$Statut_IUCN[stattopo$Statut_IUCN  %in% c("")]<-"Non connu"


stattopo <-stattopo %>% filter(SPECTAXPRIO != 'Andrena  sp.')


stattopo <- stattopo %>% group_by(type) %>% mutate(sum(N))
stattopo <- rename(stattopo,c('sum(N)'='nombre'))
stattopo <- aggregate(N ~ Statut_IUCN +type+nombre, data =stattopo, sum)

ggplot(stattopo,aes(reorder(type,nombre),N,fill = Statut_IUCN)) +
  geom_col(position = "dodge")+ 
  xlab("Type de milieu") +
  ylab("Nombre d'individu")+  
  coord_flip()+
  scale_fill_manual("Statut IUCN",values = c("#FF3030","#FFD700","purple","#76EE00"))+
  geom_text(aes(label = N,y = N+0.3), size =4,position=position_dodge(width=1))+
  scale_y_continuous(breaks = round(seq(min(0), max(240), by = 2),1))+
  labs(
    title = "Nombre d'individu en danger par milieu",
    subtitle = "",
    caption = "") +
  theme(
    plot.title = element_text(hjust = 0.5, size=18),
    plot.subtitle = element_text(face = ),
    plot.caption = element_text(face = "italic", size = 9),
    axis.title.x =element_text(size = 15 ), 
    axis.title.y =element_text(size = 15),
    axis.text.x=element_text(size = 12),
    axis.text.y=element_text(size = 13),
    legend.text =element_text(size = 12),
    legend.title =element_text(size = 14))


stattopo <- select (stationtt,Statut_IUCN,SPECTAXPRIO,type,couleur)
stattopo <-stattopo %>% filter(SPECTAXPRIO != 'Andrena  sp.')
stattopo <- distinct (stattopo)
stattopo <- count( stattopo,Statut_IUCN,SPECTAXPRIO,type )

stattopo <- stattopo %>%group_by(type)%>% mutate(sum(n)) 
stattopo$Statut_IUCN[stattopo$Statut_IUCN  %in% c("lc")]<-"LC"
stattopo <- filter (stattopo, Statut_IUCN !="LC")
stattopo <- filter (stattopo, Statut_IUCN !="DD")
stattopo <- filter (stattopo, Statut_IUCN !="")
stattopo <- filter (stattopo, Statut_IUCN !="NT")


stattopo$Statut_IUCN[stattopo$Statut_IUCN  %in% c("LC")]<-"Least Concern"
stattopo$Statut_IUCN[stattopo$Statut_IUCN %in% c("VU")]<-"Vulnerable"
stattopo$Statut_IUCN[stattopo$Statut_IUCN  %in% c("NT")]<-"Near Threatened"
stattopo$Statut_IUCN[stattopo$Statut_IUCN  %in% c("DD")]<-"Data Deficient"
stattopo$Statut_IUCN[stattopo$Statut_IUCN  %in% c("CR")]<-"Critically Endangered"
stattopo$Statut_IUCN[stattopo$Statut_IUCN  %in% c("RE")]<-"Regionally extinct"
stattopo$Statut_IUCN[stattopo$Statut_IUCN  %in% c("EN")]<-"Endangered"
stattopo$Statut_IUCN[stattopo$Statut_IUCN  %in% c("")]<-"Non connu"


stattopo <-stattopo %>% filter(SPECTAXPRIO != 'Andrena  sp.')
stattopo <- rename(stattopo,c('sum(n)'='nombre'))
stattopo <- aggregate(n ~ Statut_IUCN +type+nombre, data =stattopo, sum)

stattopo <- stattopo %>%group_by(type)%>% mutate(prop = n/ nombre *100) 
stattopo<-stattopo %>% mutate_if(is.numeric, ~round(., 1))
somme <- aggregate(prop ~type, data =stattopo, sum)
stattopo<- left_join (stattopo,somme, by ="type")

ggplot(stattopo,aes(reorder(type,prop.y),prop.x,fill = Statut_IUCN))   +
  geom_col(color="white", position = "stack")+ 
  xlab("Type de milieu") +
  ylab("Pourcentage d'espËces")+  
  coord_flip()+
  scale_fill_manual("Statut IUCN",values = c("#FF3030","#FFD700","purple","#76EE00"))+
  geom_text(aes(label = prop.x,y = prop.x-0.1), size =3.5,position=position_stack())+
  scale_y_continuous(breaks = round(seq(min(0), max(240), by = 1),1))+
  labs(
    title = "Poucentage d'espËces en danger par station",
    subtitle = "",
    caption = "Seul les catÈgories en danger ont ÈtÈ reprÈsentÈes") +
  theme(
    plot.title = element_text(hjust = 0.5, size=18),
    plot.subtitle = element_text(face = ),
    plot.caption = element_text(face = "italic", size = 9),
    axis.title.x =element_text(size = 15 ), 
    axis.title.y =element_text(size = 15),
    axis.text.x=element_text(size = 12),
    axis.text.y=element_text(size = 13),
    legend.text =element_text(size = 12),
    legend.title =element_text(size = 14))





stattopo <- select (stationtt,Statut_IUCN,SPECTAXPRIO,type,couleur)

stattopo <-stattopo %>% filter(SPECTAXPRIO != 'Andrena  sp.')
stattopo <- distinct (stattopo)
stattopo <- count( stattopo,Statut_IUCN,SPECTAXPRIO,type )
stattopo <- stattopo %>%group_by(type)%>% mutate(sum(n)) 


stattopo <- stattopo %>%group_by(type)%>% mutate(sum(n)) 
stattopo$Statut_IUCN[stattopo$Statut_IUCN  %in% c("lc")]<-"LC"
stattopo <- filter (stattopo, Statut_IUCN !="LC")
stattopo <- filter (stattopo, Statut_IUCN !="DD")
stattopo <- filter (stattopo, Statut_IUCN !="")
stattopo <- filter (stattopo, Statut_IUCN !="NT")


stattopo$Statut_IUCN[stattopo$Statut_IUCN  %in% c("LC")]<-"Least Concern"
stattopo$Statut_IUCN[stattopo$Statut_IUCN %in% c("VU")]<-"Vulnerable"
stattopo$Statut_IUCN[stattopo$Statut_IUCN  %in% c("NT")]<-"Near Threatened"
stattopo$Statut_IUCN[stattopo$Statut_IUCN  %in% c("DD")]<-"Data Deficient"
stattopo$Statut_IUCN[stattopo$Statut_IUCN  %in% c("CR")]<-"Critically Endangered"
stattopo$Statut_IUCN[stattopo$Statut_IUCN  %in% c("RE")]<-"Regionally extinct"
stattopo$Statut_IUCN[stattopo$Statut_IUCN  %in% c("EN")]<-"Endangered"
stattopo$Statut_IUCN[stattopo$Statut_IUCN  %in% c("")]<-"Non connu"

stattopo <- aggregate(n ~ Statut_IUCN +SPECTAXPRIO+type, data =stattopo, sum)
stattopo <- stattopo %>% group_by(type) %>% mutate(sum(n))
stattopo <- stattopo  %>% rename(c( "sum(n)" = "nombre"))

stattopo <- stattopo %>%group_by(type)%>% mutate(prop = n/ nombre *100) 
stattopo<-stattopo %>% mutate_if(is.numeric, ~round(., 1))
somme <- aggregate(prop ~type, data =stattopo, sum)
stattopo<- left_join (stattopo,somme, by ="type")

#stattopo <- stattopo %>% group_by(type) %>% mutate(sum(N))
stattopo <- rename(stattopo,c('sum(N)'='nombre'))



ggplot(stattopo,aes(reorder(type,prop.y),prop.x,fill = Statut_IUCN))   +
  geom_col(color="white", position = "dodge")+ 
  xlab("Type de milieu") +
  ylab("Pourcentage d'espËces")+  
  coord_flip()+
  scale_fill_manual("Statut IUCN",values = c("#FF3030","#FFD700","purple","#76EE00"))+
  geom_text(aes(label = prop.x,y = prop.x+0.1), size =3.5,position=position_dodge(width=1))+
  scale_y_continuous(breaks = round(seq(min(0), max(240), by = 1),1))+
  labs(
    title = "Poucentage d'espËces en danger par station",
    subtitle = "",
    caption = "Seul les catÈgories en danger ont ÈtÈ reprÈsentÈes") +
  theme(
    plot.title = element_text(hjust = 0.5, size=18),
    plot.subtitle = element_text(face = ),
    plot.caption = element_text(face = "italic", size = 9),
    axis.title.x =element_text(size = 15 ), 
    axis.title.y =element_text(size = 15),
    axis.text.x=element_text(size = 12),
    axis.text.y=element_text(size = 13),
    legend.text =element_text(size = 12),
    legend.title =element_text(size = 14))




type <- 
#####################################################################################lectisme 
#########par espece
#####pour tout

lecty <- select(stationtt,Lectisme,SPECTAXPRIO, type)
lecty <- distinct(lecty)
lecty <- count(lecty,Lectisme,type)

lecty <- lecty %>%group_by(type)%>%mutate(prop = n/sum(n) *100) 
lecty <-lecty %>% mutate_if(is.numeric, ~round(., 1))
lecty <- select(lecty,-n)

lecty <- xtabs(prop~type+Lectisme,lecty)
lecty <-type.convert(lecty)
write.csv2(lecty, here::here("data_output","lecty.csv"),row.names = FALSE ) 
######pour le nombre d'individu
#####pour tout

lectytt <- select(prairie,Lectisme,N)
lectytt <- aggregate(N ~ Lectisme, data =lectytt, sum)

lectytt <- lectytt %>% mutate(prop = N/ sum(N) *100) 
lectytt <-lectytt %>% mutate_if(is.numeric, ~round(., 1))


###################################################################################### socialit√© 
############par sp
######pour tout
stationtt$Sociality[stationtt$Sociality  %in% c("Socialparasite")]<-"Social parasite"
stationtt$Sociality[stationtt$Sociality  %in% c("Social parasite")]<-"Social parasite"
stationtt$Sociality[stationtt$Sociality  %in% c("Parasite social")]<-"Social parasite"

social <- select(stationtt,Sociality,SPECTAXPRIO, type)
social <- distinct(social)
social <- count(social,Sociality,type)

social <- social %>%group_by(type)%>%mutate(prop = n/sum(n) *100) 
social <-social %>% mutate_if(is.numeric, ~round(., 1))
social <- select(social,-n)

social <- xtabs(prop~type+Sociality,social)
social <-type.convert(social)
write.csv2(social, here::here("data_output","social.csv"),row.names = FALSE )
#########  par individu
####pour tout

socialind <- select(prairie,Sociality,N)
socialind <- aggregate(N ~ Sociality, data =socialind, sum)
socialind<- arrange(socialind,N)

socialind <- socialind %>% mutate(prop = N/ sum(N) *100) 
socialind <-socialind %>% mutate_if(is.numeric, ~round(., 1))


################################################################################### nidifications
###########par sp
####pour tout

nid <- select(stationtt,Nesting,SPECTAXPRIO, type)
nid <- distinct(nid)
nid <- count(nid,Nesting,type)

nid <- nid %>%group_by(type)%>%mutate(prop = n/sum(n) *100) 
nid <-nid %>% mutate_if(is.numeric, ~round(., 1))
nid <- select(nid,-n)

nid <- xtabs(prop~type+Nesting,nid)
nid <-type.convert(nid)
write.csv2(nid, here::here("data_output","nid.csv"),row.names = FALSE )
#########  par individu

#########par individu 
##### pour tout

nidi <- select(prairie,Nesting, N)
nidi <- aggregate(N ~ Nesting, data =nidi, sum)
nidi <- arrange(nidi,N)

nidi <- nidi %>% mutate(prop = N/ sum(N) *100) 
nidi <-nidi %>% mutate_if(is.numeric, ~round(., 1))

parc <- select ( stationtt, Lectisme, SPECTAXPRIO, type,TOPO)
parc <- filter (parc, type =="Parc")
parc <- distinct(parc)





prairie$Statut_IUCN[prairie$Statut_IUCN  %in% c("lc")]<-"LC"
prairie <- filter (prairie, Statut_IUCN !="LC")
prairie <- filter (prairie, Statut_IUCN !="DD")
prairie <- filter (prairie, Statut_IUCN !="")
prairie <- filter (prairie, Statut_IUCN !="NT")
