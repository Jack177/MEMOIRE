

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

cbp1 <- c( "#56B4E9", "#009E73",
           "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#E69F00","#999999")
cbp2 <- c( "#56B4E9", "#009E73",
           "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#E69F00","#8B6914","#999999",
           '#8B1C62','#32CD32',"#0000EE","#FF407A")



################################################### trait phéno, statut, fam
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
beetotal$Lectisme[beetotal$Lectisme  %in% c("Polylectic b")]<-"Polylectique"
beetotal$Lectisme[beetotal$Lectisme  %in% c("Polylectic")]<-"Polylectique"
beetotal$Lectisme[beetotal$Lectisme  %in% c("Oligolectic")]<-"Oligolectique"
beetotal$Lectisme[beetotal$Lectisme  %in% c("")]<-"Non connu"
beetotal$Sociality[beetotal$Sociality  %in% c("Primitivelyeusocial")]<-"Eusocial primitif"
beetotal$Sociality[beetotal$Sociality  %in% c("Social parasite")]<-"Parasite social"
beetotal$Sociality[beetotal$Sociality  %in% c("Solitary")]<-"Solitaire"
beetotal$Sociality[beetotal$Sociality  %in% c("Solitary+Primitivelyeusocial")]<-"Solitaire"
beetotal$Sociality[beetotal$Sociality  %in% c("Socialparasite")]<-"Parasite social"
beetotal$Nesting[beetotal$Nesting  %in% c("Carder","Renter:Existing_cavities_above_ground","Renter:Existingcavitiesaboveground", "Above")]<-"Au-dessus du sol"
beetotal$Nesting[beetotal$Nesting  %in% c("Excavator:Ground","Underground")]<-"Dans le sol"
beetotal$Nesting[beetotal$Nesting  %in% c("Wood/Existing cavities","Plant_stems/Existing_cavities","Plantstems/Existingcavities","Existingcavities")]<-"Mix"
beetotal$Nesting[beetotal$Nesting  %in% c("parasite","Social_parasite","Socialparasite","Social parasite")]<-"Parasite"


#camenbert lectisme 

lecty <- select(beetotal,Lectisme,SPECTAXPRIO)
lecty <- distinct(lecty)
lecty <- count(lecty,Lectisme)

lecty <- lecty %>% 
  arrange(desc(Lectisme)) %>%
  mutate(prop = n/ sum(lecty$n) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )

lecty <- lecty %>% mutate(prop = n/ sum(n) *100) 
lecty <-lecty %>% mutate_if(is.numeric, ~round(., 1))

ggplot(lecty,aes(x = "", y = prop, fill = Lectisme))  +
  geom_bar(width = 1,stat ="identity",color="white" )+ 
  coord_polar("y", start = 0)+
  xlab("Classe")+
  ylab("Effectifs")+
  theme_void() +
  geom_text(aes(y = ypos, label = prop), color = "black", size=6)+
  scale_fill_manual("Régime alimentaire",values = cbp2)+
  labs(
    title = "Régime alimentaire des espèces d'abeilles 
  présentes sur le site",
    subtitle = "Pourcentage d'espèces selon leur régime alimentaire",
    caption = "Calcul réalisé sur les données des années 2019,2020 et 2021, pour un total de 57 espèces +
le groupe des Terrestrisbombus (Bombus terrestris, lucorum, magnus).
  Graphique réalisé grâce à RStudio,Version 1.3.1093") +
  theme(
    plot.title = element_text(hjust = 0.5,size = 20),
    plot.subtitle = element_text(hjust = 0.5, size = 15),
    plot.caption = element_text(face = "italic", size = 11,hjust = 0.5),
    legend.title = element_text(size=18),
    legend.text = element_text(size=15))

######pour le nombre d'individu

lectytt <- select(beetotal,Lectisme,N)
lectytt <- aggregate(N ~ Lectisme, data =lectytt, sum)

lectytt <- lectytt %>% 
  arrange(desc(Lectisme)) %>%
  mutate(prop = N/ sum(lectytt$N) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )

lectytt <- lectytt %>% mutate(prop = N/ sum(N) *100) 
lectytt <-lectytt %>% mutate_if(is.numeric, ~round(., 1))

ggplot(lectytt,aes(x = "", y = prop, fill = Lectisme))  +
  geom_bar(width = 1,stat ="identity",color="white" )+ 
  coord_polar("y", start = 0)+
  xlab("Classe")+
  ylab("Effectifs")+
  theme_void() +
  geom_text(aes(y = ypos, label = prop), color = "black", size=6)+
  scale_fill_manual("Régime alimentaire",values = cbp2)+
  labs(
    title = "Régime alimentaire des abeilles 
  présentes sur le site",
    subtitle = "Pourcentage d'individus selon leur régime alimentaire",
    caption = "Calcul réalisé sur les données des années 2019,2020 et 2021, pour un total de 520 individus.
  Graphique réalisé grâce à RStudio,Version 1.3.1093") +
  theme(
    plot.title = element_text(hjust = 0.5,size = 20),
    plot.subtitle = element_text(hjust = 0.5, size = 15),
    plot.caption = element_text(face = "italic", size = 11,hjust = 0.5),
    legend.title = element_text(size=18),
    legend.text = element_text(size=15))

######Pour tout les sites 
#### par individu
stationtt <- read.csv ( "beeplanttrait.csv",header = T, sep = ";")
stationtt$Lectisme[stationtt$Lectisme  %in% c("Polylectic b")]<-"Polylectique"
stationtt$Lectisme[stationtt$Lectisme  %in% c("Polylectic")]<-"Polylectique"
stationtt$Lectisme[stationtt$Lectisme  %in% c("Oligolectic")]<-"Oligolectique"
stationtt$Lectisme[stationtt$Lectisme  %in% c("")]<-"Non connu"
lectystatt <- select(stationtt,Lectisme,N)
lectystatt[is.na(lectystatt)] <- "unknow"
lectystatt <- aggregate(N ~ Lectisme, data =lectystatt, sum)

lectystatt <- lectystatt %>% 
  arrange(desc(Lectisme)) %>%
  mutate(prop = N/ sum(lectystatt$N) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )

lectystatt <- lectystatt %>% mutate(prop = N/ sum(N) *100) 
lectystatt <-lectystatt %>% mutate_if(is.numeric, ~round(., 1))

b <-ggplot(lectystatt,aes(x = "", y = prop, fill = Lectisme))  +
  geom_bar(width = 1,stat ="identity",color="white" )+ 
  coord_polar("y", start = 0)+
  xlab("Classe")+
  ylab("Effectifs")+
  theme_void() +
  geom_text(aes(y = ypos, label = prop), color = "black", size=6)+
  scale_fill_manual("Régime alimentaire",values = c("#F0E442","#56B4E9", "#009E73"))
  labs(
    title = "Régime alimentaire des abeilles 
présentes sur tous les sites",
    subtitle = "Pourcentage d'individus selon leur régime alimentaire",
    caption = "Calcul réalisé sur les données de 2020, collectée sur 27 sites de la région de Mons, 
    pour un total de 2909 individus.Graphique réalisé grâce à RStudio,Version 1.3.1093") +
  theme(
    plot.title = element_text(hjust = 0.5,size = 20),
    plot.subtitle = element_text(hjust = 0.5, size = 15),
    plot.caption = element_text(face = "italic", size = 11,hjust = 0.5),
    legend.title = element_text(size=18),
    legend.text = element_text(size=15))

#### sp
stationtt <- read.csv ( "beeplanttrait.csv",header = T, sep = ";")
stationtt$Lectisme[stationtt$Lectisme  %in% c("Polylectic b")]<-"Polylectique"
stationtt$Lectisme[stationtt$Lectisme  %in% c("Polylectic")]<-"Polylectique"
stationtt$Lectisme[stationtt$Lectisme  %in% c("Oligolectic")]<-"Oligolectique"
stationtt$Lectisme[stationtt$Lectisme  %in% c("")]<-"Non connu"
lec <- select(stationtt,Lectisme,SPECTAXPRIO)
lec <- distinct(lec)
lec <- count(lec,Lectisme)

lec <- lec %>% 
  arrange(desc(Lectisme)) %>%
  mutate(prop = n/ sum(lec$n) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )

lec <- lec %>% mutate(prop = n/ sum(n) *100) 
lec <-lec %>% mutate_if(is.numeric, ~round(., 1))

d <-ggplot(lec,aes(x = "", y = prop, fill = Lectisme))  +
  geom_bar(width = 1,stat ="identity",color="white" )+ 
  coord_polar("y", start = 0)+
  xlab("Classe")+
  ylab("Effectifs")+
  theme_void() +
  geom_text(aes(y = ypos, label = prop), color = "black", size=6)+
  scale_fill_manual("Régime alimentaire",values = c("#F0E442","#56B4E9", "#009E73"))
  labs(
    title = "Régime alimentaire des espèces d'abeilles 
  présentes sur tous les site",
    subtitle = "Pourcentage d'espèces selon leur régime alimentaire",
    caption = "Calcul réalisé sur les données de 2020, collectée sur 27 sites de la région de Mons, 
    pour un total de 116 espèces.Graphique réalisé grâce à RStudio,Version 1.3.1093") +
  theme(
    plot.title = element_text(hjust = 0.5,size = 20),
    plot.subtitle = element_text(hjust = 0.5, size = 15),
    plot.caption = element_text(face = "italic", size = 11,hjust = 0.5),
    legend.title = element_text(size=18),
    legend.text = element_text(size=15))

###combiner les graph sur le lectisme

##ggarrange(a, b, c, d, ncol = 2, nrow = 2,common.legend = TRUE,legend = 'right')
  
##annotate_figure(figure,
 #               top.left = text_grob("Sur le site de la réserve", size = 14),
  #              top.right =text_grob( "Sur tout les site", size = 14),
  #              bottom = text_grob("Data source: \n ToothGrowth data set", color = "blue",
  #                                  hjust = 1, x = 1, face = "italic", size = 10),
  #              left = text_grob("Par individu          Par espèce", rot = 90),
  #              right = text_grob(bquote("Superscript: ("*kg~NH[3]~ha^-1~yr^-1*")"), rot = 90))

#camenbert socialité par sp

social <- select(beetotal,Sociality,SPECTAXPRIO)
social <- distinct(social)
social <- count(social,Sociality)
social<- arrange(social,n)

social <- social %>% 
  arrange(desc(Sociality)) %>%
  mutate(prop = n/ sum(social$n) *100) %>%
  mutate(ypos = cumsum(prop)-0.5*prop )

social <- social %>% mutate(prop = n/ sum(n) *100) 
social <-social %>% mutate_if(is.numeric, ~round(., 1))

ggplot(social,aes(x = "", y = prop, fill = Sociality))  +
  geom_bar(width = 1,stat ="identity",color="white" )+ 
  coord_polar("y", start = 0)+
  xlab("Classe")+
  ylab("Effectifs")+
  theme_void() +
  scale_fill_manual("Socialité",values = cbp2)+
  geom_text(aes(y = ypos, label = prop), color = "black", size=5)+
  labs(
    title = "Socialité des espèces d'abeilles 
  présentes sur le site",
    subtitle = "Pourcentage d'espèces selon leur socialité",
    caption = "Calcul réalisé sur les données des années 2019,2020 et 2021, pour un total de 57 espèces +
le groupe des Terrestrisbombus (Bombus terrestris, lucorum, magnus).
  Graphique réalisé grâce à RStudio,Version 1.3.1093") +
  theme(
    plot.title = element_text(hjust = 0.5,size = 20),
    plot.subtitle = element_text(hjust = 0.5, size = 15),
    plot.caption = element_text(face = "italic", size = 11,hjust = 0.5),
    legend.title = element_text(size=18),
    legend.text = element_text(size=15))

##### social par individu

socialind <- select(beetotal,Sociality,N)
socialind <- aggregate(N ~ Sociality, data =socialind, sum)
socialind<- arrange(socialind,N)

socialind <- socialind %>% 
  arrange(desc(Sociality)) %>%
  mutate(prop = N/ sum(socialind$N) *100) %>%
  mutate(ypos = cumsum(prop)-0.5*prop )

socialind <- socialind %>% mutate(prop = N/ sum(N) *100) 
socialind <-socialind %>% mutate_if(is.numeric, ~round(., 1))

ggplot(socialind,aes(x = "", y = prop, fill = Sociality))  +
  geom_bar(width = 1,stat ="identity",color="white" )+ 
  coord_polar("y", start = 0)+
  xlab("Classe")+
  ylab("Effectifs")+
  theme_void() +
  scale_fill_manual("Socialité",values = cbp2)+
  geom_text(aes(y = ypos, label = prop), color = "black", size=5)+
  labs(
    title = "Socialité des abeilles présentes sur le site",
    subtitle = "Pourcentage d'individus selon leur socialité",
    caption = "Calcul réalisé sur les données des années 2019,2020 et 2021, pour un total de 520 individus .
  Graphique réalisé grâce à RStudio,Version 1.3.1093") +
  theme(
    plot.title = element_text(hjust = 0.5,size = 20),
    plot.subtitle = element_text(hjust = 0.5, size = 15),
    plot.caption = element_text(face = "italic", size = 11,hjust = 0.5),
    legend.title = element_text(size=18),
    legend.text = element_text(size=15))

###### social Pour tout les sites 
#### par individu
stationtt <- read.csv ( "beeplanttrait.csv",header = T, sep = ";")
stationtt$Sociality[stationtt$Sociality  %in% c("Primitivelyeusocial")]<-"Eusocial primitif"
stationtt$Sociality[stationtt$Sociality  %in% c("Socialparasite")]<-"Parasite social"
stationtt$Sociality[stationtt$Sociality  %in% c("Social parasite")]<-"Parasite social"
stationtt$Sociality[stationtt$Sociality  %in% c("Solitary")]<-"Solitaire"
stationtt$Sociality[stationtt$Sociality  %in% c("Solitary+Primitivelyeusocial")]<-"Solitaire"
stationtt$Sociality[stationtt$Sociality  %in% c("")]<-"Non connu"
socialstatt <- select(stationtt,Sociality,N)
socialstatt <- aggregate(N ~ Sociality, data =socialstatt, sum)

socialstatt <- socialstatt %>% 
  arrange(desc(Sociality)) %>%
  mutate(prop = N/ sum(socialstatt$N) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )

socialstatt <- socialstatt %>% mutate(prop = N/ sum(N) *100) 
socialstatt <-socialstatt %>% mutate_if(is.numeric, ~round(., 1))

ggplot(socialstatt,aes(x = "", y = prop, fill = Sociality))  +
  geom_bar(width = 1,stat ="identity",color="white" )+ 
  coord_polar("y", start = 0)+
  xlab("Classe")+
  ylab("Effectifs")+
  theme_void() +
  geom_text(aes(y = ypos, label = prop), color = "black", size=5)+
  scale_fill_manual("Socialité",values = cbp2)+
  labs(
    title = "Socialité des abeilles présentes sur tous les sites",
    subtitle = "Pourcentage d'individus selon leur socialité",
    caption = "Calcul réalisé sur les données de 2020, collectée sur 27 sites de la région de Mons, 
    pour un total de 2909 individus.Graphique réalisé grâce à RStudio,Version 1.3.1093") +
  theme(
    plot.title = element_text(hjust = 0.5,size = 20),
    plot.subtitle = element_text(hjust = 0.5, size = 15),
    plot.caption = element_text(face = "italic", size = 11,hjust = 0.5),
    legend.title = element_text(size=18),
    legend.text = element_text(size=15))

#### sp
stationtt <- read.csv ( "beeplanttrait.csv",header = T, sep = ";")
stationtt$Sociality[stationtt$Sociality  %in% c("Primitivelyeusocial")]<-"Eusocial primitif"
stationtt$Sociality[stationtt$Sociality  %in% c("Socialparasite")]<-"Parasite social"
stationtt$Sociality[stationtt$Sociality  %in% c("Social parasite")]<-"Parasite social"
stationtt$Sociality[stationtt$Sociality  %in% c("Solitary")]<-"Solitaire"
stationtt$Sociality[stationtt$Sociality  %in% c("Solitary+Primitivelyeusocial")]<-"Solitaire"
stationtt$Sociality[stationtt$Sociality  %in% c("")]<-"Non connu"
socialsta <- select(stationtt,Sociality,SPECTAXPRIO)
socialsta <- distinct(socialsta)
socialsta <- count(socialsta,Sociality)

socialsta <- socialsta %>% 
  arrange(desc(Sociality)) %>%
  mutate(prop = n/ sum(socialsta$n) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )

socialsta <- socialsta %>% mutate(prop = n/ sum(n) *100) 
socialsta <-socialsta %>% mutate_if(is.numeric, ~round(., 1))

ggplot(socialsta,aes(x = "", y = prop, fill = Sociality))  +
  geom_bar(width = 1,stat ="identity",color="white" )+ 
  coord_polar("y", start = 0)+
  xlab("Classe")+
  ylab("Effectifs")+
  theme_void() +
  geom_text(aes(y = ypos, label = prop), color = "black", size=5)+
  scale_fill_manual("Socialité",values = cbp2)+
  labs(
    title = "Socialité des espèces d'abeilles présentes 
sur tous les sites",
    subtitle = "Pourcentage d'espèces selon leur socialité",
    caption = "Calcul réalisé sur les données de 2020, collectée sur 27 sites de la région de Mons, 
    pour un total de 116 espèces.Graphique réalisé grâce à RStudio,Version 1.3.1093") +
  theme(
    plot.title = element_text(hjust = 0.5,size = 20),
    plot.subtitle = element_text(hjust = 0.5, size = 15),
    plot.caption = element_text(face = "italic", size = 11,hjust = 0.5),
    legend.title = element_text(size=18),
    legend.text = element_text(size=15))

###combiner les graph sur la socialité
#combine_charts(list(a, b, c, d))


#camenbert nidifications pour site
##par sp
nid <- select(beetotal,Nesting, SPECTAXPRIO)

nid <- distinct(nid)
nid <- count(nid,Nesting)
nid <- arrange(nid,n)

nid <- nid %>% 
  arrange(desc(Nesting)) %>%
  mutate(prop = n/ sum(nid$n) *100) %>%
  mutate(ypos = cumsum(prop)-0.5*prop )

nid <- nid %>% mutate(prop = n/ sum(n) *100) 
nid <-nid %>% mutate_if(is.numeric, ~round(., 1))

ggplot(nid,aes(x = "", y = prop, fill = Nesting))  +
  geom_bar(width = 1,stat ="identity",color="white" )+ 
  coord_polar("y", start = 0)+
  xlab("Classe")+
  ylab("Effectifs")+
  theme_void() +
  scale_fill_manual("Lieu de nidification",values = cbp2)+
  geom_text(aes(y = ypos, label = prop), color = "black", size=5)+
  labs(
    title = "Lieu de nidification des espèces d'abeilles 
présentes sur le site",
    subtitle = "Pourcentage d'espèces selon leur lieu de nidification",
    caption = "Calcul réalisé sur les données des années 2019,2020 et 2021, pour un total de 57 espèces +
le groupe des Terrestrisbombus (Bombus terrestris, lucorum, magnus).
  Graphique réalisé grâce à RStudio,Version 1.3.1093")  +
  theme(
    plot.title = element_text(hjust = 0.5,size = 20),
    plot.subtitle = element_text(hjust = 0.5, size = 15),
    plot.caption = element_text(face = "italic", size = 11,hjust = 0.5),
    legend.title = element_text(size=18),
    legend.text = element_text(size=15))

#####par individu 

nidi <- select(beetotal,Nesting, N)
nidi <- aggregate(N ~ Nesting, data =nidi, sum)
nidi <- arrange(nidi,N)

nidi <- nidi %>% 
  arrange(desc(Nesting)) %>%
  mutate(prop = N/ sum(nidi$N) *100) %>%
  mutate(ypos = cumsum(prop)-0.5*prop )

nidi <- nidi %>% mutate(prop = N/ sum(N) *100) 
nidi <-nidi %>% mutate_if(is.numeric, ~round(., 1))

ggplot(nidi,aes(x = "", y = prop, fill = Nesting))+
  geom_bar(width = 1,stat ="identity",color="white" )+ 
  coord_polar("y", start = 0)+
  xlab("Classe")+
  ylab("Effectifs")+
  theme_void() +
  scale_fill_manual("Lieu de nidification",values = cbp2)+
  geom_text(aes(y = ypos, label = prop), color = "black", size=5)+
  labs(
    title = "Lieu de nidification des abeilles 
présentes sur le site",
    subtitle = "Pourcentage d'abeilles selon leur lieu de nidification",
    caption = "Calcul réalisé sur les données des années 2019,2020 et 2021, pour un total de 520 individus .
  Graphique réalisé grâce à RStudio,Version 1.3.1093")  +
  theme(
    plot.title = element_text(hjust = 0.5,size = 20),
    plot.subtitle = element_text(hjust = 0.5, size = 15),
    plot.caption = element_text(face = "italic", size = 11,hjust = 0.5),
    legend.title = element_text(size=18),
    legend.text = element_text(size=15))

###### Nid Pour tout les sites 
#### par individu
stationtt <- read.csv ( "beeplanttrait.csv",header = T, sep = ";")

nisstatt <- select(stationtt,Nesting,N)
nisstatt$Nesting[nisstatt$Nesting  %in% c("Above","Renter:Existing_cavities_above_ground","Renter:Existingcavitiesaboveground")]<-"Au-dessus du sol"
nisstatt$Nesting[nisstatt$Nesting  %in% c("Excavator:Ground","Underground")]<-"Dans le sol"
nisstatt$Nesting[nisstatt$Nesting  %in% c("Wood/Existing cavities","Plant_stems/Existing_cavities","Plantstems/Existingcavities","Existingcavities")]<-"Mix"
nisstatt$Nesting[nisstatt$Nesting  %in% c("parasite","Social_parasite","Socialparasite","Social parasite")]<-"Parasite"
nisstatt$Nesting[nisstatt$Nesting  %in% c("")]<-"Non connu"
nisstatt <- aggregate(N ~ Nesting, data =nisstatt, sum)

nisstatt <- nisstatt %>% 
  arrange(desc(Nesting)) %>%
  mutate(prop = N/ sum(nisstatt$N) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )

nisstatt <- nisstatt %>% mutate(prop = N/ sum(N) *100) 
nisstatt <-nisstatt %>% mutate_if(is.numeric, ~round(., 1))

ggplot(nisstatt,aes(x = "", y = prop, fill = Nesting))  +
  geom_bar(width = 1,stat ="identity",color="white" )+ 
  coord_polar("y", start = 0)+
  xlab("Classe")+
  ylab("Effectifs")+
  theme_void() +
  geom_text(aes(y = ypos, label = prop), color = "black", size=5)+
  scale_fill_manual("Lieu de nidification",values = cbp2)+
  labs(
    title = "Lieu de nidification des abeilles 
présentes sur tout les sites",
    subtitle = "Pourcentage d'individus selon leur lieu de nidification",
    caption = "Calcul réalisé sur les données de 2020, collectée sur 27 sites de la région de Mons, 
    pour un total de 2909 individus.Graphique réalisé grâce à RStudio,Version 1.3.1093") +
  theme(
    plot.title = element_text(hjust = 0.5,size = 20),
    plot.subtitle = element_text(hjust = 0.5, size = 15),
    plot.caption = element_text(face = "italic", size = 11,hjust = 0.5),
    legend.title = element_text(size=18),
    legend.text = element_text(size=15))

#### sp
stationtt <- read.csv ( "beeplanttrait.csv",header = T, sep = ";")

nidsta <- select(stationtt,Nesting,SPECTAXPRIO)
nidsta$Nesting[nidsta$Nesting  %in% c("Above","Renter:Existing_cavities_above_ground","Renter:Existingcavitiesaboveground")]<-"Au-dessus du sol"
nidsta$Nesting[nidsta$Nesting  %in% c("Excavator:Ground","Underground")]<-"Dans le sol"
nidsta$Nesting[nidsta$Nesting  %in% c("Wood/Existing cavities","Plant_stems/Existing_cavities","Plantstems/Existingcavities","Existingcavities")]<-"Mix"
nidsta$Nesting[nidsta$Nesting  %in% c("parasite","Social_parasite","Socialparasite","Social parasite")]<-"Parasite"
nidsta$Nesting[nidsta$Nesting  %in% c("")]<-"Non connu"

nidsta <- distinct(nidsta)
nidsta <- count(nidsta,Nesting)

nidsta <- nidsta %>% 
  arrange(desc(Nesting)) %>%
  mutate(prop = n/ sum(nidsta$n) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )

nidsta <- nidsta %>% mutate(prop = n/ sum(n) *100) 
nidsta <-nidsta %>% mutate_if(is.numeric, ~round(., 1))

ggplot(nidsta,aes(x = "", y = prop, fill = Nesting))  +
  geom_bar(width = 1,stat ="identity",color="white" )+ 
  coord_polar("y", start = 0)+
  xlab("Classe")+
  ylab("Effectifs")+
  theme_void() +
  geom_text(aes(y = ypos, label = prop), color = "black", size=5)+
  scale_fill_manual("Socialité",values = cbp2)+
  labs(
    title = "Lieu de nidification des espèces d'abeilles 
présentes sur tout les sites",
    subtitle = "Pourcentage d'espèces selon leur lieu de nidification",
    caption = "Calcul réalisé sur les données de 2020, collectée sur 27 sites de la région de Mons, 
    pour un total de 116 espèces.Graphique réalisé grâce à RStudio,Version 1.3.1093") +
  theme(
    plot.title = element_text(hjust = 0.5,size = 20),
    plot.subtitle = element_text(hjust = 0.5, size = 15),
    plot.caption = element_text(face = "italic", size = 11,hjust = 0.5),
    legend.title = element_text(size=18),
    legend.text = element_text(size=15))

#statut, ( mettre les proportion mais chiffre trop long, comment les raccourcir ? )
statut<- select(beetotal,SPECTAXPRIO,Statut_IUCN)
statut <-statut %>% filter(SPECTAXPRIO != 'Andrena  sp.')
statut <-statut %>% filter(SPECTAXPRIO != 'Nomada  sp.')
statut <- distinct(statut)
statut$Statut_IUCN[statut$Statut_IUCN  %in% c("lc")]<-"LC"
statut <- count(statut,Statut_IUCN)
statut <- arrange(statut,n)

statut <- statut %>% 
  arrange(desc(Statut_IUCN)) %>%
  mutate(prop = n/ sum(statut$n) *100) %>%
  mutate(ypos = cumsum(prop)-0.5*prop )

statut <- statut %>% mutate(prop = n/ sum(n) *100) 
statut <-statut %>% mutate_if(is.numeric, ~round(., 1))

statut$Statut_IUCN[statut$Statut_IUCN  %in% c("LC")]<-"Least Concern"
statut$Statut_IUCN[statut$Statut_IUCN  %in% c("VU")]<-"Vulnerable"
statut$Statut_IUCN[statut$Statut_IUCN  %in% c("NT")]<-"Near Threatened"
statut$Statut_IUCN[statut$Statut_IUCN  %in% c("DD")]<-"Data Deficient"
statut$Statut_IUCN[statut$Statut_IUCN  %in% c("CR")]<-"Critically Endangered,"

ggplot(statut,aes(x = "", y = prop, fill = Statut_IUCN))  +
  geom_bar(width = 1,stat ="identity",color="white" )+ 
  coord_polar("y", start = 0)+
  xlab("Classe")+
  ylab("Effectifs")+
  theme_void() +
  scale_fill_manual("Statut IUCN",values = cbp1)+
  geom_text(aes(y = ypos, label = prop), color = "black", size=3)+
  labs(
    title = "Nombre d'individu par catégorie de statut IUCN",
    subtitle = "",
    caption = "Calcul réalisé sur les données des années 2019,2020 et 2021, pour un total de 57 espèces +
le groupe des Terrestrisbombus (Bombus terrestris, lucorum, magnus).
  Graphique réalisé grâce à RStudio,Version 1.3.1093") +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(face = "bold"),
    plot.caption = element_text(face = "italic"))

ggplot(statut,aes(reorder(Statut_IUCN,n), y = n, fill = Statut_IUCN))  +
  geom_bar(stat ="identity",color="white" )+ 
  xlab("Statut IUCN")+
  ylab("Nombre d'espèces")+
  coord_flip()+
  scale_fill_manual("Statut IUCN",values = cbp1)+
  theme(legend.position = "none")+
  geom_text(aes(y = n+1, label = n), color = "black", size=3)+
  scale_y_continuous(breaks = round(seq(min(0), max(240), by = 5),1))+
  labs(
    title = "Nombre d'espèces selon leur statut IUCN",
    subtitle = "",
    caption = "Calcul réalisé sur les données des années 2019,2020 et 2021, pour un total de 57 espèces +
le groupe des Terrestrisbombus (Bombus terrestris, lucorum, magnus).
  Graphique réalisé grâce à RStudio,Version 1.3.1093") +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(face = "bold"),
    plot.caption = element_text(face = "italic"))
###combiner les graph sur la nidifications
combine_charts(list(a, b, c, d))

# statut en danger par famille. fonctionne pas pour l'instant

statfam <- select (trait,RL_CAT,Family)
statfam <- filter ( statfam, RL_CAT != "DD")
statfam <- filter ( statfam, RL_CAT != "LC")
statfam <- filter ( statfam, RL_CAT != "/")

ggplot(statfam,aes(RL_CAT,Family, fill = Family))   +
  geom_col(color="white" )+ 
  xlab("Espèce") +
  ylab("Effectifs")  +
  scale_y_continuous(breaks = round(seq(min(0), max(200), by = 4),1))+
  coord_flip()+
  scale_fill_manual("Famille d'abeille",values = cbp2)

geom_text(aes(label = nombre ,y = nombre + 1), size =3)+
  labs(
    title = "Nombre d'individu par espèces",
    subtitle = "A simple bar chart",
    caption = "Source: ImaginaryCo") +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(face = "bold"),
    plot.caption = element_text(face = "italic"))

### camembert famille 
social <- select(trait,Sociality)
social <- count(social,Sociality)
social<- arrange(social,n)

social <- social %>% 
  arrange(desc(Sociality)) %>%
  mutate(prop = n/ sum(social$n) *100) %>%
  mutate(ypos = cumsum(prop)-0.5*prop )

ggplot(social,aes(x = "", y = prop, fill = Sociality))  +
  geom_bar(width = 1,stat ="identity",color="white" )+ 
  coord_polar("y", start = 0)+
  xlab("Classe")+
  ylab("Effectifs")+
  theme_void() +
  scale_fill_manual(values = cbp1)+
  geom_text(aes(y = ypos, label = n), color = "black", size=3)+
  labs(
    title = "Nombre d'individu par espèces",
    subtitle = "A simple bar chart",
    caption = "Source: ImaginaryCo") +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(face = "bold"),
    plot.caption = element_text(face = "italic"))

Réservetot <- aggregate(N ~ Nesting, data =réserve, sum)
