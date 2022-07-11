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
library("dplyr")

cbp2 <- c( "#56B4E9", "#009E73",
           "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#E69F00","#8B6914","#999999",
           '#8B1C62','#32CD32',"#0000EE","#FF407A")

###################################################################################DATA
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


cbp2 <- c( "#56B4E9", "#009E73",
           "#0072B2", "#D55E00", "#CC79A7", "#E69F00","#999999",
           '#8B1C62','#32CD32',"#0000EE","#FF407A","#8B6914","#F0E442")



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

listeespece <- select(stationtt,N)
listeespece<- listeespece %>% mutate(sum(N))
listeespece <- distinct(listeespece)

###matrice

beestattab <- select (stationtt,SPECTAXPRIO,N,TOPO)

beestattab <- aggregate(N~SPECTAXPRIO+TOPO, data = beestattab, sum)
beestattab <- xtabs(N~TOPO+SPECTAXPRIO,beestattab)
beestattab <-type.convert(beestattab)
#####################################################matrice abeille par milieu
milieu <- select (stationtt,SPECTAXPRIO,N,type,TOPO)
milieu <-as.data.frame(milieu)

milieutab <- aggregate(N~type+SPECTAXPRIO, data = milieu, sum)
milieutab <- xtabs(N~type+SPECTAXPRIO,milieutab)
milieutab <-type.convert(milieutab)

######################################################## abeille par famille et genre
  
beefam <- select (stationtt,SPECGR2,N,SPECGEN, type )
beefam <- aggregate(N ~ SPECGEN , data =beefam, sum)
beefam <- beefam %>%group_by(SPECGEN)%>% mutate(sum(N))
beefam <- rename(beefam,c('sum(N)'='nombre'))

ggplot(beefam,aes(reorder(SPECGEN,nombre),N ))  +
  geom_col(fill = c("#56B4E9") )+ 
  xlab("Genre") +
  ylab("Nombre d'individus")  +
  coord_flip()+
  geom_text(aes(label = nombre,y = nombre + 50), size =4)+
  scale_fill_manual('Type de milieu',values = cbp2)+
  scale_y_continuous(breaks = round(seq(min(0), max(3500), by = 150),1))+
  labs(
    title = "Nombre d'individus observés par genre pour toutes les stations",
    subtitle = "",
    caption = "") +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(face = ),
    plot.caption = element_text(face = "italic") )

################################Nombre d'individus observés sur les différentes familles de plantes
plantefam <- select (stationtt,milieu,N,CONDGR2 )
plantefam$CONDGR2[plantefam$CONDGR2 %in% c("LEGUMINOSAE")]<-"Fabaceae"
plantefam$CONDGR2[plantefam$CONDGR2  %in% c("COMPOSITAE")]<-"Asteraceae"
plantefam$CONDGR2[plantefam$CONDGR2  %in% c("BETULACEAE")]<-"Betulaceae"
plantefam$CONDGR2[plantefam$CONDGR2  %in% c("BORAGINACEAE")]<-"Boraginaceae"
plantefam$CONDGR2[plantefam$CONDGR2  %in% c("LABIATAE")]<-"Lamiaceae"
plantefam$CONDGR2[plantefam$CONDGR2  %in% c("BRASSICACEAE")]<-"Brassicaceae"
plantefam$CONDGR2[plantefam$CONDGR2  %in% c("CONVOLVULACEAE")]<-"Convolvulaceae"
plantefam$CONDGR2[plantefam$CONDGR2  %in% c("SCROPHULARIACEAE")]<-"Scrophulariaceae"
plantefam$CONDGR2[plantefam$CONDGR2  %in% c("RANUNCULACEAE")]<-"Ranunculaceae"
plantefam$CONDGR2[plantefam$CONDGR2  %in% c("CYPERACEAE")]<-"Cyperaceae"
plantefam$CONDGR2[plantefam$CONDGR2  %in% c("CHENOPODIACEAE")]<-"Chenopodiaceae"
plantefam$CONDGR2[plantefam$CONDGR2  %in% c("APIACEAE")]<-"Apiaceae"
plantefam$CONDGR2[plantefam$CONDGR2  %in% c("EQUISETACEAE")]<-"Equisetaceae"
plantefam$CONDGR2[plantefam$CONDGR2  %in% c("HYPERICACEAE")]<-"Hypericaceae"
plantefam$CONDGR2[plantefam$CONDGR2  %in% c("JUNCACEAE")]<-"Juncaceae"
plantefam$CONDGR2[plantefam$CONDGR2  %in% c("ROSACEAE")]<-"Rosaceae"
plantefam$CONDGR2[plantefam$CONDGR2  %in% c("UMBELLIFERAE")]<-"Apiaceae"
plantefam$CONDGR2[plantefam$CONDGR2  %in% c("BERBERIDACEAE")]<-"Berberidaceae"
plantefam$CONDGR2[plantefam$CONDGR2  %in% c("RESEDACEAE")]<-"Resedaceae"
plantefam$CONDGR2[plantefam$CONDGR2  %in% c("GUTTIFERAE")]<-"Clusiaceae"
plantefam$CONDGR2[plantefam$CONDGR2  %in% c("GERANIACEAE")]<-"Geraniacee"
plantefam$CONDGR2[plantefam$CONDGR2  %in% c("BUDDLEJACEAE")]<-"Scrophulariaceae"
plantefam$CONDGR2[plantefam$CONDGR2  %in% c("ERICACEAE")]<-"Ericaceae"

plantefam <- plantefam %>% filter(CONDGR2 !='')


plantefam <- aggregate(N ~ CONDGR2 +milieu , data =plantefam, sum)
plantefam <- plantefam %>% group_by(CONDGR2) %>% mutate(sum(N))
plantefam <- plantefam %>% rename(c("nombre" = "sum(N)"))

plantefam <- plantefam %>% filter(nombre > 50)

ggplot(plantefam,aes(reorder(CONDGR2,nombre),N, fill = milieu)) +
  geom_col(color="white")+ 
  xlab("Famille") +
  ylab("Nombre")  +
  geom_text(aes(label = nombre,y = nombre +20), size =3)+
  coord_flip()+
  scale_fill_manual('Famille',values = cbp2)+
  scale_y_continuous(breaks = round(seq(min(0), max(2000), by = 100),1))+
  labs(title = "Nombre d'individus observés sur les différentes familles de plantes ",
       subtitle = "",
       caption = "") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(face = ),
        plot.caption = element_text(face = "italic"))

################################Nombre d'especes observés sur les différentes familles de plantes

planteind <- select (plantefam ,CONDGR2, SPECTAXPRIO )

planteind <- planteind %>% filter(CONDGR2 !='')
planteind <- distinct(planteind)
planteind <- count(planteind,SPECTAXPRIO,CONDGR2)

planteind <- aggregate(n ~ CONDGR2, data =planteind, sum)
planteind <- planteind %>% group_by(CONDGR2) %>% mutate(sum(n))
planteind <- planteind %>% rename(c("nombre" = "sum(n)"))

planteind <- planteind %>% filter(nombre > 15)

ggplot(planteind,aes(reorder(CONDGR2,nombre),n)) +
  geom_col(fill ="#009E73" )+ 
  xlab("Famille") +
  ylab("Nombre d'espèces")  +
  geom_text(aes(label = nombre,y = nombre +2), size =4)+
  coord_flip()+
  scale_fill_manual('Type de milieu',values = cbp2)+
  scale_y_continuous(breaks = round(seq(min(0), max(100), by = 5),1))+
  labs(title = "Nombre d'espèces observées sur les différentes familles de plantes ",
       subtitle = "",
       caption = "") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(face = ),
        plot.caption = element_text(face = "italic"))

#################################################Nombre d'individus observés par espèce

beettsp <- select (stationtt,SPECTAXPRIO,N )


beettspgr <- aggregate(x = beettsp$N, by = list(beettsp$SPECTAXPRIO),FUN = sum) 
beettspgr1 <- filter(beettspgr,x > 50)

ggplot(beettspgr1,aes(reorder(Group.1,x),x))   +
  geom_col( fill = "#009E73")+ 
  xlab("Espèce") +
  ylab("Nombre d'individus")  +
  theme(legend.position = "none")+
  scale_y_continuous(breaks = round(seq(min(0), max(1800), by = 100),1))+
  coord_flip()+
  geom_text(aes(label = x,y = x + 25), size =4)+
  labs(
    title = "Nombre d'individus observés par espèce",
    subtitle = "",
    caption = "") +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(face =),
    plot.caption = element_text(face = "italic"))

############################################ les 5 especes avec le plus d'individu et leur sexe
stationtt$SPECTAXPRIO[stationtt$SPECTAXPRIO%in% c("Bombus pascuorum floralis")]<-"Bombus pascuorum"
stationtt$SPECTAXPRIO[stationtt$SPECTAXPRIO%in% c("Bombus terrestris terrestris")]<-"Bombus (Bombus)  sp."
stationtt$SPECTAXPRIO[stationtt$SPECTAXPRIO%in% c("Bombus terrestris")]<-"Bombus (Bombus)  sp."
stationtt$SPECTAXPRIO[stationtt$SPECTAXPRIO%in% c("Bombus lucorum")]<-"Bombus (Bombus)  sp."


beemaxsex <- select (stationtt,SPECTAXPRIO,N,SEX )
beemaxsex$SEX[beemaxsex$SEX  %in% c("f")]<-"F"
beemaxsex$SEX[beemaxsex$SEX  %in% c("4")]<-"F"
beemaxsex <- aggregate(N ~ SPECTAXPRIO + SEX, data =beemaxsex, sum)
beemaxsex <- beemaxsex %>% group_by(SPECTAXPRIO) %>% mutate(sum(N))
beemaxsex <- beemaxsex %>% rename(c("nombre" = "sum(N)"))
beemaxsex <- filter (beemaxsex, nombre >100) 

ggplot(beemaxsex,aes(reorder(SPECTAXPRIO, nombre),N, fill = SEX))  +
  geom_col(color="white")+ 
  xlab("Espèce") +
  ylab("Nombre d'individus")  +
  coord_flip()+
  geom_text(aes(label = nombre,y = nombre + 20), size =3)+
  scale_fill_manual('Sexe',labels = c("Femelle", "Mâle", "Ouvrière"),values = cbp2)+
  scale_y_continuous(breaks = round(seq(min(0), max(1800), by = 100),1))+
  labs(
    title = "Focus sur les espèces les plus abondantes",
    subtitle = "Division des espèces par sexe",
    caption = "") +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(face = ),
    plot.caption = element_text(face = "italic"))

######################################################################diversité beta 
library(betapart)
###matrice
beestattab <- select (stationtt,SPECTAXPRIO,N,TOPO)

beestattab <- aggregate(N~SPECTAXPRIO+TOPO, data = beestattab, sum)
beestattab <- xtabs(N~TOPO+SPECTAXPRIO,beestattab)
beestattab <-type.convert(beestattab)

beestattab[beestattab > 0] <- 1

bees.beta<-beta.pair(beestattab, index.family="sor")
bees.beta



vegdist(milieutab, method="bray", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = FALSE)
dis<- vegdist(beestattab, method="bray", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = FALSE)

############################### Division par mois 

juin<- station %>% filter(DAT2 > 20200601)
juin <- juin %>% filter(DAT2 < 20200701)
juin <- juin %>%select(N,SPECTAXPRIO,TOPO,SPECGR2,SPECGEN)
juin<- aggregate(N ~ TOPO + SPECTAXPRIO+SPECGR2, data =juin, sum)

halictidae <-juin %>% filter(SPECGR2 == 'HALICTIDAE')
colletidae<-juin %>% filter(SPECGR2 == 'COLLETIDAE')
mellitidae <- juin %>% filter(SPECGR2 == 'MELITTIDAE')
andrenidae <- juin %>% filter(SPECGR2 == 'ANDRENIDAE')
apidae <- juin %>% filter(SPECGR2 == 'APIDAE')
mega <- juin %>% filter(SPECGR2 =='MEGACHILIDAE')

ggplot(mega,aes(TOPO,N, fill = SPECTAXPRIO))  +
  geom_col(color="white")+ 
  xlab("Station") +
  ylab("Nombre")  +
  coord_flip()+
  scale_y_continuous(breaks = round(seq(min(0), max(2000), by = 2),1))+
  labs(
    title = "Nombre d'individu par especes pour chaque station pour juin",
    subtitle = "",
    caption = "") +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(face = ),
    plot.caption = element_text(face = "italic"),
    legend.text = element_text(size = 8))
  legend("topright",cex = 1.5)
  #########
  juillet<- station %>% filter(DAT2 > 20200701)
  juillet <- juillet %>% filter(DAT2 < 20200801)
  juillet <- juillet %>%select(N,SPECTAXPRIO,TOPO,SPECGR2,SPECGEN)
  juillet<- aggregate(N ~ TOPO + SPECTAXPRIO+SPECGR2, data =juillet, sum)
  
  halictidae <-juillet %>% filter(SPECGR2 == 'HALICTIDAE')
  colletidae<-juillet %>% filter(SPECGR2 == 'COLLETIDAE')
  mellitidae <- juillet %>% filter(SPECGR2 == 'MELITTIDAE')
  andrenidae <- juillet %>% filter(SPECGR2 == 'ANDRENIDAE')
  apidae <- juillet %>% filter(SPECGR2 == 'APIDAE')
  mega <- juillet %>% filter(SPECGR2 =='MEGACHILIDAE')
  
  ggplot(mega,aes(TOPO,N, fill = SPECTAXPRIO))  +
    geom_col(color="white")+ 
    xlab("Station") +
    ylab("Nombre")  +
    coord_flip()+
    scale_y_continuous(breaks = round(seq(min(0), max(2000), by = 2),1))+
    labs(
      title = "Nombre d'individu par especes pour chaque station pour juillet",
      subtitle = "",
      caption = "") +
    theme(
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(face = ),
      plot.caption = element_text(face = "italic"),
      legend.text = element_text(size = 8))
  legend("topright",cex = 1.5)