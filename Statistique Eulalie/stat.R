setwd("C:/Users/lalie/Desktop/stat_prés_R_directory")

install.packages("SciViews")
SciViews::R 

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
prairie <- filter(stationtt, type == "Prairie")

####################### Rank-abundance curves

beett <- read.csv ( "combinebeeplante.csv",header = T, sep = ";")
beettsp <- select (beetotal,SPECTAXPRIO,N )

beettspgr <- aggregate(x = beettsp$N, by = list(beettsp$SPECTAXPRIO),FUN = sum) 
beett1 <- as.data.frame(t(beettspgr))
names(beett1) <- beett1[1,]
beett1 <- beett1[-1,]

beett1 <-type.convert(beett1)

RankAbun.1 <- rankabundance(beett1)
RankAbun.1
rankabunplot(RankAbun.1, scale='abundance', addit=FALSE, specnames=c(1,2,3), main="Main title",
             xlab="",
             ylab="Abondance",
             sub="Sub-title")
#scale='proportion'

######################### Rarefaction

rownamesbee<-row.names(beestattab)

test <- iNEXT(t(beestattab), q = c(0,1,2), datatype ="abundance") 
test$DataInfo
test$iNextEst
test$AsyEst
hill <- estimateD(t(beestattab), base="coverage")
hillmatrix0= hill[hill[,"order"]==0,]
hillmatrix1= hill[hill[,"order"]==1,]
hillmatrix2= hill[hill[,"order"]==2,]
hillmatrix = hillmatrix0[,c(1,2,4)]
hillmatrix$SC0 = test$DataInfo$SC
hillmatrix[,c("H0r","H0rU","H0rL")]=hillmatrix0[,c("qD","qD.UCL","qD.LCL")]
hillmatrix[,c("H1r","H1rU","H1rL")]=hillmatrix1[,c("qD","qD.UCL","qD.LCL")]
hillmatrix[,c("H2r","H2rU","H2rL")]=hillmatrix2[,c("qD","qD.UCL","qD.LCL")]
rownames(hillmatrix)= rownamesbee 

######################### histogramme multiple
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

beett <- read.csv ( "combinebeeplante.csv",header = T, sep = ";")
beett$CONDGR2[beett$CONDGR2 %in% c("LEGUMINOSAE")]<-"Fabaceae"
beett$CONDGR2[beett$CONDGR2  %in% c("COMPOSITAE")]<-"Asteraceae"
beett$CONDGR2[beett$CONDGR2  %in% c("UMBELLIFERAE")]<-"Apiaceae"
beett$CONDGR2[beett$CONDGR2  %in% c("BORAGINACEAE")]<-"Boraginaceae"
beett$CONDGR2[beett$CONDGR2  %in% c("LABIATAE")]<-"Lamiaceae"
beett$CONDGR2[beett$CONDGR2  %in% c("ROSACEAE")]<-"Rosaceae"
beett$CONDGR2[beett$CONDGR2  %in% c("CONVOLVULACEAE")]<-"Convolvulaceae"
beett$CONDGR2[beett$CONDGR2  %in% c("")]<-"Vol/Sol"

# abeille par famille et genre,( mettre le nombre total au bout ? )

beefam <- select (prairie,SPECGR2,N,SPECGEN )
beefam <- aggregate(N ~ SPECGR2 + SPECGEN , data =beefam, sum)
beefam <- beefam %>% group_by(SPECGR2) %>% mutate(sum(N))
beefam <- beefam  %>% rename(c("sum(N)" = "nombre"))

beegen <- select (prairie,SPECGR2,N,SPECGEN, SPECGR2 )
beegen <- aggregate(N ~ SPECGEN+SPECGR2 , data =beegen, sum)
beegen <- beegen %>% group_by(SPECGEN) %>% mutate(sum(N))
beegen <- beegen  %>% rename(c("sum(N)" = "nombre"))

ggplot(beefam,aes(reorder(SPECGR2, nombre),N, fill = SPECGR2))  +
  geom_col()+ 
  xlab("Famille") +
  ylab("Nombre")  +
  theme(legend.position = "none")+
  coord_flip()+
  geom_text(aes(label = nombre,y = nombre + 6), size =4)+
  scale_fill_viridis_d(name = 'Genre')+
  scale_y_continuous(breaks = round(seq(min(0), max(600), by = 20),1))+
  labs(
    title = "Nombre d'individus observ�s par famille",
    subtitle = "",
    caption = "") +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(face = ),
    plot.caption = element_text(face = "italic") )

ggplot(beegen,aes(reorder(SPECGEN, nombre),N, fill = SPECGR2))  +
  geom_col()+ 
  xlab("Genre") +
  ylab("Nombre")  +
  coord_flip()+
  geom_text(aes(label = nombre,y = nombre + 6), size =4)+
  scale_fill_viridis_d(name = 'Famille')+
  scale_y_continuous(breaks = round(seq(min(0), max(600), by = 20),1))+
  labs(
    title = "Nombre d'individus observ�s par genre",
    subtitle = "",
    caption = "") +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(face = ),
    plot.caption = element_text(face = "italic") )


#famille de plante sur lesquelle y avait bee

plantefam <- select (prairie,SPECGR2,N,CONDGR2 )
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
plantefam$CONDGR2[plantefam$CONDGR2  %in% c("SALICACEAE")]<-"Salicaceae"

plantefam <- plantefam %>% filter(CONDGR2 !='Vol/Sol')
  

plantefam <- aggregate(N ~ CONDGR2 +SPECGR2 , data =plantefam, sum)
plantefam <- plantefam %>% group_by(CONDGR2) %>% mutate(sum(N))
plantefam <- plantefam  %>% rename(c("sum(N)" = "nombre"))

ggplot(plantefam,aes(reorder(CONDGR2,nombre),N, fill = SPECGR2)) +
  geom_col(color="white")+ 
  xlab("Famille") +
  ylab("Nombre")  +
  geom_text(aes(label = nombre,y = nombre + 2), size =3)+
  coord_flip()+
  scale_fill_manual('Famille',values = cbp2)+
  scale_y_continuous(breaks = round(seq(min(0), max(200), by = 8),1))+
  labs(title = "Nombre d'individus observés sur les différentes familles de plantes ",
       subtitle = "Division des abeilles par famille",
       caption = "Calcul réalisé sur les données des années 2019, 2020 et 2021, pour un total de 520 individus.
       Les abeilles qui ont été observées au vol ou au sol ont été exclues de ce graphique.
       Graphique réalisé grâce à RStudio,Version 1.3.1093") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(face = ),
        plot.caption = element_text(face = "italic"))


#nombre d'abeille par especes

beettsp <- select (prairie,SPECTAXPRIO,N )

beettspgr <- aggregate(x = beettsp$N, by = list(beettsp$SPECTAXPRIO),FUN = sum) 
beettspgr1 <- filter(beettspgr,x > 10)

ggplot(beettspgr1,aes(reorder(Group.1,x),x))   +
  geom_col( fill = "#56B4E9")+ 
  xlab("Esp�ce") +
  ylab("Nombre")  +
  theme(legend.position = "none")+
  scale_y_continuous(breaks = round(seq(min(0), max(400), by = 10),1))+
  coord_flip()+
  geom_text(aes(label = x,y = x + 3), size =4)+
labs(
  title = "Nombre d'individus observ�s par esp�ce",
  subtitle = "",
  caption = "Seules les esp�ces pr�sentant plus de 10 individus ont �t� repr�sent�es") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18),
    plot.subtitle = element_text(face =),
    plot.caption = element_text(face = "italic", size =12),
    axis.title.x =element_text(size = 15), 
    axis.title.y =element_text(size = 15),
    axis.text.x=element_text(size = 11),
    axis.text.y=element_text(size = 11))

### les 5 especes avec le plus d'individu et leur sexe

beemaxsex <- select (prairie,SPECTAXPRIO,N,SEX )
beemaxsex$SEX[beemaxsex$SEX  %in% c("f")]<-"F"
beemaxsex <- aggregate(N ~ SPECTAXPRIO + SEX, data =beemaxsex, sum)
beemaxsex <- beemaxsex %>% group_by(SPECTAXPRIO) %>% mutate(sum(N))
beemaxsex<- beemaxsex  %>% rename (c("sum(N)" = "nombre"))
beemaxsex <- filter (beemaxsex, nombre >17) 
beemaxsex$SEX[beemaxsex$SEX  %in% c("f")]<-"F"

ggplot(beemaxsex,aes(reorder(SPECTAXPRIO, nombre),N, fill = SEX))  +
  geom_col(color="white")+ 
  xlab("Famille") +
  ylab("Nombre")  +
  coord_flip()+
  geom_text(aes(label = nombre,y = nombre + 1.5), size =3)+
  scale_fill_manual('Sexe',labels = c("Femelle", "Mâle", "Ouvrière"),values = cbp2)+
  scale_y_continuous(breaks = round(seq(min(0), max(200), by = 5),1))+
  labs(
    title = "Focus sur les 5 espèces les plus abondantes",
    subtitle = "Division des espèces par sexe",
    caption = "Calcul réalisé sur les données des années 2019, 2020 et 2021, pour un total de 520 individus.
    Réalisé grâce à RStudio,Version 1.3.1093") +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(face = ),
    plot.caption = element_text(face = "italic"))

##plante native avec bee

ggplot(prairie,aes(Indig�ne,N, fill = Mellif�re))  +
  geom_col()+ 
  xlab("Plante Indigène") +
  ylab("Nombre")  +
  scale_fill_manual('Plante mellifère',values = cbp2)+
  scale_y_continuous(breaks = round(seq(min(0), max(400), by = 10),1))+
  coord_flip()+
labs(
  title = "Nombre d'individu capturé par famille",
  subtitle = "A simple bar chart",
  caption = "Source: ImaginaryCo") +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(face = "bold"),
    plot.caption = element_text(face = "italic"))

##plante mellifere avec bee

ggplot(beett,aes(Mellifère,N, fill = Indigène))  +
  geom_col()+ 
  xlab("Plante mellifère") +
  ylab("Nombre")  +
  scale_fill_manual('Plante indigène',values = cbp2)+
  coord_flip()+
  scale_y_continuous(breaks = round(seq(min(0), max(250), by = 10),1))+
  labs(
    title = "Nombre d'individu capturé par famille",
    subtitle = "A simple bar chart",
    caption = "Calcul réalisé sur les données des années 2019, 2020 et 2021, pour un total de 520 individus.
    Réalisé grâce à RStudio,Version 1.3.1093") +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(face = "bold"),
    plot.caption = element_text(face = "italic"))

## 5sp de plante avec le plus de bee

beeplant <- select (prairie,CONDTAXPRIO,SPECGR2,N )
beeplant <- beeplant %>% filter(CONDTAXPRIO != '')
beeplant <- aggregate(N ~ CONDTAXPRIO +SPECGR2 , data =beeplant, sum)
beeplant <- beeplant %>% group_by(CONDTAXPRIO) %>% mutate(sum(N))
beeplant <- beeplant  %>% rename(c("sum(N)" = "nombre"))
beeplant <- filter(beeplant,nombre > 30)

ggplot(beeplant,aes(reorder(CONDTAXPRIO,nombre),N))+
  geom_col(fill = '#56B4E9')+ 
  xlab("Esp�ce") +
  ylab("Nombre")  +
  scale_y_continuous(breaks = round(seq(min(0), max(200), by = 4),1))+
  coord_flip()+
  scale_fill_manual("Famille d'abeille",values = cbp2)+
  geom_text(aes(label = nombre ,y = nombre + 1), size =3)+
  labs(
    title = "Les 5 esp�ces de plantes o� le plus d'abeilles ont �t� observ�es",
    subtitle = "Une division par familles d'abeilles a �t� effectu�e",
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

## 5sp de plante avec le plus sp

beepla <- select (prairie,CONDTAXPRIO,SPECGR2,SPECTAXPRIO )
beepla <- beepla %>% filter(CONDTAXPRIO != '')
beepla <- distinct(beepla)
beepla <- beepla %>% group_by(CONDTAXPRIO)%>% count(SPECTAXPRIO)
beepla <- aggregate(n ~ CONDTAXPRIO +SPECTAXPRIO , data =beepla, sum)
beepla <- beepla %>% group_by(CONDTAXPRIO) %>% mutate(sum(n))
beepla <- beepla  %>% rename(c("sum(n)" = "nombre"))
beepla <- filter(beepla,nombre > 10)

ggplot(beepla,aes(reorder(CONDTAXPRIO,nombre),n))+
  geom_col(fill = '#56B4E9')+ 
  xlab("Esp�ce") +
  ylab("Nombre d'esp�ces")  +
  scale_y_continuous(breaks = round(seq(min(0), max(200), by = 2),1))+
  coord_flip()+
  scale_fill_manual("Famille d'abeille",values = cbp2)+
  geom_text(aes(label = nombre ,y = nombre+0.2 ), size =4)+
  labs(
    title = "Les 5 esp�ces de plantes o� le plus d'abeilles ont �t� observ�es",
    subtitle = "",
    caption = "") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18),
    plot.subtitle = element_text(face = ),
    plot.caption = element_text(face = "italic"),
    axis.title.x =element_text(size = 15), 
    axis.title.y =element_text(size = 15),
    axis.text.x=element_text(size = 11),
    axis.text.y=element_text(size = 11))


##############################################################analyse paysage 

plante <-read.csv("listeplantepréstout.csv",header = T, sep = ";")
plante$Famille[plante$Famille %in% c("FABACEAE")]<-"Fabaceae"
plante$Famille[plante$Famille  %in% c("ASTERACEAE")]<-"Asteraceae"
plante$Famille[plante$Famille  %in% c("BETULACEAE")]<-"Betulaceae"
plante$Famille[plante$Famille  %in% c("BORAGINACEAE")]<-"Boraginaceae"
plante$Famille[plante$Famille  %in% c("LAMIACEAE")]<-"Lamiaceae"
plante$Famille[plante$Famille  %in% c("BRASSICACEAE")]<-"Brassicaceae"
plante$Famille[plante$Famille  %in% c("CONVOLVULACEAE")]<-"Convolvulaceae"
plante$Famille[plante$Famille  %in% c("SCROPHULARIACEAE")]<-"Scrophulariaceae"
plante$Famille[plante$Famille  %in% c("RANUNCULACEAE")]<-"Ranunculaceae"
plante$Famille[plante$Famille  %in% c("CYPERACEAE")]<-"Cyperaceae"
plante$Famille[plante$Famille  %in% c("CHENOPODIACEAE")]<-"Chenopodiaceae"
plante$Famille[plante$Famille  %in% c("APIACEAE")]<-"Apiaceae"
plante$Famille[plante$Famille  %in% c("EQUISETACEAE")]<-"Equisetaceae"
plante$Famille[plante$Famille  %in% c("HYPERICACEAE")]<-"Hypericaceae"
plante$Famille[plante$Famille  %in% c("JUNCACEAE")]<-"Juncaceae"
plante$Famille[plante$Famille  %in% c("ONAGRACEAE")]<-"Onagraceae"


##plante native 
plantenbr <- count(plante, Indigène)
ggplot(plantenbr,aes(Indigène,n, fill = Indigène))  +
  geom_col()+ 
  xlab("Plante indigène") +
  ylab("Nombre")  +
  coord_flip()+
  theme(legend.position = "none")+
  scale_y_continuous(breaks = round(seq(min(0), max(200), by = 6),1))+
  geom_text(aes(y = n +2, label = n), color = "black", size=3)+
  scale_fill_manual(values = cbp2)+
  labs(
    title = "Nombre d'espèces de plantes indigènes ou non",
    subtitle = "",
    caption = "Le caractère indigène des espèces a été évalué grâce à la nouvelle flore de la Belgique,
    du G.-D. de Luxembourg, du nord de la France et des régions voisines, sixième édition. 
    Réalisé grâce à RStudio,Version 1.3.1093 ") +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(face = "bold"),
    plot.caption = element_text(face = "italic"))

#proportion

plantenbrprop <- plantenbr %>% mutate(prop = n/ sum(n) *100) 
plantenbrprop <-plantenbrprop %>% mutate_if(is.numeric, ~round(., 1))
ggplot(plantenbrprop,aes(Indigène,prop, fill = Indigène))  +
  geom_col()+ 
  xlab("Plante indigène") +
  ylab("Pourcentage")  +
  coord_flip()+
  theme(legend.position = "none")+
  scale_y_continuous(breaks = round(seq(min(0), max(200), by = 6),1))+
  geom_text(aes(y = prop +2 , label = prop), color = "black", size=3)+
  scale_fill_manual(values = cbp2)+
  labs(
    title = "Proportions d'espèces de plante indigène ou non",
    subtitle = "",
    caption = "Le caractère indigène des espèces a été évalué grâce à la nouvelle flore de la Belgique,
    du G.-D. de Luxembourg, du nord de la France et des régions voisines, sixième édition. 
    Le nombre total d'especes de plantes est de 160, liste combinant les données d'Observations.be,
    OFFH, et les relevés de 2019 et 2020. 
    Réalisé grâce à RStudio,Version 1.3.1093 ") +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(face = "bold"),
    plot.caption = element_text(face = "italic", size =8))

###pour autre stations 
platout <-read.csv("plantetoutstation.csv",header = T, sep = ";")
planttouti <- select (platout,nom,Indigène,nom_latin )
planttouti <- planttouti %>% group_by(nom) %>% group_by(nom_latin)
planttoutinbr <- planttouti %>% group_by(nom) %>% group_by(nom_latin) %>% count(nom, Indigène)
planttoutinbr<- planttoutinbr%>% group_by(nom) %>% count(nom, Indigène)

planttoutinbr <- planttoutinbr %>% group_by(nom) %>% mutate(prop = n/ sum(n) *100) 
planttoutinbr <-planttoutinbr %>% mutate_if(is.numeric, ~round(., 1))

planttoutFA <- planttoutinbr %>% group_by(nom) %>% filter(Indigène == "FALSE")
planttoutTR <- planttoutinbr %>% group_by(nom) %>% filter(Indigène == "TRUE")%>% arrange (prop)
planttoutindi <- full_join(planttoutFA ,planttoutTR, by = "nom")
planttoutindi[is.na(planttoutindi)] <- 0

ggplot(planttoutindi,aes(reorder(nom,prop.y),prop.y, fill = Indigène.y ))+
  geom_col(position = "stack")+ 
  xlab("Station") +
  ylab("Pourcentage")  +
  coord_flip()+
  theme(legend.position = "none")+
  scale_y_continuous(breaks = round(seq(min(0), max(200), by = 5),1))+
  geom_text(aes(y = prop.y+2, label = prop.y), color = "black", size=3)+
  scale_fill_manual( values = c("#009E73"))+
  labs(
    title = "Proportion d'espèces de plantes indigènes par station",
    subtitle = "Calculé en fonction des relevés de 2020",
    caption = "Le caractère indigène des espèces a été évalué grâce à la nouvelle flore de la Belgique,du G.-D. de Luxembourg, du nord de la France et des régions voisines,
    sixième édition.Seules les plantes en fleurs pendant la période de mai-août ont été comptabilisées. Réalisé grâce à RStudio,Version 1.3.1093 ") +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(face = "italic", size = 9),
    plot.caption = element_text(face = "italic", size = 8))

###plante invasive prés

invasif <- filter(plante, Indigène == "FALSE")
invasif <- count(invasif, Famille)
invasif <- filter(invasif, n >1)

ggplot(invasif,aes(reorder(Famille,n),n, fill = Famille))  +
  geom_col(position = "dodge")+ 
  xlab("Plante non-indigène") +
  ylab("Nombre")  +
  theme(legend.position = "none")+
    scale_fill_manual(values = cbp2)+
  coord_flip()+
  scale_y_continuous(breaks = round(seq(min(0), max(200), by = 1),1))+
  labs(
    title = "Nombre d'espèces de plante non-indigène par famille",
    subtitle = "",
    caption = "Le caractère indigène des espèces a été évalué grâce à la nouvelle flore de la Belgique,
    du G.-D. de Luxembourg, du nord de la France et des régions voisines, sixième édition. 
    Les familles ne comptant qu'une seule espèce non-indigène ont été omise de ce graphique.
    Réalisé grâce à RStudio,Version 1.3.1093 ") +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(face = "bold"),
    plot.caption = element_text(face = "italic"))
  
 
###plante mellifere

plantenbrmel <- count(plante, Mellifère)


ggplot(plantenbrmel,aes(Mellifère,n, fill = Mellifère))  +
  geom_col()+ 
  xlab("Plante mellifère") +
  ylab("Pourcentage")  +
  coord_flip()+
  theme(legend.position = "none")+
  scale_y_continuous(breaks = round(seq(min(0), max(200), by = 4),1))+
  geom_text(aes(y = n +2, label = n), color = "black", size=3)+
  scale_fill_manual(values = cbp2)+
  labs(
    title = "Proportions d'espèces de plantes mellifères ou non",
    subtitle = "",
    caption = "Le caractère mellifère des espèces a été évalué grâce à la nouvelle flore de la Belgique,du G.-D. de Luxembourg, du nord de la France 
    et des régions voisines, sixième édition. Réalisé grâce à RStudio,Version 1.3.1093 ") +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(face = "bold"),
    plot.caption = element_text(face = "italic"))

###pour toute station 

platout <-read.csv("plantetoutstation.csv",header = T, sep = ";")
planttoutm <- select (platout,nom,Mellifère,nom_latin )
planttoutm <- planttoutm %>% group_by(nom) %>% group_by(nom_latin)
planttoutmelnbr <- planttoutm %>% group_by(nom) %>% group_by(nom_latin) %>% count(nom, Mellifère)
planttoutmelnbr<- planttoutmelnbr%>% group_by(nom) %>% count(nom, Mellifère)

planttoutmelnbr <- planttoutmelnbr %>% group_by(nom) %>% mutate(prop = n/ sum(n) *100) 
planttoutmelnbr <-planttoutmelnbr %>% mutate_if(is.numeric, ~round(., 1))

planttoutmel <- planttoutmelnbr %>% group_by(nom) %>% filter(Mellifère == "TRUE")%>% arrange (prop)
planttoutmel[is.na(planttoutmel)] <- 0

ggplot(planttoutmel,aes(reorder(nom,prop),prop, fill = Mellifère ))+
  geom_col(position = "stack")+ 
  xlab("Station") +
  ylab("Pourcentage")  +
  coord_flip()+
  theme(legend.position = "none")+
  scale_y_continuous(breaks = round(seq(min(0), max(200), by = 5),1))+
  geom_text(aes(y = prop+2, label = prop), color = "black", size=3)+
  scale_fill_manual( values = c("#009E73"))+
  labs(
    title = "Proportion d'espèces de plantes mellifères par station",
    subtitle = "Calculé en fonction des relevés de 2020",
    caption = "Le caractère mellifère des espèces a été évalué grâce à la nouvelle flore de la Belgique,du G.-D. de Luxembourg,du nord de la France 
    et des régions voisines, sixième édition. Seules les plantes en fleurs pendant la période de mai-août ont été comptabilisées. 
    Réalisé grâce à RStudio,Version 1.3.1093 ") +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(face = "italic", size = 9),
    plot.caption = element_text(face = "italic", size = 8))

####native et mellifère

ggplot(plante,aes(Indigène, fill = Mellifère))  +
  geom_bar()+ 
  xlab("Plante native") +
  ylab("Nombres")  +
  coord_flip()+
  theme(legend.position = "none")+
  scale_y_continuous(breaks = round(seq(min(0), max(200), by = 6),1))+
  scale_fill_brewer(palette="Accent")+
  labs(color = "Plante mellifère")





