
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

beett <- read.csv ( "combinebeeplante.csv",header = T, sep = ";")
beettsp <- select (beett,SPECTAXPRIO,N,Date )
beett19 <- beettsp %>% filter(grepl('2019',Date))
beett20 <- beettsp %>% filter(grepl('2020',Date))
beett19 <- beett19 %>% select(-'Date')
beett20 <- beett20 %>% select(-'Date')

réserve <- read.csv ( "réservetraitplante.csv",header = T, sep = ";")
beett21 <- select (réserve,SPECTAXPRIO,N)

beett19 $SPECTAXPRIO[beett19 $SPECTAXPRIO  %in% c("Bombus pascuorum ")]<-"Bombus pascuorum"
beett19 $SPECTAXPRIO[beett19 $SPECTAXPRIO  %in% c("Bombus terrestris ")]<-"Bombus terrestris"
beett19 <- aggregate(N~SPECTAXPRIO, data = beett19, sum)


beett20 <- aggregate(N~SPECTAXPRIO, data = beett20, sum)


beett21 <- aggregate(N~SPECTAXPRIO, data = beett21, sum)
réserve1920 <- full_join(beett19,beett20,by = "SPECTAXPRIO")
réserve192021 <- full_join(réserve1920,beett21,by = "SPECTAXPRIO")
names(réserve192021) <- c("Espèce","2019","2020","2021")
réserve192021[is.na(réserve192021)] = 0

réserve192021 <- as.data.frame(t(réserve192021))

names(réserve192021) <- réserve192021[1,]
réserve192021 <- réserve192021[-1,]
réserve192021 <- type.convert(réserve192021)

?rankabundance
RankAbun.1 <- rankabundance(réserve192021)
RankAbun.1
rankabunplot(RankAbun.1, scale='abundance', addit=FALSE, specnames=c(1,2,3))

#########courbe accumulations 

?specaccum
sp2 <- specaccum(réserve192021, "random", permutations = 999)
plot(sp2, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue")

plot(sp2)+
labs(
  title = "Régime alimentaire des abeilles 
présentes sur tout les sites",
  subtitle = "Pourcentage d'individus selon leur régime alimentaire",
  caption = "Calcul réalisé sur les données de 2020, collectée sur 27 sites de la région de Mons, 
    pour un total de 2909 individus.Graphique réalisé grâce à RStudio,Version 1.3.1093") +
  theme(
    plot.title = element_text(hjust = 0.5,size = 20),
    plot.subtitle = element_text(hjust = 0.5, size = 15),
    plot.caption = element_text(face = "italic", size = 11,hjust = 0.5),
    legend.title = element_text(size=18),
    legend.text = element_text(size=15))
#boxplot(sp2, col="yellow", add=TRUE, pch="+")


########courbe rarefaction

S <- specnumber(réserve192021) # observed number of species
(raremax <- min(rowSums(réserve192021)))
Srare <- rarefy(réserve192021, raremax)
plot(S, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")
abline(0, 1)
rarecurve(réserve192021, step = 5, sample = raremax, col = "blue", cex = 0.6)+

####raréfactions hill
 
abond <- iNEXT(t(réserve192021), q = c(0,1,2), datatype ="abundance") 

DataInfotabee <-type.convert(abond$DataInfo)

iNextEsttab <-type.convert(abond$iNextEst)
AsyEsttab <-type.convert(abond$AsyEst)

hill <- estimateD(t(réserve192021), base="coverage")
hillmatrix0= hill[hill[,"order"]==0,]
hillmatrix1= hill[hill[,"order"]==1,]
hillmatrix2= hill[hill[,"order"]==2,]
hillmatrix = hillmatrix0[,c(1,2,4)]
hillmatrix$SC0 = abond$DataInfo$SC
hillmatrix[,c("H0r","H0rU","H0rL")]=hillmatrix0[,c("qD","qD.UCL","qD.LCL")]
hillmatrix[,c("H1r","H1rU","H1rL")]=hillmatrix1[,c("qD","qD.UCL","qD.LCL")]
hillmatrix[,c("H2r","H2rU","H2rL")]=hillmatrix2[,c("qD","qD.UCL","qD.LCL")]
rownames(hillmatrix)=rownames(réserve192021) 

######## estimateurs diversité locale

object <- specpool(réserve192021, smallsample = TRUE)
estimateR <-estimateR(réserve192021)
#specpool2vect(object, index = c("jack1","jack2", "chao", "boot","Species"))
poolaccum <-poolaccum(réserve192021, permutations = 100, minsize = 3)
estaccumR <-estaccumR(réserve192021, permutations = 100, parallel = getOption("mc.cores"))
"summary"(object, display, alpha = 0.05)

## Accumulation model
pool <- poolaccum(réserve192021)
summary(pool, display = "chao")
plot(pool)
## Quantitative model
estimateR(réserve192021[1:3,])
