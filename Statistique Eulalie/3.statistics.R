
######################## Introductory lines

# set the working directory, which is the folder containing the files
# in my case:
#setwd("formation stat")


######################################### STATISTICS

######################## Previous: check data
# check there is no correlation between variables

library('corrplot')

database<-read.table("lineal.txt", header=T)

M <- cor(database)
corrplot(M, method = "circle")

cor.test(database$Wild_abundance,database$Flower_abundance)#two variables
ggplot(database, aes(x=Flower_abundance,y=Wild_abundance)) +
  geom_point(alpha=0.3) +
  theme_classic()

########################## PCA, fait un tets avce des varibles depensdates pour pvr chosir les mieux 
#sinon normalement fais avec des independnats , Proportion of Variance explain pour chaque axe, regrde si c'est haut
#regarde si les 2 premiers axes "explaine" bcp pour pvr l'utiliser '

database<-read.table("lineal.txt", header=T)
databasepca <- database[,6:11]

prcomp(databasepca, scale = TRUE)
plot(prcomp(databasepca))
summary(prcomp(databasepca, scale = TRUE))
biplot(prcomp(databasepca, scale = TRUE))


######################### CCA, canonical , voir comment 2 groupes sont corr?ll? 

data(varespec)
data(varechem)
## Common but bad way: use all variables you happen to have in your
## environmental data matrix
vare.cca <- cca(varespec, varechem)
vare.cca
plot(vare.cca)
## Formula interface and a better model
vare.cca <- cca(varespec ~ Al + P*(K + Baresoil), data=varechem)
vare.cca
plot(vare.cca)


database <- read.table("lineal.txt",header=T)
bees <-read.table("Composition_pollinators.txt", header=T)

vare.cca <- cca(bees ~ Flower_abundance+Flower_richness, data=database)
vare.cca
plot(vare.cca)



########################## Lineal models. Model selection + model averaging

library(MuMIn)
library(caret)

options(na.action = "na.fail")

database <- read.table("lineal.txt",header=T)
names(database)

hist(database$Flower_abundance) #check normality every variable. Not normal
hist(log(database$Flower_abundance)) #better
hist(sqrt(database$Flower_abundance))

fit <- lm(Pollinator_richness~log(Flower_abundance)+Flower_richness+Honeybee_rate+Flower_richness*Honeybee_rate, data=database) # lineal model. The dependent variable has to be normal. If not, GLM

# fit <- glm(Heterospecific_presence~log(Pollinator_richness)+Visitation_rate+Proportion_plant+Proportion_HB+Proportion_Bee+Proportion_Diptera,family=binomial, weights=Individuals_pollen, data=meandataperplotTVUF)

car::vif(fit) # check correlations between variables. There are. Remove the highest. Check

fit <- lm(Pollinator_richness~log(Flower_abundance)+Flower_richness+Honeybee_rate, data=database)

car::vif(fit) # perfect. They have to be all less than 4 in value

hist(resid(fit)) # check normality of residuals. Normal enough.

### Model selection
fit <- lm(Pollinator_richness~log(Flower_abundance)+Flower_richness+Honeybee_rate, data=database)
dd <- dredge(fit,extra="adjR^2")
ddd <- subset(dd, delta < 2) # select the ones with value of AICc less than 2 points in difference
subset(dd, delta < 2)
### Model averaging
avgmod.95delta2 <- model.avg(ddd) 
summary(avgmod.95delta2) # select the conditional average

