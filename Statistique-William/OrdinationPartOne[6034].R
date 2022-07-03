# ORDINATIONS I
library(vegan) # Un premier lot de fonctions pour ordination
library(FactoMineR) # Deux autres packages faisant la m?me chose mais en plus esth?tique
library(factoextra)
library(ape) #PCoA
library(ade4)

#On commence par d?finir le dossier de travail avec "set working directory" :
#setwd("C:/Users/User/Document/Universit?/Assistanat/MA2/Stats/Ordinations")

# On importe un jeu de donn?es test : les poissons de la rivi?re Doubs
data (doubs)
fish <- doubs$fish# Une rivi?re, 30 sites, une longue liste de poissons
env <- doubs$env # La m?me rivi?re, les m?mes sites, mais de nombreux types de descripteurs environnementaux
?doubs #Plus d'infos
# On retire le site ! car il est vide !
fish = fish[-8,]
env = env[-8,]



### PCA ####

#Fonction PCA de vegan
res.pca <- rda(env, scale = TRUE) # On rend les donn?es homog?nes avec scale = TRUE. C'est la m?me fonction que pour calculer la rda (ordination partie 2), mais avec moins de param?tres.

summary(res.pca)
# R?sum? brut de la PCA
#Toute la dispersion (inertie) repr?sent?e est non contrainte, puisqu'il s'agit d'une ordination non contrainte. La valeur totale de la dispersion est de 11.
# Le tableau qui suit donne les eigenvalues de chaque axe. La somme de toutes les eigenvalues = l'inertie totale.
# Le pourcentage du total est indiqu? en dessous : le deuxi?me axe repr?sente 20% de l'inertie totale.
# Il est suivi du pourcentage cumul?. Les deux premiers axes repr?sentent 77,7% de la dispersion initiale.
#Viennent ensuite les coordonn?es des vecteurs propres (species scores) et des donn?es dans le nouveau syst?me d'axes (site scores).

#Visualiser les eigenvalues des axes
screeplot(res.pca) #On appelle ?a un screeplot dans le m?tier

# Quels axes sont int?ressants ? Si on avait d?coup? al?atoirement l'inertie totale, quels axes auraient-on obtenus ? On peut le savoir gr?ce ? la distribution du broken stick. Comparons nos eigenvalues ? deseigenvalues al?atoires :
screeplot(res.pca, bstick = TRUE) # En rouge, la distribution du broken stick
# Les deux premiers axes repr?sentent plus de dispersion que des axes pris au hasard, ils sont donc les plus int?ressants ? repr?senter !

# Bon, et cette repr?sentation ?
biplot(res.pca, scaling = "sites")# Cette repr?sentation est celle qui conserve la distance euclidienne entre les sites. Elle permet de voir quels sites sont proches au niveau des variables environnementales consid?r?es. Les fl?ches repr?sentent les variables et permettent, par projection, de d?terminer si les sites poss?dent des valeurs hautes ou basses pour ces variables. 
biplot(res.pca, scaling = "species", display = "species") # Ce graphe permet de visualiser les relations entre les variables. Si les directions des fl?ches sont proches, les variables sont corr?l?es entre elles. La corr?lation est positive si les fl?ches sont dansle m?me sens, n?gative si les fl?ches sont en sens inverse, et nulle si les fl?ches forment un angle droit. Les coordonn?es des vecteurs repr?sentent leur corr?lation avec les axes 1 et 2. 



# La m?me chose en plus joli avec FactoMineR : 
res.pca2 = PCA(env, scale.unit = TRUE)# calculer la PCA, donne directement scaling = 2
fviz_eig(res.pca2, addlabels = TRUE, choice = "eigenvalue") # Screeplot avec eigenvalues
fviz_eig(res.pca2, addlabels = TRUE, choice = "variance") # Screeplot avec pourcentages
fviz_pca_biplot(res.pca2) # Scaling = 1

#Petit bonus : colorier la PCA en fonction de groupes
env.w <- hclust(dist(scale(env)), "ward.D") # On cr?e ici une classification des sites
gr <- as.factor(cutree(env.w, k = 4)) # Je conserve les 4 principaux groupes de cette classification
fviz_pca_ind (res.pca2,col.ind = gr, # Je colorie les points par groupes
              repel = TRUE, # ?vite le chevauchement de texte,
              legend.title = "Groups",
              palette = c("blue", "red", "black","green"))
              
####

### CA ####
res.ca = cca(fish)
summary(res.ca)# Se lit comme une PCA
screeplot(res.ca, bstick = TRUE)# L'axe 1 repr?sente beaucoup plus de dispersion que dan sun mod?le al?atoire !
plot(res.ca, scaling =1) # Des sites proches ont des communaut?s de poissons similaires. Si un site est proche d'une esp?ce, celle-ci est probablement abondante sur ce site, ou en tout cas plus abondante en comparaison des autres sites plus ?loign?s de l'esp?ce.
plot(res.ca, scaling =2) # Des poissons proches sont distribu?s de fa?on similaire parmi les sites. Si un poisson est proche d'un site, il y abonde probablement.

# L'interpr?tation des sites et esp?ces proches de l'origine est moins certaine.


### PCoA ####
bray = vegdist(fish, method ="bray") # matrice des distances de Bray-Curtis
res.pcoa = pcoa(bray)#calcul de la PCoA
res.pcoa #R?sultat
barplot(res.pcoa$values$Eigenvalues)#Screeplot, observez les valeurs n?gatives
abline(h =0, col = "red")


bray2 = sqrt(bray)# La petite astuce qui tue
res.pcoa= pcoa(bray2)
barplot(res.pcoa$values$Eigenvalues)# Screeplot, observez les valeurs positives
abline(h =0, col = "red")
biplot.pcoa(res.pcoa, fish) #A int?rpr?ter comme une CA
####


### NMDS ####
spe.nmds <- metaMDS(fish, distance = "bray", k =2, try = 999999999999)#Indiquer le nombre d'essais
spe.nmds$stress # Valeur de stress
plot(spe.nmds,
  main = paste("NMDS/Percentage difference - Stress =",round(spe.nmds$stress, 3)
  )
)

orditorp(spe.nmds,display="species",col="red",air=0.01)
orditorp(spe.nmds,display="sites",cex=0.9,air=0.01)

#Diagramme de Shepard
stressplot(spe.nmds)


# Bonus
plot(spe.nmds,
     main = paste("NMDS/Percentage difference - Stress =",round(spe.nmds$stress, 3)
     )
)

orditorp(spe.nmds,display="species",col="red",air=0.01)
colors = c("black", "blue", "yellow", "green")[gr]
orditorp(spe.nmds,display="sites",cex=0.9,air=0.01, col = colors)

####
       