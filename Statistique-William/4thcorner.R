
# 4th corner analysis (Legendre's version)

# Reference papers :
# First version of the analysis in this paper (useful to understand the second paper) :
# http://biol09.biol.umontreal.ca/numecol/Reprints/4th-corner_paper.pdf
# current version of the analysis used here (explanation of the different models):
# https://www.researchgate.net/publication/250076886_Testing_the_species_traits_environment_relationships_The_fourth-corner_problem_revisited

library(readxl)
library(ade4)
Paysage2 <- read_excel("C:/Users/User/Desktop/R.xlsx")
Paysage2$...1->a
Paysage2 = Paysage2[,-1]
rownames(Paysage2)= a
BDD_sec2 <- read_excel("C:/Users/User/Desktop/L.xlsx")
BDD_sec2$...1->a
BDD_sec2 = BDD_sec2[,-1]
rownames(BDD_sec2)= a
traits <- read_excel("C:/Users/User/Desktop/Q.xlsx")
traits$...1->a
traits = traits[,-1]
traits<-as.data.frame(apply(X = traits[,c(1:3)],FUN = function(x)as.factor(x), MARGIN = 2))
rownames(traits)= a


# The code :
library("ade4")
fourth <- fourthcorner(
  tabR = Paysage2, # Data on landscape (rows= stations, columns = environmental variables)
  tabL = BDD_sec2, # table species (columns) x stations (rows)
  tabQ = traits, # traits (species = rows, columns = traits)
  modeltype = 6, # kind of permutations of columns/rows applied to tabL. This is the best version (explained in reference paper 2)
  p.adjust.method.G = "none",
  p.adjust.method.D = "none", 
  nrepet = 999) # number of permutations. Should be really high, the higher the better

# Correction for multiple testing, here using FDR
fourth.adj <- p.adjust.4thcorner( #corriger la p-value
  fourth,
  p.adjust.method.G = "fdr",
  p.adjust.method.D = "fdr",
  p.adjust.D = "global")

# Plot
plot(fourth.adj, alpha = 0.05, stat = "D2") # seuil signf = 5%; D2 c'est le mieux, si on veut homogéneiter on met D tout seul

# Three stats can be computed :
# D2 = correlation 
# D = homogeneity of each category (for qualitative variables) 
# G is an  anova like stat for qualitative variables