# Load the required packages
library(dplyr)
library(tidyr) # fonction gather
library(vegan) # specpool () estimateR() poolaccum() estaccumR()
library(iNEXT)
library(ggplot2)
library("readxl") 
library(forcats)
library(RODBC)
library(viridis) # couleur daltonien
#library(grDevices) # pdf() dev.off()




# Créer une connexion à la base de données Access
# Pour la configuration du fichier .mdb voir README_IMPORT_DB.Rmd
con <- odbcConnect("Memoire_Jordan")

SpecCondStat <- sqlQuery(con, "
SELECT
  PROJECT.CODE,
  SPEC.*,
  COND.*,
  STAT.*
FROM
  PROJECT
INNER JOIN
  (
    STAT
    INNER JOIN
      (
        COND
        INNER JOIN SPEC ON COND.COND_ID = SPEC.COND_ID
      ) ON STAT.STAT_ID = COND.STAT_ID
  ) ON PROJECT.PROJECT_ID = STAT.PROJECT_ID
")


# Fermer la connexion à la base de données
odbcClose(con)


# Import
# Tour PC
#setwd("D:/OneDrive - UMONS/UMONS/MA2/MEMOIRE/Statistique/R/STAT-MEMOIRE/RAW")
# Laptop
# setwd("C:/Users/Jordan/OneDrive - UMONS/UMONS/MA2/MEMOIRE/Statistique/R/STAT-MEMOIRE/RAW")
# 
# # Import
# SCS <- read_excel("C:/Users/Jordan/OneDrive - UMONS/UMONS/MA2/MEMOIRE/Statistique/R/STAT-MEMOIRE/RAW/Traits abeilles")

SCS <- SpecCondStat


# # Renommer plus simplement
# rename(SCS, "sp" = "SPEC.TAXPRIO" ) -> SCS
# rename(SCS, "genus" = "SPEC.GEN" ) -> SCS
# rename(SCS, "plante" = "COND.TAXPRIO" ) -> SCS
# rename(SCS, "pl_fam" = "COND.GR2" ) -> SCS
# rename(SCS, "family" = "SPEC.GR2" ) -> SCS
# rename(SCS, "site" = "TOPO" ) -> SCS

# Renommer plus simplement
rename(SCS, "sp" = "TAXPRIO" ) -> SCS
rename(SCS, "genus" = "GEN" ) -> SCS
rename(SCS, "plante" = "TAXPRIO.1" ) -> SCS
rename(SCS, "pl_fam" = "GR2.1" ) -> SCS
rename(SCS, "family" = "GR2" ) -> SCS
rename(SCS, "site" = "TOPO" ) -> SCS
rename(SCS, "projet" = "CODE" ) -> SCS

# Remplacement des noms d'espèces désuets
SCS$sp[SCS$sp == "Bombus (Bombus)  sp."] <- "Terrestribombus  sp."
SCS$sp[SCS$sp == "Bombus lucorum"] <- "Terrestribombus  sp."
SCS$sp[SCS$sp == "Bombus terrestris"] <- "Terrestribombus  sp."
SCS$sp[SCS$sp == "Chalicodoma ericetorum"] <- "Megachile ericetorum"
SCS$sp[SCS$sp == "Halictus tumulorum"] <- "Seladonia tumulorum"

# Retirer les observations contenant l'espèce : "Apis mellifera"
SCS <- filter(SCS, sp != "Apis mellifera")

# Remplacement du noms de sites
SCS$site[SCS$site == "Les Gourmandes de la Procession"] <- "Gourmandes"

# Retirer le site de Condorcet
SCS <- filter(SCS, site != "Condorcet" )

