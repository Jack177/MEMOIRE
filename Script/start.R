# Load the required packages
library(dplyr)
library(tidyr) # fonction gather
library(vegan) # specpool () estimateR() poolaccum() estaccumR()
library(ggplot2)
library("readxl") 
library(forcats)
library(viridis) # couleur daltonien
#library(grDevices) # pdf() dev.off()


# Import
# Tour PC
#setwd("D:/OneDrive - UMONS/UMONS/MA2/MEMOIRE/Statistique/R/STAT-MEMOIRE/RAW")
# Laptop
setwd("C:/Users/Jordan/OneDrive - UMONS/UMONS/MA2/MEMOIRE/Statistique/R/STAT-MEMOIRE/RAW/")


library(RODBC)
my_database <- RODBC::odbcConnectAccess("C:/Users/Jordan/OneDrive - UMONS/UMONS/MA2/MEMOIRE/DFF/DATA/Clean/Mémoire_2022(clean).mdb")

RODBC::sqlTables(my_database)
scs <- RODBC::sqlFetch(channel.sadc,"SADCQ")


# Import
SCS <- read_excel("SpecCondStat.xls")

# Renommer plus simplement
rename(SCS, "sp" = "SPEC.TAXPRIO" ) -> SCS
rename(SCS, "genus" = "SPEC.GEN" ) -> SCS
rename(SCS, "plante" = "COND.TAXPRIO" ) -> SCS
rename(SCS, "pl_fam" = "COND.GR2" ) -> SCS
rename(SCS, "familly" = "SPEC.GR2" ) -> SCS
rename(SCS, "site" = "TOPO" ) -> SCS

# Remplacement des noms d'espèces désuets
SCS$sp[SCS$sp == "Bombus (Bombus)  sp."] <- "Terrestribombus  sp."
SCS$sp[SCS$sp == "Bombus lucorum"] <- "Terrestribombus  sp."
SCS$sp[SCS$sp == "Bombus terrestris"] <- "Terrestribombus  sp."
SCS$sp[SCS$sp == "Chalicodoma ericetorum"] <- "Megachile ericetorum"
SCS$sp[SCS$sp == "Halictus tumulorum"] <- "Seladonia tumulorum"
SCS$sp[SCS$sp == "Halictus tumulorum"] <- "Seladonia tumulorum"

# Retirer les observations contenant l'espèce : "Apis mellifera"
SCS <- filter(SCS, sp != "Apis mellifera" )

# Remplacement du noms de sites
SCS$site[SCS$site == "Les Gourmandes de la Procession"] <- "Gourmandes"

# Retirer le site de Condorcet
SCS <- filter(SCS, site != "Condorcet" )
