setwd(".")  # Repertoire du code

######## Chargement des fichiers ########
trutta <- read.table("../../Trutta/trutta_mean.txt", header = TRUE) 
gobio <- read.table("../../Gobio/gobio_mean.txt", header = TRUE)  
septimaniae <- read.table("../../Septimaniae/septimaniae_mean.txt", header = TRUE)
#_________________________________________________________________________________________________

######## Chargement des packages ########
library(plotrix)
#_________________________________________________________________________________________________

######## Préparation des données ########
#### Ajouts colonne "year" en fonction de G/sp ####
species_list <- list(
  gobio = list(data = gobio, generation_time = 2),       # Cottus gobio
  trutta = list(data = trutta, generation_time = 4),     # Salmo trutta
  septimaniae = list(data = septimaniae, generation_time = 2)  # Phoxinus septimaniae
)

start_year <- 1861  # Construction du premier barrage
start_generation <- 4000  # Fin de la période d'homogénéisation des fréquences alléliques

for (species in names(species_list)) {
  data <- species_list[[species]]$data
  generation_time <- species_list[[species]]$generation_time
  data$year <- start_year + (data$generation - start_generation) * generation_time
  species_list[[species]]$data <- data
}
gobio <- species_list$gobio$data
trutta <- species_list$trutta$data
septimaniae <- species_list$septimaniae$data



#### Données entre 1750 et 2025 ####
gobio_recent <- subset(gobio, year >= 1725 & year <= 2025)
trutta_recent <- subset(trutta, year >= 1725 & year <= 2025)
septimaniae_recent <- subset(septimaniae, year >= 1725 & year <= 2025)

#_________________________________________________________________________________________________


plot(trutta$n.adlt.rs ~ trutta$generation,
     xlab = "Year",
     ylab = "Number of Alleles per Locus",
     main = "Allelic Richness Over Time")

plot(gobio$n.adlt.rs ~ gobio$generation,
     xlab = "Year",
     ylab = "Number of Alleles per Locus",
     main = "Allelic Richness Over Time")



######## All. Richness ########

par(family = "serif")  

# Graph Trutta
plot(trutta_recent$n.adlt.rs ~ trutta_recent$year,
     xlab = "Year",
     ylab = "Number of Alleles per Locus",
     main = "Trutta Allelic Richness Over Time",
     xaxt = "n",
     pch = 16,
     col = "darkgreen",
     type="b")
axis(1, at = seq(min(trutta_recent$year) - 1, max(trutta_recent$year), by = 35))
mtext("a)", side = 3, line = 1, adj = 0, font = 2,cex=2)  # Ajoute "a)" en haut à gauche

# Graph  Septimaniae
plot(septimaniae_recent$n.adlt.rs ~ septimaniae_recent$year,
     xlab = "Year",
     ylab = "Number of Alleles per Locus",
     main = "Septimaniae Allelic Richness Over Time",
     xaxt = "n",
     pch = 17,
     col = "blue",
     type="b")
axis(1, at = seq(min(septimaniae_recent$year) - 1, max(septimaniae_recent$year), by = 25))
mtext("b)", side = 3, line = 1, adj = 0, font = 2,cex=2)  

# Graphique pour Gobio
plot(gobio_recent$n.adlt.rs ~ gobio_recent$year,
     xlab = "Year",
     ylab = "Number of Alleles per Locus",
     main = "Gobio Allelic Richness Over Time",
     xaxt = "n",
     pch = 18,
     col = "red3",
     type="b")
axis(1, at = seq(min(gobio_recent$year) - 1, max(gobio_recent$year), by = 25))
mtext("c)", side = 3, line = 1, adj = 0, font = 2,cex=2)  # Ajoute "c)" en haut à gauche


######## Fst ########

