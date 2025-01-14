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
gobio_recent <- subset(gobio, year >= 1725 & year <= 2031)
trutta_recent <- subset(trutta, year >= 1725 & year <= 2041)
septimaniae_recent <- subset(septimaniae, year >= 1725 & year <= 2031)

#_________________________________________________________________________________________________


plot(trutta$n.adlt.rs ~ trutta$generation,
     xlab = "Gen",
     ylab = "Number of Alleles per Locus",
     main = "Trutta Allelic Richness Over Time")

plot(gobio$n.adlt.rs ~ gobio$generation,
     xlab = "Gen",
     ylab = "Number of Alleles per Locus",
     main = "Gobio Allelic Richness Over Time")

plot(septimaniae$n.adlt.rs ~ septimaniae$generation,
     xlab = "Gen",
     ylab = "Number of Alleles per Locus",
     main = "Septimaniae Allelic Richness Over Time")


######## All. Richness ########

par(family = "serif")  

# Graph Trutta
plot(trutta_recent$n.adlt.rs ~ trutta_recent$year,
     xlab = "Années",
     ylab = "Nbr. All./Locus",
     cex.lab = 1.3,
     xaxt = "n",
     pch = 16,
     col = "darkgreen",
     type="b",
     bty="l",
     xlim = c(min(trutta_recent$year), 2025))
axis(1, at = seq(min(trutta_recent$year) - 1, max(trutta_recent$year), by = 35))
abline(v = c(1861, 1953), col = "red1", lty = c(5,3),lwd=1.5)
mtext("a)", side = 3, line = 1, adj = 0, font = 2,cex=2)  # Ajoute "a)" en haut à gauche

# Graph  Septimaniae
plot(septimaniae_recent$n.adlt.rs ~ septimaniae_recent$year,
     xlab = "Années",
     ylab = "Nbr. All./Locus",
     cex.lab = 1.3,
     xaxt = "n",
     pch = 17,
     col = "blue",
     type="b",
     bty="l",
     xlim = c(min(septimaniae_recent$year), 2025))
axis(1, at = seq(min(trutta_recent$year) - 1, max(trutta_recent$year), by = 35)) #Ici echelle de truite, pour avoir même axe en x partout
abline(v = c(1861, 1953), col = "red1", lty = c(5,3),lwd=1.5)
mtext("b)", side = 3, line = 1, adj = 0, font = 2,cex=2)  

# Graphique pour Gobio
plot(gobio_recent$n.adlt.rs ~ gobio_recent$year,
     xlab = "Années",
     ylab = "Nbr. All./Locus",
     cex.lab = 1.3,
     xaxt = "n",
     pch = 18,
     col = "darkgoldenrod1",
     type="b",
     bty="l",
     xlim = c(min(septimaniae_recent$year), 2025))
axis(1, at = seq(min(trutta_recent$year) - 1, max(trutta_recent$year), by = 35))
abline(v = c(1861, 1953), col = "red1", lty = c(5,3),lwd=1.5)
mtext("c)", side = 3, line = 1, adj = 0, font = 2,cex=2)  # Ajoute "c)" en haut à gauche


######## Fst ########
#### Préparation data ####

# Renommer les colonnes de chaque tableau pour éviter les conflits pour merge
colnames(gobio_recent)[colnames(gobio_recent) %in% c("n.adlt.fst.wc_p1.2", "n.adlt.fst.wc_p1.3", "n.adlt.fst.wc_p2.3")] <- 
  c("fst_gobio_1_2", "fst_gobio_1_3", "fst_gobio_2_3")
colnames(septimaniae_recent)[colnames(septimaniae_recent) %in% c("n.adlt.fst.wc_p1.2", "n.adlt.fst.wc_p1.3", "n.adlt.fst.wc_p2.3")] <- 
  c("fst_septimaniae_1_2", "fst_septimaniae_1_3", "fst_septimaniae_2_3")
colnames(trutta_recent)[colnames(trutta_recent) %in% c("n.adlt.fst.wc_p1.2", "n.adlt.fst.wc_p1.3", "n.adlt.fst.wc_p2.3")] <- 
  c("fst_trutta_1_2", "fst_trutta_1_3", "fst_trutta_2_3")

# Colonnes pertinentes
gobio_subset <- gobio_recent[, c("generation", "year", "fst_gobio_1_2", "fst_gobio_1_3", "fst_gobio_2_3")]
septimaniae_subset <- septimaniae_recent[, c("generation", "year", "fst_septimaniae_1_2", "fst_septimaniae_1_3", "fst_septimaniae_2_3")]
trutta_subset <- trutta_recent[, c("generation", "year", "fst_trutta_1_2", "fst_trutta_1_3", "fst_trutta_2_3")]

# Fusionner les tableaux par "generation" et "year"
merged_table <- Reduce(function(x, y) merge(x, y, by = "generation", all = TRUE),
                       list(gobio_subset, septimaniae_subset, trutta_subset))



















