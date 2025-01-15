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

#_________________________________________________________________________________________________

######## Fst ########
#### Visualisation data ####

# Renommer les colonnes de chaque tableau pour éviter les conflits 
colnames(gobio_recent)[colnames(gobio_recent) %in% c("n.adlt.fst.wc_p1.2", "n.adlt.fst.wc_p1.3", "n.adlt.fst.wc_p2.3")] <- 
  c("fst_gobio_1_2", "fst_gobio_1_3", "fst_gobio_2_3")
colnames(septimaniae_recent)[colnames(septimaniae_recent) %in% c("n.adlt.fst.wc_p1.2", "n.adlt.fst.wc_p1.3", "n.adlt.fst.wc_p2.3")] <- 
  c("fst_septimaniae_1_2", "fst_septimaniae_1_3", "fst_septimaniae_2_3")
colnames(trutta_recent)[colnames(trutta_recent) %in% c("n.adlt.fst.wc_p1.2", "n.adlt.fst.wc_p1.3", "n.adlt.fst.wc_p2.3")] <- 
  c("fst_trutta_1_2", "fst_trutta_1_3", "fst_trutta_2_3")

#### Barrage A: Fst entre pop 1 et 2 ####
## Graphique ##
#  Fst / (1 - Fst)
gobio_recent$fst_gobio_ratio <- gobio_recent$fst_gobio_1_2 / (1 - gobio_recent$fst_gobio_1_2)
septimaniae_recent$fst_septimaniae_ratio <- septimaniae_recent$fst_septimaniae_1_2 / (1 - septimaniae_recent$fst_septimaniae_1_2)
trutta_recent$fst_trutta_ratio <- trutta_recent$fst_trutta_1_2 / (1 - trutta_recent$fst_trutta_1_2)

# Définir les limites des axes
x_lim <- c(min(c(gobio_recent$year, septimaniae_recent$year, trutta_recent$year)), 2025)
y_lim <- c(min(c(min(gobio_recent$fst_gobio_ratio, na.rm = TRUE),
        min(septimaniae_recent$fst_septimaniae_ratio, na.rm = TRUE),
        min(trutta_recent$fst_trutta_ratio, na.rm = TRUE))),
        max(c(max(gobio_recent$fst_gobio_ratio, na.rm = TRUE),
        max(septimaniae_recent$fst_septimaniae_ratio, na.rm = TRUE),
        max(trutta_recent$fst_trutta_ratio, na.rm = TRUE)))+0.01)

# Graph vide avec bonnes limites
par(mgp = c(2.5, 1, 0))
plot(1, type = "n", 
     xlab = "Années", 
     ylab = expression(italic(F[ST]/(1 - F[ST]))), 
     xlim = x_lim, 
     ylim = y_lim, 
     cex.lab = 1.3, 
     bty = "l", 
     xaxt = "n")
axis(1, at = seq(min(trutta_recent$year)-1, 2025, by = 35))

# Points pour chaque espèce
points(gobio_recent$year, gobio_recent$fst_gobio_1_2 / (1 - gobio_recent$fst_gobio_1_2), 
       col = "darkgoldenrod1", pch = 16, type = "b", lwd = 1.5) 
points(septimaniae_recent$year, septimaniae_recent$fst_septimaniae_1_2, 
       col = "blue", pch = 17, type = "b", lwd = 1.5)
points(trutta_recent$year, trutta_recent$fst_trutta_1_2, 
       col = "darkgreen", pch = 18, type = "b", lwd = 1.5)
abline(v = 1861, col = "red1", lty = 5, lwd = 1.5)
legend("topleft", 
       legend = c(expression(italic("C. gobio")), 
                  expression(italic("P. septimaniae")), 
                  expression(italic("S. trutta"))),
       col = c("darkgoldenrod1", "blue", "darkgreen"), 
       pch = c(16, 17, 18), 
       lty = 1,
       bty="n",
       cex=1.3)
mtext("a)", side = 3, line = 1, adj = 0, font = 2,cex=2)


#### Barrage B: Fst entre pop 2 et 3 ####

## Fst / (1 - Fst)
gobio_recent$fst_gobio_ratio <- gobio_recent$fst_gobio_2_3 / (1 - gobio_recent$fst_gobio_2_3)
septimaniae_recent$fst_septimaniae_ratio <- septimaniae_recent$fst_septimaniae_2_3 / (1 - septimaniae_recent$fst_septimaniae_2_3)
trutta_recent$fst_trutta_ratio <- trutta_recent$fst_trutta_2_3 / (1 - trutta_recent$fst_trutta_2_3)

## limites des axes
x_lim <- c(min(c(gobio_recent$year, septimaniae_recent$year, trutta_recent$year)), 2025)
y_lim <- c(min(c(min(gobio_recent$fst_gobio_1_2, na.rm = TRUE),
                 min(septimaniae_recent$fst_septimaniae_1_2, na.rm = TRUE),
                 min(trutta_recent$fst_trutta_1_2, na.rm = TRUE))),
           max(c(max(gobio_recent$fst_gobio_1_2, na.rm = TRUE),
                 max(septimaniae_recent$fst_septimaniae_1_2, na.rm = TRUE),
                 max(trutta_recent$fst_trutta_1_2, na.rm = TRUE)))+0.0525) # pour avoir même echelle en Y que graphe a)

par(mgp = c(2.5, 1, 0)) 

plot(1, type = "n", 
     xlab = "Années", 
     ylab = " ", 
     xlim = x_lim, 
     ylim = y_lim, 
     cex.lab = 1.3, 
     bty = "l", 
     xaxt = "n")

axis(1, at = seq(min(trutta_recent$year) - 1, 2025, by = 35))

points(gobio_recent$year, gobio_recent$fst_gobio_ratio, 
       col = "darkgoldenrod1", pch = 16, type = "b", lwd = 1.5)
points(septimaniae_recent$year, septimaniae_recent$fst_septimaniae_ratio, 
       col = "blue", pch = 17, type = "b", lwd = 1.5)
points(trutta_recent$year, trutta_recent$fst_trutta_ratio, 
       col = "darkgreen", pch = 18, type = "b", lwd = 1.5)

abline(v = 1953, col = "red1", lty = 3, lwd = 1.5)
mtext("b)", side = 3, line = 1, adj = 0, font = 2, cex = 2)
#_________________________________________________________________________________________________

######## Chargement des fichiers pour simu. futures ########
trutta_sans_A <- read.table("../../actions_plans/Remove_bar_A/Trutta_rem_A/trutta_rem_A_mean.txt", header = TRUE) 
gobio_sans_A <- read.table("../../actions_plans/Remove_bar_A/Gobio_rem_A/gobio_rem_A_mean.txt", header = TRUE) 
septimaniae_sans_A <- read.table("../../actions_plans/Remove_bar_A/Septimaniae_rem_A/septimaniae_rem_A_mean.txt", header = TRUE)

trutta_sans_B <- read.table("../../actions_plans/Remove_bar_B/Trutta_rem_B/trutta_rem_B_mean.txt", header = TRUE) 
gobio_sans_B <- read.table("../../actions_plans/Remove_bar_B/Gobio_rem_B/gobio_rem_B_mean.txt", header = TRUE) 
septimaniae_sans_B <- read.table("../../actions_plans/Remove_bar_B/Septimaniae_rem_B/septimaniae_rem_B_mean.txt", header = TRUE)


## Dates 
species_list <- list(
  gobio_sans_A = list(data = gobio_sans_A, generation_time = 2),       # Cottus gobio_sans_A
  trutta_sans_A = list(data = trutta_sans_A, generation_time = 4),     # Salmo trutta_sans_A
  septimaniae_sans_A = list(data = septimaniae_sans_A, generation_time = 2)  # Phoxinus septimaniae_sans_A
)
start_year <- 1861  # Construction du premier barrage
start_generation <- 4000  # Fin de la période d'homogénéisation des fréquences alléliques
for (species in names(species_list)) {
  data <- species_list[[species]]$data
  generation_time <- species_list[[species]]$generation_time
  data$year <- start_year + (data$generation - start_generation) * generation_time
  species_list[[species]]$data <- data
}
gobio_sans_A <- species_list$gobio_sans_A$data
trutta_sans_A <- species_list$trutta_sans_A$data
septimaniae_sans_A <- species_list$septimaniae_sans_A$data

species_list <- list(
  gobio_sans_B = list(data = gobio_sans_B, generation_time = 2),       # Cottus gobio_sans_B
  trutta_sans_B = list(data = trutta_sans_B, generation_time = 4),     # Salmo trutta_sans_B
  septimaniae_sans_B = list(data = septimaniae_sans_B, generation_time = 2)  # Phoxinus septimaniae_sans_B
)
start_year <- 1861  # Construction du premier barrage
start_generation <- 4000  # Fin de la période d'homogénéisation des fréquences alléliques
for (species in names(species_list)) {
  data <- species_list[[species]]$data
  generation_time <- species_list[[species]]$generation_time
  data$year <- start_year + (data$generation - start_generation) * generation_time
  species_list[[species]]$data <- data
}
gobio_sans_B <- species_list$gobio_sans_B$data
trutta_sans_B <- species_list$trutta_sans_B$data
septimaniae_sans_B <- species_list$septimaniae_sans_B$data




gobio_sans_A_futur <- subset(gobio_sans_A, year >= 1861 & year <= 2541)
trutta_sans_A_futur <- subset(trutta_sans_A, year >= 1861 & year <= 2541)
septimaniae_sans_A_futur <- subset(septimaniae_sans_A, year >= 1861 & year <= 2541)

gobio_sans_B_futur <- subset(gobio_sans_B, year >= 1861 & year <= 2541)
trutta_sans_B_futur <- subset(trutta_sans_B, year >= 1861 & year <= 2541)
septimaniae_sans_B_futur <- subset(septimaniae_sans_B, year >= 1861 & year <= 2541)
#_________________________________________________________________________________________________

#### Tableau TRUTTA ####
# Ajouter des préfixes aux colonnes pour différencier les données des trois tableaux
colnames(trutta) <- c("generation", "alive.rpl", "richesse_trutta",
                      "fst_1_2_trutta", "fst_1_3_trutta", "fst_2_3_trutta", "year")

colnames(trutta_sans_A_futur) <- c("generation", "alive.rpl", "richesse_trutta_sans_A",
                                   "fst_1_2_trutta_sans_A", "fst_1_3_trutta_sans_A", 
                                   "fst_2_3_trutta_sans_A", "year")

colnames(trutta_sans_B_futur) <- c("generation", "alive.rpl", "richesse_trutta_sans_B",
                                   "fst_1_2_trutta_sans_B", "fst_1_3_trutta_sans_B", 
                                   "fst_2_3_trutta_sans_B", "year")

# Sélectionner uniquement les colonnes d'intérêt pour chaque tableau
trutta_subset <- trutta[, c("generation", "year", 
                            "richesse_trutta", "fst_1_2_trutta", "fst_2_3_trutta")]
trutta_sans_A_subset <- trutta_sans_A_futur[, c("generation", "year", 
                                                "richesse_trutta_sans_A", 
                                                "fst_1_2_trutta_sans_A", 
                                                "fst_2_3_trutta_sans_A")]
trutta_sans_B_subset <- trutta_sans_B_futur[, c("generation", "year", 
                                                "richesse_trutta_sans_B", 
                                                "fst_1_2_trutta_sans_B", 
                                                "fst_2_3_trutta_sans_B")]
# Fusionner les tableaux 
Trutta_predictions <- Reduce(function(x, y) merge(x, y, by = c("generation", "year"), all = FALSE),
                             list(trutta_subset, trutta_sans_A_subset, trutta_sans_B_subset))
#_________________________________________________________________________________________________

#### Tableau GOBIO ####
colnames(gobio) <- c("generation", "alive.rpl", "richesse_gobio",
                     "fst_1_2_gobio", "fst_1_3_gobio", "fst_2_3_gobio", "year")

colnames(gobio_sans_A_futur) <- c("generation", "alive.rpl", "richesse_gobio_sans_A",
                                  "fst_1_2_gobio_sans_A", "fst_1_3_gobio_sans_A", 
                                  "fst_2_3_gobio_sans_A", "year")

colnames(gobio_sans_B_futur) <- c("generation", "alive.rpl", "richesse_gobio_sans_B",
                                  "fst_1_2_gobio_sans_B", "fst_1_3_gobio_sans_B", 
                                  "fst_2_3_gobio_sans_B", "year")

gobio_subset <- gobio[, c("generation", "year", 
                          "richesse_gobio", "fst_1_2_gobio", "fst_2_3_gobio")]
gobio_sans_A_subset <- gobio_sans_A_futur[, c("generation", "year", 
                                              "richesse_gobio_sans_A", 
                                              "fst_1_2_gobio_sans_A", 
                                              "fst_2_3_gobio_sans_A")]
gobio_sans_B_subset <- gobio_sans_B_futur[, c("generation", "year", 
                                              "richesse_gobio_sans_B", 
                                              "fst_1_2_gobio_sans_B", 
                                              "fst_2_3_gobio_sans_B")]

Gobio_predictions <- Reduce(function(x, y) merge(x, y, by = c("generation", "year"), all = FALSE),
                            list(gobio_subset, gobio_sans_A_subset, gobio_sans_B_subset))
#_________________________________________________________________________________________________

#### Tableau SEPTIMANIAE ####

colnames(septimaniae) <- c("generation", "alive.rpl", "richesse_septimaniae",
                           "fst_1_2_septimaniae", "fst_1_3_septimaniae", "fst_2_3_septimaniae", "year")

colnames(septimaniae_sans_A_futur) <- c("generation", "alive.rpl", "richesse_septimaniae_sans_A",
                                        "fst_1_2_septimaniae_sans_A", "fst_1_3_septimaniae_sans_A", 
                                        "fst_2_3_septimaniae_sans_A", "year")

colnames(septimaniae_sans_B_futur) <- c("generation", "alive.rpl", "richesse_septimaniae_sans_B",
                                        "fst_1_2_septimaniae_sans_B", "fst_1_3_septimaniae_sans_B", 
                                        "fst_2_3_septimaniae_sans_B", "year")

septimaniae_subset <- septimaniae[, c("generation", "year", 
                                      "richesse_septimaniae", "fst_1_2_septimaniae", "fst_2_3_septimaniae")]
septimaniae_sans_A_subset <- septimaniae_sans_A_futur[, c("generation", "year", 
                                                          "richesse_septimaniae_sans_A", 
                                                          "fst_1_2_septimaniae_sans_A", 
                                                          "fst_2_3_septimaniae_sans_A")]
septimaniae_sans_B_subset <- septimaniae_sans_B_futur[, c("generation", "year", 
                                                          "richesse_septimaniae_sans_B", 
                                                          "fst_1_2_septimaniae_sans_B", 
                                                          "fst_2_3_septimaniae_sans_B")]

Septimaniae_predictions <- Reduce(function(x, y) merge(x, y, by = c("generation", "year"), all = FALSE),
                                  list(septimaniae_subset, septimaniae_sans_A_subset, septimaniae_sans_B_subset))
#_________________________________________________________________________________________________

######## Prediction Rich. All. /sp ########
#### S. Trutta ####
# Définir les limites de l'axe des Y
y_lim <- c(min(c(min(Trutta_predictions$richesse_trutta, na.rm = TRUE),
                 min(Trutta_predictions$richesse_trutta_sans_A, na.rm = TRUE),
                 min(Trutta_predictions$richesse_trutta_sans_B, na.rm = TRUE))),
           max(c(max(Trutta_predictions$richesse_trutta, na.rm = TRUE),
                 max(Trutta_predictions$richesse_trutta_sans_A, na.rm = TRUE),
                 max(Trutta_predictions$richesse_trutta_sans_B, na.rm = TRUE))))

# Tracer un graphique vide
plot(Trutta_predictions$year, Trutta_predictions$richesse_trutta, 
     type = "n",  # Ne trace rien, juste le cadre
     xlab = "Années", 
     ylab = "Nbr. All./Locus", 
     ylim = y_lim,
     cex.lab = 1.3, 
     bty = "l")

# Ajouter une zone grisée avant 2025
rect(xleft = 0, xright = 2025, 
     ybottom = par("usr")[3], ytop = par("usr")[4], 
     col = rgb(0.3, 0.3, 0.3, 0.3), border = NA)
# Ajouter les points pour chaque colonne
points(Trutta_predictions$year, Trutta_predictions$richesse_trutta, 
       col = "grey50", pch = 18, type = "b", lwd = 1.5)  # Richesse trutta (colonne 3)
points(Trutta_predictions$year, Trutta_predictions$richesse_trutta_sans_A, 
       col = "dodgerblue", pch = 17, type = "b", lwd = 1.5)  # Richesse trutta sans A (colonne 6)
points(Trutta_predictions$year, Trutta_predictions$richesse_trutta_sans_B, 
       col = "chocolate1", pch = 16, type = "b", lwd = 1.5)  # Richesse trutta sans B (colonne 9)

# Ajouter une légende
legend("topright", 
       legend = c("Aucun", "Barrage A", "Barrage B"),
       col = c("grey50", "dodgerblue", "chocolate1"), 
       pch = c(18, 17, 16), 
       lty = 1, 
       bty = "n",
       cex=1.3)
mtext("a)", side = 3, line = 1, adj = 0, font = 2,cex=2)  
#_________________________________________________________________________________________________

#### C. Gobio ####
# Définir les limites de l'axe des Y pour Gobio
y_lim_gobio <- c(
  min(c(min(Gobio_predictions$richesse_gobio, na.rm = TRUE),
        min(Gobio_predictions$richesse_gobio_sans_A, na.rm = TRUE),
        min(Gobio_predictions$richesse_gobio_sans_B, na.rm = TRUE))),
  max(c(max(Gobio_predictions$richesse_gobio, na.rm = TRUE),
        max(Gobio_predictions$richesse_gobio_sans_A, na.rm = TRUE),
        max(Gobio_predictions$richesse_gobio_sans_B, na.rm = TRUE)))
)

# Tracer un graphique vide pour Gobio
plot(Gobio_predictions$year, Gobio_predictions$richesse_gobio, 
     type = "n",  # Ne trace rien, juste le cadre
     xlab = "Années", 
     ylab = "Nbr. All./Locus", 
     ylim = y_lim_gobio,
     cex.lab = 1.3, 
     bty = "l")

rect(xleft = 0, xright = 2025, 
     ybottom = par("usr")[3], ytop = par("usr")[4], 
     col = rgb(0.3, 0.3, 0.3, 0.3), border = NA)

# Ajouter les points pour chaque colonne
points(Gobio_predictions$year, Gobio_predictions$richesse_gobio, 
       col = "grey50", pch = 18, type = "b", lwd = 1.5)  # Richesse gobio
points(Gobio_predictions$year, Gobio_predictions$richesse_gobio_sans_A, 
       col = "dodgerblue", pch = 17, type = "b", lwd = 1.5)  # Richesse gobio sans A
points(Gobio_predictions$year, Gobio_predictions$richesse_gobio_sans_B, 
       col = "chocolate1", pch = 16, type = "b", lwd = 1.5)  # Richesse gobio sans B

mtext("c)", side = 3, line = 1, adj = 0, font = 2,cex=2)  
#_________________________________________________________________________________________________

#### Septimaniae ####
# Définir les limites de l'axe des Y pour Septimaniae
y_lim_septimaniae <- c(
  min(c(min(Septimaniae_predictions$richesse_septimaniae, na.rm = TRUE),
        min(Septimaniae_predictions$richesse_septimaniae_sans_A, na.rm = TRUE),
        min(Septimaniae_predictions$richesse_septimaniae_sans_B, na.rm = TRUE))),
  max(c(max(Septimaniae_predictions$richesse_septimaniae, na.rm = TRUE),
        max(Septimaniae_predictions$richesse_septimaniae_sans_A, na.rm = TRUE),
        max(Septimaniae_predictions$richesse_septimaniae_sans_B, na.rm = TRUE)))
)

# Tracer un graphique vide pour Septimaniae
plot(Septimaniae_predictions$year, Septimaniae_predictions$richesse_septimaniae, 
     type = "n",  # Ne trace rien, juste le cadre
     xlab = "Années", 
     ylab = "Nbr. All./Locus", 
     ylim = y_lim_septimaniae,
     cex.lab = 1.3, 
     bty = "l")

rect(xleft = 0, xright = 2025, 
     ybottom = par("usr")[3], ytop = par("usr")[4], 
     col = rgb(0.3, 0.3, 0.3, 0.3), border = NA)

# Ajouter les points pour chaque colonne
points(Septimaniae_predictions$year, Septimaniae_predictions$richesse_septimaniae, 
       col = "grey50", pch = 18, type = "b", lwd = 1.5)  # Richesse septimaniae
points(Septimaniae_predictions$year, Septimaniae_predictions$richesse_septimaniae_sans_A, 
       col = "dodgerblue", pch = 17, type = "b", lwd = 1.5)  # Richesse septimaniae sans A
points(Septimaniae_predictions$year, Septimaniae_predictions$richesse_septimaniae_sans_B, 
       col = "chocolate1", pch = 16, type = "b", lwd = 1.5)  # Richesse septimaniae sans B
mtext("b)", side = 3, line = 1, adj = 0, font = 2,cex=2) 
#_________________________________________________________________________________________________

######## Prédictions Fst selon sp et traitements ########
#### Préparation data Fst barrage A ####
# Calculer Fst/(1-Fst) pour Trutta
Trutta_predictions$fst_ratio_trutta <- Trutta_predictions$fst_1_2_trutta / (1 - Trutta_predictions$fst_1_2_trutta)
Trutta_predictions$fst_ratio_trutta_sans_A <- Trutta_predictions$fst_1_2_trutta_sans_A / (1 - Trutta_predictions$fst_1_2_trutta_sans_A)

# Calculer Fst/(1-Fst) pour Gobio
Gobio_predictions$fst_ratio_gobio <- Gobio_predictions$fst_1_2_gobio / (1 - Gobio_predictions$fst_1_2_gobio)
Gobio_predictions$fst_ratio_gobio_sans_A <- Gobio_predictions$fst_1_2_gobio_sans_A / (1 - Gobio_predictions$fst_1_2_gobio_sans_A)

# Calculer Fst/(1-Fst) pour Septimaniae
Septimaniae_predictions$fst_ratio_septimaniae <- Septimaniae_predictions$fst_1_2_septimaniae / (1 - Septimaniae_predictions$fst_1_2_septimaniae)
Septimaniae_predictions$fst_ratio_septimaniae_sans_A <- Septimaniae_predictions$fst_1_2_septimaniae_sans_A / (1 - Septimaniae_predictions$fst_1_2_septimaniae_sans_A)


# Définir les limites des axes
y_lim_fst_ratio <- c(
  min(c(min(Trutta_predictions$fst_ratio_trutta, na.rm = TRUE),
        min(Trutta_predictions$fst_ratio_trutta_sans_A, na.rm = TRUE),
        min(Gobio_predictions$fst_ratio_gobio, na.rm = TRUE),
        min(Gobio_predictions$fst_ratio_gobio_sans_A, na.rm = TRUE),
        min(Septimaniae_predictions$fst_ratio_septimaniae, na.rm = TRUE),
        min(Septimaniae_predictions$fst_ratio_septimaniae_sans_A, na.rm = TRUE))),
  max(c(max(Trutta_predictions$fst_ratio_trutta, na.rm = TRUE),
        max(Trutta_predictions$fst_ratio_trutta_sans_A, na.rm = TRUE),
        max(Gobio_predictions$fst_ratio_gobio, na.rm = TRUE),
        max(Gobio_predictions$fst_ratio_gobio_sans_A, na.rm = TRUE),
        max(Septimaniae_predictions$fst_ratio_septimaniae, na.rm = TRUE),
        max(Septimaniae_predictions$fst_ratio_septimaniae_sans_A, na.rm = TRUE)))
)

#### Graphique Fst barrage A ####
par(mgp = c(2.5, 1, 0)) 
plot(1, type = "n", 
     xlab = "Années", 
     ylab = expression(italic(F[ST]/(1 - F[ST]))), 
     xlim = c(min(Trutta_predictions$year), max(Trutta_predictions$year)), 
     ylim = y_lim_fst_ratio, 
     cex.lab = 1.3, 
     bty = "l")

# Définir les couleurs avec transparence pour les tirets
transparent_darkgreen <- adjustcolor("darkgreen", alpha.f = 0.7)       # Trutta
transparent_darkgoldenrod <- adjustcolor("darkgoldenrod1", alpha.f = 0.7) # Gobio
transparent_blue <- adjustcolor("blue", alpha.f = 0.7)                 # Septimaniae

# Ajouter les courbes pour Trutta
lines(Trutta_predictions$year, Trutta_predictions$fst_ratio_trutta, 
      col = transparent_darkgreen, lty = 2, lwd = 1.75)  # Barrage intact (tirets, transparent)
lines(Trutta_predictions$year, Trutta_predictions$fst_ratio_trutta_sans_A, 
      col = "darkgreen", lty = 1, lwd = 1.75)  # Sans barrage A (ligne continue)

# Ajouter les courbes pour Gobio
lines(Gobio_predictions$year, Gobio_predictions$fst_ratio_gobio, 
      col = transparent_darkgoldenrod, lty = 2, lwd = 1.75)  # Barrage intact (tirets, transparent)
lines(Gobio_predictions$year, Gobio_predictions$fst_ratio_gobio_sans_A, 
      col = "darkgoldenrod1", lty = 1, lwd = 1.75)  # Sans barrage A (ligne continue)

# Ajouter les courbes pour Septimaniae
lines(Septimaniae_predictions$year, Septimaniae_predictions$fst_ratio_septimaniae, 
      col = transparent_blue, lty = 2, lwd = 1.75)  # Barrage intact (tirets, transparent)
lines(Septimaniae_predictions$year, Septimaniae_predictions$fst_ratio_septimaniae_sans_A, 
      col = "blue", lty = 1, lwd = 1.75)  # Sans barrage A (ligne continue)

# Ajouter la zone grisée avant 2025
rect(xleft = 0, xright = 2025, 
     ybottom = par("usr")[3], ytop = par("usr")[4], 
     col = rgb(0.3, 0.3, 0.3, 0.3), border = NA)

# Ajouter une légende
legend("topleft", 
       legend = c(expression(italic("C. gobio")), 
                  expression(italic("P. septimaniae")), 
                  expression(italic("S. trutta"))),
       col = c("darkgoldenrod1", "blue", "darkgreen"),  # Couleurs des espèces
       pch = 15,  # Rectangle coloré pour chaque espèce
       bty = "n",  # Pas de bordure autour de la légende
       pt.cex = 1.5,  # Taille des rectangles de couleur
       cex = 1.3)  # Taille du texte
mtext("a)", side = 3, line = 1, adj = 0, font = 2,cex=2)
#_________________________________________________________________________________________________


#### Préparation data Fst barrage B ####


Trutta_predictions$fst_ratio_trutta <- Trutta_predictions$fst_2_3_trutta / (1 - Trutta_predictions$fst_2_3_trutta)
Trutta_predictions$fst_ratio_trutta_sans_B <- Trutta_predictions$fst_2_3_trutta_sans_B / (1 - Trutta_predictions$fst_2_3_trutta_sans_B)

Gobio_predictions$fst_ratio_gobio <- Gobio_predictions$fst_2_3_gobio / (1 - Gobio_predictions$fst_2_3_gobio)
Gobio_predictions$fst_ratio_gobio_sans_B <- Gobio_predictions$fst_2_3_gobio_sans_B / (1 - Gobio_predictions$fst_2_3_gobio_sans_B)

Septimaniae_predictions$fst_ratio_septimaniae <- Septimaniae_predictions$fst_2_3_septimaniae / (1 - Septimaniae_predictions$fst_2_3_septimaniae)
Septimaniae_predictions$fst_ratio_septimaniae_sans_B <- Septimaniae_predictions$fst_2_3_septimaniae_sans_B / (1 - Septimaniae_predictions$fst_2_3_septimaniae_sans_B)


y_lim_fst_ratio <- c(
  min(c(min(Trutta_predictions$fst_ratio_trutta, na.rm = TRUE),
        min(Trutta_predictions$fst_ratio_trutta_sans_B, na.rm = TRUE),
        min(Gobio_predictions$fst_ratio_gobio, na.rm = TRUE),
        min(Gobio_predictions$fst_ratio_gobio_sans_B, na.rm = TRUE),
        min(Septimaniae_predictions$fst_ratio_septimaniae, na.rm = TRUE),
        min(Septimaniae_predictions$fst_ratio_septimaniae_sans_B, na.rm = TRUE))),
  max(c(max(Trutta_predictions$fst_ratio_trutta, na.rm = TRUE),
        max(Trutta_predictions$fst_ratio_trutta_sans_B, na.rm = TRUE),
        max(Gobio_predictions$fst_ratio_gobio, na.rm = TRUE),
        max(Gobio_predictions$fst_ratio_gobio_sans_B, na.rm = TRUE),
        max(Septimaniae_predictions$fst_ratio_septimaniae, na.rm = TRUE),
        max(Septimaniae_predictions$fst_ratio_septimaniae_sans_B, na.rm = TRUE)))
)

#### Graphique Fst barrage B ####
par(mgp = c(2.5, 1, 0)) 
plot(1, type = "n", 
     xlab = "Années", 
     ylab= "  ",
     xlim = c(min(Trutta_predictions$year), max(Trutta_predictions$year)), 
     ylim = y_lim_fst_ratio, 
     cex.lab = 1.3, 
     bty = "l")

transparent_darkgreen <- adjustcolor("darkgreen", alpha.f = 0.7)       # Trutta
transparent_darkgoldenrod <- adjustcolor("darkgoldenrod1", alpha.f = 0.7) # Gobio
transparent_blue <- adjustcolor("blue", alpha.f = 0.7)                 # Septimaniae

lines(Trutta_predictions$year, Trutta_predictions$fst_ratio_trutta, 
      col = transparent_darkgreen, lty = 2, lwd = 1.75)  
lines(Trutta_predictions$year, Trutta_predictions$fst_ratio_trutta_sans_B, 
      col = "darkgreen", lty = 1, lwd = 1.75)  

lines(Gobio_predictions$year, Gobio_predictions$fst_ratio_gobio, 
      col = transparent_darkgoldenrod, lty = 2, lwd = 1.75)  
lines(Gobio_predictions$year, Gobio_predictions$fst_ratio_gobio_sans_B, 
      col = "darkgoldenrod1", lty = 1, lwd = 1.75)  

lines(Septimaniae_predictions$year, Septimaniae_predictions$fst_ratio_septimaniae, 
      col = transparent_blue, lty = 2, lwd = 1.75) 
lines(Septimaniae_predictions$year, Septimaniae_predictions$fst_ratio_septimaniae_sans_B, 
      col = "blue", lty = 1, lwd = 1.75)  

rect(xleft = 0, xright = 2025, 
     ybottom = par("usr")[3], ytop = par("usr")[4], 
     col = rgb(0.3, 0.3, 0.3, 0.3), border = NA)
mtext("b)", side = 3, line = 1, adj = 0, font = 2, cex = 2)
#_________________________________________________________________________________________________

#### hugo
#Stats####
gobio_recent$esp<-"gobio"
septimaniae_recent$esp<-"septimaniae"
trutta_recent$esp<-"trutta"
names(gobio_recent)<-c("generation", "alive","Na","Fst_12","Fst_13","Fst_23","year","Fst_ratio","esp")
names(septimaniae_recent)<-names(gobio_recent)
names(trutta_recent)<-names(gobio_recent)

recent<-rbind(gobio_recent, septimaniae_recent)
recent<-rbind(recent, trutta_recent)

##tout les trucs en hastags c'est pour les conditions d'applications mais frère la je m'en blc un peu
#hist(recent$Na)
mod1=lm(Na~year*esp,data=recent)
anova(mod1)  
#hist(residuals(mod1))
summary(mod1) 
#summary(mod1)$r.squared #pour savoir varaition expliquée
#lillie.test(mod1$residuals)
#bptest(mod1)

mod2=lm(Fst_12~year*esp,data=recent)
anova(mod2)  
summary(mod2) 

mod3=lm(Fst_23~year*esp,data=recent)
anova(mod3)  
summary(mod3) 

###faut faire le futur

#F-index####