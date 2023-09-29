################################################################################

##########                ECOLOGICAL DYNAMIC REGIMES:                 ##########
####            IDENTIFICATION, CHARACTERIZATION, AND COMPARISON            ####

################################################################################


# Install and load ecoregime
install.packages("ecoregime")
library(ecoregime)

# Load other useful packages
library(vegan)
library(ecotraj)
library(dbscan)

## Inventory data --------------------------------------------------------------

# Species abundances
abun <- rbind(EDR_data$EDR1$abundance,
              EDR_data$EDR2$abundance,
              EDR_data$EDR3$abundance)

# ID trajectory
abun$ID <- paste0(abun$EDR, "_", abun$traj)

## State space -----------------------------------------------------------------

# State dissimilarities (e.g., Bray Curtis)
dStates <- vegdist(x = abun[, paste0("sp", 1:12)], 
                   method = "bray")

# State space (PCoA)
pcoa_states <- cmdscale(dStates, k = nrow(as.matrix(dStates)) - 1, add = T)
state_coord <- pcoa_states$points

# Plot the states
plot(x = state_coord[, 1], y = state_coord[, 2],
     xlab = "Axis 1", ylab = "Axis 2",
     main = "State space")

# Plot trajectories in the state space
trajectoryPCoA(d = dStates, 
               sites = abun$ID, 
               surveys = abun$state)
title("State space")

## Trajectory space ------------------------------------------------------------

# Trajectory dissimilarities
dTraj <- trajectoryDistances(d = dStates, 
                             sites = abun$ID, 
                             surveys = abun$state, 
                             distance.type = "DSPD")

# Trajectory space (PCoA)
pcoa_traj <- cmdscale(dTraj, k = nrow(as.matrix(dTraj)) - 1, add = T)
traj_coord <- pcoa_traj$points

# Plot the trajectory space
plot(x = traj_coord[, 1], y = traj_coord[, 2],
     xlab = "Axis 1", ylab = "Axis 2",
     main = "Trajectory space")

## Identify EDRs ---------------------------------------------------------------

# Clustering analysis (e.g., HDBSCAN)
EDR <- hdbscan(x = dTraj, minPts = 10)
EDR_cluster <- data.frame(ID = unique(abun$ID),
                          EDR_cluster = EDR$cluster)

# Plot the trajectory space identifying EDRs
plot(x = traj_coord[, 1], y = traj_coord[, 2],
     xlab = "Axis 1", ylab = "Axis 2",
     main = "Trajectory space", 
     col = EDR_cluster$EDR_cluster + 1)

# Plot trajectories in the state space identifying EDRs
trajectoryPCoA(d = dStates, 
               sites = abun$ID, 
               surveys = abun$state,
               traj.colors = EDR_cluster$EDR_cluster + 1)
title("State space")

## Representative trajectories -------------------------------------------------

# Select the EDR
ID_EDR <- which(abun$EDR == 1)

# Apply RETRA-EDR
RT <- retra_edr(d = as.matrix(dStates)[ID_EDR, ID_EDR],
                trajectories = abun[ID_EDR]$traj,
                states = abun[ID_EDR]$state, 
                minSegs = 5)

# Explore representative trajectories
RT$T1

# Summarize representative trajectories
summary(RT)

# Plot representative trajectories
plot(x = RT, d = as.matrix(dStates)[ID_EDR, ID_EDR], 
     trajectories = abun[ID_EDR]$traj,
     states = abun[ID_EDR]$state,
     select_RT = "T2",
     main = "Representative trajectories")

# Extract field data for representative trajectories
seg_components <- strsplit(gsub("\\]", "", gsub("\\[", "-", RT$T2$Segments)), "-")
RT_data <- do.call(rbind, lapply(seg_components, function(iseg){
  data.frame(traj = rep(iseg[[1]], 2),
             state = c(iseg[[2]], iseg[[3]]))
}))
RT_data <- merge(RT_data, abun[EDR == 1], all.x = T, sort = F)

# Plot changes in species abundances
plot(x = 1:nrow(RT_data), y = RT_data$sp1, type = "l",
     xlab = "RT state", ylab = "Species abundance",
     main = "Species abundances in RT")
for (i in 2:4) {
  lines(x = 1:nrow(RT_data), y = RT_data[, 3 + i], col = i)
}
legend("topleft", paste0("sp", 1:4), lty = 1, col = 1:4)

## Distribution of trajectories in the EDR -------------------------------------

# Dynamic dispersion
dDis_28 <- dDis(d = as.matrix(dStates)[ID_EDR, ID_EDR],
             d.type = "dStates",
             trajectories = abun[ID_EDR]$traj,
             states = abun[ID_EDR]$state, 
             reference = 28)
dDis_4 <- dDis(d = as.matrix(dStates)[ID_EDR, ID_EDR],
             d.type = "dStates",
             trajectories = abun[ID_EDR]$traj,
             states = abun[ID_EDR]$state, 
             reference = 4)

# Dynamic evenness
dEve <- dEve(d = as.matrix(dStates)[ID_EDR, ID_EDR],
           d.type = "dStates",
           trajectories = abun[ID_EDR]$traj,
           states = abun[ID_EDR]$state)

# Dynamic beta diversity
dBD <- dBD(d = as.matrix(dStates)[ID_EDR, ID_EDR],
           d.type = "dStates",
           trajectories = abun[ID_EDR]$traj,
           states = abun[ID_EDR]$state)

## Compare EDRs ----------------------------------------------------------------

# Identify EDRs in the abundance matrix
abun <- merge(abun, EDR_cluster, by = "ID", all.x = T)

# EDR dissimilarity
dDR <- dist_edr(d = dStates, d.type = "dStates", 
                trajectories = abun$ID,
                states = abun$state, 
                edr = abun$EDR_cluster)

################################################################################