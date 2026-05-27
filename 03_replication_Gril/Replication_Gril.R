##### Replication of Gril (2024) ###############################################
# Martin Möhler, 2026
# based on:
# L. Gril, L. Steinkemper, M. Gross, U. Rendtel: Kernel Heaping - Kernel Density
#   Estimation from regional aggregates via measurement error model (2024)
#   https://journal.r-project.org/articles/RJ-2024-026/
############################################################################## #

library(Kernelheaping)  # apply kernel heaping algorithm
library(ptable)         # compute CKM transition probabilities
library(sf)             # working with geo data
library(ggplot2)        # plotting
library(patchwork)      # plotting
library(dplyr)          # data wrangling
library(tidyr)          # data wrangling


# adapt file path to local setting:
fpath <- "~/ckm_regio_2026/03_replication_Gril/data/"
spath <- "~/ckm_regio_2026/01_replication_Prenzel/data/"


# ----- (0) Helper functions -----

get_perturbation_value <- function(x_ck, pt) {
  
  if(is.na(x_ck[1])) { 
    NA 
  } else {
    max_x <- max(pt$i)
    pt$v[pt$i == min(max_x, x_ck[1]) & data.table::between(x_ck[2], pt$p_int_lb, pt$p_int_ub)]
  }
}


# ----- (1) Prepare Data -----

# ShapeFile
Berlin <- sf::st_read(file.path(fpath, "RBS_OD_LOR_2015_12/RBS_OD_LOR_2015_12.shp")) %>%
  dplyr::select(PLR, PLRNAME)
Rheinl <- sf::st_read(file.path(spath, "vg250_01-01.utm32s.shape.ebenen/vg250_ebenen_0101/VG250_GEM.shp")) %>%
  dplyr::filter(GF == 4 & SN_L == "07") %>% dplyr::select(ARS, GEN)

# Inhabitant data - BE
HK <- c("HK_Arab", "HK_EheSU", "HK_Polen")
EWRMigra <- read.csv2(file.path(fpath, "EWRMIGRA201512H_Matrix.csv"),
                      colClasses = list(RAUMID = "character")) %>%
  dplyr::select(RAUMID, all_of(HK))

# Inhabitant data - RP
GL <- c("Serbia", "Ukraine", "Russian Federation")
EWRGebLa <- read.csv2(file.path(spath, "1000A-1016_en_flat.csv")) %>%
  dplyr::filter(substr(X1_variable_attribute_code, 1, 2) == "07" &
                X2_variable_attribute_label %in% GL) %>%
  dplyr::rename(ARS = X1_variable_attribute_code, GL = X2_variable_attribute_label) %>%
  dplyr::select(ARS, GL, value)
EWRGebLa$value[EWRGebLa$value == "-"] <- 0
EWRGebLa$value <- as.numeric(EWRGebLa$value)

Berlin <- merge(Berlin, EWRMigra, by.x = "PLR", by.y = "RAUMID")
Rheinl <- merge(Rheinl, EWRGebLa, by.x = "ARS", by.y = "ARS")

nPLR <- nrow(Berlin)
nGEM <- nrow(Rheinl) / length(GL)

rm(EWRMigra, EWRGebLa)


## Prepare perturbation tables & set up data sets to be perturbed

# make cell keys
set.seed(20260224)
for(hk in HK) { Berlin[, paste0("ck_", hk)] <- runif(nPLR) }
Rheinl$ck <- runif(nGEM * length(GL))

# prepare p-tables
ptabs <- vector("list", 3)

D   <- c(5, 10, 20)
V   <- c(2.5, 5, 10)
js  <- c(2, 2, 2)
lbl <- c("original", "CKM - mittel", "CKM - stark", "CKM - sehr stark")

for(p in 1:3) { 
  
  # make p-tables
  ptabs[[p]] <- ptable::create_cnt_ptable(D = D[p], V = V[p], js = js[p])
  ptabs[[p]]@pTable$variant <- factor(p, levels = 1:3, labels = lbl[2:4])
  
  ## apply perturbation
  
  # BE
  for(hk in HK) {
    x <- st_drop_geometry(Berlin[, c(hk, paste0("ck_", hk))])
    y <- x[, 1] + apply(x[, 1:2], 1, get_perturbation_value, pt = ptabs[[p]]@pTable)
    Berlin[, paste0(hk, "_pt", p)] <- y
    Berlin[, paste0(hk, "_pt", p, "_diff")] <- y - x[, 1]
  }
  # RP
  x <- cbind(Rheinl$value, Rheinl$ck)
  y <- x[, 1] + apply(x[, 1:2], 1, get_perturbation_value, pt = ptabs[[p]]@pTable)
  Rheinl[, paste0("value_pt", p)] <- y
  Rheinl[, paste0("value_pt", p, "_diff")] <- y - x[, 1]
}


# ----- (2) Estimate Densities -----

# centroids for initital estimate
centr_BE <- st_coordinates(st_centroid(Berlin))
centr_RP <- st_coordinates(st_centroid(Rheinl[Rheinl$GL == GL[1], ]))

# set common parameters

burn <- 5
samp <- 10
grsz <- 325
runs <- 10

# quantile thresholds for Jaccard similarity
visquantile <- c(.95, .75)

# What part of the evaluation rectangle lies within the shape?
bbox_Berlin <- as_Spatial(Berlin)@bbox
bbox_Rheinl <- as_Spatial(Rheinl)@bbox

BE_outl <- st_sf(st_union(Berlin))
RP_outl <- st_sf(st_union(Rheinl))

gridx_BE <- seq(bbox_Berlin[1, 1], bbox_Berlin[1, 2], length = grsz)
gridy_BE <- seq(bbox_Berlin[2, 1], bbox_Berlin[2, 2], length = grsz)

gridx_RP <- seq(bbox_Rheinl[1, 1], bbox_Rheinl[1, 2], length = grsz)
gridy_RP <- seq(bbox_Rheinl[2, 1], bbox_Rheinl[2, 2], length = grsz)

grpts_BE <- expand.grid(x = gridx_BE, y = gridy_BE) %>%
  st_as_sf(coords = c("x", "y"), remove = FALSE, crs = st_crs(Berlin))
grpts_RP <- expand.grid(x = gridx_RP, y = gridy_RP) %>%
  st_as_sf(coords = c("x", "y"), remove = FALSE, crs = st_crs(Rheinl))

grpts_BE$in_shape <- as.logical(st_intersects(grpts_BE, BE_outl, sparse = FALSE))
grpts_RP$in_shape <- as.logical(st_intersects(grpts_RP, RP_outl, sparse = FALSE))

grpts_BE <- st_drop_geometry(grpts_BE)
grpts_RP <- st_drop_geometry(grpts_RP)

# results data.frames
results_BE <- expand.grid(x     = gridx_BE,
                          y     = gridy_BE, 
                          HK    = HK,
                          ptab  = 0:3,
                          dens  = NA,
                          hdr95 = NA,
                          hdr75 = NA)
results_RP <- expand.grid(x     = gridx_RP,
                          y     = gridy_RP, 
                          GL    = GL,
                          ptab  = 0:3,
                          dens  = NA,
                          hdr95 = NA,
                          hdr75 = NA)


## function for computing HDRs

get_hdr_kh <- function(x, runs, Q, g, grpts, centr, shp) {
  
  nq  <- length(Q)
  hdr <- matrix(nrow = g^2, ncol = nq)
  dat <- as.data.frame(cbind(centr, x))
  
  # run kernelheaping several times
  est <- dshapebivr(data      = dat,
                    burnin    = burn, 
                    samples   = samp,
                    adaptive  = FALSE,
                    shapefile = as_Spatial(shp),
                    gridsize  = g, 
                    boundary  = TRUE,
                    numChains = runs)
  
  # extract density estimate 
  dens <- as.vector(est$Mestimates$estimate)
  
  # subset to evaluation points above quantile
  for(j in 1:nq) {
    qj <- quantile(dens[grpts$in_shape], Q[j])
    hdr[, j] <- ifelse(dens > qj, 1, 0)
  }
  
  cbind(dens, hdr)
}

## Run kernel heaping and calculate HDRs

for(p in 0:3) {
  
  # BE
  for(hk in HK) {
    
    val <- ifelse(p == 0, hk, paste0(hk, "_pt", p))
    results_BE[results_BE$HK == hk & results_BE$ptab == p, 
               c("dens", "hdr95", "hdr75")] <-
      get_hdr_kh(st_drop_geometry(Berlin[, val]), 
                 runs  = runs, 
                 Q     = visquantile, 
                 g     = grsz, 
                 grpts = grpts_BE, 
                 centr = centr_BE, 
                 shp   = Berlin)
  }
  
  # RP
  for(gl in GL) {
    
    val <- ifelse(p == 0, "value", paste0("value_pt", p))
    results_RP[results_RP$GL == gl & results_RP$ptab == p, 
               c("dens", "hdr95", "hdr75")] <-
      get_hdr_kh(st_drop_geometry(Rheinl[Rheinl$GL == gl, val]), 
                 runs  = runs, 
                 Q     = visquantile, 
                 g     = grsz, 
                 grpts = grpts_RP, 
                 centr = centr_RP, 
                 shp   = Rheinl[Rheinl$GL == gl, ])
  }
}


# derive further metrics

results_BE$variant    <- factor(results_BE$ptab, levels = 0:3, labels = lbl)
results_RP$variant    <- factor(results_RP$ptab, levels = 0:3, labels = lbl)

results_BE$dens_diff  <- results_BE$dens - results_BE$dens[results_BE$ptab == 0]
results_RP$dens_diff  <- results_RP$dens - results_RP$dens[results_RP$ptab == 0]

results_BE$dens_diff[abs(results_BE$dens_diff) == 0] <- NA
results_RP$dens_diff[abs(results_RP$dens_diff) == 0] <- NA

hs_type <- c("FN", "WAHR", "FP") # pixel-wise hot spot class

results_BE$hdr95_diff <- results_BE$hdr95 - results_BE$hdr95[results_BE$ptab == 0]
results_RP$hdr95_diff <- results_RP$hdr95 - results_RP$hdr95[results_RP$ptab == 0]
results_BE$hdr75_diff <- results_BE$hdr75 - results_BE$hdr75[results_BE$ptab == 0]
results_RP$hdr75_diff <- results_RP$hdr75 - results_RP$hdr75[results_RP$ptab == 0]

results_BE$hdr95_diff <- factor(results_BE$hdr95_diff, labels = hs_type)
results_RP$hdr95_diff <- factor(results_RP$hdr95_diff, labels = hs_type)
results_BE$hdr75_diff <- factor(results_BE$hdr75_diff, labels = hs_type)
results_RP$hdr75_diff <- factor(results_RP$hdr75_diff, labels = hs_type)

results_BE$hdr95_any <- as.integer(results_BE$hdr95 | results_BE$hdr95[results_BE$ptab == 0])
results_RP$hdr95_any <- as.integer(results_RP$hdr95 | results_RP$hdr95[results_RP$ptab == 0])
results_BE$hdr75_any <- as.integer(results_BE$hdr75 | results_BE$hdr75[results_BE$ptab == 0])
results_RP$hdr75_any <- as.integer(results_RP$hdr75 | results_RP$hdr75[results_RP$ptab == 0])

#save(results_BE, results_RP, file = paste0(fpath, "results_KH.RData"))


# ----- (3) Assess HDR similarity ----- 

#load(paste0(fpath, "results_KH.RData"))

# mean population in areas

sort(c(mean_arab  = mean(Berlin$HK_Arab), 
       mean_ehesu = mean(Berlin$HK_EheSU), 
       mean_polen = mean(Berlin$HK_Polen)))

sort(c(mean_serb = mean(Rheinl$value[Rheinl$GL == "Serbia"]),
       mean_ukra = mean(Rheinl$value[Rheinl$GL == "Ukraine"]),
       mean_russ = mean(Rheinl$value[Rheinl$GL == "Russian Federation"])))


## calculate Jaccard coefficients and MISE / RMISE

mise_scale_BE <- (bbox_Berlin[1, 2] - bbox_Berlin[1, 1]) * (bbox_Berlin[2, 2] - bbox_Berlin[2, 1])
mise_scale_RP <- (bbox_Rheinl[1, 2] - bbox_Rheinl[1, 1]) * (bbox_Rheinl[2, 2] - bbox_Rheinl[2, 1])

get_jaccard <- function(x, y) { sum(x & y) / sum(x | y) }
get_mise    <- function(x, y, scale) { mean((y - x)^2) * scale }

results_jacc <- expand.grid(var  = c(HK, GL),
                            hdr  = visquantile * 100,
                            ptab = 1:3,
                            jacc = NA,
                            mise = NA)

for(pt in 1:3) {
  for(hdr in visquantile * 100) {
    
    for(hk in HK) {
      
      results_jacc$jacc[results_jacc$var  == hk & 
                        results_jacc$ptab == pt & 
                        results_jacc$hdr  == hdr] <- 
        round(
          get_jaccard(x = results_BE[results_BE$HK == hk & results_BE$ptab == 0,  paste0("hdr", hdr)],
                      y = results_BE[results_BE$HK == hk & results_BE$ptab == pt, paste0("hdr", hdr)]) * 100, 1)
      
      results_jacc$mise[results_jacc$var  == hk & 
                        results_jacc$ptab == pt & 
                        results_jacc$hdr  == hdr] <- 
        get_mise(x = results_BE[results_BE$HK == hk & results_BE$ptab == 0,  "dens"],
                 y = results_BE[results_BE$HK == hk & results_BE$ptab == pt, "dens"],
                 scale = mise_scale_BE)
    }
    
    for(gl in GL) {
      
      results_jacc$jacc[results_jacc$var  == gl & 
                        results_jacc$ptab == pt & 
                        results_jacc$hdr  == hdr] <- 
        round(
          get_jaccard(x = results_RP[results_RP$GL == gl & results_RP$ptab == 0,  paste0("hdr", hdr)],
                      y = results_RP[results_RP$GL == gl & results_RP$ptab == pt, paste0("hdr", hdr)]) * 100, 1)
      
      results_jacc$mise[results_jacc$var  == gl & 
                        results_jacc$ptab == pt & 
                        results_jacc$hdr  == hdr] <- 
        get_mise(x = results_RP[results_RP$GL == gl & results_RP$ptab == 0,  "dens"],
                 y = results_RP[results_RP$GL == gl & results_RP$ptab == pt, "dens"],
                 scale = mise_scale_BE)
    }
  }
}

results_jacc$rmise <- round(sqrt(results_jacc$mise), 7)


# ----- (4) Plot -----

map_BE <- Berlin %>% pivot_longer(cols = ends_with("_diff")) %>%
  dplyr::select(-starts_with("ck_"), -starts_with("HK_"), -name)
map_BE$HK <- rep(HK, nPLR * 3)

map_RP <- Rheinl[, names(Rheinl) != "value"] %>% 
  pivot_longer(cols = ends_with("_diff")) %>%
  dplyr::select(-starts_with("value_"), -name)

map_BE$pt <- factor(rep(rep(1:3, each = 3), nPLR), labels = lbl[2:4])
map_RP$pt <- factor(rep(1:3, 3 * nGEM), labels = lbl[2:4])


# absolute changes in PLR
m1a <- ggplot(map_BE[map_BE$HK == HK[1], ]) +
  geom_sf(aes(fill = value), color = alpha("grey", .3)) +
  scale_fill_gradient2(low = "orange", high = "blue", 
                       name = "Diff. Anzahl") +
  geom_sf(data = BE_outl, fill = NA) +
  facet_grid(~ pt) +
  theme_void() +
  theme(legend.key.width= unit(.25, "cm")) +
  ggtitle("a)")

m1b <- ggplot(map_RP[map_RP$GL == GL[1], ]) +
  geom_sf(aes(fill = value), color = alpha("grey", .3)) +
  scale_fill_gradient2(low = "orange", high = "blue", 
                       name = "Diff. Anzahl") +
  geom_sf(data = RP_outl, fill = NA) +
  facet_grid(~ pt) +
  theme_void() +
  theme(legend.key.width= unit(.25, "cm")) +
  ggtitle("a)")

# pixel-wise density change
m2a <- ggplot() +
  geom_tile(data = results_BE[results_BE$HK == HK[1] & 
                             results_BE$ptab != 0 & 
                             !is.na(results_BE$dens_diff), ], 
              aes(x, y, fill = dens_diff)) +
  scale_fill_gradient2(low = "orange", high = "blue", 
                       name = "Diff. Dichte") +
  geom_sf(data = BE_outl, fill = NA) +
  facet_wrap(~variant) +
  theme_void() +
  theme(legend.key.width = unit(.25, "cm")) +
  ggtitle("b)")

m2b <- ggplot() +
  geom_tile(data = results_RP[results_RP$GL == GL[3] & 
                              results_RP$ptab != 0 & 
                              !is.na(results_RP$dens_diff), ], 
            aes(x, y, fill = dens_diff)) +
  scale_fill_gradient2(low = "orange", high = "blue", 
                       name = "Diff. Dichte") +
  geom_sf(data = RP_outl, fill = NA) +
  facet_wrap(~variant) +
  theme_void() +
  theme(legend.key.width = unit(.25, "cm")) +
  ggtitle("b)")

# change in HDRs

closeup_BE <- rbind(x = c(382000, 398000), y = c(5815000, 5825000))
closeup_RP <- rbind(x = c(400000, 460000), y = c(5480000, 5540000))

m3a <- ggplot() +
  geom_sf(data = Berlin, fill = NA, color = alpha("lightgrey", .3)) +
  geom_sf(data = BE_outl, fill = NA) +
  geom_tile(data = results_BE[results_BE$HK == HK[1] & 
                              results_BE$ptab != 0 & 
                              results_BE$hdr95 == 1, ], 
              aes(x, y), fill = "darkgrey") +
  #annotate(geom = "rect", fill = NA, color = "darkred",
  #         xmin = closeup_BE[1, 1], xmax = closeup_BE[1, 2],
  #         ymin = closeup_BE[2, 1], ymax = closeup_BE[2, 2]) +
  facet_wrap(~variant) +
  theme_void() +
  theme(legend.key.size = unit(.25, "cm")) +
  ggtitle("c)")

m3b <- ggplot() +
  geom_sf(data = Rheinl, fill = NA, color = alpha("lightgrey", .3)) +
  geom_sf(data = RP_outl, fill = NA) +
  geom_tile(data = results_RP[results_RP$GL == GL[1] & 
                              results_RP$ptab != 0 & 
                              results_RP$hdr75 == 1, ], 
            aes(x, y), fill = "darkgrey") +
  #annotate(geom = "rect", fill = NA, color = "darkred",
  #         xmin = closeup_RP[1, 1], xmax = closeup_RP[1, 2],
  #         ymin = closeup_RP[2, 1], ymax = closeup_RP[2, 2]) +
  facet_wrap(~variant) +
  theme_void() +
  theme(legend.key.size = unit(.25, "cm")) +
  ggtitle("c)")

m4a <- ggplot() +
  geom_sf(data = Berlin, fill = NA, color = alpha("lightgrey", .3)) +
  geom_sf(data = BE_outl, fill = NA) +
  coord_sf(xlim = closeup_BE[1, ], 
           ylim = closeup_BE[2, ]) +
  geom_tile(data = results_BE[results_BE$HK == HK[1] & 
                              results_BE$ptab != 0 & 
                              results_BE$hdr95_any != 0, ], 
              aes(x, y, fill = as.factor(hdr95_diff), 
                  alpha = as.factor(hdr95_diff))) +
  scale_fill_manual(values = c("orange", "darkgrey", "blue"), 
                    name = "HDR") +
  scale_alpha_manual(values = c(1, 0.4, 1), 
                     name = "HDR") +
  facet_wrap(~variant) +
  theme_void() +
  theme(legend.key.size = unit(.25, "cm")) +
  ggtitle("d)")

m4b <- ggplot() +
  geom_sf(data = Rheinl, fill = NA, color = alpha("lightgrey", .3)) +
  geom_sf(data = RP_outl, fill = NA) +
  coord_sf(xlim = closeup_RP[1, ], 
           ylim = closeup_RP[2, ]) +
  geom_tile(data = results_RP[results_RP$GL == GL[1] & 
                              results_RP$ptab != 0 & 
                              results_RP$hdr75_any != 0, ], 
            aes(x, y, fill = as.factor(hdr75_diff), 
                alpha = as.factor(hdr75_diff))) +
  scale_fill_manual(values = c("orange", "darkgrey", "blue"), 
                    name = "HDR") +
  scale_alpha_manual(values = c(1, 0.4, 1), 
                     name = "HDR") +
  facet_wrap(~variant) +
  theme_void() +
  theme(legend.key.size = unit(.25, "cm")) +
  ggtitle("d)")


m1a / m2a / m3a / m4a
#ggsave("maps_BE.png", width = 2100, height = 3200, units = "px")
m1b / m2b / m3b / m4b
#ggsave("maps_RP.png", width = 2100, height = 3200, units = "px")

