##### Grid data accessibility zones ############################################
# Martin Möhler, 2026
# earlier version: 
# https://cros.ec.europa.eu/system/files/2025-02/Book%20of%20Abstracts%20-%20NTTS%202025%20-%2027.02.25.pdf
# needs preparation step:
# 00_Simulate_data.R
############################################################################## #

library(sf)         # working with geo data
library(dplyr)      # data wrangling
library(ptable)     # compute CKM transition probabilities
library(ggplot2)    # plotting
library(latex2exp)  # plotting
library(patchwork)  # plotting

# adapt file path to local setting:
fpath <- "~/ckm_regio_2026/02_accessibility_zones/data/"
spath <- "~/ckm_regio_2026/01_replication_Prenzel/data/vg250_01-01.utm32s.shape.ebenen/vg250_ebenen_0101/"


# ----- (0) Helper functions -----

get_perturbation_value <- function(x_ck, pt) {
  
  if(is.na(x_ck[1])) { 
    NA 
  } else {
    max_x <- max(pt$i)
    pt$v[pt$i == min(max_x, x_ck[1]) & data.table::between(x_ck[2], pt$p_int_lb, pt$p_int_ub)]
  }
}


# ----- (1) Prepare data -----

# read in routing results
pharm_ped_gzdb <- read.csv(paste0(fpath, "apotheken_20230801_131025_ped_100_2024-08-29/apotheken_20230801_131025_pedestrian_100.csv"))
pharm_ped_gzdb <- pharm_ped_gzdb %>% select(id, dist_min)

# read in shapes of LAU
shp <- read_sf(paste0(spath, "VG250_GEM.shp")) %>% filter(GF == 4)
shp <- shp %>%  mutate(area = as.numeric(st_area(shp))) %>%
  select(ARS, GEN, area) %>% st_transform(crs = "EPSG:3035")

# read in population data
load(paste0(fpath, "pop_data.RData"))

# draw record keys for population
set.seed(820240829)
pop$rk <- runif(nrow(pop), 0, 1)

# aggregate record-level pop. data by grid cell
pop <- pop %>% group_by(Gitter_ID_100m) %>%
  summarise(npop = n(), ck = sum(rk) %% 1,
            E = unique(x_mp_100m), N = unique(y_mp_100m),
            GEN = unique(GEN), ARS = unique(ARS))
# unify ID for merging
pop$id <- paste0(substr(pop$Gitter_ID_100m, 11, 20), 
                 substr(pop$Gitter_ID_100m, 23, 28))

# subset to relevant LAU
bl  <- "05"
shp <- shp %>% filter(substr(ARS, 1, 2) == bl)

# make ptable
ptab <- create_cnt_ptable(D = 6, V = 2.2, js = 2)

# apply CKM
pop$npop_mod <- apply(pop[, c("npop", "ck")], 1, get_perturbation_value, pt = ptab@pTable)
pop$npop_ckm <- pop$npop + pop$npop_mod


## Merge cell-level pop. counts and cell-level routed distance measures

# merge pop. counts (original & perturbed) to routing results
pharm_ped <- pharm_ped_gzdb[pharm_ped_gzdb$id %in% pop$id, ] %>%
  merge(pop, by.x = "id", by.y = "id")
rm(pharm_ped_gzdb, pop)

# distance categories
dist_labs <- c("bis 200m", "200m - 500m", "500m - 1km", "1km - 2km", "2km - 4km", "mehr als 4km")
pharm_ped$dist_min <- factor(pharm_ped$dist_min, levels = 1:6, labels = dist_labs)

## plot

lau <- "Solingen"

map1 <- ggplot() +
  geom_tile(data = pharm_ped[pharm_ped$GEN %in% lau, ], 
            aes(E, N, fill = npop)) +
  geom_sf(data = shp[shp$GEN %in% lau, ], 
          fill = NA) +
  scale_fill_viridis_c(trans = "log10", 
                       name = "Bev.\n65+") +
  theme_void() +
  theme(legend.key.width = unit(.25, "cm"))

map2 <- ggplot() +
  geom_tile(data = pharm_ped[pharm_ped$GEN %in% lau, ], 
            aes(E, N), fill = "blue",
            show.legend = FALSE) +
  geom_sf(data = shp[shp$GEN %in% lau, ], 
          fill = NA) +
  facet_wrap(~dist_min) +
  theme_void()

map1 + map2
#ggsave("map_distmin.png", width = 2000, height = 700, units = "px")


# ------ (2) Simulation -----

gen  <- unique(pharm_ped$GEN) # LAU for which simulated data is available
nsim <- length(gen)           # no. of LAUs
nCI  <- 1000                  # no. of iterations for Monte-Carlo-CIs

# data.frame to collect results
res_ped <- expand.grid(dist_min   = dist_labs, 
                       GEN        = gen,
                       ncell      = 0, 
                       npop       = 0, 
                       npop_ckm   = 0,
                       npop_uq99  = 0,
                       npop_lq99  = 0,
                       npop_uq95  = 0,
                       npop_lq95  = 0,
                       npop_med   = 0,
                       var_npop   = 0)

# intermediate result table for conf. intervals
CI_ped <- expand.grid(dist_min = levels(pharm_ped$dist_min), re = 1:nCI, npop_ci = 0)

set.seed(20251125)
for(i in 1:nsim) {
  
  # subset to LAU
  ped <- pharm_ped[pharm_ped$GEN == gen[i], ]
  
  # count pop. (before and after CKM) in focus area by distance category
  nc_ped     <- table(ped$dist_min)
  np_ped     <- xtabs(npop     ~ dist_min, data = ped)
  np_ckm_ped <- xtabs(npop_ckm ~ dist_min, data = ped)
  # store single-run results
  res_ped$npop[res_ped$GEN == gen[i]]       <- np_ped
  res_ped$npop_ckm[res_ped$GEN == gen[i]]   <- np_ckm_ped
  res_ped$ncell[res_ped$GEN == gen[i]]      <- nc_ped
  
  # get MC-based CIs
  for(j in 1:nCI) {
    # nCI times ...
    # ... draw new set of CKs, ...
    ped$ck_ci <- runif(nrow(ped))
    # ... get alternative perturbation, ...
    ped$npop_mod_ci <- apply(ped[, c("npop", "ck_ci")], 1, get_perturbation_value, pt = ptab@pTable)
    # ... calculate alternative aggregates, ...
    ped$npop_ci <- ped$npop + ped$npop_mod_ci
    # ... store for CI calculation.
    CI_ped$npop_ci[CI_ped$re == j] <- xtabs(npop_ci ~ dist_min, data = ped)
  }
  for(j in dist_labs) {
    ped_ci <- CI_ped$npop_ci[CI_ped$dist_min == j]
    # CIs and empirical variance for frequencies
    res_ped$npop_med[res_ped$GEN  == gen[i] & res_ped$dist_min == j] <- median(ped_ci,          na.rm = TRUE)
    res_ped$npop_lq95[res_ped$GEN == gen[i] & res_ped$dist_min == j] <- quantile(ped_ci, 0.025, na.rm = TRUE)
    res_ped$npop_uq95[res_ped$GEN == gen[i] & res_ped$dist_min == j] <- quantile(ped_ci, 0.975, na.rm = TRUE)
    res_ped$npop_lq99[res_ped$GEN == gen[i] & res_ped$dist_min == j] <- quantile(ped_ci, 0.005, na.rm = TRUE)
    res_ped$npop_uq99[res_ped$GEN == gen[i] & res_ped$dist_min == j] <- quantile(ped_ci, 0.995, na.rm = TRUE)
    res_ped$var_npop[res_ped$GEN  == gen[i] & res_ped$dist_min == j] <- var(ped_ci,             na.rm = TRUE)
  }
  
  print(paste("Done:", i, "To Do:", nsim - i))
}

# difference between pre-CKM and post-CKM
res_ped$npop_mod   <- res_ped$npop_ckm - res_ped$npop
# theoretical variance of composite
res_ped$vtheo_npop <- res_ped$ncell * ptab@pParams@V
# pop. in composite (true and estimated)
res_ped$avg_pop     <- ifelse(res_ped$ncell == 0, 0, res_ped$npop     / res_ped$ncell)
res_ped$avg_pop_ckm <- ifelse(res_ped$ncell == 0, 0, res_ped$npop_ckm / res_ped$ncell)

#save(res_ped, file = "02_accessibility_zones/data/res_ped.R")


# ----- (3) Analyse results -----

#load("02_accessibility_zones/data/res_ped.R")

## Mean absolute errors
res_ped %>% group_by(dist_min) %>%
  summarise(npop_mod = mean(abs(npop_mod)))

# average cell count per distance category
avg_cc <- res_ped[res_ped$ncell > 10, ] %>%
  group_by(dist_min) %>% summarise(ncell    = sum(ncell), 
                                   npop     = sum(npop), 
                                   npop_ckm = sum(npop_ckm))
avg_cc$avg_pop     <- round(avg_cc$npop     / avg_cc$ncell, 1)
avg_cc$avg_pop_ckm <- round(avg_cc$npop_ckm / avg_cc$ncell, 1)

# theoretical confidence bands
CI_theo <- expand.grid(dist_min = dist_labs, 
                       ncell = 10:max(res_ped$ncell)) %>%
  merge(select(avg_cc, dist_min, avg_pop)) %>%
  arrange(dist_min, ncell) %>%
  mutate(npop_calc  = ncell * avg_pop, 
         vtheo_npop = npop_calc * ptab@pParams@V)

CI_theo$npop_lq95 <- qnorm(.025, mean = CI_theo$npop_calc, sd = sqrt(CI_theo$vtheo_npop))
CI_theo$npop_uq95 <- qnorm(.975, mean = CI_theo$npop_calc, sd = sqrt(CI_theo$vtheo_npop))
CI_theo$npop_lq99 <- qnorm(.005, mean = CI_theo$npop_calc, sd = sqrt(CI_theo$vtheo_npop))
CI_theo$npop_uq99 <- qnorm(.995, mean = CI_theo$npop_calc, sd = sqrt(CI_theo$vtheo_npop))


## graphs

c1 <- ggplot(res_ped[res_ped$npop_ckm > 0 & res_ped$npop_ckm < 4e+4, ], aes(npop/ 1000, npop_ckm/1000)) +
  geom_abline(slope = 1, lty = "dashed") +
  geom_point(alpha = .6, color = "grey20") +
  xlab("Bev. 65+ (aus original) in Tausend") +
  ylab("Bev. 65+ (aus überlagert) in Tausend") +
  theme_bw() +
  theme(axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 10))) +
  ggtitle("a)")

c2 <- ggplot(res_ped[res_ped$npop_ckm < 500, ],  aes(npop)) +
  geom_abline(slope = 1, lty = "dashed") +
  geom_point(aes(y = npop_ckm), show.legend = FALSE, color = "grey20", alpha = .6) +
  facet_wrap(~dist_min, ncol = 3) +
  xlab("Bev. 65+ (aus original)") +
  ylab("Bev. 65+ (aus überlagert)") +
  scale_x_continuous(breaks = (0:5)*100) +
  scale_y_continuous(breaks = (0:5)*100) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 10))) +
  ggtitle("b)")

c3 <- ggplot(res_ped[res_ped$ncell > 10, ]) +
  geom_errorbar(aes(ncell, ymin = (npop_lq95 - npop)/npop * 100,  
                           ymax = (npop_uq95 - npop)/npop * 100), 
                color = "grey20", alpha = .4, show.legend = FALSE) +
  geom_line(data = CI_theo, aes(ncell, (npop_lq95 - npop_calc)/npop_calc * 100),
            lty = "dashed", color = "blue") +
  geom_line(data = CI_theo, aes(ncell, (npop_uq95 - npop_calc)/npop_calc * 100),
            lty = "dashed", color = "blue") +
  geom_line(data = CI_theo, aes(ncell, (npop_lq99 - npop_calc)/npop_calc * 100),
            lty = "dashed", color = "red") +
  geom_line(data = CI_theo, aes(ncell, (npop_uq99 - npop_calc)/npop_calc * 100),
            lty = "dashed", color = "red") +
  geom_label(data = avg_cc, aes(1000, 40, label = paste("N/K =", format(avg_pop_ckm, decimal.mark = ","))), 
             color = "black", show.legend = FALSE) +
  scale_x_log10() +
  scale_y_continuous(breaks = c(-60, -40, -20, 0, 20, 40, 60)) +
  facet_wrap(~dist_min, ncol = 3) +
  theme_bw() +
  xlab("Anz. Rasterzellen (K)") +
  ylab("relativer Fehler (%)") +
  theme(panel.grid.minor = element_blank(),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 10))) +
  ggtitle("c)")

layout <- "AABBB\nAABBB\nCCCCC\nCCCCC"
c1 + c2 + c3 + plot_layout(design = layout)
#ggsave("dist_counts.png", width = 2800, height = 2600, units = "px")


## 'Bach-style' plots

cvals <- c(0.1, 0.05, 0.025) # accuracy targets

# bounds for plotting
ncell_ub   <- 1500 
avg_pop_ub <- 100  

# range of noise variances considered
sig_rng <- expand.grid(ncell = 0:ncell_ub, avg_pop_ckm = seq(0.1, avg_pop_ub, 0.1), c = cvals)
sig_rng$sig2 <- sig_rng$ncell * sig_rng$avg_pop_ckm^2 * sig_rng$c^2
sig_rng$c_lvl <- factor(sig_rng$c, levels = cvals, 
                        labels = paste0("c = ", format(cvals*100, decimal.mark = ","), "%"))

# true noise variance used
sig_real <- expand.grid(ncell = seq(0.1, ncell_ub, 0.1), c = cvals)
sig_real$avg_pop_ckm <- sqrt(ptab@pParams@V / (sig_real$ncell * sig_real$c^2))
sig_real <- sig_real[sig_real$avg_pop_ckm <= avg_pop_ub, ]
sig_real$c_lvl <- factor(sig_real$c, levels = cvals, 
                         labels = paste0("c = ", format(cvals*100, decimal.mark = ","), "%"))

# levels of noise variance for robustness
breaks <- c(0, 1:5, 10, max(sig_rng$sig2))
labs   <- c(paste(breaks[c(-7, -8)], "–", breaks[c(-1, -8)]), "10+")

# classify results into robust vs. non-robust
comp_pts <- rbind(res_ped, res_ped, res_ped)
comp_pts$c <- rep(cvals, each = nrow(res_ped))
comp_pts$c_lvl <- factor(comp_pts$c, levels = cvals, 
                         labels = paste0("c = ", format(cvals*100, decimal.mark = ","), "%"))
comp_pts$maxvar <- comp_pts$ncell * comp_pts$avg_pop_ckm^2 * comp_pts$c^2
comp_pts$robust <- comp_pts$maxvar >= ptab@pParams@V

ggplot(sig_rng[sig_rng$ncell >= 1, ], aes(ncell, avg_pop_ckm)) +
  geom_contour_filled(aes(z = sig2), breaks = breaks, alpha = .5) +
  geom_line(data = sig_real, lty = "dashed") +
  geom_point(data = comp_pts[comp_pts$npop_ckm > 0, ], 
             color = "grey20",
             size = .3,
             alpha = .4) +
  scale_x_log10(limits = c(1, ncell_ub)) +
  scale_y_log10(limits = c(0.5, avg_pop_ub)) +
  scale_fill_viridis_d(option = "plasma", 
                       name = "level\nVar",
                       labels = labs) +
  theme_bw() +
  ylab(TeX("\\bar{N}")) +
  xlab("K") +
  facet_wrap(~ c_lvl)
#ggsave("bachplot.png", width = 1900, height = 700, units = "px")

# aggregate results on robustness
table(comp_pts$robust, comp_pts$c_lvl)
table(comp_pts$robust[comp_pts$dist_min %in% dist_labs[5:6]],
      comp_pts$c_lvl[comp_pts$dist_min  %in% dist_labs[5:6]])
table(comp_pts$robust[comp_pts$dist_min %in% dist_labs[1:3]],
      comp_pts$c_lvl[comp_pts$dist_min  %in% dist_labs[1:3]])

