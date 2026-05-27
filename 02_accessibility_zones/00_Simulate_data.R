##### Grid data accessibility zones (preparation step) #########################
# Martin Möhler, 2026
# main script: 01_Accessibility_zones.R
############################################################################## #

library(data.table)
library(dplyr)
library(sf)

# adapt file path to local setting:
fpath <- "~/ckm_regio_2026/02_accessibility_zones/data/"
spath <- "~/ckm_regio_2026/01_replication_Prenzel/data/vg250_01-01.utm32s.shape.ebenen/vg250_ebenen_0101/"


# ----- (1) Prepare Data -----

## load data

# population at 100m grid cell level
# https://www.destatis.de/DE/Themen/Gesellschaft-Umwelt/Bevoelkerung/Zensus2022/_inhalt.html#1403950
pop_grd <- fread(paste0(fpath, "Zensus2022_Bevoelkerungszahl_100m-Gitter.csv"))

# LAU-level shape file
# https://gdz.bkg.bund.de/index.php/default/digitale-geodaten/verwaltungsgebiete/verwaltungsgebiete-1-250-000-stand-01-01-vg250-01-01.html
shp_gem <- st_read(paste0(spath, "VG250_GEM.shp")) %>% filter(GF == 4) %>%
  select(ARS, GEN, SN_L) %>%
  st_transform(crs = st_crs("EPSG:3035"))

# LAU-level pop. aged 65+
lau_65p <- fread(paste0(fpath, "1000A-1003_en_flat.csv"))
lau_65p$value[lau_65p$value == "-"] <- "0"
lau_65p$value <- as.numeric(lau_65p$value)
names(lau_65p)[8] <- "ARS"


## subset to a single state

bl <- "05" # select federal state to simulate data for (05 = NW)

# subset LAU info to bl
lau_65p <- lau_65p %>% select(ARS, value) %>% filter(substr(ARS, 1, 2) == bl)

# join LAU info to grid cells
pop_grd <- pop_grd %>%
  st_as_sf(coords = c("x_mp_100m", "y_mp_100m"), crs = st_crs(shp_gem)) %>%
  st_join(filter(shp_gem[, "ARS"], ARS %in% lau_65p$ARS)) %>%
  filter(!is.na(ARS))


# ----- (2) Simulate Data -----

# unique LAU in selected federal state
ars <- unique(pop_grd$ARS)
N <- sum(lau_65p$value)

pop <- data.frame(person_id = 1:N, 
                  Gitter_ID_100m = NA, 
                  GEN = NA, ARS = NA)

set.seed(20260302)
j <- 1
for(i in ars) {
  
  # no. of people aged 65+ in LAU to allocate 
  n <- lau_65p$value[lau_65p$ARS == i]
  # reserve space in unit-level data.frame
  j_next <- j + n
  rng <- j:(j_next - 1)

  # make sample frame 
  cells  <- pop_grd$GITTER_ID_100m[pop_grd$ARS == i]
  sframe <- rep(seq(cells), pop_grd$Einwohner[pop_grd$ARS == i])
  # draw sample
  samp <- sample(sframe, n, replace = FALSE)
  # store grid cells drawn
  pop$Gitter_ID_100m[rng] <- cells[samp]
  # add LAU ID and name
  pop$ARS[rng] <- i
  pop$GEN[rng] <- shp_gem$GEN[shp_gem$ARS == i]
  
  j <- j_next
}
# add grid center point info
pop$x_mp_100m <- as.numeric(substr(pop$Gitter_ID_100m, 24, 30)) + 50
pop$y_mp_100m <- as.numeric(substr(pop$Gitter_ID_100m, 16, 22)) + 50
# simulate exact locations (optional)
pop$x <- floor(pop$x_mp_100m + runif(N, -50, 50))
pop$y <- floor(pop$y_mp_100m + runif(N, -50, 50))

save(pop, file = paste0(fpath, "pop_data.RData"))

