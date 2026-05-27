##### Extension of Prenzel (2024) ##############################################
# Martin Möhler, 2026
# based on:
# P. Prenzel: 'Kann man das überhaupt messen?' - Der Bedarf an detaillierten 
#              räumlichen Bevölkerungsdaten in der Wirtschaftsgeographie, 2024
# earlier version: 
# https://statshorts.blogspot.com/2025/04/herfindahl-hirschman-index-als-ma-fur.html
############################################################################## #


library(dplyr)        # data wrangling
library(readxl)       # reading in xlsx files
library(sf)           # working with geo data
library(ggplot2)      # plotting
library(patchwork)    # plotting

# adapt file path to local setting:
fpath <- "~/ckm_regio_2026/01_replication_Prenzel/data/"


# ----- Helper functions -----

getV_hhi <- function(Nj, N, K = 1) { 
  
  (4 * K / N^4) * sum(Nj^2) 
}

getV_ent <- function(Nj, N, K = 1) { 
  sum((log(N) - log(Nj) - 1)^2) * (K / N^2) 
}

getV_div1 <- function(Nj, N, K = 1, M = length(Nj), normfact = 1) { 
  
  normfact^2 * getV_hhi(Nj, N, K) / (1 - 1/M)^2 
}

getV_div2 <- function(Nj, N, K = 1, M = length(Nj), normfact = 1) { 
  
  normfact^2 * getV_ent(Nj, N, K) / log(M)^2 
}


# ----- (1) Prepare data -----

# data freely available from:
# https://ergebnisse.zensus2022.de/datenbank/online/url/b2afd8ab
cob1_path <- paste0(fpath, "1000A-1016_en_flat.csv")
cob2_path <- paste0(fpath, "1000A-1016_en_flat_rb.csv")
cob3_path <- paste0(fpath, "1000A-1016_en_flat_bl.csv")

# shapefile of Germany for plotting from:
# https://gdz.bkg.bund.de/index.php/default/digitale-geodaten/verwaltungsgebiete/verwaltungsgebiete-1-250-000-stand-01-01-vg250-01-01.html
shp1_path <- paste0(fpath, "vg250_01-01.utm32s.shape.ebenen/vg250_ebenen_0101/VG250_GEM.shp")
shp2_path <- paste0(fpath, "vg250_01-01.utm32s.shape.ebenen/vg250_ebenen_0101/VG250_LAN.shp")

# additional classifications of municipalities by BBSR
# https://www.bbsr.bund.de/BBSR/DE/forschung/raumbeobachtung/Raumabgrenzungen/downloads/download-referenzen.html
bbsr_path <- paste0(fpath, "raumgliederungen-referenzen-2022.xlsx")

# classification of municipalities to election districts & district-level election results
# https://www.bundeswahlleiterin.de/dam/jcr/aa868597-0e60-476c-bd2b-279c1e9a142a/btw25_wkr_gemeinden_20241130_utf8.csv
# https://www.bundeswahlleiterin.de/dam/jcr/f49a47a1-735b-4e9b-b4e1-4c73cad2292e/btw25_kerg2.csv
elec1_path <- paste0(fpath, "btw25_wkr_gemeinden_20241130_utf8.csv")
elec2_path <- paste0(fpath, "btw25_kerg2.csv")


## read in data

cob_data  <- read.csv2(cob1_path) # LAU level
cob_data2 <- read.csv2(cob2_path) # RB level needed for NUTS-2 
cob_data3 <- read.csv2(cob3_path) # BL level needed for NUTS-2

BL_NUTS2 <- c("01", "02", "04", "07", "10", "11", "12", "13", "15", "16")
cob_data3 <- cob_data3 %>% filter(X1_variable_attribute_code %in% BL_NUTS2)

# mark data by level of regional aggregation
cob_data$lvl  <- 0
cob_data2$lvl <- cob_data3$lvl <- 1
cob_data$lvl[cob_data$X1_variable_attribute_code == "DG"] <- 2

# bind together
cob_data <- rbind(cob_data, 
                  cob_data2[cob_data2$X1_variable_attribute_code != "DG", ],
                  cob_data3[cob_data3$X1_variable_attribute_code != "DG", ]) %>%
  filter(value != "-") %>%
  select(X1_variable_attribute_code, 
         X1_variable_attribute_label,
         X2_variable_attribute_code,
         X2_variable_attribute_label,
         value, lvl) %>%
  arrange(desc(lvl), X1_variable_attribute_code, X2_variable_attribute_code)

rm(cob_data2, cob_data3)

cob_data$value <- as.numeric(cob_data$value)
names(cob_data)[1:4] <- c("regio_code", "regio_label", "cob_code", "cob_label")

# discard category 'unknown / stateless' (see Prenzel, 2024, p.79)

regio_cob999 <- cob_data$regio_code[cob_data$cob_code == "LAND999"]

cob_data$value[cob_data$regio_code %in% regio_cob999 & cob_data$cob_code == ""] <-
  cob_data$value[cob_data$regio_code %in% regio_cob999 & cob_data$cob_code == ""] -
  cob_data$value[cob_data$regio_code %in% regio_cob999 & cob_data$cob_code == "LAND999"]

cob_data <- cob_data %>% filter(cob_code != "LAND999")

# no. of COB-categories
M <- length(unique(cob_data$cob_code[cob_data$cob_code != ""])) 

lau_codes   <- sort(unique(cob_data$regio_code[cob_data$lvl == 0]))
nuts2_codes <- sort(unique(cob_data$regio_code[cob_data$lvl == 1]))
n_lau   <- length(lau_codes) 
n_nuts2 <- length(nuts2_codes)

# list of LAU and NUTS-2 regions
frac_data <- data.frame(regio = c("DG", nuts2_codes, lau_codes),
                        region_typ = "admin.",
                        lvl = rep(c(2, 1, 0), c(1, n_nuts2, n_lau)),
                        K = 1, N = NA,
                        hhi  = NA, div1  = NA, ent  = NA, div2  = NA,
                        Vhhi = NA, Vdiv1 = NA, Vent = NA, Vdiv2 = NA)

## Compute regional indices

for(i in 1:nrow(frac_data)) {
  
  regio_totl <- cob_data$value[cob_data$regio_code == frac_data$regio[i] & cob_data$cob_code == ""]
  regio_cobs <- cob_data$value[cob_data$regio_code == frac_data$regio[i] & cob_data$cob_code != ""]
  
  cobs_shares <- regio_cobs / regio_totl
  frac_data$N[i] <- regio_totl
  
  # compute indices
  frac_data$hhi[i]  <- sum((cobs_shares)^2)
  frac_data$ent[i]  <- -sum(cobs_shares * log(cobs_shares))
  
  # compute variance multiplier
  frac_data$Vhhi[i]  <- getV_hhi(Nj  = regio_cobs, N = regio_totl, K = 1)
  frac_data$Vent[i]  <- getV_ent(Nj  = regio_cobs, N = regio_totl, K = 1)
  frac_data$Vdiv1[i] <- getV_div1(Nj = regio_cobs, N = regio_totl, K = 1, M = M, normfact = 100)
  frac_data$Vdiv2[i] <- getV_div2(Nj = regio_cobs, N = regio_totl, K = 1, M = M, normfact = 100)
}

boxplot(frac_data$hhi)
frac_data[which.min(frac_data$hhi), ]
# discard outlier for interpretability & enforce logical bounds
frac_data[which.min(frac_data$hhi), c("hhi", "ent", "Vhhi", "Vent", "Vdiv1", "Vdiv2")] <- NA
frac_data$hhi[frac_data$hhi > 1] <- 1
frac_data$ent[frac_data$ent < 0] <- 0

# standardize to DIV1, DIV2
frac_data$div1 <- (1 - frac_data$hhi) / (1 - 1/M) * 100
frac_data$div2 <- frac_data$ent / log(M) * 100

# index over national benchmark
frac_data$div1_diff <- frac_data$div1 - frac_data$div1[frac_data$regio == "DG"]
frac_data$div2_diff <- frac_data$div2 - frac_data$div2[frac_data$regio == "DG"]

maxDiv1 <- max(frac_data$div1, na.rm = TRUE)
maxDiv2 <- max(frac_data$div2, na.rm = TRUE)

# thresholds for indicator quality
max(frac_data[frac_data$lvl == 0 & sqrt(frac_data$Vdiv1) >= 1,   "N"], na.rm = TRUE)
min(frac_data[frac_data$lvl == 0 & sqrt(frac_data$Vdiv1) <= 0.1, "N"], na.rm = TRUE)
table(sqrt(frac_data[frac_data$lvl == 0, ]$Vdiv1) >= 1.0)
table(sqrt(frac_data[frac_data$lvl == 0, ]$Vdiv1) <= 0.1)


# ----- (2) Analyse composite regions ----

# ----- (2.1) BBSR categories -----

# read in BBSR municipality classifications
bbsr_raumg <- readxl::read_xlsx(bbsr_path, sheet = "Gemeindereferenz (inkl. Kreise)") %>%
  select(GEM2022_RS, GTU2022, GTU_NAME, RLG2022, RLG_NAME, GWS2022, GWS_NAME)
bbsr_raumg <- bbsr_raumg[-1, ]

# correcting leading zeros in regional code
bbsr_raumg$GEM2022_RS <- ifelse(nchar(bbsr_raumg$GEM2022_RS) == 11,
                                paste0("0", bbsr_raumg$GEM2022_RS), 
                                bbsr_raumg$GEM2022_RS)

# categories to factors
GTU_cat <- unique(bbsr_raumg$GTU_NAME)[order(unique(bbsr_raumg$GTU2022))]
RLG_cat <- unique(bbsr_raumg$RLG_NAME)[order(unique(bbsr_raumg$RLG2022))]
GWS_cat <- unique(bbsr_raumg$GWS_NAME)[order(unique(bbsr_raumg$GWS2022), decreasing = TRUE)]

bbsr_raumg$GTU2022 <- factor(bbsr_raumg$GTU2022, 
                             levels = sort(unique(bbsr_raumg$GTU2022)),
                             labels = GTU_cat)
bbsr_raumg$RLG2022 <- factor(bbsr_raumg$RLG2022, 
                             levels = sort(unique(bbsr_raumg$RLG2022)),
                             labels = RLG_cat)
bbsr_raumg$GWS2022 <- factor(bbsr_raumg$GWS2022, 
                             levels = sort(unique(bbsr_raumg$GWS2022), decreasing = TRUE),
                             labels = GWS_cat)

# check for full coverage
all(cob_data$regio_code[cob_data$lvl == 0] %in% bbsr_raumg$GEM2022_RS)

# merge BBSR categories to COB-data
cob_data <- merge(x = cob_data, y = select(bbsr_raumg, GEM2022_RS, GTU2022, RLG2022, GWS2022), 
                  by.x = "regio_code", by.y = "GEM2022_RS",
                  all.x = TRUE, all.y = FALSE)

# aggregate by BBSR categories
bbsr_data <- rbind(cob_data[cob_data$lvl == 0, ] %>% group_by(GTU2022, cob_code) %>%
                     summarise(value = sum(value), K = n()) %>% rename(regio = GTU2022),
                   cob_data[cob_data$lvl == 0, ] %>% group_by(RLG2022, cob_code) %>%
                     summarise(value = sum(value), K = n()) %>% rename(regio = RLG2022),
                   cob_data[cob_data$lvl == 0, ] %>% group_by(GWS2022, cob_code) %>%
                     summarise(value = sum(value), K = n()) %>% rename(regio = GWS2022))

# list composite regions
frac_bbsr <- data.frame(regio = unique(bbsr_data$regio),
                        region_typ = "BBSR",
                        lvl = rep(c("GTU", "RLG", "GWS"), 
                                  c(length(GTU_cat), 
                                    length(RLG_cat), 
                                    length(GWS_cat) - 1)),
                        K = bbsr_data$K[bbsr_data$cob_code == ""],
                        N = bbsr_data$value[bbsr_data$cob_code == ""],
                        hhi  = NA, div1  = NA, ent  = NA, div2  = NA,
                        Vhhi = NA, Vdiv1 = NA, Vent = NA, Vdiv2 = NA)

# Compute indices for composite regions

for(i in 1:nrow(frac_bbsr)) {
  
  regio_totl <- bbsr_data$value[bbsr_data$regio == frac_bbsr$regio[i] & bbsr_data$cob_code == ""]
  regio_cobs <- bbsr_data$value[bbsr_data$regio == frac_bbsr$regio[i] & bbsr_data$cob_code != ""]
  
  cobs_shares <- regio_cobs / regio_totl
  K <- frac_bbsr$K[i]
  
  # compute indices
  frac_bbsr$hhi[i] <- sum((cobs_shares)^2)
  frac_bbsr$ent[i] <- -sum(cobs_shares * log(cobs_shares))
  
  # compute variance multiplier
  frac_bbsr$Vhhi[i] <- getV_hhi(Nj = regio_cobs, N = regio_totl, K = K)
  frac_bbsr$Vent[i] <- getV_ent(Nj = regio_cobs, N = regio_totl, K = K)
  frac_bbsr$Vdiv1[i] <- getV_div1(Nj = regio_cobs, N = regio_totl, K = K, M = M, normfact = 100)
  frac_bbsr$Vdiv2[i] <- getV_div2(Nj = regio_cobs, N = regio_totl, K = K, M = M, normfact = 100)
}

# standardize to DIV1, DIV2
frac_bbsr$div1 <- (1 - frac_bbsr$hhi) / (1 - 1/M) * 100
frac_bbsr$div2 <- frac_bbsr$ent / log(M) * 100
# index over national benchmark
frac_bbsr$div1_diff <- frac_bbsr$div1 - frac_data$div1[frac_data$regio == "DG"]
frac_bbsr$div2_diff <- frac_bbsr$div2 - frac_data$div2[frac_data$regio == "DG"]
# multiplier for the standard error
frac_bbsr$SEdiv1 <- sqrt(frac_bbsr$Vdiv1)
frac_bbsr$SEdiv2 <- sqrt(frac_bbsr$Vdiv2)

rbind(index_ratio = frac_bbsr$div2/ frac_bbsr$div1, 
      SE_multiplier_ratio = frac_bbsr$SEdiv2 / frac_bbsr$SEdiv1)


# ----- (2.2) Election results -----

# read in electoral district classification
elec_class <- read.csv2(elec1_path, skip = 7, colClasses = "character")

# create region key for electoral district data + match key in census data

elec_class$regio_code8 <- paste0(elec_class$RGS_Land, elec_class$RGS_RegBez,
                                 elec_class$RGS_Kreis, elec_class$RGS_Gemeinde)

cob_data$regio_code8 <- ifelse(nchar(cob_data$regio_code) < 12, NA, 
                               paste0(substr(cob_data$regio_code, 1, 5), 
                                      substr(cob_data$regio_code, 10, 12)))

# region keys are coarser in electoral district data and some need to be corrected;
# some redistricting needs to also be corrected
regio8_nf <- unique(cob_data$regio_code8[(!cob_data$regio_code8 %in% elec_class$regio_code8) 
                                               & cob_data$lvl == 0])
nf_pos <- match(cob_data$regio_code8[cob_data$regio_code8 %in% regio8_nf], regio8_nf)

regio8_new <- c("01059126", "01059126", "06635001", "13071043", "13071006",
                "14522275", "14522275", "16061045", "16061045", "16061119",
                "16061119", "16061119", "16061119", "16062062", "16061119",
                "16061119", "16061119", "16061119", "16061119", "16061119",
                "16063104", "16064074", "16061118", "16064071", "16064071",
                "16061118", "16064076", "16066042", "16067092", "16068064",
                "16071004", "16076094", "16076094", "16076003", "16076039")

cob_data$regio_code8[cob_data$regio_code8 %in% regio8_nf] <- regio8_new[nf_pos]

# merge electoral districts to COB-data
cob_data <- merge(x = cob_data, 
                  y = elec_class[!duplicated(elec_class$regio_code8), 
                                 c("Wahlkreis.Nr", "regio_code8")],
                  by = "regio_code8", all.x = TRUE, all.y = FALSE)

# aggregate by electoral district
elec_data <- cob_data[!is.na(cob_data$Wahlkreis.Nr), ] %>% 
  group_by(Wahlkreis.Nr, cob_code) %>%
  summarise(value = sum(value), K = n())

# list composite regions
frac_elec <- data.frame(regio = unique(elec_data$Wahlkreis.Nr),
                        region_typ = "Wahlkreise",
                        lvl = "elec.dist.",
                        K = elec_data$K[elec_data$cob_code == ""],
                        N = elec_data$value[elec_data$cob_code == ""],
                        hhi  = NA, div1  = NA, ent  = NA, div2  = NA,
                        Vhhi = NA, Vdiv1 = NA, Vent = NA, Vdiv2 = NA)

# Compute indices for composite regions

for(i in 1:nrow(frac_elec)) {
  
  regio_totl <- elec_data$value[elec_data$Wahlkreis.Nr == frac_elec$regio[i] & elec_data$cob_code == ""]
  regio_cobs <- elec_data$value[elec_data$Wahlkreis.Nr == frac_elec$regio[i] & elec_data$cob_code != ""]
  
  cobs_shares <- regio_cobs / regio_totl
  K <- frac_elec$K[i]
  
  # compuhte indices
  frac_elec$hhi[i] <- sum((cobs_shares)^2)
  frac_elec$ent[i] <- -sum(cobs_shares * log(cobs_shares))
  
  # compute variance multiplier
  frac_elec$Vhhi[i] <- getV_hhi(Nj = regio_cobs, N = regio_totl, K = K)
  frac_elec$Vent[i] <- getV_ent(Nj = regio_cobs, N = regio_totl, K = K)
  frac_elec$Vdiv1[i] <- getV_div1(Nj = regio_cobs, N = regio_totl, K = K, M = M, normfact = 100)
  frac_elec$Vdiv2[i] <- getV_div2(Nj = regio_cobs, N = regio_totl, K = K, M = M, normfact = 100)
}

# standardize to DIV1, DIV2
frac_elec$div1 <- (1 - frac_elec$hhi) / (1 - 1/M) * 100
frac_elec$div2 <- frac_elec$ent / log(M) * 100
# index over national benchmark
frac_elec$div1_diff <- frac_elec$div1 - frac_data$div1[frac_data$regio == "DG"]
frac_elec$div2_diff <- frac_elec$div2 - frac_data$div2[frac_data$regio == "DG"]
# multiplier for the standard error
frac_elec$SEdiv1 <- sqrt(frac_elec$Vdiv1)
frac_elec$SEdiv2 <- sqrt(frac_elec$Vdiv2)
# to order of magnitude for plotting
frac_elec$SEOOMdiv1 <- floor(log10(frac_elec$SEdiv1))
frac_elec$SEOOMdiv2 <- floor(log10(frac_elec$SEdiv2))
maxOOMdiv1 <- max(frac_elec$SEOOMdiv1)
maxOOMdiv2 <- max(frac_elec$SEOOMdiv2)
frac_elec$div1_weak <- factor(frac_elec$SEOOMdiv1 == maxOOMdiv1, levels = c(TRUE, FALSE),
                              labels = c(paste(maxOOMdiv1), paste("<", maxOOMdiv1)))
frac_elec$div2_weak <- factor(frac_elec$SEOOMdiv2 == maxOOMdiv2, levels = c(TRUE, FALSE),
                              labels = c(paste(maxOOMdiv2), paste("<", maxOOMdiv2)))
table(frac_elec$div1_weak)
table(frac_elec$div2_weak)


## add election results

# read in results data
cClass <- c(rep("character", 9), rep("numeric", 8), rep("character", 2))
elec_res <- read.csv2(elec2_path, skip = 9, colClasses = cClass) %>%
  filter(Stimme == 2 & Gruppenart == "Partei") %>%
  select(Gebietsart, Gebietsnummer, Gruppenart, Gruppenname, Anzahl)

elec_full <- elec_res %>% filter(Gebietsart == "Bund") %>% select(-Gebietsart)
elec_res  <- elec_res %>% filter(Gebietsart == "Wahlkreis") %>% select(-Gebietsart)

# municipalities with several electoral districts need to have theirs combined
# in the results data (since we only have HHI at municipality level)
double_dist <- unique(elec_class$regio_code8[duplicated(elec_class$regio_code8)])

for(i in seq(double_dist)) {
  # re-name electoral districts to reflect combination
  ed_nrs <- elec_class$Wahlkreis.Nr[elec_class$regio_code8 == double_dist[i]]
  elec_res$Gebietsnummer[elec_res$Gebietsnummer %in% ed_nrs] <- ed_nrs[1]
}

# re-compute percentages
elec_res <- elec_res %>% group_by(Gebietsnummer, Gruppenname) %>%
  summarise(Anzahl = sum(Anzahl, na.rm = TRUE))
elec_sums <- elec_res %>% group_by(Gebietsnummer) %>%
  summarise(Total = sum(Anzahl))
elec_res <- merge(elec_res, elec_sums, by = "Gebietsnummer")
elec_res$Prozent <- (elec_res$Anzahl / elec_res$Total) * 100
# compute country-level percentage benchmark
elec_full$Total = sum(elec_full$Anzahl)
elec_full$Prozent_Bund <- (elec_full$Anzahl / elec_full$Total) * 100 

# merge country results baseline
elec_res <- merge(x = elec_res, y = elec_full[, c("Gruppenname", "Prozent_Bund")],
                  by = "Gruppenname", all.x = TRUE, all.y = FALSE)
elec_res$Prozent_diff <- elec_res$Prozent - elec_res$Prozent_Bund

# merge election results and electoral districts HHIs
elec_res <- merge(x = elec_res, y = frac_elec,
                  by.x = "Gebietsnummer", by.y = "regio",
                  all.x = TRUE, all.y = FALSE)

# subset to relevant parties
elec_res <- elec_res %>% filter(Gruppenname %in% c("AfD", "CDU", "FDP", "SPD", "GRÜNE", "BSW"))


# ----- (3) Plots -----

# ----- (3.1) Choropleth maps -----

## prepare shape files
shp1 <- read_sf(shp1_path) %>% filter(GF == 4) # LAU shapes
shp2 <- read_sf(shp2_path) %>% filter(GF == 4) # BL shapes

# derive NUTS 2 shapes
shp3 <- shp1
shp3$ARS <- ifelse(substr(shp1$ARS_0, 1, 2) %in% BL_NUTS2, 
                   substr(shp1$ARS_0, 1, 2), substr(shp1$ARS_0, 1, 3))
shp3 <- shp3 %>% group_by(ARS) %>% summarise() %>% st_cast()

# derive composite shapes
shp1 <- merge(shp1, bbsr_raumg, by.x = "ARS", by.y = "GEM2022_RS", all.x = TRUE)
shp1 <- merge(shp1, select(elec_class, Wahlkreis.Nr, Wahlkreis.Bez, regio_code8),
              by.x = "AGS", by.y = "regio_code8", all.x = TRUE)

shp_gtu  <- shp1 %>% group_by(GTU2022)      %>% summarise() %>% st_cast()
shp_rlg  <- shp1 %>% group_by(RLG2022)      %>% summarise() %>% st_cast()
shp_gws  <- shp1 %>% group_by(GWS2022)      %>% summarise() %>% st_cast()
shp_elec <- shp1 %>% group_by(Wahlkreis.Nr) %>% summarise() %>% st_cast()

## prepare data for plotting

# elementary regions

shp1$lvl <- 0
shp2$lvl <- shp3$lvl <- 1

map_data <- rbind(select(shp1, ARS, lvl), shp3)
map_data <- merge(x = map_data, y = select(frac_data, -lvl),
                  by.x = "ARS", by.y = "regio",
                  all.x = TRUE, all.y = FALSE)

n_lau_map   <- nrow(map_data[map_data$lvl == 0, ])
n_nuts2_map <- nrow(map_data[map_data$lvl == 1, ])

map_div <- rbind(map_data[map_data$lvl == 0, ], 
                 map_data[map_data$lvl == 0, ],
                 map_data[map_data$lvl == 1, ], 
                 map_data[map_data$lvl == 1, ])

map_div$Indexwert <- c(map_data$div1[map_data$lvl == 0], 
                       map_data$div2[map_data$lvl == 0],
                       map_data$div1[map_data$lvl == 1], 
                       map_data$div2[map_data$lvl == 1])

map_div$lvl <- factor(map_div$lvl, levels = c(1, 0), 
                      labels = c("NUTS 2", "Gemeinden"))

map_div$Index <- factor(c(rep(c(1, 2), each = n_lau_map), 
                          rep(c(1, 2), each = n_nuts2_map)),
                        labels = c("DIV1", "DIV2"))

map_lnd <- rbind(shp2, shp2) %>% select(ARS)
map_lnd$lvl   <- factor(rep(c(0, 1), each = 16), levels = c(1, 0), 
                        labels = c("NUTS 2", "Gemeinden"))

# composite regions

names(shp_gtu)[1] <- names(shp_rlg)[1] <- names(shp_gws)[1] <- names(shp_elec)[1] <- "regio"

shp_gtu  <- merge(shp_gtu,  frac_bbsr,         by = "regio", all.x = TRUE)
shp_rlg  <- merge(shp_rlg,  frac_bbsr,         by = "regio", all.x = TRUE)
shp_gws  <- merge(shp_gws,  frac_bbsr,         by = "regio", all.x = TRUE)
shp_elec <- merge(shp_elec, frac_elec[, 1:17], by = "regio", all.x = TRUE)

labsc <- c("Stadt- und Gemeindetyp", "Lagetyp", "Wachstumstendenz", "Wahlkreise")

shp_gtu$lvl  <- 1
shp_rlg$lvl  <- 2
shp_gws$lvl  <- 3
shp_elec$lvl <- 4
shp_gtu$region_typ <- shp_rlg$region_typ <- shp_gws$region_typ <- "BBSR"
shp_elec$region_typ <- "Wahlkreise"

mapc <- rbind(shp_gtu, shp_rlg, shp_gws, shp_elec) %>%
  filter(!is.na(K))
mapc$lvl <- factor(mapc$lvl, levels = 1:4, labels = labsc)
mapc$lbl <- factor(mapc$regio, labels = paste0(mapc$regio, "\nK = ", mapc$K))

mapc_lnd <- rbind(shp2, shp2, shp2, shp2) %>% select(ARS)
mapc_lnd$lvl <- factor(rep(1:4, each = 16), levels = 1:4, labels = labsc)


## maps

# DIV1 maps (NUTS 2 and LAU)
ggplot(map_div[map_div$Index == "DIV1", ]) +
  geom_sf(aes(fill = Indexwert), color = NA) +
  scale_fill_viridis_c(direction = -1, na.value = "grey50", limits = c(0, maxDiv1),
                       name = "Index-\nwert") +
  geom_sf(data = map_lnd, fill = NA, color = "black") +
  theme_void() +
  theme(legend.key.width = unit(.25, "cm"),
        plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~lvl) +
  ggtitle("DIV1")
#ggsave("DIV1_map.png", width = 1600, height = 1000, units = "px")

# DIV2 maps (NUTS 2 and LAU)
ggplot(map_div[map_div$Index == "DIV2", ]) +
  geom_sf(aes(fill = Indexwert), color = NA) +
  scale_fill_viridis_c(direction = -1, na.value = "grey50", limits = c(0, maxDiv2),
                       name = "Index-\nwert") +
  geom_sf(data = map_lnd, fill = NA, color = "black") +
  theme_void() +
  theme(legend.key.width = unit(.25, "cm"),
        plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~lvl) +
  ggtitle("DIV2")
#ggsave("DIV2_map.png", width = 1600, height = 1000, units = "px")

# DIV1 maps (selected composite regions)
c_select <- labsc[c(2, 4)]
ggplot(mapc[mapc$lvl %in% c_select, ]) +
  geom_sf(aes(fill = div1), color = NA) +
  scale_fill_viridis_c(direction = -1, na.value = "grey50", limits = c(0, maxDiv1),
                       name = "Index-\nwert") +
  geom_sf(data = mapc_lnd[mapc_lnd$lvl %in% c_select, ], fill = NA, color = "black") +
  theme_void() +
  theme(legend.key.width = unit(.25, "cm"),
        plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~lvl) +
  ggtitle("DIV1")
#ggsave("DIV1_map_c.png", width = 1600, height = 1000, units = "px")

# Map: Lagetyp

ggplot(mapc[mapc$lvl == labsc[2] & !is.na(mapc$regio), ]) +
  geom_sf(data = map_lnd) +
  geom_sf(fill = "blue", color = NA) +
  facet_wrap(~lbl, nrow = 1) +
  theme_void() +
  theme(strip.text = element_text(size = 10))
#ggsave("map_lagetyp.png", height = 900, width = 2100, units = "px")


# ----- (3.2) Election scatter plot -----

parties <- sort(unique(elec_res$Gruppenname))
elec_res$pseudonym <- paste("Partei", match(elec_res$Gruppenname, parties))

ggplot(elec_res) +
  geom_hline(yintercept = 0, lty = "dashed") +
  geom_vline(xintercept = 0, lty = "dashed") +
  geom_point(aes(div1_diff, Prozent_diff, color = SEdiv1), 
             size = .5, alpha = .6) +
  scale_color_gradient(low = "blue", high = "orange", trans = "log10",
                       name = "Multiplikator des Standardfehlers") +
  facet_wrap(~pseudonym, nrow = 2) +
  xlab("DIV1 (Diff. zu Deutschland insg.)") +
  ylab("Zweitstimmen (Diff. zum Bundesergebnis, %pkt.)") +
  theme_bw() +
  theme(legend.position = "bottom", 
        legend.key.height = unit(.25, "cm"),
        legend.key.width = unit(.7, "cm"))
#ggsave("DIV1_elec.png", width = 2000, height = 1300, units = "px")


# ----- (3.3) Error plots -----

plotvars <- c("regio", "region_typ", "lvl", "K", "N", "div1", "div2", "Vdiv1", "Vdiv2")

frac_full <- rbind(frac_bbsr[, plotvars], frac_elec[, plotvars],
                   frac_data[frac_data$lvl == 0 & !is.na(frac_data$div1), plotvars])

frac_full$region_typ <- factor(frac_full$region_typ, 
                               levels = c("admin.", "BBSR", "Wahlkreise"),
                               labels = c("Gemeinden", "BBSR-Gliederung", "Wahlkreise"))

p1 <- ggplot(frac_full[frac_full$region_typ == "Gemeinden", ], aes(N/K, sqrt(Vdiv1))) +
  geom_point(size = .5, alpha = .6, color = "grey30") +
  xlab("N") +
  ylab("Multiplikator des Standardfehlers") +
  scale_x_continuous(trans = "log10") +
  facet_wrap(~region_typ) +
  theme_bw() +
  ggtitle("a)") +
  theme(legend.key.width = unit(.25, "cm"), legend.position = "inside", 
        legend.position.inside = c(0.85, 0.6), legend.background = element_rect(color = "grey50"))

p2 <- ggplot(frac_full[frac_full$region_typ != "Gemeinden", ], aes(N/K, sqrt(Vdiv1))) +
  geom_point(size = .5, alpha = .6, color = "grey30") +
  xlab("N / K") +
  ylab(NULL) +
  scale_x_continuous(trans = "log10") +
  facet_wrap(~region_typ) +
  theme_bw() +
  ggtitle("b)") +
  theme(legend.key.width = unit(.25, "cm"))

layout <- c(area(1, 1, 1, 1), area(1, 2, 1, 3))
p1 + p2 + plot_layout(design = layout)
#ggsave("error_plot_div1.png", width = 2400, height = 1100, units = "px")

