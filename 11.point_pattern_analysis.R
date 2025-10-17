#
# Spatial clustering of ranked university (& population) with Gbif observations
# Drosophila Paper - Hjelman & Curnow
# Point Pattern Analysis - Olson 07/2025
#
# # # # # # # # # # # # #

# libraries
library(terra);library(sf)
library(dplyr);library(tidyr)
library(ggplot2)
library(spatstat)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#     LOAD DATASETS
###########################################################

# GBIF DATA
# pth_gbif = "data/output/all.locale.gbif.csv" # SUBSET (11022 obs)
pth_gbif <- "data/output/all.drosophila.gbif.csv" # ALL DATA (63320 obs)
dfgbif <- read.csv(pth_gbif)
gbif <- vect(dfgbif, geom=c("long", "lat"), crs="EPSG:4326")

# POPULATION DATA
# file exceeds 100 MB github filesize limit
# download from https://data.worldpop.org/GIS/Population/Global_2021_2022_1km_UNadj/
popath <- 'path_to_file/global_ppp_2022_1km_UNadj_constrained.tif'
popul <- rast(popath)
pop <- as.points(popul)
# pop density > 10K
pop2 = pop[pop$global_ppp_2022_1km_UNadj_constrained>10e3,]

# UNIVERSITY DATA
pth_u <-  "data/output/university_loc_rankings.csv"
dfuniv <-  read.csv(pth_u)
# Top 500 ranked universities
uv1500 <- vect(dfuniv[dfuniv$Rank<=500, ], geom=c("long", "lat"), crs="EPSG:4326")

# BIOCLIMATIC REGIONS
eco_realm2 = vect("data/dros_spatial/eco-realm2.geojson")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#     HELPER FUNCTIONS
###########################################################

# FUNC - reprojected global window
create_window_global <- function(global=TRUE){
  r = rast("data/dros_spatial/wc2_10m_bio_land.tif")
  ry = patches(r)
  ry[ry==1 | ry>=1634] = NA
  rz = zonal(cellSize(ry,unit="km"), ry, sum, as.raster=TRUE)
  # remove small patches
  s = ifel(rz < 1.2e5, NA, 1)
  sp2 = project(as.polygons(s), "EPSG:6933")
  return(sp2)
}

# FUNC - data prep for point pattern object
prep_k_data <- function(univ_vect, gbif_vect, sp2){
  #
  # required
  #   uv - university vector object
  #   gv - gbif spatial object
  #   sp2 - window object
  # 
  
  # project and crop to extent
  univ = project(univ_vect, sp2)[sp2,]
  gbif = project(gbif_vect, sp2)[sp2,]
  
  # create dataframes
  obs_coords <- as.data.frame(geom(gbif)[, c("x", "y")])
  school_coords <- as.data.frame(geom(univ)[, c("x", "y")])
  school_coords = school_coords %>% na.omit()
  
  # remove duplicates
  if(sum(duplicated(obs_coords))>100){
    print(paste("removing", sum(duplicated(obs_coords)), "duplicates..."))
    obs_coords = obs_coords[!duplicated(obs_coords[, c("x","y")]),]
  }
  
  # create ppp
  my_window <- as.owin(st_as_sf(sp2))  
  obs_ppp <- ppp(obs_coords$x, obs_coords$y, window = my_window)
  school_ppp <- ppp(school_coords$x, school_coords$y, window = my_window)
  
  # rescale
  ppp_obs <- rescale(obs_ppp, 1000, "km")
  ppp_univ <- rescale(school_ppp, 1000, "km")
  
  # combine
  # Ensure both ppp objects share the same window
  common_window <- union.owin(ppp_obs$window, ppp_univ$window)
  # Set marks as you create the objects (not afterward)
  ppp_obs_marked <- ppp(x = ppp_obs$x, y = ppp_obs$y, window = common_window,
                        marks = factor(rep("obs", ppp_obs$n)))
  ppp_univ_marked <- ppp(x = ppp_univ$x, y = ppp_univ$y, window = common_window,
                         marks = factor(rep("univ", ppp_univ$n)))
  # Combine them using superimpose()
  ppp_all <- superimpose(ppp_obs_marked, ppp_univ_marked)
  
  return(ppp_all)
}


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#     ANALYSIS (PCF inhom varying)
###########################################################

# LOAD window
sp2 <- create_window_global()

# # # # # # # # # # # # # # # # # # # #
# CLUSTERING FROM UNIVERSITIES
univ_lists = c(10, 50, 100, 500, 1000, 1500)  # ranked subsets

# LOOP
for (rank_level in univ_lists){
  # iterate through subsets
  print(paste("...running rank:", rank_level,"or", match(rank_level, univ_lists), "of", length(univ_lists)))
  uv0 <- vect(dfuniv[dfuniv$Rank<=rank_level, ], geom=c("long", "lat"), crs="EPSG:4326")
  
  # prep data
  cross_all = prep_k_data(uv0, gv, sp2)
  
  # run cross (PCF inhom)
  L_obs_cross <- pcfcross.inhom(cross_all, i = "univ", j = "obs")
  
  # convert to df and save
  L_obs_cross = as.data.frame(L_obs_cross)
  L_obs_cross$univ_rank = rank_level
  
  if(match(rank_level, univ_lists)==1){
    df_cross = L_obs_cross
  }else{
    # c("r", "theo", "border", "trans", "iso")
    df_cross = rbind(df_cross, L_obs_cross)
  }
}
names(df_cross) = c("r", "theo", "trans", "iso", "univ_rank")
head(df_cross);dim(df_cross)

# PLOT UNIVERSITY
df_long_univ <- df_cross %>%
  pivot_longer(cols = !c(r, univ_rank),
               names_to = "variable",
               values_to = "value") %>%
  mutate(univ_rank = as.factor(univ_rank))

df_long_univ %>% 
  filter(variable=='trans') %>% # only show one border correction method
  # filter(univ_rank%in%c(100, 500, 1000)) %>%
  ggplot( aes(x = r, y = value, color = univ_rank)) +
  geom_line(linewidth=1.2) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.5) +
  scale_color_grey(start = 0.1, end = 0.7) +
  labs(x = "r (km)", y = expression(hat(g)[obs/univ]), color = "Rank") +
  xlim(c(0,100)) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 14), 
    axis.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.key.width = unit(1,"cm"),
    legend.position = c(0.7,0.75),
    legend.background = element_rect(fill = "white", color = "black")) +
  guides(linetype = guide_legend(override.aes = list(linewidth = 4)))


# # # # # # # # # # # # # # # # # # # # # # #
# CLUSTERING FROM POPULATION CENTERS
pop_lists = c(10, 20, 30, 40)  *1e3
#
# LOOP
for (pop_level in pop_lists){
  print(paste("...running rank:", pop_level,"or", match(pop_level, pop_lists), "of", length(pop_lists)))
  pop0 <- pop[pop$global_ppp_2022_1km_UNadj_constrained>pop_level,]
  print(dim(pop0))
  # prep data
  cross_all = prep_k_data(pop0, gbif, sp2)
  
  # run cross (PCF inhom)
  L_obs_cross <- pcfcross.inhom(cross_all, i = "univ", j = "obs")
  
  L_obs_cross = as.data.frame(L_obs_cross)
  L_obs_cross$univ_rank = pop_level
  
  if(match(pop_level, pop_lists)==1){
    df_cross_pop = L_obs_cross
  }else{
    # c("r", "theo", "border", "trans", "iso")
    df_cross_pop = rbind(df_cross_pop, L_obs_cross)
  }
}
names(df_cross_pop) = c("r", "theo", "trans", "iso", "pop")
head(df_cross_pop);dim(df_cross_pop)

#PLOT
df_long_pop <- df_cross_pop %>%
  pivot_longer(cols = !c(r, pop),
               names_to = "variable",
               values_to = "value") %>%
  mutate(pop = as.factor(pop))

df_long_pop %>% 
  filter(variable=='trans') %>%
  # filter(pop%in%c(10000, 20000)) %>%
  ggplot( aes(x = r, y = value, color = pop)) +
  geom_line(linewidth=1.2) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.5) +
  scale_color_grey(start = 0.1, end = 0.7) +
  labs(x = "r (km)", y = expression(hat(g)[obs/pop]),color = "Population") +
  xlim(c(0,100)) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 14), 
    axis.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.key.width = unit(1,"cm"),
    legend.position = c(0.7,0.75),
    legend.background = element_rect(fill = "white", color = "black")) +
  guides(linetype = guide_legend(override.aes = list(linewidth = 4)))


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#     ANALYSIS (PCF inhom ECO REALMS)
###########################################################

# classify points in eco
gbif$eco = NA
for(ue in unique(eco_realm2$WWF_REALM2)){
  # gv[eco_realm2[eco_realm2$WWF_REALM2=="Australasia"],]
  pts_nam = gbif[eco_realm2[eco_realm2$WWF_REALM2==ue],]
  gv$eco[pts_nam$X] = ue
}

# UNIVERSITY ECO
# FIXED (vary univ rank here)
uv <- vect(dfuniv[dfuniv$Rank<=500, ], geom=c("long", "lat"), crs="EPSG:4326")
ue_names = setdiff(unique(eco_realm2$WWF_REALM2), c("Antarctic", "Oceania"))
#
# LOOP
for(ue in ue_names){
  print(paste("...running realm:", ue,"or", match(ue, unique(eco_realm2$WWF_REALM2)), "of", length(unique(eco_realm2$WWF_REALM2))))
  
  # subsets
  uv_sub = uv[eco_simple2[eco_simple2$WWF_REALM2==ue],]
  gv_sub = gbif[eco_simple2[eco_simple2$WWF_REALM2==ue],]
  
  # prep data
  cross_eco_all = prep_k_data(uv_sub, gv_sub, sp2)
  
  # run cross
  L_eco_cross <- pcfcross.inhom(cross_eco_all, i = "univ", j = "obs")
  
  L_eco_cross = as.data.frame(L_eco_cross)
  L_eco_cross$realm = ue
  
  if(match(ue, unique(eco_simple2$WWF_REALM2))==1){
    df_eco_cross = L_eco_cross
  }else{
    # c("r", "theo", "border", "trans", "iso")
    df_eco_cross = rbind(df_eco_cross, L_eco_cross)
  }
  
}
head(df_eco_cross)

# Plot
df_eco_cross %>%
  pivot_longer(cols = !c(r, realm),
               names_to = "variable",
               values_to = "value") %>% 
  filter(variable=='trans') %>%
  ggplot( aes(x = r, y = value,
              color = realm)) +
  geom_line(linewidth=1.2) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.5) +
  scale_color_grey(start = 0.1, end = 0.7) +
  labs(x = "r (km)", y = expression(hat(g)[obs/univ]), color = "Realm") +
  xlim(c(0,100)) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 14), 
    axis.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.key.width = unit(1,"cm"),
    legend.position = c(0.7,0.75),
    legend.background = element_rect(fill = "white", color = "black")) +
  guides(linetype = guide_legend(override.aes = list(linewidth = 4)))


# POPULATION
# FIXED pop size
pop_sub = pop[pop$global_ppp_2022_1km_UNadj_constrained>10000,]
ue_names = setdiff(unique(eco_realm2$WWF_REALM2), c("Antarctic", "Oceania"))
#
# LOOP
for(ue in ue_names){
  
  print(paste("...running realm:", ue,"or", match(ue, unique(eco_realm2$WWF_REALM2)), "of", length(unique(eco_realm2$WWF_REALM2))))
  
  # subsets
  pop_sub = pop_sub[eco_simple2[eco_simple2$WWF_REALM2==ue],]
  gv_sub = gbif[eco_simple2[eco_simple2$WWF_REALM2==ue],]
  
  # prep data
  cross_eco_all = prep_k_data(pop_sub, gv_sub, sp2)
  
  # run cross
  L_eco_cross <- pcfcross.inhom(cross_eco_all, i = "univ", j = "obs")
  
  L_eco_cross = as.data.frame(L_eco_cross)
  L_eco_cross$realm = ue
  
  if(match(ue, unique(eco_simple2$WWF_REALM2))==1){
    df_eco_pop = L_eco_cross
  }else{
    # c("r", "theo", "border", "trans", "iso")
    df_eco_pop = rbind(df_eco_cross, L_eco_cross)
  }
  
}
head(df_eco_pop)

# PLOT
df_eco_pop %>%
  pivot_longer(cols = !c(r, realm),
               names_to = "variable",
               values_to = "value") %>% 
  filter(variable=='trans') %>%
  ggplot( aes(x = r, y = value,
              color = realm)) +
  geom_line(linewidth=1.2) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.5) +
  # Option 1 (automatic greyscale)
  scale_color_grey(start = 0.1, end = 0.7) +
  labs(x = "r (km)",
       y = expression(hat(g)[obs/univ]),
       color = "Realm") +
  xlim(c(0,100)) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 14), 
    axis.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.key.width = unit(1,"cm"),
    legend.position = c(0.7,0.75),
    legend.background = element_rect(fill = "white", color = "black")) +
  guides(linetype = guide_legend(override.aes = list(linewidth = 4)))


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#     SIMULATION ENVELOPES
###########################################################

# determine distances where g is significant compared with random labeling

# ASSESS UNIVERSITY STRENGTH
uv0 <- vect(dfuniv[dfuniv$Rank<=500,], 
            geom=c("long", "lat"), crs="EPSG:4326")
# prep data
cross_all = prep_k_data(uv0, gbif, sp2)
# run cross
r_vals <- seq(0, 100, by = 2.5)
pcf_cross <- pcfcross.inhom(cross_all, i = "univ", j = "obs", r = r_vals)

set.seed(123)  # reproducibility
envelope_cross <- envelope(cross_all,
                           fun = pcfcross.inhom,
                           i = "univ",
                           j = "obs",
                           nsim = 99,
                           simulate = expression(rlabel(cross_all)),  
                           correction = "translation",
                           savefuns = TRUE,
                           verbose = FALSE,
                           r = r_vals)

# PLOT
plot(envelope_cross, main = "Cross-PCF with Simulation Envelopes")
head(as.data.frame(envelope_cross)[,c('r','obs','hi')], 15)

# Universities
# Below top 100 not significant
# Top 100 - 0-25 km was significant (up to 300K times more likely )
# Top 500 - 0-25 km significant with x1328 times likelehood of observing
#.    again 45-100 km up to 96x
# Top 1000 & 2000 not significant

# # # # # # # # # # # # # # # # # # #

# ASSESS POPULATION STRENGTH
pop0 <- pop[pop$global_ppp_2022_1km_UNadj_constrained>10000,]
# prep data
cross_all_p = prep_k_data(pop0, gbif, sp2)
# run cross
r_vals <- seq(0, 100, by = 2.5)
pcf_cross_p <- pcfcross.inhom(cross_all_p, i = "univ", j = "obs", r = r_vals)

set.seed(123)  
envelope_cross_p <- envelope(cross_all_p,
                             fun = pcfcross.inhom,
                             i = "univ",
                             j = "obs",
                             nsim = 99,
                             simulate = expression(rlabel(cross_all_p)),
                             correction = "translation",
                             savefuns = TRUE,
                             verbose = FALSE,
                             r = r_vals)

plot(envelope_cross_p, main = "Cross-PCF with Simulation Envelopes")
head(as.data.frame(envelope_cross_p)[,c('r','obs','hi')], 15)

# POPULATION
# Repulsion observed on global scale for 10k-50k population density
