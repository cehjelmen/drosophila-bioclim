#
# Spatial distribution based on species temperature tolerance from WorldClim v2.1
# Drosophila Project - Hjelman & Curnow
# Figures - Olson 10/2024
#
# # # # # # # # # # # # # # # # # # # # # # # # # #

# load world political data
library(sf);library(terra);library(dplyr)
library(ggplot2);library(gridExtra);library(tidyterra)

# define folder with results (if saved)
results     <- "data/dros_spatial" 

## load output data for visualization
# OUTPUT files from 05.WorldClim_dros_spatial
df          <- read.csv(file.path(results, "climate_data_drosophila_Sept17.csv"))
dfgs        <- read.csv(file.path(results, "dros_gs_occur_june.csv"))
world       <- read_sf(file.path(results, "world-simple.geojson"))
dros_heat   <- rast(file.path(results, "drosophila_heatmap_n67_10m_1.grd"))
df_comb     <- read.csv(file.path(results, "drosophila_area_table.csv"))
df_cont     <- read.csv(file.path(results, "drosophila_continent_table.csv"))
df_eco      <- read.csv(file.path(results, "drosophila_eco_table.csv"))
eco_regions <- st_read(file.path(results, "eco-realm-summary.geojson"))

# # # # # # # # # # # # # # # # # #
# Figure 1. global species heatmap
world %>% ggplot() + 
  geom_spatvector(fill = alpha("grey40", 0.5) ) +
  geom_spatraster(data = dros_heat) +
  scale_fill_whitebox_c(
    palette = "muted", # can use "viridi"
    n.breaks = 10,
    guide = guide_legend(reverse = TRUE)) +
  labs(
    fill = "Counts",
    title = "Drosophila Distribution (n = 67)",
    subtitle = "based on temperature tolerance") +
  ylim(c(-60,85)) + theme_minimal()

# # # # # # # # # # # # # # # # # #
# Figure 1.1 global species w/ heatmap heatmap
# point map
p1 <- world %>% ggplot() + 
  geom_spatvector(fill = alpha("grey40", 0.5)) +
  geom_point(data = dfgs, aes(x = long, y = lat, colour = GS), size = 0.5) +
  labs(
    tag = "A",
    fill = "Species Occurences",
    title = "GBIF observations") +
  scale_color_viridis_c() + 
  ylim(c(-60,85)) + theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 10, face = "bold"))
# heat map
p2 <- world %>% ggplot() +
  geom_spatvector(fill = alpha("grey40", 0.5) ) +
  geom_spatraster(data = dros_heat) +
  scale_fill_whitebox_c(
    palette = "muted", 
    n.breaks = 8,
    guide = guide_legend(reverse = TRUE)) +
  labs(
    tag = "B", 
    fill = "Species\nOccurences",
    title = "Temperature tolerance") +
  ylim(c(-60,85)) + theme_minimal() +
  theme(plot.title = element_text(size = 10, face = "bold"))
# adjust margins to align
p1_adjusted <- p1 + theme(plot.margin = margin(t = 0, r = 18, b = 10, l = 0))  
p2_adjusted <- p2 + theme(plot.margin = margin(t = 0, r = 0, b = 10, l = 5)) 
## plot combined with layout adjustment
grid.arrange(p1_adjusted, p2_adjusted, ncol = 1, heights = c(1, 1) )

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Figure 2. MbDna female vs temperature-based area distribution (WorldClim)
df_comb %>% ggplot(aes(x = area_mkm2, y = MbDNA_Female, colour = MbDNA_Male)) +
  labs(
    title = "Species MbDNA as a function of spatial distribution",
    subtitle = "According to WorldClim v2.1",
    colour = "MbDNA (Male)",
    y = "MbDNA (Female)",
    x = bquote("Species distribution (million km"^2*")")
  ) +
  geom_point() + theme_minimal() +
  scale_color_viridis_c()

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Figure 2.1 MbDna female vs temperature-based area distribution (WorldClim)
## plot with trend & metrics
lm_eqn <- function(df, x, y){ # function - plot standard metrics
  m <- lm(y ~ x, df)
  pval <- summary(m)$coefficients[2, 4]
  r2 <- summary(m)$r.squared
  eq <- substitute( italic(y) == a ~italic(x) + b * "," ~ italic(pvalue) ~ "=" ~ pval * "," ~ italic(r)^2 ~ "=" ~ r2, 
                    list(a = format(unname(coef(m)[2]), digits = 2),
                         b = format(unname(coef(m)[1]), digits = 2),
                         pval = format(pval, digits = 3),
                         r2 = format(r2, digits = 2) ) )
  as.character(as.expression(eq))
}
data.label <- data.frame( # create label
  x = 25, y = 300,
  label = c(lm_eqn(df_comb, x = df_comb$area_mkm2, y = df_comb$MbDNA_Female))
)
# plot
df_comb %>% ggplot(aes(x = area_mkm2, y = MbDNA_Female)) +
  labs(
    title = "Species MbDNA as a function of spatial distribution",
    subtitle = "based on temperature threshold and WorldClim v2.1",
    colour = "MbDNA (Male)",
    y = "MbDNA (Female)",
    x = bquote("Species distribution (million km"^2*")")
  ) + 
  geom_smooth(method = "lm", color="red", linetype = "dashed", size = 0.5, formula = y ~ x) +
  geom_text(data = data.label, aes(x = x , y = y , label = label), size=4, parse = TRUE) +
  geom_point() + theme_minimal() +
  scale_color_viridis_c()


# # # # # # # # # # # # # # # # # # # # # #
# Figure 3. total distribution by continent
df_cont %>% 
  pivot_longer(cols = !Species,
               names_to = "Continent",
               values_to = "values") %>% 
  filter(Continent != "Antarctica") %>% 
  mutate(values = values/1e6) %>% 
  group_by(Continent) %>% summarize(values = sum(values)) %>% 
  ggplot(aes(x = Continent, y = values)) +
  labs(
    x = "", y = bquote("Area distribution (million"~km^2~")"),
    title = "Drosophila spatial distribution by continent"
  ) +
  geom_bar(stat="identity") + 
  theme_classic()

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Figure 4. species distribution by continent (excluding Antarctica and Oceania)
left_join(df_cont, df %>% select(Species, Subgenus, MbDNA_Male, MbDNA_Female), by='Species') %>%
  select(-c(Antarctica, Oceania)) %>% 
  pivot_longer(cols = !c(Species,Subgenus, MbDNA_Male, MbDNA_Female),
               names_to = "Continent",
               values_to = "values") %>%  
  mutate(values = values/1e6) %>%
  ggplot(aes(x = values, y = MbDNA_Female, colour=Subgenus)) +
  labs(
    y = "MbDNA (Female)", x = bquote("Area distribution (million"~km^2~")"),
    title = "Drosophila genome size and continental spatial distribution"
  ) +
  guides(fill = guide_legend(show=F)) +
  geom_point() + facet_wrap(.~Continent) +
  theme_minimal() + scale_color_brewer(palette = "Set2")


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Figure 5. species distribution by ecoregions (excluding Antarctic and Oceania)
df_eco %>% 
  pivot_longer(cols = !Species,
               names_to = "Biogeographic",
               values_to = "values") %>% 
  mutate(values = values/1e6) %>% 
  filter(Biogeographic != "Antarctic", Biogeographic != "Oceania" ) %>% 
  group_by(Biogeographic) %>% summarize(values = sum(values)) %>% 
  ggplot(aes(x = Biogeographic, y = values)) +
  labs(
    x = "", y = bquote("Area distribution (million"~km^2~")"),
    title = "Drosophila biogeographic distributions"
  ) +
  geom_bar(stat="identity") + 
  theme_classic()

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Figure 6. biogeographic distribution count (excluding Antarctic and Oceania)
m1 = left_join(df_eco, df_comb %>% select(Species, Subgenus, MbDNA_Male, MbDNA_Female), by='Species') %>% 
  pivot_longer(cols = !c(Species,Subgenus, MbDNA_Male, MbDNA_Female),
               names_to = "Biogeographic",
               values_to = "values") %>% 
  filter(Biogeographic != "Antarctic", Biogeographic != "Oceania", Biogeographic != "X"  ) %>% 
  mutate(Biogeographic = as.factor(Biogeographic)) %>% 
  filter(values > 0) %>% 
  mutate(values = values/1e6) %>% 
  ggplot(aes(x = values, y = MbDNA_Female, colour=Biogeographic)) +
  labs(title = "Drosophila distribution and genome size",
       subtitle = "by biogeographic region", y = "MbDNA (Female)",
       x = bquote("Area by temperature tolerance (million"~km^2~")")) +
  scale_colour_brewer(palette = "Dark2", guide = 'none') +
  geom_point() + facet_wrap(.~Biogeographic) +
  theme(legend.text=element_text(size=rel(0.5))) +
  theme_minimal()
m2 = eco_regions %>%
  filter(WWF_REALM2 != "Antarctic", WWF_REALM2 != "Oceania" ) %>% 
  st_transform(6933) %>% 
  ggplot() + 
  labs(x = "",y="") +
  geom_sf(aes(fill = WWF_REALM2)) + 
  guides(fill=guide_legend(title="Species count (n=67)")) +
  scale_fill_brewer(palette = "Dark2") +
  geom_sf_label(aes(label = species_max, alpha = 0.7))  +
  scale_alpha(guide = 'none') +
  theme_minimal() 
## plot combined
grid.arrange(m1,m2)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Figure 6.1 biogeographic distribution count (excluding Antarctic and Oceania)
## plot with trend
# join and pivot data
df_eco2 <- left_join(df_eco, df_comb %>% select(Species, Subgenus, MbDNA_Male, MbDNA_Female), by='Species') %>% 
  pivot_longer(cols = !c(Species,Subgenus, MbDNA_Male, MbDNA_Female),
               names_to = "Biogeographic",
               values_to = "values") %>% 
  filter(Biogeographic == "Afrotropic" | Biogeographic == "Neotropic") %>% 
  mutate(Biogeographic = as.factor(Biogeographic)) %>% 
  filter(values > 0) %>% 
  mutate(values = values/1e6)
# create labels
data.label <- df_eco2 %>%
  group_by(Biogeographic) %>%
  summarise(
    x = 10,  
    y = 310,  
    label = lm_eqn(., values, MbDNA_Female)
  )
# plot
df_eco2 %>% 
  ggplot(aes(x = values, y = MbDNA_Female, colour=Biogeographic)) +
  labs(title = "Drosophila distribution and genome size",
       subtitle = "by biogeographic region", y = "MbDNA (Female)",
       x = bquote("Area by temperature tolerance (million"~km^2~")")) +
  scale_colour_brewer(palette = "Dark2", guide = 'none') +
  geom_point() + facet_wrap(.~Biogeographic) +
  geom_smooth(method = "lm", color="red", linetype = "dashed", size = 0.5, formula = y ~ x) +
  geom_text(data = data.label, aes(x = x , y = y , label = label), size=3.2, parse = TRUE) +
  theme(legend.text=element_text(size=rel(0.5))) +
  theme_minimal()

## END