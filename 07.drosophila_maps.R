#
# Area distribution based on species temperature tolerance from WorldClim v2.1
# Drosophila Project - Hjelman & Curnow
# Map figures only - Olson 10/2024
#
# # # # # # # # # # # # # # # # # # # # # # # # # #

# load world political data
library(sf);library(terra);library(dplyr)
library(ggplot2);library(gridExtra);library(tidyterra)

# load output data for plotting
results     <- "./data/mapping"
dfgs        <- read.csv("data/output/dros_gs_occur.csv") # 2-4-25 DATA
world       <- read_sf(file.path(results, "world-simple.geojson"))
dros_heat   <- rast(file.path(results, "drosophila_heatmap_n67_10m_1.grd"))

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
    tag = "a",
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
    tag = "b", 
    fill = "Species\nOccurences",
    title = "Temperature tolerance") +
  ylim(c(-60,85)) + theme_minimal() +
  theme(plot.title = element_text(size = 10, face = "bold"))
# adjust margins to align
p1_adjusted <- p1 + theme(plot.margin = margin(t = 0, r = 18, b = 10, l = 0))  
p2_adjusted <- p2 + theme(plot.margin = margin(t = 0, r = 0, b = 10, l = 5)) 
## plot combined with layout adjustment
grid.arrange(p1_adjusted, p2_adjusted, ncol = 1, heights = c(1, 1) )
