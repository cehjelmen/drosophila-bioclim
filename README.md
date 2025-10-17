# drosophila-bioclim
#### Scripts and data necessary to complete analysis of genome size versus bioclimatic variables in Drosophila species


## 01.consensus_analyses.R 

Purpose:  PIC and linear model of consensus tree and mined data

## 02.PIC.Figure.R

Purpose:  PIC of consensus tree and figure developement for panelled PIC figure

## 03.phylo.plot.dros.R

Purpose: Gets geographic information from GBIF, Makes phylomaps, Plots locales on worldmap, Makes scatterplot of GS by latitude, Makes plot of R1 institutions vs. localities

## 04.WorldClimdata-dros-code.R

Purpose: Gets bioclimatic variables for each location for species

## 05.looping_climate_analyses.R

Purpose: PICs for bioclimatic variables vs GS in Drosophila, Uses 100 trees randomly selected from bayes distribution, Loops through each tree 100 times with data randomly selected for species with more than one record

## 06.WorldClim_dros_spatial.R

Purpose: Spatial overlay for Drosophila distribution based on species temperature tolerance in GFIB and WorldClim v2.1 data.

## 07.WorldClim_dros_figures.R

Purpose: Visualizations of spatial overlay outputs including: 
  1. Global species heatmap
  2. Species MbDNA as a function of spatial distribution
  3. Total distribution by continent
  4. Species distribution by continent (excluding Antarctica and Oceania)
  5. Species distribution by eco-regions (excluding Antarctic and Oceania)
  6. Biogeographic distribution count (excluding Antarctic and Oceania)

## 08.drosophila_maps.R

Purpose: Visualizations for plot 1 in `07.WorldClim_dros_figures.R`. Uses thermal tolerance map figure and locality info from GBIF to compare distributions (figure 1)


## 09.all_gbif_dros_map_.R

Purpose: download all gbif occurence records for the search term "Drosophila". Saves csv with species, lat, long, reference, and datasetkey in output folder of data folder


## 10.research_university_mapping.R

Purpose: uses list of top research universities to get lat-long information and then plots locations of universities on world map

## 11. point_pattern_analysis

Purpose: uses multi type pair correlation function (pdf) to explore the spatial relationship between universities, population density, and GBIF observations