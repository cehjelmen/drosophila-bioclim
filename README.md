# drosophila-bioclim
#### Scripts and data necessary to complete analysis of genome size versus bioclimatic variables in Drosophila species


## 01.Script is for Sam's PIC 

Purpose:  PIC of consensus tree and figure developement for panelled PIC figure

## 02.phylo.plot.dros.R

Purpose: Gets geographic information from GBIF, Makes phylomaps, Plots locales on worldmap, Makes scatterplot of GS by latitude, Makes plot of R1 institutions vs. localities

## 03.WorldClimdata-dros-code.R

Purpose: Gets bioclimatic variables for each location for species

## 04.looping_climate_analyses.R

Purpose: PICs for bioclimatic variables vs GS in Drosophila, Uses 100 trees randomly selected from bayes distribution, Loops through each tree 100 times with data randomly selected for species with more than one record

## 05.WorldClim_dros_spatial.R

Purpose: Spatial overlay for Drosophila distribution based on species temperature tolerance in GFIB and WorldClim v2.1 data.

## 06.WorldClim_dros_figures.R

Purpose: Visualizations of spatial overlay outputs including: 
  1. Global species heatmap
  2. Species MbDNA as a function of spatial distribution
  3. Total distribution by continent
  4. Species distribution by continent (excluding Antarctica and Oceania)
  5. Species distribution by eco-regions (excluding Antarctic and Oceania)
  6. Biogeographic distribution count (excluding Antarctic and Oceania)

## 07.drosophila_maps.R

Purpose: Visualizations for plot 1 in `06.WorldClim_dros_figures.R`. Uses thermal tolerance map figure and locality info from GBIF to compare distributions (figure 1)


