# drosophila-bioclim
#### Scripts necessary to complete analysis of genome size versus bioclimatic variables in Drosophila species
###01.Script is for Sam's PIC of conensus tree and figure developement for panelled PIC figure

###02.phylo.plot.dros.R
•	Purpose:
o	Gets geographic information from gbif
o	Makes phylomaps
	By subgenus
	By whole phylogeny
	Plots locales on world map
	Makes scatterplot of GS by latitude
	Makes plot of R1 institutions vs. localities
•	Packages:
o	spocc
o	ape
o	phytools
o	viridis
o	ggplot2
o	maps
o	ggExtra
•	Uses:	
o	“mytree2.tre”
o	“climate_data_drosophila_Setp17”
	Used for the GS and the “groups” of species
o	“R1_Univ.csv”
	Used for R1 vs. locality figure.  GPS coordinates of R1 university cities
•	Generates: 
o	all.locale.gbif.csv (from output.latlong, re-read in as phy.dat)
o	dros_gs_occur.csv (to be used in next code)
o	“lat_long_ref.csv” (references for all the locations pulled off of GBIF)
o	Phyloplots
o	gs vs lat scatterplot
o	R1 university vs. locality figure
![image](https://github.com/user-attachments/assets/bb4ec4b6-26da-46ee-ba45-6d7be93ecd6f)
