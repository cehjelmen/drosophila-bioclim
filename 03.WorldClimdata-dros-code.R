#This script Gets bioclimatic variables for each location for species
####load appropriate libraries####
library(raster)
library(dismo)
library(geodata)
library(dplyr)
resolution <- 2.5  # choose the desired resolution.
#variables that you can use to pull data
#temperatures are celsius *10
# tmean
# tmin
# tmax
# prec
# bio
# alt
#a calculation to make bio from the tmean, tmax, pre values

#summary of bioclimatic variables
# BIO1 = Annual Mean Temperature
# BIO2 = Mean Diurnal Range (Mean of monthly (max temp - min temp))
# BIO3 = Isothermality (BIO2/BIO7) (×100)
# BIO4 = Temperature Seasonality (standard deviation ×100)
# BIO5 = Max Temperature of Warmest Month
# BIO6 = Min Temperature of Coldest Month
# BIO7 = Temperature Annual Range (BIO5-BIO6)
# BIO8 = Mean Temperature of Wettest Quarter
# BIO9 = Mean Temperature of Driest Quarter
# BIO10 = Mean Temperature of Warmest Quarter
# BIO11 = Mean Temperature of Coldest Quarter
# BIO12 = Annual Precipitation
# BIO13 = Precipitation of Wettest Month
# BIO14 = Precipitation of Driest Month
# BIO15 = Precipitation Seasonality (Coefficient of Variation)
# BIO16 = Precipitation of Wettest Quarter
# BIO17 = Precipitation of Driest Quarter
# BIO18 = Precipitation of Warmest Quarter
# BIO19 = Precipitation of Coldest Quarter
#pull in worldclim_data at the resolution you want and the variable you want
worldclim_data <- worldclim_global(var = 'bio', path=tempdir(), res = resolution)

#read in dros locality data
dros<-read.csv("data/output/dros_gs_occur.csv")
#remove redundant number column
dros$X<-NULL
#make empty matrix to fill with data
drosdata<-matrix(,nrow=length(dros$Species), ncol=23)
#name columns for matrix
colnames(drosdata)<-c("Species", "GS", "Latitude", "Longitude", "Bio1","Bio2",
                      "Bio3","Bio4","Bio5","Bio6","Bio7","Bio8","Bio9","Bio10","Bio11",
                      "Bio12","Bio13","Bio14","Bio15","Bio16","Bio17","Bio18","Bio19")
#loop to pull in all the bio climatic data for each locality
i<-1
for(i in 1:length(dros$Species)){
  latitude<-dros$lat[i] #gets the latitude for the row
  longitude<-dros$long[i] #gets the longitude for the row
  location<-cbind(longitude,latitude) #location info for the worldclim data
  bioclim_values<-extract(worldclim_data, location) #pulls bioclimatic data
  bioclim_values<-as.matrix(bioclim_values)
  drosdata[i,1]<-dros$Species[i] #puts name of species in matrix
  drosdata[i,2]<-dros$GS[i] #puts GS for each species in matrix
  drosdata[i,3]<-dros$lat[i] #puts latitude in
  drosdata[i,4]<-dros$long[i] #puts longitude in
  drosdata[i,5:23]<-bioclim_values[1:19] #pastes all 19 bioclimatic variables into the matrix
}

#make dataframe from matrix
drosdata.frame<-as.data.frame(drosdata)

#make certain columns numeric
drosdata.frame <- drosdata.frame %>%
  mutate(across(c(GS, Latitude, Longitude, Bio1, Bio2, Bio3, Bio4, Bio5, Bio6,
                  Bio7, Bio8, Bio9, Bio10, Bio11, Bio12, Bio13, Bio14, Bio15,
                  Bio16, Bio17, Bio18, Bio19), as.numeric))
#make species a factor
drosdata.frame$Species<-as.factor(drosdata.frame$Species)
write.csv(drosdata.frame, "data/output/drosophila_bioclim.csv")
