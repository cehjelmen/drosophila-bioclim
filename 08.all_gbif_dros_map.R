#This script is to download information for all GBIF localities that include the 
#term "Drosophila"
####Load packages####
library(spocc)
library(maps)
library(ggplot2)
library(viridis)
library(rgbif)
library(ggExtra)
library(dismo)
library(rgbif)

####do the search on gbif####
#this looks for a max of 100,000 records, but limits the data to those that include
#coordinates
drosophila_data <- occ_search(
  taxonKey = 1522683,
  limit = 100000, 
  hasCoordinate = TRUE
)
#has 968 unique species names
length(unique(drosophila_data$data$scientificName))
#and 63,320 records
length(drosophila_data$data$scientificName)

####making a dataset with all the relevant information
#long
#lat
#species
#ref

drosophila_data$data[1]

dros.mat<-matrix(,nrow=length(drosophila_data$data$scientificName), ncol=6)
colnames(dros.mat)<-c("species", "lat", "long", "ref", "gbifID","datasetKey")

for(i in 1:length(drosophila_data$data$scientificName)){
  dros.mat[i,1]<-drosophila_data$data$species[i]
  dros.mat[i,2]<-drosophila_data$data$decimalLatitude[i]
  dros.mat[i,3]<-drosophila_data$data$decimalLongitude[i]
  dros.mat[i,4]<-drosophila_data$data$references[i]
  dros.mat[i,5]<-drosophila_data$data$gbifID[i]
  dros.mat[i,6]<-drosophila_data$data$datasetKey[i]
}

write.csv(dros.mat, "data/output/all.drosophila.gbif.csv")

####Making map of locality by GS with no tree
#make that matrix into a datafram with labels
dros.df<-as.data.frame(dros.mat)
dros.df$lat<-as.numeric(dros.df$lat)
dros.df$long<-as.numeric(dros.df$long)

#load worldmap
worldmap<-map_data("world")
#plot of occurrences on map 
Dros.map<-ggplot()+
  geom_polygon(data=worldmap, aes(x=long, y=lat, group=group),
               fill="grey45",alpha=0.2,
               color="grey40",
               linewidth=0.12)+
  geom_point(data=dros.df, aes(x=long, y=lat),alpha=0.4, size=0.8, 
             fill="maroon2", pch=21)+
  #coord_fixed(ratio=1.3, xlim=c(-160,-70),
  #            ylim=c(20,50))+
  scale_x_continuous(breaks=c(-120, -60, 0, 60, 120))+
  scale_y_continuous(breaks=c(-60,-40,-20,0,20,40,60,80))+
  #scale_color_viridis_c(end=0.9, name="GS (Mbp)")+
  theme_minimal()+
  ylab("Latitude")+
  xlab("Longitude")+
  ggtitle(expression("All GBIF observations"))+
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        panel.grid.minor = element_blank(),
        legend.key.size = unit(0.75, 'cm'))

#####Save a pdf of map####
pdf("data/output/all_gbif_map.pdf", width=8.81, height=5.54)
Dros.map
dev.off()
