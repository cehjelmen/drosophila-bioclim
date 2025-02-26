#This script pulls in data from rankings of research universities from SCImago
#this searches up the geographic locations of universities through a search
#this also generates a plot comparing distributions of university locations
####Load Packages####
library(maps)
library(ggplot2)
library(viridis)
#install.packages("tidygeocoder")
library(tidygeocoder)
library(tidyverse)
library(ggpubr)

####Read in university data####
#this script and search was previously run, so you have the option to pull in the data
#from the output data folder and not have to run this whole search again
#see step further below before running the script!

#you have to read as csv2 to account for semicolon data
#this data was downloaded as a csv on Feb 25, 2025 at
#https://www.scimagoir.com/rankings.php?sector=Higher+educ.&ranking=Research 
univ<-read.csv2("data/university_info/ScimagoIR 2024 - Research Rank - Universities.csv")

#make a dataframe of university names
#this pulls the first/top 1000 research universities
univ2<-data.frame(univ$Institution)
#this pulls the lat and long data for each name when searched
universities_geo <- univ2 %>%
  geocode(address = univ.Institution, method = "osm")

#this adds a ranking component to the university, as the original document
#is in order of research university rank
univ.rank<-data.frame(universities_geo, c(1:length(universities_geo$univ.Institution)))
colnames(univ.rank)[4]<-"Rank"

#save output of data for future use. If you don't want to run this whole thing again
#pull in the file from output data
write.csv(univ.rank, "data/output/university_loc_rankings.csv")
#read this file as univ.rank in if you want to not run the whole search
#univ.rank<-read.csv("data/output/university_loc_rankings.csv")

#separate out ranking clusters
#ranks 1-500
univ.rank500<-univ.rank[1:500,]
univ.rank1000<-univ.rank[1:1000,]
univ.rank1500<-univ.rank[1:1500,]
univ.rank2000<-univ.rank[1:2000,]

####Plotting University Data####
#load worldmap
worldmap<-map_data("world")
#plot of occurrences on map 

####function to map in ggplot map based on your dataset####
rank.mapping<-function(univ.dat){
  ggplot()+
  geom_polygon(data=worldmap, 
               aes(x=long, y=lat, group=group),
               fill="grey45",alpha=0.2,
               color="grey40",
               linewidth=0.12)+
  geom_point(data=univ.dat, 
             aes(x=long, y=lat,color=Rank),
             alpha=0.7, size=1, 
             )+
  scale_x_continuous(breaks=c(-120, -60, 0, 60, 120))+
  scale_y_continuous(breaks=c(-60,-40,-20,0,20,40,60,80))+
  scale_color_viridis_c(end=0.9,option="B",
    breaks = c(1, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500),  # Manually define breaks
    labels = c("1", "500", "1000", "1500", "2000", "2500", "3000", "3500", "4000", "4500")  # Ensure "1" is labeled
  )+ 
  theme_minimal()+
  ylab("Latitude")+
  xlab("Longitude")+
  ggtitle(expression("Top Ranked Research Universities"))+
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        panel.grid.minor = element_blank(),
        legend.key.size = unit(0.75, 'cm'))+
  guides(color = guide_colorbar(reverse = TRUE))
}


####plotting each set of rankings####
top500<-rank.mapping(univ.rank500)
top1000<-rank.mapping(univ.rank1000)
top1500<-rank.mapping(univ.rank1500)
all<-rank.mapping(univ.rank)

#multiple arrangement of plots
ggarrange(top500, top1000, top1500, all, nrow=2, ncol=2)

#saving as a pdf
pdf(file="data/output/research_univ_map.pdf", width=12.19, height=6.83)
ggarrange(top500, top1000, top1500, all, nrow=2, ncol=2)
dev.off()
