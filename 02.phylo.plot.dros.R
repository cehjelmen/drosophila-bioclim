library(spocc)
library(maps)
library(ggplot2)
library(ape)
library(phytools)
library(viridis)
library(rgbif)
library(ggExtra)

#read in list of species
gsdat<-read.csv("data/climate_data_drosophila_Sept17.csv", as.is=T)
gs<-gsdat$Species
####Getting geographic info from gbif####
#we will make a large table of geographic info
#starting with first species out of the loop, then loop through the rest to complete table
i<-1
test<-occ(query=(gs[i]), from='gbif')

#get dataframe of lat and longitude
latlong<-data.frame(test$gbif$data[[1]]$latitude, test$gbif$data[[1]]$longitude, test$gbif$data[[1]]$datasetKey)
#name lat long columns
colnames(latlong)<-c("lat", "long", "dataset_key")
#make a matrix
latlong.mat<-as.matrix(latlong)
#get species name for each
specname<-rep((gs[i]), length.out=length(test$gbif$data[[1]]$latitude))
#make species name a matrix
specname<-as.data.frame(specname)
#name rows species names
rownames(latlong.mat)<-specname$specname
output.latlong<-latlong.mat
#loop of species
for(i in 2:length(gs)){
  print(i)
  test<-occ(query=(gs[i]), from='gbif')
  #get dataframe of lat and longitude
  if(test$gbif$meta$found == 0){
    (i<-i+1)
  }
  else if(is.na(test$gbif$data[[1]]$latitude[1])>0){
    (i<-i+1)
    }

else{
    latlong<-data.frame(test$gbif$data[[1]]$latitude, test$gbif$data[[1]]$longitude, test$gbif$data[[1]]$datasetKey)
  #name lat long columns
  colnames(latlong)<-c("lat", "long", "dataset_key")
  #make a matrix
  latlong.mat<-as.matrix(latlong)
  #get species name for each
  specname<-rep((gs[i]), length.out=length(test$gbif$data[[1]]$latitude))
  #make species name a matrix
  specname<-as.data.frame(specname)
  #name rows species names
  rownames(latlong.mat)<-specname$specname
  #append matrix
  output.latlong<-rbind(output.latlong,latlong.mat)
}
}
str(output.latlong)

#getting citation information for records
i<-1
output.latlong.ref<-matrix(, nrow=length(output.latlong[,1]), ncol=4)
colnames(output.latlong.ref)<-c("lat", "long", "datasetKey", "citation")
rownames(output.latlong.ref)<-rownames(output.latlong)
for(i in 5062:length(output.latlong[,1])){
  print(i)
  output.latlong.ref[i,1]<-output.latlong[i,1]
  output.latlong.ref[i,2]<-output.latlong[i,2]
  output.latlong.ref[i,3]<-output.latlong[i,3]
  ref.info<-gbif_citation(x=output.latlong[i,3])
  output.latlong.ref[i,4]<-ref.info$citation$citation
}
#11022 records  
write.csv(output.latlong.ref, "data/output/lat_long_ref.csv")
locale.output<-data.frame(row.names(output.latlong), output.latlong[,1], output.latlong[,2])
colnames(locale.output)<-c("Species", "lat", "long")

#write it as a csv, but were going to read it as phydat first!
phy.dat<-locale.output
write.csv(phy.dat, "data/output/all.locale.gbif.csv")
#phy.dat<-read.csv("data/output/all.locale.gbif.csv")

####pulling in tree for plot####
mytree<-read.tree("data/mytree2.tre")

#remove na's
phy.dat.na<-na.omit(phy.dat)
#phy.dat$X<-NULL
colnames(phy.dat.na)<-c("Species", "lat", "long")
#species names as factor
phy.dat.na$Species<-as.factor(phy.dat.na$Species)
length(unique(phy.dat.na$Species))
#keep only species with data in the tree
mytree2<-keep.tip(mytree, tip=levels(phy.dat.na$Species))
#make a dataframe of just lat and long
foo<-data.frame(phy.dat.na$lat, phy.dat.na$long)
colnames(foo)<-c("lat", "long")
foo$lat<-as.numeric(foo$lat)
foo$long<-as.numeric(foo$long)
#maek into a matrix
testagain2<-as.matrix(foo)
#name rows after species
rownames(testagain2)<-phy.dat.na$Species

#####All Species Phylo to Map####
#make map object
dros.phymap<-phylo.to.map(mytree2, testagain2, plot=FALSE)
#general plot is crazy
#plot(dros.phymap, ftype="i", cols="red")
#modified later to be less so

#####Read in file (if you haven't already) to get subgenus info and color info####
#read in file with subgenus breakup
#gsdat<-read.csv("data/alldata.csv")
str(gsdat)
gsdat$Subgenus<-as.factor(gsdat$Subgenus)
#make list of speices which are sophophora
soph<-gsdat$Species[gsdat$Subgenus=="Sophophora"]

#list of species which are Drosophila subgenus
dros<-gsdat$Species[gsdat$Subgenus=="Drosophila"]
#interesection of species which are sophophora and in occurence data
soph2<-intersect(phy.dat.na$Species, soph)
#subset by sophophora species
soph.phy.da<-phy.dat.na[phy.dat.na$Species %in% soph2,]

#do the same for Drosophila subgenus
dros2<-intersect(phy.dat.na$Species, dros)
dros.phy.da<-phy.dat.na[phy.dat.na$Species %in% dros2,]

#make matrix files for phymap Sophophora
soph.da<-data.frame(soph.phy.da$lat, soph.phy.da$long)
colnames(soph.da)<-c("lat", "long")
soph.da$lat<-as.numeric(soph.da$lat)
soph.da$long<-as.numeric(soph.da$long)
#maek into a matrix
soph.da2<-as.matrix(soph.da)
#name rows after species
rownames(soph.da2)<-soph.phy.da$Species

drossub.da<-data.frame(dros.phy.da$lat, dros.phy.da$long)
colnames(drossub.da)<-c("lat", "long")
str(drossub.da)
drossub.da$lat<-as.numeric(drossub.da$lat)
drossub.da$long<-as.numeric(drossub.da$long)
#maek into a matrix
drossub.da2<-as.matrix(drossub.da)
#name rows after species
rownames(drossub.da2)<-dros.phy.da$Species


#sophtree
sophtree<-keep.tip(mytree2, intersect(phy.dat.na$Species, soph))
#drostree
drostree<-keep.tip(mytree2, intersect(phy.dat.na$Species, dros))
#set colors for plot soph
#sophcols<-setNames(sample(viridis(n=Ntip(sophtree))),
#               sophtree$tip.label)

#####Colors for Phylomaps####

####Colors by subgenus and GS####
#dros gs data for colors
dros.gs.da<-gsdat[gsdat$Species %in% dros2,]
g9<-ggplot(dros.gs.da, aes(x=Haploid_Number, y=MbDNA_Female, col=MbDNA_Female))+
  geom_point()+
  scale_color_viridis()

p<-ggplot_build(g9)
color.list<-p$data[[1]]["colour"]
#put the list of colors which correspond here
droscoltest<-setNames(color.list$colour, dros.gs.da$Species)

soph.gs.da<-gsdat[gsdat$Species %in% soph2,]
g10<-ggplot(soph.gs.da, aes(x=Haploid_Number, y=MbDNA_Female, col=MbDNA_Female))+
  geom_point()+
  scale_color_viridis()

p2<-ggplot_build(g10)
color.list2<-p2$data[[1]]["colour"]
#put the list of colors which correspond here
sophcoltest<-setNames(color.list2$colour, soph.gs.da$Species)

#colors for full tree
all.spec.da<-gsdat[gsdat$Species %in% mytree2$tip.label,]
g11<-ggplot(all.spec.da, aes(x=Haploid_Number, y=MbDNA_Female, col=MbDNA_Female))+
  geom_point()+
  scale_color_viridis()
p3<-ggplot_build(g11)
color.list3<-p3$data[[1]]["colour"]

#put the list of colors which correspond here
allcoltest<-setNames(color.list3$colour, all.spec.da$Species)

####Plot of trees with colors by GS####
drossub.phy.map<-phylo.to.map(drostree, drossub.da2, plot=FALSE)
sophsub.phy.map<-phylo.to.map(sophtree, soph.da2, plot=FALSE)
#dros subgenus tree
plot(drossub.phy.map, direction="rightwards",
     fsize=0.5, ftype="i", asp=1, lty="solid", map.bg="grey",
     lwd=0.4, pts=FALSE, colors=droscoltest, cex.points=0.4, delimit_map=TRUE)
#soph subgenus tree
plot(sophsub.phy.map, direction="rightwards",
     fsize=0.5, ftype="i", asp=1, lty="solid", map.bg="grey",
     lwd=0.4, pts=FALSE, colors=sophcoltest, cex.points=0.4, delimit_map=TRUE)
#all species tree
plot(dros.phymap, fsize=0.5,direction="rightwards",
     ftype="i", asp=1, lty="solid", map.bg="grey",
     lwd=0.4, pts=FALSE, colors=allcoltest, cex.points=0.4, delimit_map=TRUE)

####Making map of locality by GS with no tree
####This uses info from the all species tree color plot
newinfo<-matrix(,nrow=length(testagain2[,1]), ncol=4)
phy.ma.gs<-as.matrix(phy.dat.na)
#head(phy.ma)
for(i in 1:length(testagain2[,1])){
  newinfo[i,1]<-phy.ma.gs[i,1]
  newinfo[i,2]<-all.spec.da$MbDNA_Female[all.spec.da$Species==phy.ma.gs[i,1]]
  newinfo[i,3]<-phy.ma.gs[i,2]
  newinfo[i,4]<-phy.ma.gs[i,3]
}
#make that matrix into a datafram with labels
newinfo2<-as.data.frame(newinfo)
names(newinfo2)<-c("Species", "GS", "lat", "long")
#speices as factor and other information as numeric
newinfo2$Species<-as.factor(newinfo2$Species)
newinfo2$GS<-as.numeric(newinfo2$GS)
newinfo2$lat<-as.numeric(newinfo2$lat)
newinfo2$long<-as.numeric(newinfo2$long)
str(newinfo2)
write.csv(newinfo2, "data/output/dros_gs_occur.csv")
#ggplot of gs by latitude
cor.gs<-ggplot(newinfo2, aes(x=abs(lat), y=GS, color=GS))+
  geom_point(alpha=0.4)+
  scale_color_viridis(end=0.9, name="Genome Size \n(Mbp)")+
  theme_minimal()+
  geom_smooth(method = "lm", col="red")+
  xlab("Distance from Equator (latitude)")+
  ylab("Genome Size (Mbp)")+
  theme(axis.title=element_text(face="bold", size=14),
        axis.text=element_text(size=12, color="black"),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10),
        legend.key.width = unit(3, "cm"),
        legend.position = "bottom")
ggMarginal(cor.gs, type="histogram", fill="grey90")
summary(lm(GS~abs(lat), data=newinfo2))
#load worldmap
worldmap<-map_data("world")

#plot of occurence on map colored by genome size
#gs.map<-
  
  ggplot()+
  geom_polygon(data=worldmap, aes(x=long, y=lat, group=group),
               fill="grey45",alpha=0.5,
               color="grey40",
               linewidth=0.12)+
  geom_point(data=newinfo2, aes(x=long, y=lat, color=GS),alpha=0.6, size=0.8)+
  coord_fixed(ratio=1.3, xlim=c(-160,-70),
              ylim=c(20,50))+
  scale_x_continuous(breaks=c(-120, -60, 0, 60, 120))+
  scale_y_continuous(breaks=c(-60,-40,-20,0,20,40,60,80))+
  scale_color_viridis_c(end=0.9, name="GS (Mbp)")+
  theme_minimal()+
  ylab("Latitude")+
  xlab("Longitude")+
  ggtitle(expression("GBIF observations"))+
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        panel.grid.minor = element_blank(),
        legend.key.size = unit(0.75, 'cm'))


r1<-read.csv("data/R1_Univ.csv")
us.dat<-map_data("state")
ggplot()+
  geom_polygon(data=us.dat, aes(x=long, y=lat, group=group),
               fill="grey70", alpha=0.5,
               color="grey20",
               linewidth=0.12)+
  coord_fixed(ratio=1.3)+
  ylim(25,50)+xlim(-125,-65)+
  geom_point(data=newinfo2, aes(x=long, y=lat, fill="GBIF Record"), 
             shape=21, size=2, alpha=0.7)+
  geom_point(data=r1, aes(x=Long, y=Lat, fill="R1 University"),
             shape=24,  size=2, alpha=0.7)+
  theme_minimal()+
  scale_fill_manual(values = c("slateblue", "forestgreen"), 
                     labels = c("GBIF Record", "R1 University"))+
  ylab("Latitude")+xlab("Longitude")+
  ggtitle("GBIF Occurrences vs. R1 Universities")+
  theme(axis.text=element_text(size=12, face="bold"),
        axis.title=element_text(size=14, face="bold"),
        plot.title=element_text(size=18, face="bold"), 
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size=14, face="bold"))+
  guides(fill = guide_legend(override.aes = list(size = 5)))  # Adjust the size here



