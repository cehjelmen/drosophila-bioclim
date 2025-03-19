#---Script to use consensus tree to perform PICs for data mined variables--#
####load packages####
library(ape)
library(geiger)
library(caper)
library(ggplot2)
library(ggpubr)

####make sure you working directory is set####
#I have this set to work in our github Repo with the data folder#
#setwd("~/Coding/R/Drosophila Climate Research/Data Drosophila")
dat<-read.csv("data/climate_data_drosophila_Sept17.csv")
tree<-read.tree("data/mytree2.tre")

####Trim datasets to those that have values for species####
#tmin
gstmin<-as.numeric(dat$MbDNA_Female[is.na(dat$TMin)==0])
names(gstmin)<-dat$Species[is.na(dat$TMin)==0]
tmin<-dat$TMin[is.na(dat$TMin)==0]
names(tmin)<-dat$Species[is.na(dat$TMin)==0]
#names to keep
namestokeep<-dat$Species[is.na(dat$TMin)==0]
#drop extra tips
tmintree<-keep.tip(tree, namestokeep)
#plot(tmintree)

#tmax
gstmax<-as.numeric(dat$MbDNA_Female[is.na(dat$TMax)==0])
names(gstmax)<-dat$Species[is.na(dat$TMax)==0]
tmax2<-dat$TMax[is.na(dat$TMax)==0]
names(tmax2)<-dat$Species[is.na(dat$TMax)==0]
#names to keep
namestokeeptmax<-dat$Species[is.na(dat$TMax)==0]
#drop extra tips
tmaxtree<-keep.tip(tree, namestokeeptmax)
#plot(tmaxtree)

#ctmin
gsctmin<-as.numeric(dat$MbDNA_Female[is.na(dat$CTMin_Overall)==0])
names(gsctmin)<-dat$Species[is.na(dat$CTMin_Overall)==0]
ctmin<-dat$CTMin_Overall[is.na(dat$CTMin_Overall)==0]
names(ctmin)<-dat$Species[is.na(dat$CTMin_Overall)==0]
#names to keep
namestokeepctmin<-dat$Species[is.na(dat$CTMin_Overall)==0]
#drop extra tips
ctmintree<-keep.tip(tree, namestokeepctmin)
#plot(tmintree)

#ctmax
gsctmax<-as.numeric(dat$MbDNA_Female[is.na(dat$CTMax_Overall)==0])
names(gsctmax)<-dat$Species[is.na(dat$CTMax_Overall)==0]
ctmax2<-dat$CTMax_Overall[is.na(dat$CTMax_Overall)==0]
names(ctmax2)<-dat$Species[is.na(dat$CTMax_Overall)==0]
#names to keep
namestokeepctmax<-dat$Species[is.na(dat$CTMax_Overall)==0]
#drop extra tips
ctmaxtree<-keep.tip(tree, namestokeepctmax)
#plot(tmaxtree)

#dev time
gsdevtime<-as.numeric(dat$MbDNA_Female[is.na(dat$Dev_Time)==0])
names(gsdevtime)<-dat$Species[is.na(dat$Dev_Time)==0]
devtime2<-dat$Dev_Time[is.na(dat$Dev_Time)==0]
names(devtime2)<-dat$Species[is.na(dat$Dev_Time)==0]
#names to keep
namestokeepdevtime<-dat$Species[is.na(dat$Dev_Time)==0]
#drop extra tips
devtimetree<-keep.tip(tree, namestokeepdevtime)
#plot(devtimetree)

#haploid
gshaploid<-as.numeric(dat$MbDNA_Female[is.na(dat$Haploid_Number)==0])
names(gshaploid)<-dat$Species[is.na(dat$Haploid_Number)==0]
haploid2<-dat$Haploid_Number[is.na(dat$Haploid_Number)==0]
names(haploid2)<-dat$Species[is.na(dat$Haploid_Number)==0]
#names to keep
namestokeephaploid<-dat$Species[is.na(dat$Haploid_Number)==0]
#drop extra tips
haploidtree<-keep.tip(tree, namestokeephaploid)
#plot(haploidtree)

#precipitation
gsprecip<-as.numeric(dat$MbDNA_Female[is.na(dat$Precip)==0])
names(gsprecip)<-dat$Species[is.na(dat$Precip)==0]
precip2<-dat$Precip[is.na(dat$Precip)==0]
names(precip2)<-dat$Species[is.na(dat$Precip)==0]
#names to keep
namestokeepprecip<-dat$Species[is.na(dat$Precip)==0]
#drop extra tips
preciptree<-keep.tip(tree, namestokeepprecip)
#plot(preciptree)



####Performing PICs and LMs for each variable####
#pic for tmin
pic.gstmin<-pic(gstmin, tmintree)
pic.tmin<-pic(tmin, tmintree)

#fit linear model without intercept
fit.pic.tmin<-lm(pic.gstmin ~ pic.tmin +0)
summary(fit.pic.tmin)

#tmax
pic.gstmax<-pic(gstmax, tmaxtree)
pic.tmax<-pic(tmax2, tmaxtree)

fit.pic.tmax<-lm(pic.gstmax ~ pic.tmax +0)
summary(fit.pic.tmax)

#ctmin
pic.gsctmin<-pic(gsctmin, ctmintree)
pic.ctmin<-pic(ctmin, ctmintree)

fit.pic.ctmin<-lm(pic.gsctmin ~ pic.ctmin +0)
summary(fit.pic.ctmin)

#ctmax
pic.gsctmax<-pic(gsctmax, ctmaxtree)
pic.ctmax<-pic(ctmax2, ctmaxtree)

fit.pic.ctmax<-lm(pic.gsctmax ~ pic.ctmax +0)
summary(fit.pic.ctmax)

#dev time
pic.gsdevtime<-pic(gsdevtime, devtimetree)
pic.devtime<-pic(devtime2, devtimetree)

fit.pic.devtime<-lm(pic.gsdevtime ~ pic.devtime +0)
summary(fit.pic.devtime)

#haploid
pic.gshaploid<-pic(gshaploid, haploidtree)
pic.haploid<-pic(haploid2, haploidtree)

fit.pic.haploid<-lm(pic.gshaploid ~ pic.haploid +0)
summary(fit.pic.haploid)

#precipitation
pic.gsprecip<-pic(gsprecip, preciptree)
pic.precip<-pic(precip2, preciptree)

fit.pic.precip<-lm(pic.gsprecip ~ pic.precip +0)
summary(fit.pic.precip)

#pic ctmax
pic.gs<-pic(gsctmax, ctmaxtree)
pic.ctmax<-pic(ctmax2, ctmaxtree)

#fit linear model without intercept ctmax
fit.pic.ctmax<-lm(pic.gs ~ pic.ctmax +0)
summary(fit.pic.ctmax)

#pic time ctmin
pic.gsctmin<-pic(gsctmin, ctmintree)
pic.ctmin<-pic(ctmin, ctmintree)

#fit linear model without intercept ctmin
fit.pic.ctmax<-lm(pic.gs ~ pic.ctmax +0)
summary(fit.pic.ctmax)

#fit linearmodel 
fit.pic.ctmin<-lm(pic.gsctmin ~pic.ctmin+0)
summary(fit.pic.ctmin)



####Making individual plots####
#ggplot (make individual plots)
#PIC for TMin
correctedtmin<-ggplot(, aes(y=pic.gstmin, x=pic.tmin))+
  geom_point(size=3)+
  geom_smooth(method=lm, level=0.95, colour="red")+
  ylab("PICs for Female GS (Mbp)")+
  xlab("PICs for T Min")+
  ggtitle("PICs for T Min")+
  theme(axis.title = element_text(face="bold", size=14),
        axis.text = element_text(size=12),
        title = element_text(size=16, face="bold"))

#Uncorrected TMin
lineartmin<-ggplot(dat, aes(y=MbDNA_Female, x=TMin))+
    geom_point(size=3)+
    geom_smooth(method=lm, level=0.95)+
  ylab("Female GS (Mbp)")+
    xlab("T Min")+
  ggtitle("Uncorrected T Min")+
    theme(axis.title = element_text(face="bold", size=14),
          axis.text = element_text(size=12),
          title = element_text(size=16, face="bold"))
  
#PIC for Tmax
correctedtmax<-ggplot(, aes(y=pic.gstmax, x=pic.tmax))+
  geom_point(size=3)+
  geom_smooth(method=lm, level=0.95, colour="red")+
  ylab("PICs for Female GS (Mbp)")+
  xlab("PICs for T Max")+
  ggtitle("PICs for T Max")+
  theme(axis.title = element_text(face="bold", size=14),
        axis.text = element_text(size=12),
        title = element_text(size=16, face="bold"))

#Uncorrected for Tmax
lineartmax<-ggplot(dat, aes(y=MbDNA_Female, x=TMax))+
  geom_point(size=3)+
  geom_smooth(method=lm, level=0.95)+
  ylab("Female GS (Mbp)")+
  xlab("T Max")+
  ggtitle("Uncorrected T Max")+
  theme(axis.title = element_text(face="bold", size=14),
        axis.text = element_text(size=12),
        title = element_text(size=16, face="bold"))


#PIC for CTMax
correctedctmax<-ggplot(, aes(y=pic.gsctmax, x=pic.ctmax))+
  geom_point(size=3)+
  geom_smooth(method=lm, level=0.95, colour="red")+
  ylab("PICs for Female GS (Mbp)")+
  xlab("PICs for CT Max")+
  ggtitle("PICs for CT Max")+
  theme(axis.title = element_text(face="bold", size=14),
        axis.text = element_text(size=12),
        title = element_text(size=16, face="bold"))

#Uncorrected for CTmax
linearctmax<-ggplot(dat, aes(y=MbDNA_Female, x=CTMax_Overall))+
  geom_point(size=3)+
  geom_smooth(method=lm, level=0.95)+
  ylab("Female GS (Mbp)")+
  xlab("CT Max")+
  ggtitle("Uncorrected CT Max")+
  theme(axis.title = element_text(face="bold", size=14),
        axis.text = element_text(size=12),
        title = element_text(size=16, face="bold"))
  
#PIC for CTMin
correctedctmin<-ggplot(, aes(y=pic.gsctmin, x=pic.ctmin))+
  geom_point(size=3)+
  geom_smooth(method=lm, level=0.95, colour="red")+
  ylab("PICs for Female GS (Mbp)")+
  xlab("PICs for CT Min")+
  ggtitle("PICs for CT Min")+
  theme(axis.title = element_text(face="bold", size=14),
        axis.text = element_text(size=12),
        title = element_text(size=16, face="bold"))

#Uncorrected for CTmin
linearctmin<-ggplot(dat, aes(y=MbDNA_Female, x=CTMin_Overall))+
  geom_point(size=3)+
  geom_smooth(method=lm, level=0.95)+
  ylab("Female GS (Mbp)")+
  xlab("CT Min")+
  ggtitle("Uncorrected CT Min")+
  theme(axis.title = element_text(face="bold", size=14),
        axis.text = element_text(size=12),
        title = element_text(size=16, face="bold"))
  
#PIC for Dev Time
correcteddevtime<-ggplot(, aes(y=pic.gsdevtime, x=pic.devtime))+
  geom_point(size=3)+
  geom_smooth(method=lm, level=0.95, colour="red")+
  ylab("PICs for Female GS (Mbp)")+
  xlab("PICs for Dev Time")+
  ggtitle("PICs for Dev Time")+
  theme(axis.title = element_text(face="bold", size=14),
        axis.text = element_text(size=12),
        title = element_text(size=16, face="bold"))

#Uncorrected for Dev Time
lineardevtime<-ggplot(dat, aes(y=MbDNA_Female, x=Dev_Time))+
  geom_point(size=3)+
  geom_smooth(method=lm, level=0.95)+
  ylab("Female GS (Mbp)")+
  xlab("Development Time")+
  ggtitle("Uncorrected CT Min")+
  theme(axis.title = element_text(face="bold", size=14),
        axis.text = element_text(size=12),
        title = element_text(size=16, face="bold"))
  
#PIC for Haploid
correctedhaploid<-ggplot(, aes(y=pic.gshaploid, x=pic.haploid))+
  geom_point(size=3)+
  geom_smooth(method=lm, level=0.95, colour="red")+
  ylab("PICs for Female GS (Mbp)")+
  xlab("PICs for Haploid Number")+
  ggtitle("PICs for Dev Time")+
  theme(axis.title = element_text(face="bold", size=14),
        axis.text = element_text(size=12),
        title = element_text(size=16, face="bold"))

#Uncorrected for Haploid
linearhaploid<-ggplot(dat, aes(y=MbDNA_Female, x=Haploid_Number))+
  geom_point(size=3)+
  geom_smooth(method=lm, level=0.95)+
  ylab("Female GS (Mbp)")+
  xlab("Haploid Number")+
  ggtitle("Uncorrected Haploid Number")+
  theme(axis.title = element_text(face="bold", size=14),
        axis.text = element_text(size=12),
        title = element_text(size=16, face="bold"))

#PIC for Precipitation
correctedprecipitation<-ggplot(, aes(y=pic.gsprecip, x=pic.precip))+
  geom_point(size=3)+
  geom_smooth(method=lm, level=0.95, colour="red")+
  ylab("PICs for Female GS (Mbp)")+
  xlab("PICs for Precipitation")+
  ggtitle("PICs for Precipitation")+
  theme(axis.title = element_text(face="bold", size=14),
        axis.text = element_text(size=12),
        title = element_text(size=16, face="bold"))

#Uncorrected for Precipitation
linearprecipitation<-ggplot(dat, aes(y=MbDNA_Female, x=Precip))+
  geom_point(size=3)+
  geom_smooth(method=lm, level=0.95)+
  ylab("Female GS (Mbp)")+
  xlab("Precipitation")+
  ggtitle("Uncorrected Precipitation")+
  theme(axis.title = element_text(face="bold", size=14),
        axis.text = element_text(size=12),
        title = element_text(size=16, face="bold"))

#ggplot
corrected<-ggplot(, aes(y=pic.gs, x=pic.ctmax))+
  geom_point(size=3)+
  geom_smooth(method=lm, level=0.95, colour="red")+
  ylab("PICs for Female GS (Mbp)")+
  xlab("PICs for CT Max")+
  ggtitle("PICs for CT Max")+
  theme(axis.title = element_text(face="bold", size=14),
        axis.text = element_text(size=12),
        title = element_text(size=16, face="bold"))


linear<-ggplot(dat, aes(y=MbDNA_Female, x=CTMax_Overall))+
  geom_point(size=3)+
  geom_smooth(method=lm, level=0.95)+
  ylab("Female GS (Mbp)")+
  xlab("CT Max")+
  ggtitle("Uncorrected CT Max")+
  theme(axis.title = element_text(face="bold", size=14),
        axis.text = element_text(size=12),
        title = element_text(size=16, face="bold"))

correctedmin<-ggplot(, aes(y=pic.gsctmin, x=pic.ctmin))+
  geom_point(size=3)+
  geom_smooth(method=lm, level=0.95, colour="red")+
  ylab("PICs for Female GS (Mbp)")+
  xlab("PICs for CT Min")+
  ggtitle("PICs for CT Min")+
  theme(axis.title = element_text(face="bold", size=14),
        axis.text = element_text(size=12),
        title = element_text(size=16, face="bold"))


linearmin<-ggplot(dat, aes(y=MbDNA_Female, x=CTMin_Overall))+
  geom_point(size=3)+
  geom_smooth(method=lm, level=0.95)+
  ylab("Female Genome Size (Mbp)")+
  xlab("CT Min")+
  ggtitle("Uncorrected CT Min")+
  theme(axis.title = element_text(face="bold", size=14),
        axis.text = element_text(size=12),
        title = element_text(size=16, face="bold"))

#arrange plots

ggarrange(linearctmin, correctedctmin, linearctmax, correctedctmax, lineartmin, correctedtmin, lineartmax, correctedtmax,  nrow=4, ncol=2, labels='AUTO', font.label=list(size=18, face="bold"))

pdf("data/output/figures3.pdf", height=14, width=8.5)
ggarrange(linearctmin, correctedctmin, linearctmax, correctedctmax, lineartmin, correctedtmin, lineartmax, correctedtmax,  nrow=4, ncol=2, labels='AUTO', font.label=list(size=18, face="bold"))
dev.off()
