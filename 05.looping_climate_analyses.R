#PICs for bioclimatic variables vs GS in drosophila
#Uses 100 trees randomly selected from bayes distribution
#Loops through each tree 100 
####times with data randomly selected for species with more than one record
####Read in Packages####
library(ape)
library(phytools)
library(ggplot2)
library(viridis)
library(dplyr)
library(ggpubr)
library(coda)

####read in climate data generated from 03.WorldClimdata-dros-code.R####
dros<-read.csv("data/output/drosophila_bioclim.csv")
#remove data which is missing values
dros<-na.omit(dros)
#remove integer column
dros$X<-NULL
#make species a factor
dros$Species<-as.factor(dros$Species)
####read in trees file####
#100 trees sampled from prior
dros.trees<-read.nexus("data/drostrees.nex")
#list of species we have in our dataset
to.keep<-levels(dros$Species)

####Make a matrix to put in outputs####
#rows are for 100 trees with 100 datasets generated from random selection of 
#bioclim variables from possible lat/longs
output.table.pic<-matrix(,nrow=19*10000, ncol=12)
colnames(output.table.pic)<-c("Tree", "Data_set", "Variable", "Estimate", "St. Error", "t-value", "p-value", "r-square",
                              "Adj. r-square", "F-stat", "Corr.P", "FDR")

#loop to do it all
for(i in 1:length(dros.trees)){
#data table for analysis
  print(i)
output<-matrix(,nrow=length(levels(dros$Species)),ncol=(ncol(dros)))
colnames(output)<-colnames(dros)
#loop to pull values randomly when we have more than one
#l<-1
for(l in 1:100){
for(j in 1:length(levels(dros$Species))){
  currentspecies<-levels(dros$Species)[j]
  #make subset for species
  foo<-subset(dros, dros$Species == currentspecies)
  #get one value for species
  if(length(foo$Species)>1){
    # Get a random row index
    random_row_index <- sample(nrow(foo), 1)
    # Subset the dataframe to get the random row
    output[j,] <- as.matrix(foo[random_row_index, ])
  }
  else{
    output[j,]<-as.matrix(foo)
  }
}
#output as data.frame.
output.df<-as.data.frame(output)
#make numeric columns numeric rather than characters
output.df <- output.df %>%
  mutate(across(c(GS, Latitude, Longitude, Bio1, Bio2, Bio3, Bio4, Bio5, Bio6,
                  Bio7, Bio8, Bio9, Bio10, Bio11, Bio12, Bio13, Bio14, Bio15,
                  Bio16, Bio17, Bio18, Bio19), as.numeric))

#keep only relevant species
tree2<-keep.tip(dros.trees[[i]], to.keep)
#fix polytomies if they are there
tree2<-multi2di(tree2)
###Phylogenetic independent contrast portion
#making GS pic
gspic<-pic(log(output.df$GS), phy=tree2)
#k<-5
#we can't do logs of 0 or less, so we need to modify by adding 1+abs(min(output.df$Bio6))
#filling the table with values for PIC analyses
for(k in 5:23){
  output.table.pic[((i-1)*1900+(l-1)*19+(k-4)),1]<-i
  output.table.pic[((i-1)*1900+(l-1)*19+(k-4)),2]<-l
  output.table.pic[((i-1)*1900+(l-1)*19+(k-4)),3]<-names(output.df[k])
  #make pic of variable
  if(min(na.omit(output.df[[k]]))<1){
    var.dat<-output.df[[k]]+(1+abs(min(na.omit(output.df[[k]]))))
  }
  else{
  var.dat<-output.df[[k]]
  }
  varpic<-pic(log(var.dat), phy=tree2)
  picmod<-lm(gspic~varpic+0)
  testing<-summary(picmod)
  #estimate
  output.table.pic[((i-1)*1900+(l-1)*19+(k-4)),4]<-testing$coefficients[1]
  #st error
  output.table.pic[((i-1)*1900+(l-1)*19+(k-4)),5]<-testing$coefficients[2]
  #t-value
  output.table.pic[((i-1)*1900+(l-1)*19+(k-4)),6]<-testing$coefficients[3]
  #p-value
  output.table.pic[((i-1)*1900+(l-1)*19+(k-4)),7]<-testing$coefficients[4]
  #r-square
  output.table.pic[((i-1)*1900+(l-1)*19+(k-4)),8]<-testing$r.squared
  #adj. R square
  output.table.pic[((i-1)*1900+(l-1)*19+(k-4)),9]<-testing$adj.r.squared[1]
  #f stat
  output.table.pic[((i-1)*1900+(l-1)*19+(k-4)),10]<-testing$fstatistic[1]
}
#here is where the correction goes
#p values are column 7 in matrix
#need to find range for nineteen variables
startrow<-((i-1)*1900+(l-1)*19+(1))
endrow<-((i-1)*1900+(l-1)*19+(19))
output.table.pic[startrow:endrow, 11]<-p.adjust(output.table.pic[startrow:endrow,7], method = "bonferroni")
output.table.pic[startrow:endrow, 12]<-p.adjust(output.table.pic[startrow:endrow,7], method = "BH")

}
}

#write the output
write.csv(output.table.pic, "data/output/looped_output/pic.output.everything.csv")
#make a dataframe
#reading in old file
#output.table.pic.df<-read.csv("data/output/looped_output/pic.output.everything.csv")
#if you're reading in the old data, remove the first "X" variable
#output.table.pic.df$X<-NULL
output.table.pic.df<-as.data.frame(output.table.pic)
str(output.table.pic.df)
#fix structure of data
for(i in 1:2){
  output.table.pic.df[[i]]<-as.numeric(output.table.pic.df[[i]])
}
for(i in 4:12){
  output.table.pic.df[[i]]<-as.numeric(output.table.pic.df[[i]])
}
#order variables
output.table.pic.df$Variable<-factor(output.table.pic.df$Variable, ordered=TRUE,
                                     levels=c("Bio1", "Bio2", "Bio3", "Bio4",
                                              "Bio5", "Bio6","Bio7", "Bio8",
                                              "Bio9", "Bio10", "Bio11",
                                              "Bio12", "Bio13", "Bio14",
                                              "Bio15","Bio16", "Bio17",
                                              "Bio18", "Bio19"))
#output.table.pic.df$Variable<-as.factor(output.table.pic.df$Variable)
str(output.table.pic.df)
output.table.pic.df$X<-NULL
####plots of PIC####
##facet plot of p-values
pic.facet<-ggplot(output.table.pic.df, aes(x=`p-value`, fill=Variable))+
  geom_density()+
  facet_wrap(~Variable, scales="free")+
  geom_vline(xintercept = 0.05, color="red", linewidth=1.1, linetype="dotdash")+
  scale_fill_viridis_d()+
  theme_bw()+
  theme(legend.position = "none",
        strip.text = element_text(face="bold", size=14),
        axis.text=element_text(color="black", size=12),
        axis.title.y=element_blank(),
        axis.title.x=element_text(size=14, face="bold", color="black"))+
  xlim(c(-0.2, 1.1))
pic.corr.facet<-ggplot(output.table.pic.df, aes(x=`Corr.P`, fill=Variable))+
  geom_density()+
  facet_wrap(~Variable, scales="free")+
  geom_vline(xintercept = 0.05, color="red", linewidth=1.1, linetype="dotdash")+
  scale_fill_viridis_d()+
  theme_bw()+
  theme(legend.position = "none",
        strip.text = element_text(face="bold", size=14),
        axis.text=element_text(color="black", size=12),
        axis.title.y=element_blank(),
        axis.title.x=element_text(size=14, face="bold", color="black"))+
  xlim(c(-0.2, 1.1))
pic.corrbh.facet<-ggplot(output.table.pic.df, aes(x=`FDR`, fill=Variable))+
  geom_density()+
  facet_wrap(~Variable, scales="free")+
  geom_vline(xintercept = 0.05, color="red", linewidth=1.1, linetype="dotdash")+
  scale_fill_viridis_d()+
  theme_bw()+
  theme(legend.position = "none",
        strip.text = element_text(face="bold", size=14),
        axis.text=element_text(color="black", size=12),
        axis.title.y=element_blank(),
        axis.title.x=element_text(size=14, face="bold", color="black"))+
  xlim(c(-0.2, 1.1))

#to print pdf
pdf("data/output/supp_facet_pic_fig.pdf", width=9.84, height=7.77)
pic.facet
dev.off()
pdf("data/output/supp_facet_pic_corr_fig.pdf", width=9.84, height=7.77)
pic.corr.facet
dev.off()
pdf("data/output/supp_facet_pic_corrbh_fig.pdf", width=9.84, height=7.77)
pic.corrbh.facet
dev.off()

#hpd calculations
#needs to set up a table with 19 bioclimatic variables
#and lower and upper bounds
hpd_dat<-matrix(,nrow=19, ncol=3)
colnames(hpd_dat)<-c("Variable", "Lower", "Upper")

i<-1
for(i in 1:19){
  levels(output.table.pic.df$Variable)
  foo.dat<-subset(output.table.pic.df, Variable==levels(output.table.pic.df$Variable)[i])
  runs_beta <- as.mcmc(foo.dat$Estimate)
  hpd_results <- HPDinterval(runs_beta)
  hpd_dat[i,1]<-levels(output.table.pic.df$Variable)[i]
  hpd_dat[i,2]<-hpd_results[1]
  hpd_dat[i,3]<-hpd_results[2]
}
hpd_dat<-as.data.frame(hpd_dat)
hpd_dat$Lower<-as.numeric(hpd_dat$Lower)
hpd_dat$Upper<-as.numeric(hpd_dat$Upper)
hpd_dat$Variable<-factor(hpd_dat$Variable, ordered=TRUE,
                                     levels=rev(c("Bio1", "Bio2", "Bio3", "Bio4",
                                              "Bio5", "Bio6","Bio7", "Bio8",
                                              "Bio9", "Bio10", "Bio11",
                                              "Bio12", "Bio13", "Bio14",
                                              "Bio15","Bio16", "Bio17",
                                              "Bio18", "Bio19")))

str(hpd_dat)
#making a plot of higher posterior density
hpd.plot<-ggplot(hpd_dat, aes(x=Variable, ymin=Lower, ymax=Upper))+
  geom_linerange(color="slateblue", size=3)+
  geom_point(aes(y=(Lower+Upper)/2), color="maroon", size=3)+
  coord_flip()+
  ylim(-0.5,0.5)+
  ylab("Beta Coefficient")+
  geom_hline(yintercept = 0, color="red", size=1.5, linetype="dashed")+
  theme_minimal()+
  ggtitle("95% HPD of Beta Coefficent for PIC Analyses")+
  theme(axis.text=element_text(size=14, face="bold", color="black"),
        axis.title.y=element_blank(),
        axis.title.x=element_text(size=16, face="bold"),
        title=element_text(size=18, face="bold"))

#save pdf of figure
pdf("data/output/hpd_bioclim_figure.pdf", width=9.2, height=9.5)
hpd.plot
dev.off()

#save table
write.csv(hpd_dat, "data/output/looped_output/hpd_range.csv")

#break it up by bio5
bio.5<-subset(output.table.pic.df, Variable=="Bio5")
#the below plot facet by all 100 trees
# ggplot(bio.5, aes(x=`p-value`, fill=Tree))+
#   geom_density()+
#   facet_wrap(~Tree, scales="free", ncol=3)+
#   geom_vline(xintercept = 0.05, color="red", linewidth=1.1, linetype="dotdash")+
#   theme_bw()+
#   theme(legend.position = "none",
#         strip.text = element_text(face="bold", size=14),
#         axis.text=element_text(color="black", size=12),
#         axis.title.y=element_blank(),
#         axis.title.x=element_text(size=14, face="bold", color="black"))+
#   xlim(c(-0.2, 1.1))

#does hte data vary across trees?
kruskal.test(`p-value` ~ Tree, data=bio.5)
bio5.comp<-pairwise.wilcox.test(bio.5$`p-value`, bio.5$Tree,
                     p.adjust.method = "bonferroni")
bio5.comp.mat<-bio5.comp$p.value
#it varies in 0.07414 comparisons
sum(bio5.comp.mat<0.05, na.rm=TRUE)/sum(bio5.comp.mat>0, na.rm=TRUE)


#pairwise wilcox test for bio10
bio.10<-subset(output.table.pic.df, Variable=="Bio10")
# ggplot(bio.10, aes(x=`p-value`, fill=Tree))+
#   geom_density()+
#   facet_wrap(~Tree, scales="free", ncol=3)+
#   geom_vline(xintercept = 0.05, color="red", linewidth=1.1, linetype="dotdash")+
#   theme_bw()+
#   theme(legend.position = "none",
#         strip.text = element_text(face="bold", size=14),
#         axis.text=element_text(color="black", size=12),
#         axis.title.y=element_blank(),
#         axis.title.x=element_text(size=14, face="bold", color="black"))+
#   xlim(c(-0.2, 1.1))

kruskal.test(`p-value` ~ Tree, data=bio.10)
bio10.comp<-pairwise.wilcox.test(bio.10$`p-value`, bio.10$Tree,
                                p.adjust.method = "bonferroni")
bio10.comp.mat<-bio10.comp$p.value

#varies in 0.1109091 comparisons
sum(bio10.comp.mat<0.05, na.rm=TRUE)/sum(bio10.comp.mat>0, na.rm=TRUE)


#now for variable bio9#now for TRUEvariable bio9
bio.9<-subset(output.table.pic.df, Variable=="Bio9")
# ggplot(bio.9, aes(x=`p-value`, fill=Tree))+
#   geom_density()+
#   facet_wrap(~Tree, scales="free", ncol=3)+
#   geom_vline(xintercept = 0.05, color="red", linewidth=1.1, linetype="dotdash")+
#   theme_bw()+
#   theme(legend.position = "none",
#         strip.text = element_text(face="bold", size=14),
#         axis.text=element_text(color="black", size=12),
#         axis.title.y=element_blank(),
#         axis.title.x=element_text(size=14, face="bold", color="black"))+
#   xlim(c(-0.2, 1.1))

kruskal.test(`p-value` ~ Tree, data=bio.9)
bio9.comp<-pairwise.wilcox.test(bio.9$`p-value`, bio.9$Tree,
                                 p.adjust.method = "bonferroni")
bio9.comp.mat<-bio9.comp$p.value
#varies in 0.1872727 comparisons
sum(bio9.comp.mat<0.05, na.rm=TRUE)/sum(bio9.comp.mat>0, na.rm=TRUE)


#now for bio19
bio.19<-subset(output.table.pic.df, Variable=="Bio19")
# ggplot(bio.19, aes(x=`p-value`, fill=Tree))+
#   geom_density()+
#   facet_wrap(~Tree, scales="free", ncol=3)+
#   geom_vline(xintercept = 0.05, color="red", linewidth=1.1, linetype="dotdash")+
#   theme_bw()+
#   theme(legend.position = "none",
#         strip.text = element_text(face="bold", size=14),
#         axis.text=element_text(color="black", size=12),
#         axis.title.y=element_blank(),
#         axis.title.x=element_text(size=14, face="bold", color="black"))+
#   xlim(c(-0.2, 1.1))

kruskal.test(`p-value` ~ Tree, data=bio.19)
bio19.comp<-pairwise.wilcox.test(bio.19$`p-value`, bio.19$Tree,
                                p.adjust.method = "bonferroni")
bio19.comp.mat<-bio19.comp$p.value
#varies in 0.01454545 comparisons
sum(bio19.comp.mat<0.05, na.rm=TRUE)/sum(bio19.comp.mat>0, na.rm=TRUE)


##Facet plot of coefficents "estimates"
ggplot(output.table.pic.df, aes(x=Estimate, fill=Variable))+
  geom_density()+
  facet_wrap(~Variable)+geom_vline(xintercept=0, colour="red", linewidth=1.1, linetype="dotdash")+
  theme_bw()+
  theme(legend.position = "none",
        strip.text = element_text(face="bold", size=14),
        axis.text=element_text(color="black", size=12),
        axis.title.y=element_blank(),
        axis.title.x=element_text(size=14, face="bold", color="black"))+
  xlim(c((min(output.table.pic.df$Estimate)-0.2), (abs(min((output.table.pic.df$Estimate)-0.2)))))
#facet plots of adj r squared
ggplot(output.table.pic.df, aes(x=`Adj. r-square`, fill=Variable))+
  geom_density()+
  facet_wrap(~Variable, scales="free_y")+
  geom_vline(xintercept=0, colour="red", linewidth=1.1, linetype="dotdash")+
  theme_bw()+
  theme(legend.position = "none",
        strip.text = element_text(face="bold", size=14),
        axis.text=element_text(color="black", size=12),
        axis.title.y=element_blank(),
        axis.title.x=element_text(size=14, face="bold", color="black"))+
  xlim(c(-0.2,1))

####signifiance of PIC models####
signf.table.pic<-matrix(, nrow=19, ncol=16)
colnames(signf.table.pic)<-c("Variable", "<0.05", "<0.01", "<0.001",
                             "minimum.p", "maximum.p", "mean.p", "median.p",
                             "minimim.adj.r2", "maximum.adj.r2", "mean.adj.r2", "median.adj.r2",
                             "minimum.est", "maximum.est", "mean.est", "median.est"
                             )
#fill table with proportions
for(i in 1:length(levels(output.table.pic.df$Variable))){
  #input varialbe
  signf.table.pic[i,1]<-levels(output.table.pic.df$Variable)[i]
  #subset data by variable
  dat<-subset(output.table.pic.df, output.table.pic.df$Variable==levels(output.table.pic.df$Variable)[i])
  #percent significant at 0.05 level
  signf.table.pic[i,2]<-length(dat[[7]][dat[[7]]<0.05])/length(dat[[7]])
  #percent significnat at 0.01 level
  signf.table.pic[i,3]<-length(dat[[7]][dat[[7]]<0.01])/length(dat[[7]])
  #percent signifiant at 0.001 level
  signf.table.pic[i,4]<-length(dat[[7]][dat[[7]]<0.001])/length(dat[[7]])
  #minimum
  signf.table.pic[i,5]<-min(dat[[7]])
  #maximum
  signf.table.pic[i,6]<-max(dat[[7]])
  #mean
  signf.table.pic[i,7]<-mean(dat[[7]])
  #median
  signf.table.pic[i,8]<-median(dat[[7]])
  #min adj R2
  signf.table.pic[i,9]<-min(dat[[9]])
  #max adj r2
  signf.table.pic[i,10]<-max(dat[[9]])
  #mean adj r2
  signf.table.pic[i,11]<-mean(dat[[9]])
  #median adj r2
  signf.table.pic[i,12]<-median(dat[[9]])
  #min est.
  signf.table.pic[i,13]<-min(dat[[4]])
  #max est
  signf.table.pic[i,14]<-max(dat[[4]])
  #mean est
  signf.table.pic[i,15]<-mean(dat[[4]])
  #median est
  signf.table.pic[i,16]<-median(dat[[4]])
}
#make matrix and fix structure
signf.table.pic.df<-as.data.frame(signf.table.pic)
for(i in 2:8){
  signf.table.pic.df[[i]]<-as.numeric(signf.table.pic.df[[i]])
}
signf.table.pic.df$Variable<-as.factor(signf.table.pic.df$Variable)
# # Apply multiple comparison corrections across the 19 bioclimatic variables
# # using median p-value as the per-variable summary statistic
# signf.table.pic.df$p.bonferroni <- p.adjust(signf.table.pic.df$median.p, 
#                                             method = "bonferroni")
# signf.table.pic.df$p.fdr        <- p.adjust(signf.table.pic.df$median.p, 
#                                             method = "BH")
# # Flag which variables survive each correction
# signf.table.pic.df$sig.bonferroni <- signf.table.pic.df$p.bonferroni < 0.05
# signf.table.pic.df$sig.fdr        <- signf.table.pic.df$p.fdr < 0.05


write.csv(signf.table.pic.df, "data/output/looped_output/significance_proportions_pic.csv")

####signifiance of PIC corrected bonferroni models####
signf.table.pic.bon<-matrix(, nrow=19, ncol=16)
colnames(signf.table.pic.bon)<-c("Variable", "<0.05", "<0.01", "<0.001",
                             "minimum.p", "maximum.p", "mean.p", "median.p",
                             "minimim.adj.r2", "maximum.adj.r2", "mean.adj.r2", "median.adj.r2",
                             "minimum.est", "maximum.est", "mean.est", "median.est"
)
#fill table with proportions
for(i in 1:length(levels(output.table.pic.df$Variable))){
  #input varialbe
  signf.table.pic.bon[i,1]<-levels(output.table.pic.df$Variable)[i]
  #subset data by variable
  dat<-subset(output.table.pic.df, output.table.pic.df$Variable==levels(output.table.pic.df$Variable)[i])
  #percent significant at 0.05 level
  signf.table.pic.bon[i,2]<-length(dat[[11]][dat[[11]]<0.05])/length(dat[[11]])
  #percent significnat at 0.01 level
  signf.table.pic.bon[i,3]<-length(dat[[11]][dat[[11]]<0.01])/length(dat[[11]])
  #percent signifiant at 0.001 level
  signf.table.pic.bon[i,4]<-length(dat[[11]][dat[[11]]<0.001])/length(dat[[11]])
  #minimum
  signf.table.pic.bon[i,5]<-min(dat[[11]])
  #maximum
  signf.table.pic.bon[i,6]<-max(dat[[11]])
  #mean
  signf.table.pic.bon[i,7]<-mean(dat[[11]])
  #median
  signf.table.pic.bon[i,8]<-median(dat[[11]])
  #min adj R2
  signf.table.pic.bon[i,9]<-min(dat[[9]])
  #max adj r2
  signf.table.pic.bon[i,10]<-max(dat[[9]])
  #mean adj r2
  signf.table.pic.bon[i,11]<-mean(dat[[9]])
  #median adj r2
  signf.table.pic.bon[i,12]<-median(dat[[9]])
  #min est.
  signf.table.pic.bon[i,13]<-min(dat[[4]])
  #max est
  signf.table.pic.bon[i,14]<-max(dat[[4]])
  #mean est
  signf.table.pic.bon[i,15]<-mean(dat[[4]])
  #median est
  signf.table.pic.bon[i,16]<-median(dat[[4]])
}
#make matrix and fix structure
signf.table.pic.bon.df<-as.data.frame(signf.table.pic.bon)
for(i in 2:8){
  signf.table.pic.bon.df[[i]]<-as.numeric(signf.table.pic.bon.df[[i]])
}
signf.table.pic.bon.df$Variable<-as.factor(signf.table.pic.bon.df$Variable)

write.csv(signf.table.pic.bon.df, "data/output/looped_output/significance_proportions_pic_bonferroni.csv")

####signifiance of PIC corrected bh models####
signf.table.pic.bh<-matrix(, nrow=19, ncol=16)
colnames(signf.table.pic.bh)<-c("Variable", "<0.05", "<0.01", "<0.001",
                                 "minimum.p", "maximum.p", "mean.p", "median.p",
                                 "minimim.adj.r2", "maximum.adj.r2", "mean.adj.r2", "median.adj.r2",
                                 "minimum.est", "maximum.est", "mean.est", "median.est"
)
#fill table with proportions
for(i in 1:length(levels(output.table.pic.df$Variable))){
  #input varialbe
  signf.table.pic.bh[i,1]<-levels(output.table.pic.df$Variable)[i]
  #subset data by variable
  dat<-subset(output.table.pic.df, output.table.pic.df$Variable==levels(output.table.pic.df$Variable)[i])
  #percent significant at 0.05 level
  signf.table.pic.bh[i,2]<-length(dat[[12]][dat[[12]]<0.05])/length(dat[[12]])
  #percent significnat at 0.01 level
  signf.table.pic.bh[i,3]<-length(dat[[12]][dat[[12]]<0.01])/length(dat[[12]])
  #percent signifiant at 0.001 level
  signf.table.pic.bh[i,4]<-length(dat[[12]][dat[[12]]<0.001])/length(dat[[12]])
  #minimum
  signf.table.pic.bh[i,5]<-min(dat[[12]])
  #maximum
  signf.table.pic.bh[i,6]<-max(dat[[12]])
  #mean
  signf.table.pic.bh[i,7]<-mean(dat[[12]])
  #median
  signf.table.pic.bh[i,8]<-median(dat[[12]])
  #min adj R2
  signf.table.pic.bh[i,9]<-min(dat[[9]])
  #max adj r2
  signf.table.pic.bh[i,10]<-max(dat[[9]])
  #mean adj r2
  signf.table.pic.bh[i,11]<-mean(dat[[9]])
  #median adj r2
  signf.table.pic.bh[i,12]<-median(dat[[9]])
  #min est.
  signf.table.pic.bh[i,13]<-min(dat[[4]])
  #max est
  signf.table.pic.bh[i,14]<-max(dat[[4]])
  #mean est
  signf.table.pic.bh[i,15]<-mean(dat[[4]])
  #median est
  signf.table.pic.bh[i,16]<-median(dat[[4]])
}
#make matrix and fix structure
signf.table.pic.bh.df<-as.data.frame(signf.table.pic.bh)
for(i in 2:8){
  signf.table.pic.bh.df[[i]]<-as.numeric(signf.table.pic.bh.df[[i]])
}
signf.table.pic.bh.df$Variable<-as.factor(signf.table.pic.bh.df$Variable)

write.csv(signf.table.pic.bh.df, "data/output/looped_output/significance_proportions_pic_bh.csv")


####reading things in for some updated plots
#signf.table.pic.df<-read.csv("data/output/looped_output/pic.output.everything.csv")
#output.table.pic.df<-signf.table.pic.df
str(signf.table.pic.df)
#which variable has highest proportion of significance, it's now BIO19
signf.table.pic.df$Variable[signf.table.pic.df$`<0.05`==max(signf.table.pic.df$`<0.05`)]
signf.table.pic.bon.df$Variable[signf.table.pic.bon.df$`<0.05`==max(signf.table.pic.bon.df$`<0.05`)]
signf.table.pic.bh.df$Variable[signf.table.pic.bh.df$`<0.05`==max(signf.table.pic.bh.df$`<0.05`)]

#bio5, bio9, bio10, bio19
#sig.dat<-subset(output.table.pic.df, output.table.pic.df=="Bio5")
sig.dat <- output.table.pic.df[output.table.pic.df$Variable %in% c("Bio5", "Bio9",
                                                                   "Bio10", "Bio19"), ]
sig.dat$Variable<-as.character(sig.dat$Variable)
sig.dat$Variable<-factor(sig.dat$Variable, ordered=TRUE,
                         levels=c("Bio5", "Bio9", "Bio10", "Bio19"))

str(sig.dat)
colnames(sig.dat)[7]<-"p.value"
colnames(sig.dat)[9]<-'Adj..r.square'
library(plyr)
#sometimes these have to have modified names in order to match
#depending on if you re-read in data, it might be `p-value` vs. p.value
mu <- ddply(sig.dat, "Variable", summarise, grp.median=median(p.value))
mu.r <- ddply(sig.dat, "Variable", summarise, grp.median=median(Adj..r.square))
mu.e <- ddply(sig.dat, "Variable", summarise, grp.median=median(Estimate))
# me.ev.<-ddply(output.table.pic.df, "Variable", summarise, grp.median=median(`p-value`))

# ggplot(output.table.pic.df, aes(x=`p-value`, fill=Variable))+
#   geom_density(alpha=0.4)+
#   xlim(-0.3,1)+
#   geom_vline(xintercept = 0.05, color="black", linetype=
#                "dashed", linewidth=2)+
#   geom_vline(data=me.ev., aes(xintercept=grp.median, color=Variable), linewidth=1)+
#   theme_bw()#+
#  coord_flip()


str(sig.dat)
str(mu)
mu$grp.median
p.val<-ggplot(sig.dat, aes(x=p.value, fill=Variable))+
  geom_density(alpha=0.4)+
  xlim(-0.3,1)+
  scale_fill_viridis_d(end=0.8)+
  scale_color_viridis_d(end=0.8)+
  geom_vline(xintercept = 0.05, color="red",
             linewidth=2, linetype="dashed")+
  geom_vline(xintercept = 0, color="black",
             linetype="dashed", linewidth=1.5)+
  
  geom_vline(data=mu, aes(xintercept=grp.median, color=Variable),
             linewidth=1)+
  annotate("text", x=1, y=7, 
           label=("Median p-value\nBio5 = 0.074\nBio9 = 0.097\nBio10 = 0.058\nBio19 = 0.061"),
           size=4, hjust=1)+
  theme_bw()+
  theme(
    axis.title.x = element_text(size=16, color="black", face="bold"),
    axis.title.y = element_blank(),
    axis.text=element_text(size=14, color="black"),
    legend.title = element_blank(),
    legend.text = element_text(size=14, color="black"),
    legend.position = "bottom"
  )
mu.r$grp.median
adj.r<-ggplot(sig.dat, aes(x=Adj..r.square, fill=Variable))+
  geom_density(alpha=0.4)+
  xlim(-0.3,0.35)+
  ylim(0,14)+
  scale_fill_viridis_d(end=0.8)+
  scale_color_viridis_d(end=0.8)+
  geom_vline(xintercept = 0, color="black", 
             linetype="dashed", linewidth=1.5)+
  geom_vline(data=mu.r, aes(xintercept=grp.median, color=Variable),
             linewidth=1)+
  annotate("text", x=0.35, y=11.75, 
           label=("Median Adj. R-squared\nBio5 = 0.031\nBio9 = 0.0325\nBio10 = 0.036\nBio19 = 0.035"),
           size=4, hjust=1)+
  theme_bw()+
  theme(
    axis.title.x = element_text(size=16, color="black", face="bold"),
    axis.title.y = element_blank(),
    axis.text=element_text(size=14, color="black"),
    legend.title = element_blank(),
    legend.text = element_text(size=14, color="black"),
    legend.position = "bottom"
  )
mu.e$grp.median 
est<-ggplot(sig.dat, aes(x=`Estimate`, fill=Variable))+
    geom_density(alpha=0.4)+
    xlim(-0.7,0.7)+
  ylim(0,22)+
  scale_fill_viridis_d(end=0.8)+
  scale_color_viridis_d(end=0.8)+
  geom_vline(xintercept = 0, color="black",
             linetype="dashed", linewidth=1.5)+
  geom_vline(data=mu.e, aes(xintercept=grp.median, color=Variable),
             linewidth=1)+
  annotate("text", x=0.7, y=18.5, 
           label=("Median Estimate\nBio5 = -0.245\nBio9 = -0.058\nBio10 = -0.210\nBio19 = -0.041"),
           size=4, hjust=1)+
  theme_bw()+
  theme(
    axis.title.x = element_text(size=16, color="black", face="bold"),
    axis.title.y = element_blank(),
    axis.text=element_text(size=14, color="black"),
    legend.title = element_blank(),
    legend.text = element_text(size=14, color="black"),
    legend.position = "bottom"
  )
  

unique.bioclim<-ggarrange(p.val, adj.r, est, common.legend = TRUE, legend="bottom",
          nrow=1, labels="auto",
          hjust=-0.1, align="hv",
          font.label = list(size=20))
pdf("data/output/figure3.pdf", width=14.21, height=4.25)
unique.bioclim
dev.off()

#doing the same for BH correction
####reading things in for some updated plots
#signf.table.pic.df<-read.csv("data/output/looped_output/pic.output.everything.csv")
#output.table.pic.df<-signf.table.pic.df
str(signf.table.pic.df)
#which variable has highest proportion of significance, it's now BIO19
signf.table.pic.df$Variable[signf.table.pic.df$`<0.05`==max(signf.table.pic.df$`<0.05`)]
signf.table.pic.bon.df$Variable[signf.table.pic.bon.df$`<0.05`==max(signf.table.pic.bon.df$`<0.05`)]
signf.table.pic.bh.df$Variable[signf.table.pic.bh.df$`<0.05`==max(signf.table.pic.bh.df$`<0.05`)]

#bio5, bio9, bio10, bio19
#sig.dat<-subset(output.table.pic.df, output.table.pic.df=="Bio5")
sig.dat.bh <- output.table.pic.df[output.table.pic.df$Variable %in% c("Bio5", "Bio9",
                                                                   "Bio10", "Bio19"), ]
sig.dat.bh$Variable<-as.character(sig.dat.bh$Variable)
sig.dat.bh$Variable<-factor(sig.dat.bh$Variable, ordered=TRUE,
                         levels=c("Bio5", "Bio9", "Bio10", "Bio19"))

str(sig.dat.bh)
colnames(sig.dat.bh)[7]<-"p.value"
colnames(sig.dat.bh)[9]<-'Adj..r.square'
library(plyr)
#sometimes these have to have modified names in order to match
#depending on if you re-read in data, it might be `p-value` vs. p.value
mu <- ddply(sig.dat.bh, "Variable", summarise, grp.median=median(FDR))
mu.r <- ddply(sig.dat.bh, "Variable", summarise, grp.median=median(Adj..r.square))
mu.e <- ddply(sig.dat.bh, "Variable", summarise, grp.median=median(Estimate))
# me.ev.<-ddply(output.table.pic.df, "Variable", summarise, grp.median=median(`p-value`))

# ggplot(output.table.pic.df, aes(x=`p-value`, fill=Variable))+
#   geom_density(alpha=0.4)+
#   xlim(-0.3,1)+
#   geom_vline(xintercept = 0.05, color="black", linetype=
#                "dashed", linewidth=2)+
#   geom_vline(data=me.ev., aes(xintercept=grp.median, color=Variable), linewidth=1)+
#   theme_bw()#+
#  coord_flip()


str(sig.dat.bh)
str(mu)
mu$grp.median
p.val<-ggplot(sig.dat.bh, aes(x=FDR, fill=Variable))+
  geom_density(alpha=0.4)+
  xlim(-0.3,1)+
  scale_fill_viridis_d(end=0.8)+
  scale_color_viridis_d(end=0.8)+
  geom_vline(xintercept = 0.05, color="red",
             linewidth=2, linetype="dashed")+
  geom_vline(xintercept = 0, color="black",
             linetype="dashed", linewidth=1.5)+
  
  geom_vline(data=mu, aes(xintercept=grp.median, color=Variable),
             linewidth=1)+
  annotate("text", x=1, y=7, 
           label=("Median p-value\nBio5 = 0.271\nBio9 = 0.320\nBio10 = 0.246\nBio19 = 0.289"),
           size=4, hjust=1)+
  theme_bw()+
  xlab("Benjamini-Hochberg\nCorrection")+
  theme(
    axis.title.x = element_text(size=16, color="black", face="bold"),
    axis.title.y = element_blank(),
    axis.text=element_text(size=14, color="black"),
    legend.title = element_blank(),
    legend.text = element_text(size=14, color="black"),
    legend.position = "bottom"
  )
mu.r$grp.median
adj.r<-ggplot(sig.dat.bh, aes(x=Adj..r.square, fill=Variable))+
  geom_density(alpha=0.4)+
  xlim(-0.3,0.35)+
  ylim(0,14)+
  scale_fill_viridis_d(end=0.8)+
  scale_color_viridis_d(end=0.8)+
  geom_vline(xintercept = 0, color="black", 
             linetype="dashed", linewidth=1.5)+
  geom_vline(data=mu.r, aes(xintercept=grp.median, color=Variable),
             linewidth=1)+
  annotate("text", x=0.35, y=11.75, 
           label=("Median Adj. R-squared\nBio5 = 0.031\nBio9 = 0.0325\nBio10 = 0.036\nBio19 = 0.035"),
           size=4, hjust=1)+
  theme_bw()+
  xlab("Adj. r-square")+
  theme(
    axis.title.x = element_text(size=16, color="black", face="bold"),
    axis.title.y = element_blank(),
    axis.text=element_text(size=14, color="black"),
    legend.title = element_blank(),
    legend.text = element_text(size=14, color="black"),
    legend.position = "bottom"
  )
mu.e$grp.median 
est<-ggplot(sig.dat.bh, aes(x=`Estimate`, fill=Variable))+
  geom_density(alpha=0.4)+
  xlim(-0.7,0.7)+
  ylim(0,22)+
  scale_fill_viridis_d(end=0.8)+
  scale_color_viridis_d(end=0.8)+
  geom_vline(xintercept = 0, color="black",
             linetype="dashed", linewidth=1.5)+
  geom_vline(data=mu.e, aes(xintercept=grp.median, color=Variable),
             linewidth=1)+
  annotate("text", x=0.7, y=18.5, 
           label=("Median Estimate\nBio5 = -0.245\nBio9 = -0.058\nBio10 = -0.210\nBio19 = -0.041"),
           size=4, hjust=1)+
  theme_bw()+
  theme(
    axis.title.x = element_text(size=16, color="black", face="bold"),
    axis.title.y = element_blank(),
    axis.text=element_text(size=14, color="black"),
    legend.title = element_blank(),
    legend.text = element_text(size=14, color="black"),
    legend.position = "bottom"
  )


unique.bioclim<-ggarrange(p.val, adj.r, est, common.legend = TRUE, legend="bottom",
                          nrow=1, labels="auto",
                          hjust=-0.1, align="hv",
                          font.label = list(size=20))
pdf("data/output/figure3.bh.pdf", width=14.21, height=4.25)
unique.bioclim
dev.off()

# ?ggarrange
# ggplot(sig.dat, aes(x=`p-value`, y=Variable,fill=Variable))+
#   geom_violin(alpha=0.6, draw_quantiles = TRUE)+
#   geom_vline(xintercept = 0.05, color="red", linetype=
#                "dashed", linewidth=2)+
#   theme_bw()
# 
# ggplot(sig.dat, aes(x=`p-value`, y=Variable,fill=Variable))+
#   geom_boxplot(alpha=0.6)+
#   geom_vline(xintercept = 0.05, color="red", linetype=
#                "dashed", linewidth=2)



#####linear Models####
output.table.lm<-matrix(,nrow=19*1000, ncol=11)
colnames(output.table.lm)<-c("Data_set",
                             "Variable", "Estimate", "St. Error",
                             "t-value", "p-value", "r-square",
                              "Adj. r-square", "F-stat", "Corr.P", "FDR")

#loop to do 1000 sets of data randomly selecting geographical records for species
#with more than one record
i<-1
for(i in 1:1000){
  output<-matrix(,nrow=length(levels(dros$Species)),ncol=(ncol(dros)))
  colnames(output)<-colnames(dros)
  
  for(j in 1:length(levels(dros$Species))){
    currentspecies<-levels(dros$Species)[j]
    #make subset for species
    foo<-subset(dros, dros$Species == currentspecies)
    #get one value for species
    if(length(foo$Species)>1){
      # Get a random row index
      random_row_index <- sample(nrow(foo), 1)
      # Subset the dataframe to get the random row
      output[j,] <- as.matrix(foo[random_row_index, ])
    }
    else{
      output[j,]<-as.matrix(foo)
    }
  }
  #output as data.frame.
  output.df<-as.data.frame(output)
  #make numeric columns numeric rather than characters
  output.df <- output.df %>%
    mutate(across(c(GS, Latitude, Longitude, Bio1, Bio2, Bio3, Bio4, Bio5, Bio6,
                    Bio7, Bio8, Bio9, Bio10, Bio11, Bio12, Bio13, Bio14, Bio15,
                    Bio16, Bio17, Bio18, Bio19), as.numeric))

  for(k in 5:23){
    output.table.lm[((i-1)*19+(k-4)),1]<-i
    output.table.lm[((i-1)*19+(k-4)),2]<-names(output.df[k])
    
    if(min(as.numeric(na.omit(output.df[[k]])))<1){
       var.dat<-as.numeric(output.df[[k]])+(1+abs(min(as.numeric(na.omit(output.df[[k]])))))
     }
     else{
      var.dat<-as.numeric(output.df[[k]])
     }
    lmmod<-lm(log(as.numeric(output.df$GS))~log(var.dat))
    testing<-summary(lmmod)
    #estimate
    output.table.lm[((i-1)*19+(k-4)),3]<-testing$coefficients[2]
    #st error
    output.table.lm[((i-1)*19+(k-4)),4]<-testing$coefficients[4]
    #t-value
    output.table.lm[((i-1)*19+(k-4)),5]<-testing$coefficients[6]
    #p-value
    output.table.lm[((i-1)*19+(k-4)),6]<-testing$coefficients[8]
    #r-square
    output.table.lm[((i-1)*19+(k-4)),7]<-testing$r.squared
    #adj. R square
    output.table.lm[((i-1)*19+(k-4)),8]<-testing$adj.r.squared[1]
    #f stat
    output.table.lm[((i-1)*19+(k-4)),9]<-testing$fstatistic[1]
  }
  #here is where the correction goes
  #p values are column 6 in matrix
  #need to find range for nineteen variables
  startrow<-((i-1)*19+(1))
  endrow<-((i-1)*19+(19))
  output.table.lm[startrow:endrow, 10]<-p.adjust(output.table.lm[startrow:endrow,6], method = "bonferroni")
  output.table.lm[startrow:endrow, 11]<-p.adjust(output.table.lm[startrow:endrow,6], method = "BH")
}
#write output
write.csv(output.table.lm, "data/output/looped_output/lm.output.everything.csv")
#output.table.lm.df<-read.csv("data/output/looped_output/lm.output.everything.csv")
#output.table.lm.df$X<-NULL
#make it a dataframe
output.table.lm.df<-as.data.frame(output.table.lm)
#fix structure issues
for(i in 1:1){
  output.table.lm.df[[i]]<-as.numeric(output.table.lm.df[[i]])
}
for(i in 3:11){
  output.table.lm.df[[i]]<-as.numeric(output.table.lm.df[[i]])
}
output.table.lm.df$Variable<-factor(output.table.lm.df$Variable, ordered=TRUE, 
                                    levels=c("Bio1", "Bio2", "Bio3", "Bio4", "Bio5", "Bio6",
                                             "Bio7", "Bio8", "Bio9", "Bio10", "Bio11",
                                             "Bio12", "Bio13", "Bio14", "Bio15",
                                             "Bio16", "Bio17", "Bio18", "Bio19"))
####Linear model plots####
str(output.table.lm.df)
#facet plots of p-values
lm.facet<-ggplot(output.table.lm.df, aes(x=`p-value`, fill=Variable))+
  geom_density()+
  facet_wrap(~Variable, scales="free")+
  geom_vline(xintercept = 0.05, color="red", linewidth=1.1, linetype="dotdash")+
  theme_bw()+
  theme(legend.position = "none",
        strip.text = element_text(face="bold", size=14),
        axis.text=element_text(color="black", size=8),
        axis.title.y=element_blank(),
        axis.title.x=element_text(size=14, face="bold", color="black"))+
  xlim(c(-0.2, 1))

lm.facet.bon<-ggplot(output.table.lm.df, aes(x=`Corr.P`, fill=Variable))+
  geom_density()+
  facet_wrap(~Variable, scales="free")+
  geom_vline(xintercept = 0.05, color="red", linewidth=1.1, linetype="dotdash")+
  theme_bw()+
  theme(legend.position = "none",
        strip.text = element_text(face="bold", size=14),
        axis.text=element_text(color="black", size=8),
        axis.title.y=element_blank(),
        axis.title.x=element_text(size=14, face="bold", color="black"))+
  xlim(c(-0.2, 1))


lm.facet.bh<-ggplot(output.table.lm.df, aes(x=`FDR`, fill=Variable))+
  geom_density()+
  facet_wrap(~Variable, scales="free")+
  geom_vline(xintercept = 0.05, color="red", linewidth=1.1, linetype="dotdash")+
  theme_bw()+
  theme(legend.position = "none",
        strip.text = element_text(face="bold", size=14),
        axis.text=element_text(color="black", size=8),
        axis.title.y=element_blank(),
        axis.title.x=element_text(size=14, face="bold", color="black"))+
  xlim(c(-0.2, 1))

pdf("data/output/looped_output/lm.model.facet.pvalue.pdf", width=9.83, height=7.76)
lm.facet
dev.off()
pdf("data/output/looped_output/lm.model.facet.pvalue.bon.pdf", width=9.83, height=7.76)
lm.facet.bon
dev.off()
pdf("data/output/looped_output/lm.model.facet.pvalue.bh.pdf", width=9.83, height=7.76)
lm.facet.bh
dev.off()

##Facet plot of "estimates" coefficients 
ggplot(output.table.lm.df, aes(x=Estimate, fill=Variable))+
  geom_density()+
  facet_wrap(~Variable)+
  theme_bw()+
  geom_vline(xintercept=0, colour="red", linewidth=1.1, linetype="dotdash")+
  theme(legend.position = "none",
        strip.text = element_text(face="bold", size=14),
        axis.text=element_text(color="black", size=12),
        axis.title.y=element_blank(),
        axis.title.x=element_text(size=14, face="bold", color="black"))
  
##Facet plot of adj. r-square
ggplot(output.table.lm.df, aes(x=`Adj. r-square`, fill=Variable))+
  geom_density()+
  facet_wrap(~Variable)+
  geom_vline(xintercept=0, colour="red", linewidth=1.1, linetype="dotdash")+
  theme_bw()+
  theme(legend.position = "none",
        strip.text = element_text(face="bold", size=14),
        axis.text=element_text(color="black", size=12),
        axis.title.y=element_blank(),
        axis.title.x=element_text(size=14, face="bold", color="black"))

####signifiance of models####
#emptpy matrix of values
signf.table.lm<-matrix(, nrow=19, ncol=4)
colnames(signf.table.lm)<-c("Variable", "%p<0.05", "%p<0.01", "%p<0.001")

#loop through to get proportions of significant values
for(i in 1:length(levels(output.table.lm.df$Variable))){
  #input varialbe
  signf.table.lm[i,1]<-levels(output.table.lm.df$Variable)[i]
  #subset data by variable
  dat<-subset(output.table.lm.df, output.table.lm.df$Variable==levels(output.table.lm.df$Variable)[i])
  #percent significant at 0.05 level
  signf.table.lm[i,2]<-length(dat[[6]][dat[[6]]<0.05])/length(dat[[6]])
  #percent significnat at 0.01 level
  signf.table.lm[i,3]<-length(dat[[6]][dat[[6]]<0.01])/length(dat[[6]])
  #percent signifiant at 0.001 level
  signf.table.lm[i,4]<-length(dat[[6]][dat[[6]]<0.001])/length(dat[[6]])
}
#make matrix a dataframe
signf.table.lm.df<-as.data.frame(signf.table.lm)
for(i in 2:4){
  signf.table.lm.df[[i]]<-as.numeric(signf.table.lm.df[[i]])
}
signf.table.lm.df$Variable<-as.factor(signf.table.lm.df$Variable)
str(signf.table.lm.df)
write.csv(signf.table.lm.df, "data/output/looped_output/significance_proportions_lm.csv")

#emptpy matrix of values
signf.table.lm.bon<-matrix(, nrow=19, ncol=4)
colnames(signf.table.lm.bon)<-c("Variable", "%p<0.05", "%p<0.01", "%p<0.001")
colnames(output.table.lm.df)
#loop through to get proportions of significant values
for(i in 1:length(levels(output.table.lm.df$Variable))){
  #input varialbe
  signf.table.lm.bon[i,1]<-levels(output.table.lm.df$Variable)[i]
  #subset data by variable
  dat<-subset(output.table.lm.df, output.table.lm.df$Variable==levels(output.table.lm.df$Variable)[i])
  #percent significant at 0.05 level
  signf.table.lm.bon[i,2]<-length(dat[[10]][dat[[10]]<0.05])/length(dat[[10]])
  #percent significnat at 0.01 level
  signf.table.lm.bon[i,3]<-length(dat[[10]][dat[[10]]<0.01])/length(dat[[10]])
  #percent signifiant at 0.001 level
  signf.table.lm.bon[i,4]<-length(dat[[10]][dat[[10]]<0.001])/length(dat[[10]])
}
#make matrix a dataframe
signf.table.lm.bon.df<-as.data.frame(signf.table.lm.bon)
for(i in 2:4){
  signf.table.lm.bon.df[[i]]<-as.numeric(signf.table.lm.bon.df[[i]])
}
signf.table.lm.bon.df$Variable<-as.factor(signf.table.lm.bon.df$Variable)
str(signf.table.lm.bon.df)
write.csv(signf.table.lm.bon.df, "data/output/looped_output/significance_proportions_lm_bon.csv")

#emptpy matrix of values
signf.table.lm.bh<-matrix(, nrow=19, ncol=4)
colnames(signf.table.lm.bh)<-c("Variable", "%p<0.05", "%p<0.01", "%p<0.001")
colnames(output.table.lm.df)
#loop through to get proportions of significant values
for(i in 1:length(levels(output.table.lm.df$Variable))){
  #input varialbe
  signf.table.lm.bh[i,1]<-levels(output.table.lm.df$Variable)[i]
  #subset data by variable
  dat<-subset(output.table.lm.df, output.table.lm.df$Variable==levels(output.table.lm.df$Variable)[i])
  #percent significant at 0.05 level
  signf.table.lm.bh[i,2]<-length(dat[[11]][dat[[11]]<0.05])/length(dat[[11]])
  #percent significnat at 0.01 level
  signf.table.lm.bh[i,3]<-length(dat[[11]][dat[[11]]<0.01])/length(dat[[11]])
  #percent signifiant at 0.001 level
  signf.table.lm.bh[i,4]<-length(dat[[11]][dat[[11]]<0.001])/length(dat[[11]])
}
#make matrix a dataframe
signf.table.lm.bh.df<-as.data.frame(signf.table.lm.bh)
for(i in 2:4){
  signf.table.lm.bh.df[[i]]<-as.numeric(signf.table.lm.bh.df[[i]])
}
signf.table.lm.bh.df$Variable<-as.factor(signf.table.lm.bh.df$Variable)
str(signf.table.lm.bh.df)
write.csv(signf.table.lm.bh.df, "data/output/looped_output/significance_proportions_lm_bh.csv")

