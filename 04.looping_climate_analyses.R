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
output.table.pic<-matrix(,nrow=19*10000, ncol=10)
colnames(output.table.pic)<-c("Tree", "Data_set", "Variable", "Estimate", "St. Error", "t-value", "p-value", "r-square",
                              "Adj. r-square", "F-stat")

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
}
}

#write the output
write.csv(output.table.pic, "data/output/looped_output/pic.output.everything.csv")
#make a dataframe
#reading in old file
#output.table.pic<-read.csv("data/output/looped_output/pic.output.everything.csv")
output.table.pic.df<-as.data.frame(output.table.pic)
str(output.table.pic.df)
#fix structure of data
for(i in 1:2){
  output.table.pic.df[[i]]<-as.numeric(output.table.pic.df[[i]])
}
for(i in 4:10){
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
ggplot(output.table.pic.df, aes(x=`p.value`, fill=Variable))+
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
#it varies in 0.044 comparisons
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

#varies in 0.06828 comparisons
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
#varies in 0.1365657 comparisons
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
#varies in 0.01393939 comparisons
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
write.csv(signf.table.pic.df, "data/output/looped_output/significance_proportions_pic.csv")

####reading things in for some updated plots
signf.table.pic.df<-read.csv("data/output/looped_output/pic.output.everything.csv")
output.table.pic.df<-signf.table.pic.df
str(signf.table.pic.df)
#which variable has highest proportion of significance, it's now BIO19
signf.table.pic.df$Variable[signf.table.pic.df$`<0.05`==max(signf.table.pic.df$`<0.05`)]

#bio5, bio9, bio10, bio19
#sig.dat<-subset(output.table.pic.df, output.table.pic.df=="Bio5")
sig.dat <- output.table.pic.df[output.table.pic.df$Variable %in% c("Bio5", "Bio9",
                                                                   "Bio10", "Bio19"), ]
sig.dat$Variable<-as.character(sig.dat$Variable)
sig.dat$Variable<-factor(sig.dat$Variable, ordered=TRUE,
                         levels=c("Bio5", "Bio9", "Bio10", "Bio19"))

str(sig.dat)
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
           label=("Median p-value\nBio5 = 0.087\nBio9 = 0.078\nBio10 = 0.065\nBio19 = 0.060"),
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
           label=("Median Adj. R-squared\nBio5 = 0.028\nBio9 = 0.030\nBio10 = 0.034\nBio19 = 0.036"),
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
           label=("Median Estimate\nBio5 = -0.241\nBio9 = -0.059\nBio10 = -0.212\nBio19 = -0.0420"),
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
  

ggarrange(p.val, adj.r, est, common.legend = TRUE, legend="bottom",
          nrow=1, labels="auto",
          hjust=-0.1, align="hv",
          font.label = list(size=20))
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
output.table.lm<-matrix(,nrow=19*1000, ncol=9)
colnames(output.table.lm)<-c("Data_set",
                             "Variable", "Estimate", "St. Error",
                             "t-value", "p-value", "r-square",
                              "Adj. r-square", "F-stat")

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
}
#write output
write.csv(output.table.lm, "data/output/looped_output/lm.output.everything.csv")
#make it a dataframe
output.table.lm.df<-as.data.frame(output.table.lm)
#fix structure issues
for(i in 1:1){
  output.table.lm.df[[i]]<-as.numeric(output.table.lm.df[[i]])
}
for(i in 3:9){
  output.table.lm.df[[i]]<-as.numeric(output.table.lm.df[[i]])
}
output.table.lm.df$Variable<-factor(output.table.lm.df$Variable, ordered=TRUE, 
                                    levels=c("Bio1", "Bio2", "Bio3", "Bio4", "Bio5", "Bio6",
                                             "Bio7", "Bio8", "Bio9", "Bio10", "Bio11",
                                             "Bio12", "Bio13", "Bio14", "Bio15",
                                             "Bio16", "Bio17", "Bio18", "Bio19"))
####Linear model plots####
#facet plots of p-values
ggplot(output.table.lm.df, aes(x=`p-value`, fill=Variable))+
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



