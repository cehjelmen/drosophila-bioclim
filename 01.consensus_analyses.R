####Script to Run Analyses on Consensus tree
####Load in packages####
library(ggplot2)
library(phytools)
library(dplyr)
library(viridis)
####Read in data####
#pulling in data
dat<-read.csv("data/climate_data_drosophila_Sept17.csv")
#remove data we don't need here
dat$MbDNA_Male<-NULL
length(names(dat))
dat[14:17]<-NULL
dat[8:10]<-NULL
dat$Subgenus<-NULL
#read in tree file
#consensus tree
con.tree<-read.tree("data/mytree2.tre")

to.compare<-colnames(dat)
####Running Linear Analyses
#we want to run the analysis for every tree and every variable
#variables we want to look at
to.compare<-colnames(dat)
#make an output table
#change number of rows to correspond to value
output.table.lm<-matrix(,nrow=(length(to.compare)-2), ncol=9)
#name columns
colnames(output.table.lm)<-c("Variable", 
                              "Estimate", "StError", "tvalue", 
                              "pvalue", "rsquare",
                              "Adjrsquare", "Fstat", "Corr.P")

#loop through linear model and make table
for(i in 3:length(to.compare)){
  #subset data for variable
  #spec.data<-dat$Species[is.na(dat[i])==FALSE]
  foo<-subset(dat, subset=(is.na(dat[i])==FALSE))
  #str(foo)
  #make linear model
  
  linearmod<-lm(log(foo$MbDNA_Female)~foo[[i]])
  testing<-summary(linearmod)

  #info we want
  #variable
  output.table.lm[(i-2),1]<-names(dat[i])
  #estimate
  output.table.lm[(i-2),2]<-testing$coefficients[2]
  #st error
  output.table.lm[(i-2),3]<-testing$coefficients[4]
  #t-value
  output.table.lm[(i-2),4]<-testing$coefficients[6]
  #p-value
  output.table.lm[(i-2),5]<-testing$coefficients[8]
  #r-square
  output.table.lm[(i-2),6]<-testing$r.squared
  #adj. R square
  output.table.lm[(i-2),7]<-testing$adj.r.squared[1]
  #f stat
  output.table.lm[(i-2),8]<-testing$fstatistic[1]
}

lm.out<-as.data.frame(output.table.lm)
lm.out <- lm.out %>%
  mutate(across(c(Estimate,StError,tvalue,pvalue,rsquare, Adjrsquare, Fstat), as.numeric))
str(lm.out)
lm.out$Corr.P<-p.adjust(lm.out$pvalue, method="bonferroni")

#writes output as csv file to your working directory
write.csv(lm.out, "data/output/consensus_analyses/lm_out_data.csv")

#make an output table
#change number of rows to correspond to value
output.table.pic.con<-matrix(,nrow=(length(to.compare)-2), ncol=9)
#name columns
colnames(output.table.pic.con)<-c("Variable", 
                             "Estimate", "StError", "tvalue", 
                             "pvalue", "rsquare",
                             "Adjrsquare", "Fstat", "Corr.P")


for(i in 3:length(to.compare)){
  #get data for species names with the variable of interest
  spec.data<-dat$Species[is.na(dat[i])==FALSE]
  foo<-subset(dat, subset=(is.na(dat[i])==FALSE))
    #keep only relevant species
    tree2<-keep.tip(con.tree, spec.data)
    #fix polytomies if they are there
    tree2<-multi2di(tree2)
    #get gs PIC
    gspic<-pic(log(foo$MbDNA_Female), phy=tree2)
    if(min(na.omit(foo[[i]]))<1){
      var.dat<-foo[[i]]+(1+abs(min(na.omit(foo[[i]]))))
    }
    else{
      var.dat<-foo[[i]]
    }
    #get log pic value for variable
    varpic<-pic(log(var.dat), phy=tree2)
    picmod<-lm(gspic~varpic+0)
    testing<-summary(picmod)
  #info we want
  #variable
  output.table.pic.con[(i-2),1]<-names(dat[i])
  #estimate
  output.table.pic.con[(i-2),2]<-testing$coefficients[1]
  #st error
  output.table.pic.con[(i-2),3]<-testing$coefficients[2]
  #t-value
  output.table.pic.con[(i-2),4]<-testing$coefficients[3]
  #p-value
  output.table.pic.con[(i-2),5]<-testing$coefficients[4]
  #r-square
  output.table.pic.con[(i-2),6]<-testing$r.squared
  #adj. R square
  output.table.pic.con[(i-2),7]<-testing$adj.r.squared[1]
  #f stat
  output.table.pic.con[(i-2),8]<-testing$fstatistic[1]
}

pic.con.out<-as.data.frame(output.table.pic.con)
pic.con.out <- pic.con.out %>%
  mutate(across(c(Estimate,StError,tvalue,pvalue,rsquare, Adjrsquare, Fstat), as.numeric))
str(pic.con.out)
pic.con.out$Corr.P<-p.adjust(pic.con.out$pvalue, method="bonferroni")

#writes output as csv file to your working directory
write.csv(pic.con.out, "data/output/consensus_analyses/pic_con_out_data.csv")
