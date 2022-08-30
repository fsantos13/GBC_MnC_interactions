#Fernanda Santos and Elizabeth Herndon, 2022
#Special thanks to Dr. Steven Hall in the Department of Ecology, Evolution, and Organismal Biology at Iowa State University, for sharing this code
# Original code from Dr. Hall can be found in EDI (https://portal.edirepository.org/nis/mapbrowse?packageid=edi.575.1), NEON_NMR_heatmap_fig3.R 

#Figure 3
#remove old data
rm(list=ls())

setwd("C:/Users/nanda/Documents/Postdoc_2020/ORNL_Mnprojects/NEONdata/Manuscript_NEON/Code_GlobalBiogeochemicalCycles")

library(viridis)
library(Hmisc)
library(GGally)
library(reshape)
library(scales)

# load NEON data
data1<-read.csv("neon_correlationmatrix.csv",header=TRUE) #predictors which are most NEON variables
data2<-read.csv("neon_correlationmatrix_stocks.csv",header=TRUE) #response variables which are the C and N stocks

#select NEON variables that make sense to add to the correlation matrix
data1_sel <-data1[,c("temp",
                      "precip",
                      "mnMjelm",
                      "sulf",
                      "alMjelm",
                      "caMjelm",
                      "feMjelm",
                      "kMjelm",
                      "mgMjelm",
                      "naMjelm",
                      "pMjelm",
                      "siMjelm",
                      "srMjelm",
                      "tiMjelm",
                      "zrMjelm",
                      "phCacl2",
                      "phH2o",
                      "ec12pre",
                      "mnKcl",
                      "caNh4d",
                      "kNh4d",
                      "mgNh4d",
                      "cecdNh4",
                      "baseSumCecd10",
                      "bsesatCecd10"
                      )]

#Original variable names can be ugly sometimes so you want a clear but short name for the correlation matrix
names(data1)<-c("Temp",
                      "Precip",                 
                      "Mn",
                      "S",
                      "Al",
                      "Ca",
                      "Fe",
                      "K",
                      "Mg",
                      "Na",
                      "P",
                      "Si",
                      "Sr",
                      "Ti",
                      "Zr",
                      "pH (Exchang)",
                      "pH (Water)",
                      "EC",
                      "Mn (Exchang)",
                      "Ca (Exchang)",
                      "K (Exchang)",
                      "Mg (Exchang)",
                      "CEC",
                      "Extractable bases",
                      "Base saturation")

names(data2)<-c("N stock",
                      "C stock"             
                      )

#create a correlation vector
correlation<-as.matrix(cor(data2,data1_sel,use= "pairwise.complete.obs", method="spearman"))
correlation_mat <- melt(correlation)
names(correlation_mat) [3] <- "r"

#make matrix of significance of correlations, using rcorr function in Hmisc
#extract the P-value matrix
Pmat<-as.data.frame(rcorr(as.matrix(data2),as.matrix(data1_sel), type="spearman")$P)
str(Pmat)
#get the same format as for the matrices above
Pmat_red<-Pmat[1:2,3:length(Pmat)]
str(Pmat_red)
#replace non sig values with NA; threshold is p = 0.1/(6*20)
Pmat_red[Pmat_red>0.000833 & Pmat_red<1]<-NA
Pmat_red<-Pmat_red*60
#use Bonferroni correction; replace sig values with 1 or 2 or 3
Pmat_red[Pmat_red< 0.05 & Pmat_red> 0.01]<-4
Pmat_red[Pmat_red< 0.01 & Pmat_red> 0.001]<-3
Pmat_red[Pmat_red< 0.001 & Pmat_red> 0.0001]<-2
Pmat_red[Pmat_red< 0.0001 ]<-1

#now replace with stars
Pmat_red[Pmat_red==1]<-"****"
Pmat_red[Pmat_red==2]<-"***"
Pmat_red[Pmat_red==3]<-"**"
Pmat_red[Pmat_red==4]<-"*"

#melt the matrix to get the same format as above
Pmat_red_melt<-melt(as.matrix(Pmat_red))

#paste the correlation vector to the p value data frame
#for color labels
Pmat_red_melt$Corr<-correlation_mat$r
Pmat_red_melt$Corr[Pmat_red_melt$Corr<0]<- -1
Pmat_red_melt$Corr[Pmat_red_melt$Corr>0]<- 1
Pmat_red_melt$Corr[Pmat_red_melt$Corr==1]<-"black"
Pmat_red_melt$Corr[Pmat_red_melt$Corr==-1]<-"white"

#plot the corr matrix 
#quartz(,10,3)
ggplot(data=correlation_mat)+
  geom_tile(aes(x=X2,y=X1,fill=r),color = "white", lwd = 0.5, linetype = 1)+
  scale_fill_viridis(name=(expression(rho)), option="viridis",limits=c(-0.95,0.95))+
  scale_y_discrete("Carbon and Nitrogen")+
  scale_x_discrete("", labels=c("Temp",
                      "Precip",                 
                      "Mn",
                      "S",
                      "Al",
                      "Ca",
                      "Fe",
                      "K",
                      "Mg",
                      "Na",
                      "P",
                      "Si",
                      "Sr",
                      "Ti",
                      "Zr",
                      "pH (Exchang)",
                      "pH (Water)",
                      "EC",
                      "Mn (Exchang)",
                      "Ca (Exchang)",
                      "K (Exchang)",
                      "Mg (Exchang)",
                      "CEC",
                      "Extractable bases",
                      "Base saturation"))+
  theme(axis.text.x=element_text(angle = +45, vjust=1, hjust = 1,size=12, colour="black", face="bold"),axis.title.x=element_blank(),
        axis.text.y=element_text(size=12, colour="black", face="bold"),
        axis.title.y=element_blank(), axis.ticks.x = element_line (colour="black"))+ coord_fixed()+
  geom_text(data=Pmat_red_melt,
            mapping=aes(x=X2,y=X1,label=value), 
            color=Pmat_red_melt$Corr,
            size=7)



