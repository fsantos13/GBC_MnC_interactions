#Fernanda Santos and Elizabeth Herndon, 2022

setwd("C:/Users/nanda/Documents/Postdoc_2020/ORNL_Mnprojects/NEONdata/Manuscript_NEON/Code_GlobalBiogeochemicalCycles")
#Remove the # if you need to install the packages below
#install.packages("ggplot2")
#install.packages("easyGgplot2")
#install.packages("ggthemes")
#install.packages("scales")
#install.packages("Rmisc")
#install.packages("RColorBrewer")
#install.packages("gridExtra")
#install.packages("readr")
#install.packages("stats")
#install.packages("minpack.lm")
#install.packages("ggpubr")
#install.packages("ppcor") 
#install.packages("ggsci") 

#open libraries
library(ggplot2)
library(easyGgplot2)
library(ggthemes)
library(scales)
library(Rmisc)
library(RColorBrewer)
library(gridExtra)
library(readr)
library(stats)
library(minpack.lm)
library(ggsci) # uses jco palette

#FIGURE 2
corr = read.csv("allsites_ind.csv", header=TRUE)


model1 <- nls(cstock ~ SSasymp (mnMjelm, yf, y0, log_alpha), data = corr)
model2 <- nls(nstock ~ SSasymp (mnMjelm, yf, y0, log_alpha), data = corr)

#Add predicted data to the dataset 

predlma = predict(model1, interval="confidence")
predlmb = predict(model2, interval="confidence")

#Figure 2a
plot1_hor <- ggplot(corr, aes(x= mnMjelm, y= cstock)) + ggtitle ("(a)") + geom_point(data=corr, aes (fill=horizonName),size=5, pch=21) + labs(y=(expression(paste("C stock (kg m"^"-2",")"))), x=(expression(paste("Manganese concentration (mg Mn kg"^" -1",")")))) + theme_bw() + xlim(0,2500) + ylim(-0.5,20) + theme(legend.text=element_text(size=14), axis.title.x=element_text(colour="black", face="bold", size=14), axis.title.y = element_text(colour="black", face="bold", size=14), axis.text = element_text(colour="black", size=14)) + geom_line(aes(y=predlma), size=1, linetype = "dashed", color = "black")  + scale_fill_brewer (palette = "Greys", name= "Soil horizons", breaks=c("Oi","Oe","Oa"), labels=c("Oi","Oe","Oa"))

#Figure 2c
plot2_hor <- ggplot(corr, aes(x= mnMjelm, y= nstock)) + ggtitle ("(c)") + geom_point(data=corr, aes (fill=horizonName),size=5, pch=21) + labs(y=(expression(paste("N stock (kg m"^"-2",")"))), x=(expression(paste("Manganese concentration (mg Mn kg"^" -1",")")))) + theme_bw() + xlim(0,2500) + ylim(-0.1,1.5) + theme(legend.text=element_text(size=14), axis.title.x=element_text(colour="black", face="bold", size=14), axis.title.y = element_text(colour="black", face="bold", size=14), axis.text = element_text(colour="black", size=14)) + geom_line(aes(y=predlmb), size=1, linetype = "dashed", color = "black")  + scale_fill_brewer (palette = "Blues", name= "Soil horizons", breaks=c("Oi","Oe","Oa"), labels=c("Oi","Oe","Oa"))


#Figure 2b
corr3 <- ggplot(corr, aes(y=oneovercstock, x=mnMjelm)) + geom_point(size=2, shape=19, color ="black", fill ="black") + labs(y=(expression(paste("1/C stock"))), x=(expression(paste("Manganese concentration (mg Mn kg"^"  -1",")")))) + theme_bw() + ylim(0,3) + xlim(0,2500) + theme(legend.text=element_text(size=14), axis.title.x=element_text(colour="black", face="bold", size=16), axis.title.y = element_text(colour="black", face="bold", size=16), axis.text = element_text(colour="black", size=16))
#Fit regression line
require(stats)
reg2 <- lm (oneovercstock ~ mnMjelm, data = corr)
reg2
coeff=coefficients(reg2)
rsquared=summary(reg2)$r.squared
# Equation of the line : 
#Change the second number after coeff to change the number 
#of decimals of the slope and intercept values
eq = paste0("(b)  y = ", round(coeff[2],2), "*x + ", round(coeff[1],2), ",", " ","R squared = ", round(rsquared,3)) 
# Plot and fit model to data and add line type, color and size
corr4 <- corr3 + stat_smooth(method="lm", se=FALSE, color="black", linetype="dashed", size=1.5) + geom_point(aes(fill=horizonName), pch=21, size= 5) + scale_fill_brewer (palette = "Greys", name= "Soil horizons", breaks=c("Oi","Oe","Oa"), labels=c("Oi","Oe","Oa")) + ggtitle(eq) + stat_cor(p.accuracy = 0.001, label.x = 3)

#Figure 2d
corrn <- ggplot(corr, aes(y=oneovernstock, x=mnMjelm)) + geom_point(size=2, shape=19, color ="black", fill ="black") + labs(y=(expression(paste("1/N stock"))), x=(expression(paste("Manganese concentration (mg Mn kg"^"  -1",")")))) + theme_bw() + ylim(0,90) + xlim(0,2500) + theme(legend.text=element_text(size=14), axis.title.x=element_text(colour="black", face="bold", size=16), axis.title.y = element_text(colour="black", face="bold", size=16), axis.text = element_text(colour="black", size=16))
corrn
#Fit regression line
#First variable after lm is your y variable and the second variable is your x variable
require(stats)
reg2a <- lm (oneovernstock ~ mnMjelm, data = corr)
reg2a
coeff=coefficients(reg2a)
rsquared=summary(reg2a)$r.squared
# Equation of the line : 
#Change the second number after coeff to change the number 
#of decimals of the slope and intercept values
eq = paste0("(d)  y = ", round(coeff[2],2), "*x + ", round(coeff[1],2), ",", " ","R squared = ", round(rsquared,3)) 
# Plot and fit model to data and add line type, color and size
corrN <- corrn + stat_smooth(method="lm", se=FALSE, color="black", linetype="dashed", size=1.5) + geom_point(aes(fill=horizonName), pch=21, size= 5) + scale_fill_brewer (palette = "Blues", name= "Soil horizons", breaks=c("Oi","Oe","Oa"), labels=c("Oi","Oe","Oa")) + ggtitle(eq) + stat_cor(p.accuracy = 0.001, label.x = 3)


#Figure 2
#now they do not have a common.legend
ggarrange(plot1_hor, plot2_hor, corr4, corrN, nrow=2, ncol=2)








