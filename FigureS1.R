#Fernanda Santos and Elizabeth Herndon, 2022


#Figure S1

# Run Code_GBC2022_Figure1 first


#Set the colors so that there are 16 different colors instead of 8
mycolors = c(colorRampPalette(brewer.pal(name="Dark2", n = 8))(16))

plot1_sites <- ggplot(cor, aes(x= Mn, y= C_stock)) + ggtitle("(a)") + geom_point(data=cor, aes (color=site), size=4, shape=19) + labs(y=(expression(paste("C stock (kg m"^"-2",")"))), x=(expression(paste("Manganese concentration (mg Mn kg"^" -1",")")))) + theme_bw() + xlim(0,2500) + ylim(0,20) + theme(legend.text=element_text(size=14), axis.title.x=element_text(colour="black", face="bold", size=16), axis.title.y = element_text(colour="black", face="bold", size=14), axis.text = element_text(colour="black", size=14)) + geom_line(aes(y=predlmC), size=1, linetype = "dashed", color = "black")   + scale_color_manual(name= "NEON sites", values = mycolors)


plot2_sites <- ggplot(cor, aes(x= Mn, y= N_stock)) + ggtitle("(b)") + geom_point(data=cor, aes (color=site), size=4, shape=19) + labs(y=(expression(paste("N stock (kg m"^"-2",")"))), x=(expression(paste("Manganese concentration (mg Mn kg"^" -1",")")))) + theme_bw() + xlim(0,2500) + ylim(0,1.1) + theme(legend.position = "none", legend.text=element_text(size=14), axis.title.x=element_text(colour="black", face="bold", size=14), axis.title.y = element_text(colour="black", face="bold", size=16), axis.text = element_text(colour="black", size=14)) + geom_line(aes(y=predlmN), size=1, linetype = "dashed", color = "black") + scale_color_manual(name= "NEON sites", values = mycolors)


#Figure S1

ggarrange(plot1_sites, plot2_sites, nrow=2, common.legend = TRUE)