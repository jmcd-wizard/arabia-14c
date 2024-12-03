#load packages
library(tidyverse)

#brute force stops scientific numbers in favour of full digits.
options(scipen=999)

##read the curves which originate in the 'rcarbon spd + kde' r-script
All_Curves = readRDS("Master Paper/rds/All_Curves.rds")


##this is the taphonomic correction (equation 4) from Bluhm and Surovell 2019. Define it for use below.##
taph.loss = (21149.57*(All_Curves$calBP+1788.03)^-1.26)

##this is the old taphonomic correction from Surovell et al. 2009##
#taph.loss = (5.726442*(10^6*(All_Curves$calBP+2176.4)^-1.3925309))


##Here we create a new column with 'taphonomically corrected' values for each curve.##
##Each chunk also normalises the corrected values by dividing by their sum to unity (1)##
All_Curves = All_Curves%>% 
  mutate(corrected.SPD = (All_Curves$spd_dens/taph.loss),
         sum1 = sum(corrected.SPD),
         corrected.SPD = corrected.SPD/sum1,
         sum1=NULL,
         
         corrected.NSPD = (All_Curves$nspd_dens/taph.loss),
         sum2=sum(corrected.NSPD),
         corrected.NSPD = corrected.NSPD/sum2,
         sum2=NULL,
         
         corrected.GMM = (All_Curves$GaussMix_50/taph.loss),
         sum3=sum(corrected.GMM),
         corrected.GMM = corrected.GMM/sum3,
         sum3=NULL,
         
         corrected.cKDE = (All_Curves$cKDE.mean/taph.loss),
         sum4=sum(corrected.cKDE),
         corrected.cKDE = corrected.cKDE/sum4,
         sum4=NULL,
         
         corrected.Oxcal = (All_Curves$probability/taph.loss),
         sum5=sum(corrected.Oxcal),
         corrected.Oxcal = corrected.Oxcal/sum5,
         sum5=NULL,
  )

###Now, we take these and compare the corrected values with their original ###
### counterparts using pearson's correlation coefficient. Also a test plot is generated. ###

##unnormalised spd##
cor.spd = cor.test(All_Curves$spd_dens, All_Curves$corrected.SPD, method = "pearson")
##plot
Figure_S4a = ggplot(data = All_Curves, aes(x=calBP))+
  geom_line(aes(y=spd_dens, colour = "SPD"))+
  geom_line(aes(y=corrected.SPD, colour = "Corrected SPD"))+
  annotate("text", x=-4000,y=.0004, label = paste0("rho = ",round(cor.spd$estimate, 3)), color = "black", 
           angle = 0, hjust = 0, size = 2.5, fontface = "bold")+
  annotate("text", x=-4000,y=.000385, label = paste0("p-value = <0.05"), color = "black", 
           angle = 0, hjust = 0, size = 2.5, fontface = "bold")+
  scale_y_continuous(limits = c(0,0.0008), breaks=seq(0,0.0008, .0004))+
  scale_x_reverse(limits = c(5200,3200), breaks=seq(5200,3200, -200), sec.axis = sec_axis(~ . *-1+1950, name = "Years (calBC)", breaks = seq(-3300,-1300, 200))) +
  labs(x='Years (calBP)', y= 'Summed Probability')+
  guides(color = guide_legend(title = NULL), fill = guide_legend(title = NULL))+
  theme_bw()+
  theme(plot.margin = unit(c(0.1,0.25,0.1,0.25), "cm"), legend.key.size = unit(.5, 'cm'),legend.position = c(.1315, .925),legend.justification = c("right", "top"),legend.box.just = "left",legend.margin = margin(1, 1, 1, 1),legend.key.width = unit(.65,"cm"), axis.text.y = element_text(angle = 90,hjust=0.5))



##normalised spd##
cor.nspd = cor.test(All_Curves$nspd_dens, All_Curves$corrected.NSPD, method = "pearson")
##plot
Figure_S4b = ggplot(data = All_Curves, aes(x=calBP))+
  geom_line(aes(y=corrected.NSPD, colour = "Corrected Normalised SPD"))+
  geom_line(aes(y=nspd_dens, colour = "Normalised SPD"))+
  annotate("text", x=-4000,y=.0005, label = paste0("rho = ",round(cor.nspd$estimate, 3)), color = "black", 
           angle = 0, hjust = 0, size = 2.5, fontface = "bold")+
  annotate("text", x=-4000,y=.000475, label = paste0("p-value = <0.05"), color = "black", 
           angle = 0, hjust = 0, size = 2.5, fontface = "bold")+
  scale_y_continuous(limits = c(0,0.0008), breaks=seq(0,0.0008, .0004))+
  scale_x_reverse(limits = c(5200,3200), breaks=seq(5200,3200, -200), sec.axis = sec_axis(~ . *-1+1950, name = "Years (calBC)", breaks = seq(-3300,-1300, 200))) +
  labs(x='Years (calBP)', y= 'Summed Probability')+
  guides(color = guide_legend(title = NULL), fill = guide_legend(title = NULL))+
  theme_bw()+
  theme(plot.margin = unit(c(0.1,0.25,0.1,0.25), "cm"), legend.key.size = unit(.5, 'cm'),legend.position = c(.175, .925),legend.justification = c("right", "top"),legend.box.just = "left",legend.margin = margin(1, 1, 1, 1),legend.key.width = unit(.65,"cm"), axis.text.y = element_text(angle = 90,hjust=0.5))



##Gaussian Mixture Model##
cor.GMM = cor.test(All_Curves$GaussMix_50, All_Curves$corrected.GMM, method = "pearson")
##plot
Figure_S4c = ggplot(data = All_Curves, aes(x=calBP))+
  geom_line(aes(y=GaussMix_50, colour = "Gaussian Mixture Model"))+
  geom_line(aes(y=corrected.GMM, colour = "Corrected GMM"))+
  annotate("text", x=-4000,y=.0005, label = paste0("rho = ",round(cor.GMM$estimate, 3)), color = "black", 
           angle = 0, hjust = 0, size = 2.5, fontface = "bold")+
  annotate("text", x=-4000,y=.000475, label = paste0("p-value = <0.05"), color = "black", 
           angle = 0, hjust = 0, size = 2.5, fontface = "bold")+
  scale_y_continuous(limits = c(0,0.0008), breaks=seq(0,0.0008, .0004))+
  scale_x_reverse(limits = c(5200,3200), breaks=seq(5200,3200, -200), sec.axis = sec_axis(~ . *-1+1950, name = "Years (calBC)", breaks = seq(-3300,-1300, 200))) +
  labs(x='Years (calBP)', y= 'Probability')+
  guides(color = guide_legend(title = NULL), fill = guide_legend(title = NULL))+
  theme_bw()+
  theme(plot.margin = unit(c(0.1,0.25,0.1,0.25), "cm"), legend.key.size = unit(.5, 'cm'),legend.position = c(.1615, .925),legend.justification = c("right", "top"),legend.box.just = "left",legend.margin = margin(1, 1, 1, 1),legend.key.width = unit(.65,"cm"), axis.text.y = element_text(angle = 90,hjust=0.5))


####connect figure 3a & 3b together into one.
Figure_S4 = ggarrange(Figure_S4a,Figure_S4b,Figure_S4c,
                      labels = c("A", "B", "C"),
                      ncol = 1, nrow = 3)
#plot figure 3a and 3b
Figure_S4

#save figure S3 plot to file.
ggsave("Figure_S4.pdf",  plot = last_plot(), device = "pdf",
       path = "Master Paper/paper pdf/", width = 3508, height = 1750,
       units = 'px', dpi = 300)




###extra plots for general interest.
##cKDE##
cor.ckde = cor.test(All_Curves$cKDE.mean, All_Curves$corrected.cKDE, method = "pearson")
##plot
ggplot(data = All_Curves, aes(x=time))+
  geom_line(aes(y=corrected.cKDE, colour = "Corrected cKDE"))+
  geom_line(aes(y=cKDE.mean, colour = "cKDE"))+
  annotate("text", x=-4000,y=.0005, label = paste0("rho = ",round(cor.ckde$estimate, 3)), color = "black", 
           angle = 0, hjust = 0, size = 2.5, fontface = "bold")+
  annotate("text", x=-4000,y=.000485, label = paste0("p-value = <0.05"), color = "black", 
           angle = 0, hjust = 0, size = 2.5, fontface = "bold")



##Bayesian cKDE##
cor.oxcal = cor.test(All_Curves$probability, All_Curves$corrected.Oxcal, method = "pearson")
##plot
ggplot(data = All_Curves, aes(x=time))+
  geom_line(aes(y=corrected.Oxcal, colour = "Corrected Oxcal cKDE"))+
  geom_line(aes(y=probability, colour = "Oxcal cKDE"))+
  annotate("text", x=-2700,y=.0005, label = paste0("rho = ",round(cor.oxcal$estimate, 3)), color = "black", 
           angle = 0, hjust = 0, size = 2.5, fontface = "bold")+
  annotate("text", x=-2700,y=.000475, label = paste0("p-value = <0.05"), color = "black", 
           angle = 0, hjust = 0, size = 2.5, fontface = "bold")


#plot all values, in the same plot as the normal values
ggplot(All_Curves, aes(x = time)) +
  geom_line(aes(y = corrected.cKDE, colour = "Composite KDE"), size=1)+
  geom_line(aes(y = corrected.GMM, colour = "Gaussian Mixture"), size=1, linetype = 2)+
  geom_line(aes(y = corrected.NSPD, colour = "Normalised SPD"), size=1)+
  geom_line(aes(y = corrected.SPD, colour = "Unnormalised SPD"), size=1)+
  geom_line(aes(y = corrected.Oxcal, colour = "Oxcal cKDE"), size=1)+
  annotate("rect", xmin = -4500, xmax = -3200, ymin = 0, ymax = 0.0012, alpha = 0.2, fill= "gold") + #Neolithic
  annotate("rect", xmin = -3200, xmax = -2700, ymin = 0, ymax = 0.0012, alpha = 0.2, fill= "blue") + #wadi suq
  annotate("rect", xmin = -2700, xmax = -2000, ymin = 0, ymax = 0.0012, alpha = 0.2, fill= "darkgreen") + #umm an-nar
  annotate("rect", xmin = -2000, xmax = -1600, ymin = 0, ymax = 0.0012, alpha = 0.2, fill= "red") + #middle bronze
  annotate("rect", xmin = -1600, xmax = -1300, ymin = 0, ymax = 0.0012, alpha = 0.2, fill= "black") + #late bronze
  annotate("rect", xmin = -1300, xmax = -500, ymin = 0, ymax = 0.0012, alpha = 0.2, fill= "purple") + #early iron age
  scale_y_continuous(limits = c(0,0.0012), breaks=seq(0,0.0012, .0002))+
  scale_x_continuous(limits = c(-4500,-500), breaks=seq(-4500,-500, 250), sec.axis = sec_axis(~ . *-1+1950, name = "Years (calBP)", breaks = seq(1000,8000, 250))) +
  labs(x='Years (calBC)', y= 'Probability')+
  guides(color = guide_legend(title = "With Taphonomic Correction"), fill = guide_legend(title = 'Mixture Envelope'))+
  theme(plot.margin = unit(c(0,0,0,0), "cm"))+
  theme_bw()+
  theme(legend.position = c(.185, .95),legend.justification = c("right", "top"),legend.box.just = "left",legend.margin = margin(1, 1, 1, 1),legend.key.width = unit(.6,"cm"))
