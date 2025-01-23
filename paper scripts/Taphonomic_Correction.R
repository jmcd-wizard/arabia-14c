#load packages
library(tidyverse)

#brute force stops scientific numbers in favour of full digits.
options(scipen=999)

##read the curves written in the 'Arabia Analysis' script.
All_Curves = readRDS("rds/All_Curves.rds")


##this is the taphonomic correction (equation 4) from Bluhm and Surovell 2019. Define it for use below.##
taph.loss = (21149.57*(All_Curves$calBP+1788.03)^-1.26)

##this is the old taphonomic correction from Surovell et al. 2009##
#taph.loss = (5.726442*(10^6*(All_Curves$calBP+2176.4)^-1.3925309))


##Here we create a new column with 'taphonomically corrected' values for each curve.##
##Each chunk also normalises the corrected values by dividing by their sum to unity (1)##
All_Curves = All_Curves%>% 
  mutate(corrected.SPD = (All_Curves$normspd_dens/taph.loss),
         sum1 = sum(corrected.SPD),
         corrected.SPD = corrected.SPD/sum1,
         sum1=NULL,
         
         corrected.GMM = (All_Curves$GaussMix_50/taph.loss),
         sum3=sum(corrected.GMM),
         corrected.GMM = corrected.GMM/sum3,
         sum3=NULL,
  )

###Now, we take these and compare the corrected values with their original ###
### counterparts using pearson's correlation coefficient. Also a test plot is generated. ###

##unnormalised spd##
cor.spd = cor.test(All_Curves$normspd_dens, All_Curves$corrected.SPD, method = "pearson")
##plot
Figure_S3a = ggplot(data = All_Curves, aes(x=calBP))+
  geom_line(aes(y=normspd_dens, colour = "SPD"))+
  geom_line(aes(y=corrected.SPD, colour = "Corrected SPD"))+
  annotate("text", x=-4000,y=.0004, label = paste0("rho = ",round(cor.spd$estimate, 3)), color = "black", 
           angle = 0, hjust = 0, size = 2.5, fontface = "bold")+
  annotate("text", x=-4000,y=.000385, label = paste0("p-value = <0.05"), color = "black", 
           angle = 0, hjust = 0, size = 2.5, fontface = "bold")+
  scale_y_continuous(limits = c(0,0.00145), breaks=seq(0,0.0014, .0004))+
  scale_x_reverse(limits = c(5200,3200), breaks=seq(5200,3200, -200), sec.axis = sec_axis(~ . *-1+1950, name = "Years (calBC)", breaks = seq(-3300,-1300, 200))) +
  labs(x='Years (calBP)', y= 'Summed Probability')+
  guides(color = guide_legend(title = NULL), fill = guide_legend(title = NULL))+
  theme_bw()+
  theme(plot.margin = unit(c(0.1,0.25,0.1,0.25), "cm"), legend.key.size = unit(.5, 'cm'),legend.position = c(.1315, .925),legend.justification = c("right", "top"),legend.box.just = "left",legend.margin = margin(1, 1, 1, 1),legend.key.width = unit(.65,"cm"), axis.text.y = element_text(angle = 90,hjust=0.5))


##Gaussian Mixture Model##
cor.GMM = cor.test(All_Curves$GaussMix_50, All_Curves$corrected.GMM, method = "pearson")
##plot
Figure_S3b = ggplot(data = All_Curves, aes(x=calBP))+
  geom_line(aes(y=GaussMix_50, colour = "Gaussian Mixture Model 50%"))+
  geom_line(aes(y=corrected.GMM, colour = "Corrected GMM 50%"))+
  annotate("text", x=-4000,y=.0005, label = paste0("rho = ",round(cor.GMM$estimate, 3)), color = "black", 
           angle = 0, hjust = 0, size = 2.5, fontface = "bold")+
  annotate("text", x=-4000,y=.000475, label = paste0("p-value = <0.05"), color = "black", 
           angle = 0, hjust = 0, size = 2.5, fontface = "bold")+
  scale_y_continuous(limits = c(0,0.001), breaks=seq(0,0.001, .0004))+
  scale_x_reverse(limits = c(5200,3200), breaks=seq(5200,3200, -200), sec.axis = sec_axis(~ . *-1+1950, name = "Years (calBC)", breaks = seq(-3300,-1300, 200))) +
  labs(x='Years (calBP)', y= 'Probability')+
  guides(color = guide_legend(title = NULL), fill = guide_legend(title = NULL))+
  theme_bw()+
  theme(plot.margin = unit(c(0.1,0.25,0.1,0.25), "cm"), legend.key.size = unit(.35, 'cm'),legend.position = c(.175, .925),legend.justification = c("right", "top"),legend.box.just = "left",legend.margin = margin(1, 1, 1, 1),legend.key.width = unit(.45,"cm"), axis.text.y = element_text(angle = 90,hjust=0.5))


####connect figure 3a & 3b together into one.
Figure_S3 = ggarrange(Figure_S3a,Figure_S3b,
                      labels = c("A", "B"),
                      ncol = 1, nrow = 2)
#plot figure 3a and 3b
Figure_S3

#save figure S3 plot to file.
ggsave("Figure_S3.pdf",  plot = last_plot(), device = "pdf",
       path = "paper pdf/", width = 3508, height = 1750,
       units = 'px', dpi = 300)






