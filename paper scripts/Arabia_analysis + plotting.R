#load packages
library(rcarbon)
library(tidyverse)
library(ggpubr)


#brute force stops scientific numbers in favour of full digits.
options(scipen=999)

#####IMPORT THE GAUSSIAN MIXTURE MODEL TEST RESULTS#####
#####The script for these can be found in elsewhere in the supplementary material#####
#read the Gaussian Mixture test results
GMM = read_rds("GMM/outputs/K_final/arabia_dens_Kbest.rds")

#####IMPORT THE normalised date Posterior Predictive Model TEST RESULTS#####
Post.Pred.Check = read_rds("rds/Double_Exp_post.pred.check.rds") %>% 
  mutate(time = -3250:-1250)

#####IMPORT THE unnormalised date Posterior Predictive Model TEST RESULTS#####
Post.Pred.Check.unnorm = read_rds("rds/Double_Exp_post.pred.check.unnorm.rds") %>% 
  mutate(time = -3250:-1250)

###and finally, the MCMC posteriors####
MCMC_posteriors = read_rds("rds/MCMC_Double_Exp_posteriors.rds")

#calculating the mean growth rates 1 and 2 from the MCMC posterior distributions for comparison with GMM Growth rates
r1_mean = mean(c(MCMC_posteriors$r1)*100)
r2_mean = mean(c(MCMC_posteriors$r2)*100)
chp_mean = mean(c(MCMC_posteriors$chp))
BC_chp_mean = mean(c(MCMC_posteriors$chp)*-1+1950)

#Extract the 95% confidence envelope for r1, with median.
mcmc_r1_conf = quantile(MCMC_posteriors$r1, probs=c(.1,.5,.9))
#Extract the 95% confidence envelope for r2, with median.
mcmc_r2_conf = quantile(MCMC_posteriors$r2, probs=c(.1,.5,.9))
#Extract the 95% confidence envelope for chp, with median.
mcmc_chp_conf = quantile(MCMC_posteriors$chp, probs=c(.1,.5,.9))


#join the GMM model to the normalised and unnormalised spds.
joined = GMM %>% 
  left_join(Post.Pred.Check, join_by(time == time)) %>%
  left_join(Post.Pred.Check.unnorm, join_by(time == time, calBP == calBP)) %>% 
  rename(normspd_dens = PrDens.x, unnormspd_dens = PrDens.y) %>% 
  filter(calBP <= 5300 & calBP >=3100) %>% 
  mutate(Rate_50 = Rate_50*100,
         Rate_2.5 = Rate_2.5*100,
         Rate_97.5 = Rate_97.5*100)


####save the resulting joined dataframe to file####
write_rds(joined,"rds/All_Curves.rds")


##### GMM GROWTH RATES
####Table 2: This uses the GMM 50% quantiles to create approximations of average####
####growth for periods, age, and the different millenium. ####
#periods (by the southeastern chronology)
GMM.Growth_Rates.Periods = joined %>% 
  mutate(period = ifelse((time>= -3200&time< -2700), "Hafit",
                         ifelse((time>= -2700&time< -2000), "Umm an-Nar",
                                ifelse((time>= -2000&time< -1300), "Wadi Suq", NA_real_)))) %>% 
  group_by(period) %>% 
  summarise(mean = mean(c(Rate_50, Rate_2.5, Rate_97.5))) %>% 
  filter(period != "NA") 
####save  dataframe to file####
write_csv(GMM.Growth_Rates.Periods,"csv/Results/Growth_Rates_by_Periods.csv")


#archaeological ages, stone, iron, bronze etc
GMM.Growth_Rates.Ages = joined %>% 
  mutate(period = ifelse((time>= -3200&time< -2000), "Early Bronze",
                         ifelse((time>= -2000&time< -1600), "Middle Bronze",
                                ifelse((time>= -1600&time< -1300), "Late Bronze", NA_real_)))) %>% 
  group_by(period) %>% 
  summarise(mean = mean(c(Rate_50, Rate_2.5, Rate_97.5))) %>% 
  filter(period != "NA") 
####save  dataframe to file####
write_csv(GMM.Growth_Rates.Ages,"csv/Results/Growth_Rates_by_Age.csv")

#by millennium
GMM.Growth_Rates.Millennium = joined %>% 
  mutate(period = ifelse((time>= -3200&time< -3000), "4thBC",
                         ifelse((time>= -3000&time< -2000), "3rdBC",
                                ifelse((time>= -2000&time< -1300), "2ndBC", NA_real_)))) %>% 
  group_by(period) %>% 
  summarise(mean = mean(c(Rate_50, Rate_2.5, Rate_97.5))) %>% 
  filter(period != "NA")  
####save  dataframe to file####
write_csv(GMM.Growth_Rates.Millennium,"csv/Results/Growth_Rates_by_Millennium.csv")

#Before and after Change point (SPD DOUBLE EXPONENTIAL MODEL CHANGEPOINT MEAN, SEE POSTERIOR IN FIGURE 2E!!)
GMM.Growth_Rates.Model_rates = joined %>% 
  mutate(period = ifelse((time>= -3200&time< (mean(c(-3593,-3922)*+1)+1950)), "GMM_r1",
                         ifelse((time>= (mean(c(-3593,-3922)*+1)+1950)&time< -1300), "GMM_r2", NA_real_))) %>% 
  group_by(period) %>% 
  summarise(mean = mean(c(Rate_50, Rate_2.5, Rate_97.5))) %>% 
  filter(period != "NA")  
####save  dataframe to file####
write_csv(GMM.Growth_Rates.Model_rates,"csv/Results/Growth_Rates_for_model_comparison.csv")
####


#Calculate the Pearson correlation coefficient between the unnormalised and#
#normalised SPDS, BayDem SPDs, and the midpoint of the Gaussian Mixture model#
#obtained calibrating radiocarbon samples#
pears.cor.table = cor(joined[c(3:5,15,21)], method = "pearson")
#now round to three decimal places
pears.cor.table = round(pears.cor.table, 3)

####Table 1--- save the correlation matrix as a .csv file, which will be####
#####wrangled in EXCEL for the product seen in the final table 1.####
write_csv(as.data.frame(pears.cor.table), "csv/Results/All_Models_Correlation.csv")


#Calculate the Pearson correlation coefficient between the unnormalised and normalised SPDS obtained calibrating radiocarbon samples
normalisation.cor = cor.test(joined$normspd_dens, joined$unnormspd_dens, method = "pearson")
#the two curves are highly correlated. rho=0.8, p-value <0.00001


####plotting!#####

######FIGURE2!!!!###
#
#subset only those being used in correlation
Post.Pred.Cor = joined[c(9:13,15)]
###This provides the correlation results for all proxies compared with 
cor.res = as.numeric(cor(Post.Pred.Cor, method = "pearson")[1,-1])
#calculate the min and max values to give the range of the correlation matrix.
cor.res.min = round(min(cor.res), 3)
cor.res.max = round(max(cor.res), 3)



###figure2 prep and plot####
##### Figure 2a -- SPD curves.####
Figure_2a = ggplot(joined, aes(x = calBP)) +
  geom_line(aes(y = normspd_dens, colour = "Normalised Date SPD"), size=1)+
  geom_line(aes(y = unnormspd_dens, colour = "Unnormalised Date SPD"), size=1)+
  scale_colour_manual(values=c("Unnormalised Date SPD"="#50786B","Normalised Date SPD"="#F4C431"))+
  #this is to show the pearsons correlation between normalised and unnormalised SPDs
  annotate("text", x=5200,y=0.000775, label = paste0("rho = ",round(normalisation.cor$estimate, 3)), color = "grey25", 
           angle = 0, hjust = 0, vjust = 0.5, size = 3.25, fontface = "bold", family = "sans")+
  annotate("text", x=5200,y=0.0007, label = "p-value = <0.05", color = "grey25", 
           angle = 0, hjust = 0, vjust = 0.5, size = 3.25, fontface = "bold", family = "sans")+
  scale_y_continuous(limits = c(0,0.0013), breaks=seq(0,0.0013, .0004))+
  scale_x_reverse(limits = c(5200,3200), breaks=seq(5200,3200, -200), sec.axis = sec_axis(~ . *-1+1950, name = "Years (calBC)", breaks = seq(-3300,-1300, 200))) +
  labs(x='Years (calBP)', y= 'Summed Probability')+
  guides(color = guide_legend(title = NULL), fill = guide_legend(title = NULL))+
  theme_bw()+
  theme(plot.margin = unit(c(0.1,0.25,0.1,0.25), "cm"), legend.key.size = unit(.5, 'cm'),legend.position = c(.165, .925),legend.justification = c("right", "top"),legend.box.just = "left",legend.margin = margin(1, 1, 1, 1),legend.key.width = unit(.65,"cm"), axis.text.y = element_text(angle = 90,hjust=0.5))
#
####a bit of wrangling to make sure the ribbon behaves properly for red and blue colouration in firgure 2b
Post.Pred.Check_joined = joined %>% 
  mutate(above = ifelse((normspd_dens>=`95%.x`), normspd_dens, NA_real_), below = ifelse((normspd_dens<=`5%.x`), normspd_dens, NA_real_))
##### Figure 2a -- SPD Model####
Figure_2b = ggplot(Post.Pred.Check_joined, aes(x = calBP)) +
  geom_ribbon(aes(ymin=`5%.x`, ymax=`95%.x`, fill = "Double Exponential Model Envelope"), alpha = 0.5)+
  geom_ribbon(aes(ymin=`95%.x`, ymax=above, fill = ">95%"), alpha = 0.5)+
  geom_ribbon(aes(ymin=`5%.x`, ymax=below, fill = "<5%"), alpha = 0.5)+
  geom_line(aes(y = normspd_dens, colour = "SPD"), size=1, linetype = "solid")+
  scale_colour_manual(values=c("SPD"="black"))+
  scale_fill_manual(values=c("Double Exponential Model Envelope"="grey75",">95%" = "blue", "<5%" = "red"))+
  ##this is for the extra regional periods
  #Southeastern chronology
  annotate("segment", x = 3200, y = 0.0013, xend = 5150, yend = 0.0013, color = "navyblue", linetype = "solid",linewidth = .5)+
  annotate("segment", x = 5150, y = 0.0013, xend = 5200, yend = 0.0013, color = "navyblue", linetype = "dotted",linewidth = .5)+
  #Hafit start
  annotate("text", x=4900,y=0.00125, label = "Hafit", color = "navyblue", 
           angle = 0, hjust = 0.5, size = 3, fontface = "bold", family = "sans")+
  annotate("segment", x = 5150, xend = 5150, y = 0.0013, yend = 0.001225,  color = "navyblue",linetype="solid",linewidth = .5) +
  #Hafit end
  annotate("segment", x = 4650, xend = 4650, y = 0.0013, yend = 0.001225, color = "navyblue",linetype="solid",linewidth = .5) +
  #Umm an-Nar start
  annotate("text", x=4300,y=0.00125, label = "Umm an-Nar", color = "navyblue", 
           angle = 0, hjust = 0.5, size = 3, fontface = "bold", family = "sans")+
  annotate("segment", x = 4650, xend = 4650, y = 0.0013, yend = 0.001225,  color = "navyblue",linetype="solid",linewidth = .5) +
  #Umm an-Nar end
  annotate("segment", x = 3950, xend = 3950, y = 0.0013, yend = 0.001225,  color = "navyblue",linetype="solid",linewidth = .5) +
  #Wadi Suq start
  annotate("text", x=3775,y=0.00125, label = "Wadi Suq", color = "navyblue", 
           angle = 0, hjust = 0.5, size = 3, fontface = "bold", family = "sans")+
  annotate("segment", x = 4650, xend = 4650, y = 0.0013, yend = 0.001225,  color = "navyblue",linetype="solid",linewidth = .5) +
  #Wadi Suq end
  annotate("segment", x = 3600, xend = 3600, y = 0.0013, yend = 0.001225,  color = "navyblue",linetype="solid",linewidth = .5) +
  #Late Bronze start
  annotate("text", x=3400,y=0.00125, label = "Late Bronze", color = "navyblue", 
           angle = 0, hjust = 0.5, size = 3, fontface = "bold", family = "sans")+
  annotate("segment", x = 3600, xend = 3600, y = 0.0013, yend = 0.001225,  color = "navyblue",linetype="solid",linewidth = .5) +
  #Late Bronze start
  annotate("segment", x = 3200, xend = 3200, y = 0.0013, yend = 0.001225,  color = "navyblue",linetype="solid",linewidth = .5) +
  ##
  #end of captioning
  #this is to show the pearsons correlation between spd and model
  annotate("text", x=5200,y=0.00065, label = paste0("rho = ~",cor.res.min,"-",cor.res.max), color = "grey25", 
           angle = 0, hjust = 0, vjust = 0.5, size = 3.25, fontface = "bold", family = "sans")+
  annotate("text", x=5200,y=0.00056, label = "p-value = <0.05", color = "grey25", 
           angle = 0, hjust = 0, vjust = 0.5, size = 3.25, fontface = "bold", family = "sans")+
  scale_y_continuous(limits = c(0,0.0014), breaks=seq(0,0.0014, .0004))+
  scale_x_reverse(limits = c(5200,3200), breaks=seq(5200,3200, -200), sec.axis = sec_axis(~ . *-1+1950, name = "Years (calBC)", breaks = seq(-3300,-1300, 200))) +
  labs(x='Years (calBP)', y= 'Summed Probability')+
  guides(fill = guide_legend(title = NULL, order = 3),color = guide_legend(title = NULL))+
  theme_bw()+
  theme(plot.margin = unit(c(0.1,0.25,0.1,0.25), "cm"), legend.key.size = unit(.1, 'cm'),legend.position = c(.2025, .4695),legend.justification = c("right", "bottom"),legend.box.just = "left",legend.margin = margin(1, 1, 1, 1),legend.key.width = unit(.25,"cm"), axis.text.y = element_text(angle = 90,hjust=0.5))

###### wrangle for the next three sub-figures####
MCMC_posteriors=MCMC_posteriors %>% 
  mutate(r1_quant = ifelse((r1>=quantile(r1, probs=.15)&r1<=quantile(r1, probs=.95)), r1,NA_real_),
         r2_quant = ifelse((r2>=quantile(r2, probs=.15)&r2<=quantile(r2, probs=.95)), r2,NA_real_),
         chp_quant = ifelse((chp>=quantile(chp, probs=.15)&chp<=quantile(chp, probs=.95)), chp,NA_real_))
##### Figure 2b -- MCMC posterior r1####
Figure_2c = ggplot(MCMC_posteriors)+
  geom_density(aes(r1),colour = "grey25", fill = "lightblue")+
  labs(x='Growth Rate 1', y= 'Frequency')+
  scale_y_continuous(limits = c(0,2750), breaks=seq(0,2500, 500), guide = guide_axis(angle = 90))+
  scale_x_continuous(limits = c(0.00005,0.00075), breaks=seq(0.00005,0.00075, 0.0002))+
  annotate("text", x=0.0004,y=2600, label = paste0("95%HDPI: ",round(mcmc_r1_conf[1],5),"~",round(mcmc_r1_conf[3],5)), color = "grey25", 
           angle = 0, hjust = 0.5, vjust = 0.5, size = 3.25, fontface = "bold", family = "sans")+
  theme_bw()+
  theme(plot.margin = unit(c(0.1,0.25,0.25,0.25), "cm"))

##### Figure 2d -- MCMC posterior r2####
Figure_2d = ggplot(MCMC_posteriors)+
  geom_density(aes(r2),colour = "grey25", fill= "lightblue")+
  labs(x='Growth Rate 2', y= 'Frequency')+
  scale_y_continuous(limits = c(0,1500), breaks=seq(0,1500, 500), guide = guide_axis(angle = 90))+
  scale_x_continuous(limits = c(-0.0011,0.0004), breaks=seq(-0.001,0.0004, 0.0004))+
  annotate("text", x=-0.00035,y=1400, label = paste0("95%HDPI: ",round(mcmc_r2_conf[1],5),"~",round(mcmc_r2_conf[3],5)), color = "grey25", 
           angle = 0, hjust = 0.5, vjust = 0.5, size = 3.25, fontface = "bold", family = "sans")+
  theme_bw()+
  theme(plot.margin = unit(c(0.1,0.25,0.25,0.25), "cm"))

##### Figure 2e -- MCMC posterior chp####
Figure_2e = ggplot(MCMC_posteriors)+
  geom_density(aes(chp), colour = "grey25", fill= "lightblue")+
  labs(x='Changepoint (calBP)', y= 'Density')+
  scale_y_continuous(limits = c(0,0.00325), breaks=seq(0,0.003, 0.001), guide = guide_axis(angle = 90))+
  scale_x_continuous(limits = c(3300,4200), breaks=seq(3300,4200, 100))+
  annotate("text", x=3750,y=0.003, label = paste0("95%HDPI: ",round(mcmc_chp_conf[1],0),"~",round(mcmc_chp_conf[3],0)), color = "grey25", 
           angle = 0, hjust = 0.5, vjust = 0.5, size = 3.25, fontface = "bold", family = "sans")+
  theme_bw()+
  theme(plot.margin = unit(c(0.1,0.25,0.25,0.25), "cm"))
#Figure 2 arrangement#
Figure_2 = ggarrange(Figure_2a, Figure_2b,
                     ggarrange(Figure_2c,Figure_2d,Figure_2e,
                               ncol = 3, labels = c("C", "D", "E")),
                     nrow = 3,labels = c("A", "B"))

##plot figure 2 in its entirety                    
Figure_2
#save figure 2 plot to file.
ggsave("Figure_2.pdf",  plot = last_plot(), device = "pdf",
       path = "paper pdf/", width = 3508, height = 2480,
       units = 'px', dpi = 300)



#####FIGURE 3######
##### Figure 3a -- Gaussian Mixture Model####
Fig3a = ggplot(joined, aes(x = calBP)) +
  geom_line(aes(y = GaussMix_50, colour = "GMM 50%"), size=1, linetype = "dashed")+
  geom_ribbon(aes(ymin=GaussMix_2.5, ymax=GaussMix_97.5, fill = "GMM 2.5%-97.5% Quantiles"), alpha = 0.5)+
  scale_colour_manual(values=c("GMM 50%"="#A91C72"))+
  #Southeastern chronology
  annotate("segment", x = 3200, y = 0.00095, xend = 5150, yend = 0.00095, color = "navyblue", linetype = "solid",linewidth = .5)+
  annotate("segment", x = 5150, y = 0.00095, xend = 5200, yend = 0.00095, color = "navyblue", linetype = "dotted",linewidth = .5)+
  #Hafit start
  annotate("text", x=4900,y=0.0009, label = "Hafit", color = "navyblue", 
           angle = 0, hjust = 0.5, size = 3, fontface = "bold", family = "sans")+
  annotate("segment", x = 5150, xend = 5150, y = 0.00095, yend = 0.000875,  color = "navyblue",linetype="solid",linewidth = .5) +
  #Hafit end
  annotate("segment", x = 4650, xend = 4650, y = 0.00095, yend = 0.000875, color = "navyblue",linetype="solid",linewidth = .5) +
  #Umm an-Nar start
  annotate("text", x=4300,y=0.0009, label = "Umm an-Nar", color = "navyblue", 
           angle = 0, hjust = 0.5, size = 3, fontface = "bold", family = "sans")+
  annotate("segment", x = 4650, xend = 4650, y = 0.00095, yend = 0.000875,  color = "navyblue",linetype="solid",linewidth = .5) +
  #Umm an-Nar end
  annotate("segment", x = 3950, xend = 3950, y = 0.00095, yend = 0.000875,  color = "navyblue",linetype="solid",linewidth = .5) +
  #Wadi Suq start
  annotate("text", x=3775,y=0.0009, label = "Wadi Suq", color = "navyblue", 
           angle = 0, hjust = 0.5, size = 3, fontface = "bold", family = "sans")+
  annotate("segment", x = 4650, xend = 4650, y = 0.00095, yend = 0.000875,  color = "navyblue",linetype="solid",linewidth = .5) +
  #Wadi Suq end
  annotate("segment", x = 3600, xend = 3600, y = 0.00095, yend = 0.000875,  color = "navyblue",linetype="solid",linewidth = .5) +
  #Late Bronze start
  annotate("text", x=3400,y=0.0009, label = "Late Bronze", color = "navyblue", 
           angle = 0, hjust = 0.5, size = 3, fontface = "bold", family = "sans")+
  annotate("segment", x = 3600, xend = 3600, y = 0.00095, yend = 0.000875,  color = "navyblue",linetype="solid",linewidth = .5) +
  #Late Bronze start
  annotate("segment", x = 3200, xend = 3200, y = 0.00095, yend = 0.000875,  color = "navyblue",linetype="solid",linewidth = .5) +
  ##
  #end of captioning
  scale_y_continuous(limits = c(0,0.001), breaks=seq(0,0.001, .0004))+
  scale_x_reverse(limits = c(5200,3200), breaks=seq(5200,3200, -200), sec.axis = sec_axis(~ . *-1+1950, name = "Years (calBC)", breaks = seq(-3300,-1300, 200))) +
  labs(x='Years (calBP)', y= 'Probability')+
  guides(color = guide_legend(title = NULL), fill = guide_legend(title = NULL))+
  theme_bw()+
  theme(plot.margin = unit(c(0.1,0.25,0.1,0.25), "cm"), legend.key.size = unit(.5, 'cm'),legend.position = c(.185, .835),legend.justification = c("right", "top"),legend.box.just = "left",legend.margin = margin(1, 1, 1, 1),legend.key.width = unit(.65,"cm"), axis.text.y = element_text(angle = 90,hjust=0.5))


##
####a bit of wrangling to make sure the ribbon behaves properly for red and blue colouration
joined_rates = joined %>% 
  mutate(above1 = ifelse((Rate_2.5>=0), Rate_2.5,NA_real_), below1 = ifelse((above1!=0), Rate_97.5, NA_real_),
         above2 = ifelse((Rate_97.5<=0), Rate_97.5,NA_real_), below2 = ifelse((above2!=0), Rate_2.5, NA_real_))

##### Figure 3b -- Gaussian Mixture Model####
Fig3b = ggplot(joined_rates, aes(x = calBP)) +
  geom_ribbon(aes(ymin=Rate_2.5, ymax=Rate_97.5, fill = "Growth Rate 2.5%-97.5% Quantiles"), alpha = 0.5)+
  geom_ribbon(aes(ymin=above1, ymax=below1, fill = "2.5% above 0"), alpha = 0.5)+
  geom_ribbon(aes(ymin=above2, ymax=below2, fill = "97.5% below 0"), alpha = 0.5)+
  geom_line(aes(y = Rate_50, colour = "Median Growth Rate"), size=.75, linetype = "solid")+
  geom_hline(yintercept = 0, alpha = 0.3, linetype = 2)+
  scale_colour_manual(values=c("Median Growth Rate"="grey10"))+
  scale_fill_manual(values=c("Growth Rate 2.5%-97.5% Quantiles"="grey75", "2.5% above 0" = "green", "97.5% below 0" = "red"))+
  labs(x='Years (calBP)', y= 'Growth Rate')+
  scale_x_reverse(limits = c(5200,3200), breaks=seq(5200,3200, -200), sec.axis = sec_axis(~ . *-1+1950, name = "Years (calBC)", breaks = seq(-3300,-1300, 200))) +
  guides(color = guide_legend(title = NULL), fill = guide_legend(title = NULL))+
  theme_bw()+
  theme(plot.margin = unit(c(0.1,0.25,0.1,0.25), "cm"), legend.key.size = unit(.35, 'cm'),legend.position = c(.625, .95),legend.justification = c("right", "top"),legend.box.just = "left",legend.margin = margin(1, 1, 1, 1),legend.key.width = unit(.65,"cm"), axis.text.y = element_text(angle = 90,hjust=0.5))

####connect figure 3a & 3b together into one.
Figure_3 = ggarrange(Fig3a,Fig3b,
                     labels = c("A", "B"),
                     ncol = 1, nrow = 2)
#plot figure 3a and 3b
Figure_3
#save figure 3  to file.
ggsave("Figure_3.pdf",  plot = last_plot(), device = "pdf",
       path = "paper pdf/", width = 3508, height = 1750,
       units = 'px', dpi = 300)




####a bit of wrangling to make sure the ribbon behaves properly for red and green colouration
Post.Pred.Check_joined = joined %>% 
  mutate(above = ifelse((normspd_dens>=`95%.x`), normspd_dens, NA_real_), below = ifelse((normspd_dens<=`5%.x`), normspd_dens, NA_real_))
##### Figure S1a  -- SPD Model####
Figure_S1a = ggplot(Post.Pred.Check_joined, aes(x = calBP)) +
  geom_ribbon(aes(ymin=`5%.x`, ymax=`95%.x`, fill = "Double Exponential Model Envelope"), alpha = 0.5)+
  geom_ribbon(aes(ymin=`95%.x`, ymax=above, fill = ">95%"), alpha = 0.5)+
  geom_ribbon(aes(ymin=`5%.x`, ymax=below, fill = "<5%"), alpha = 0.5)+
  geom_line(aes(y = normspd_dens, colour = "SPD"), size=1, linetype = "solid")+
  scale_colour_manual(values=c("SPD"="black"))+
  scale_fill_manual(values=c("Double Exponential Model Envelope"="grey75", ">95%" = "blue", "<5%" = "red"))+
  #this is to show the pearsons correlation between spd and model
  annotate("text", x=5200,y=0.000625, label = paste0("rho = ~",cor.res.min,"-",cor.res.max), color = "grey25", 
           angle = 0, hjust = 0, vjust = 0.5, size = 3.25, fontface = "bold", family = "sans")+
  annotate("text", x=5200,y=0.000555, label = "p-value = <0.05", color = "grey25", 
           angle = 0, hjust = 0, vjust = 0.5, size = 3.25, fontface = "bold", family = "sans")+
  scale_y_continuous(limits = c(0,0.0013), breaks=seq(0,0.0012, .0004))+
  scale_x_reverse(limits = c(5200,3200), breaks=seq(5200,3200, -200), sec.axis = sec_axis(~ . *-1+1950, name = "Years (calBC)", breaks = seq(-3300,-1300, 200))) +
  labs(x='Years (calBP)', y= 'Summed Probability')+
  guides(fill = guide_legend(title = NULL),color = guide_legend(title = NULL))+
  theme_bw()+
  theme(plot.margin = unit(c(0.1,0.25,0.1,0.25), "cm"), legend.key.size = unit(.1, 'cm'),legend.position = c(.20, .60),legend.justification = c("right", "bottom"),legend.box.just = "left",legend.margin = margin(1, 1, 1, 1),legend.key.width = unit(.25,"cm"), axis.text.y = element_text(angle = 90,hjust=0.5))
####a bit of wrangling to make sure the ribbon behaves properly for red and green colouration
Log_Post.Pred.Check = read_rds("rds/Logistic_post_pred_check.rds")
###and the MCMC posteriors
Log_MCMC_posteriors = read_rds("rds/MCMC_Logistic_posteriors.rds")
#
#subset only those being used in correlation
Log_Post.Pred_cor = Log_Post.Pred.Check[c(1:5,7)]
###This provides the correlation results for all proxies compared with 
log.cor.res = as.numeric(cor(Log_Post.Pred_cor, method = "pearson")[1,-1])
#calculate the min and max values to give the range of the correlation matrix.
log_cor.res.min = round(min(log.cor.res), 3)
log_cor.res.max = round(max(log.cor.res), 3)
#
Log_Post.Pred.Check = Log_Post.Pred.Check %>% 
  mutate(above = ifelse((PrDens>=`95%`), PrDens, NA_real_), below = ifelse((PrDens<=`5%`), PrDens, NA_real_))
##### Figure S1b -- SPD Model####
Figure_S1b = ggplot(Log_Post.Pred.Check, aes(x = calBP)) +
  geom_ribbon(aes(ymin=`5%`, ymax=`95%`, fill = "Logistic Model Envelope"), alpha = 0.5)+
  geom_ribbon(aes(ymin=`95%`, ymax=above, fill = ">95%"), alpha = 0.5)+
  geom_ribbon(aes(ymin=`5%`, ymax=below, fill = "<5%"), alpha = 0.5)+
  geom_line(aes(y = PrDens, colour = "SPD"), size=1, linetype = "solid")+
  scale_colour_manual(values=c("SPD"="black"))+
  scale_fill_manual(values=c("Logistic Model Envelope"="grey75", ">95%" = "blue", "<5%" = "red"))+
  #this is to show the pearsons correlation between spd and model
  annotate("text", x=5200,y=0.000635, label = paste0("rho = ~",log_cor.res.min,"-",log_cor.res.max), color = "grey25", 
           angle = 0, hjust = 0, vjust = 0.5, size = 3.25, fontface = "bold", family = "sans")+
  annotate("text", x=5200,y=0.000555, label = "p-value = <0.05", color = "grey25", 
           angle = 0, hjust = 0, vjust = 0.5, size = 3.25, fontface = "bold", family = "sans")+
  scale_y_continuous(limits = c(0,0.0013), breaks=seq(0,0.0012, .0004))+
  scale_x_reverse(limits = c(5200,3200), breaks=seq(5200,3200, -200), sec.axis = sec_axis(~ . *-1+1950, name = "Years (calBC)", breaks = seq(-3300,-1300, 200))) +
  labs(x='Years (calBP)', y= 'Summed Probability')+
  guides(fill = guide_legend(title = NULL),color = guide_legend(title = NULL))+
  theme_bw()+
  theme(plot.margin = unit(c(0.1,0.25,0.1,0.25), "cm"), legend.key.size = unit(.225, 'cm'),legend.position = c(.17, .565),legend.justification = c("right", "bottom"),legend.box.just = "left",legend.margin = margin(1, 1, 1, 1),legend.key.width = unit(.4,"cm"), axis.text.y = element_text(angle = 90,hjust=0.5))

####connect figure 2a & 2b together into one.
Figure_S1 = ggarrange(Figure_S1a,Figure_S1b,
                      labels = c("A", "B"),
                      ncol = 1, nrow = 2)
#plot figure 2a and 2b
Figure_S1
#save figure S1 plot to file.
ggsave("Figure_S1.pdf",  plot = last_plot(), device = "pdf",
       path = "paper pdf/", width = 3508, height = 1750,
       units = 'px', dpi = 300)



#wrangle for S2
#calculating the mean growth rates 1 and 2 from the MCMC posterior distributions for comparison with GMM Growth rates
r_mean = mean(c(Log_MCMC_posteriors$r)*100)
k_mean = mean(c(Log_MCMC_posteriors$k)*100)

#Extract the 95% confidence envelope for r1, with median.
mcmc_r_conf = quantile(Log_MCMC_posteriors$r, probs=c(.1,.5,.9))
#Extract the 95% confidence envelope for r2, with median.
mcmc_k_conf = quantile(Log_MCMC_posteriors$k, probs=c(.1,.5,.9))
##### Figure S2a -- MCMC posterior r####
Figure_S2a = ggplot(Log_MCMC_posteriors)+
  geom_density(aes(r),colour = "grey25", fill= "lightblue")+
  labs(x='Growth Rate', y= 'Frequency')+
  scale_y_continuous(limits = c(0,4600), breaks=seq(0,4500, 1000), guide = guide_axis(angle = 90))+
  scale_x_continuous(limits = c(-0.0002,0.0005), breaks=seq(-0.0002,0.0006, 0.0002))+
  annotate("text", x=0.000155,y=4350, label = paste0("95%HDPI: ",round(mcmc_r_conf[1],5),"~",round(mcmc_r_conf[3],5)), color = "grey25", 
           angle = 0, hjust = 0.5, vjust = 0.5, size = 3.25, fontface = "bold", family = "sans")+
  theme_bw()+
  theme(plot.margin = unit(c(0.1,0.25,0.25,0.25), "cm"))

##### Figure S2b -- MCMC posterior k####
Figure_S2b = ggplot(Log_MCMC_posteriors)+
  geom_density(aes(as.numeric(k)),colour = "grey25", fill= "lightblue")+
  labs(x='Carrying Capacity (k)', y= 'Frequency')+
  scale_y_continuous(limits = c(0,11.5), breaks=seq(1,13, 2), guide = guide_axis(angle = 90))+
  scale_x_continuous(limits = c(-0.02,0.21), breaks=seq(0,0.2, 0.06))+
  annotate("text", x=0.1,y=11, label = paste0("95%HDPI: ",round(mcmc_k_conf[1],5),"~",round(mcmc_k_conf[3],5)), color = "grey25", 
           angle = 0, hjust = 0.5, vjust = 0.5, size = 3.25, fontface = "bold", family = "sans")+
  theme_bw()+
  theme(plot.margin = unit(c(0.1,0.25,0.25,0.25), "cm"))


####connect figure 3a & 3b together into one.
Figure_S2 = ggarrange(Figure_S2a,Figure_S2b,
                      labels = c("A", "B"),
                      ncol = 2, nrow = 1)
#plot figure 3a and 3b
Figure_S2

#save figure S2 plot to file.
ggsave("Figure_S2.pdf",  plot = last_plot(), device = "pdf",
       path = "paper pdf/", width = 3508, height = 1750,
       units = 'px', dpi = 300)




