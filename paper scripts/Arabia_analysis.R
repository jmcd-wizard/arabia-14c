#load packages
library(rcarbon)
library(tidyverse)
library(ggpubr)


#brute force stops scientific numbers in favour of full digits.
options(scipen=999)

#load the c14 dates
mydates = read.csv("Master Paper/csv/arabia_dates_20.09.24.csv")

##do a bin sensitivity test before moving on. Adjust as seems fit. 
#test calibration
testcaldates=calibrate(x=mydates$CRA, errors=mydates$Error, calCurves='intcal20', normalised=FALSE, spdnormalised = FALSE, type = 'full')
#actual bin sensitivity test
binsense(x=testcaldates,y=mydates$SiteName,h=seq(0,200,50),timeRange=c(6500,2500)) 
#define bin clustering for rcarbon
bins.n = 50
#define the bin object for use below
bins = binPrep(sites=mydates$SiteName, ages=mydates$CRA, h = bins.n)
#smoothing of SPDs
runm <- 50

####normalised SPDS
#calibrate the normalised dates
ncaldates=calibrate(x=mydates$CRA, errors=mydates$Error, calCurves='intcal20', normalised=TRUE, spdnormalised = TRUE, type = 'full')
#run the normalised spd
nspd = spd(ncaldates,bins = bins, runm=runm, timeRange=c(6500,2500), spdnormalised = TRUE)

#extract the values into a tidyverse dataframe and add the BC dates
extracted_nspd = as.data.frame(nspd$grid) %>% 
  mutate(time = (calBP*-1)+1950)


####UNnormalised SPDS
#calibrate the unnormalised dates
caldates=calibrate(x=mydates$CRA, errors=mydates$Error, calCurves='intcal20', normalised=FALSE, spdnormalised = FALSE, type = 'full')
#run the normalised spd
spd = spd(caldates, bins = bins, runm=runm,timeRange=c(6500,2500), spdnormalised = TRUE)

#extract the values into a tidyverse dataframe and add the BC dates
extracted_spd = as.data.frame(spd$grid) %>% 
  mutate(time = (calBP*-1)+1950)




####RUN THE RCARBON cKDE######
#recalibrate some dates & make sure they are normalised
ckde.caldates=calibrate(x=mydates$CRA, errors=mydates$Error, calCurves='intcal20', normalised=TRUE, datenormalised=FALSE, spdnormalised=TRUE, type = 'full')

#sample random dates from the distribution of each calibrated date
randates = sampleDates(ckde.caldates,bins=bins,nsim=1000)

#run the ckde generation.
ckde = ckde(randates,timeRange=c(6500,2500),bw=50, normalised = TRUE)

#par(mfrow=c(1, 1))
plot(ckde, type = 'envelope')

#extract the posterior-values
cKDE.means = as.data.frame(t(ckde$res.matrix))
#calculate their means & sum to unity to normalise.
cKDE.means = cKDE.means %>%
  summarise(cKDE.mean = colMeans(cKDE.means))%>%
  mutate(calBP = 6500:2500,
         sum = sum(cKDE.mean),
         cKDE.mean = cKDE.mean/sum,
         sum=NULL)
#test to see if the ckde sums to 1
sum(cKDE.means$cKDE.mean)


#####Add and wrangle the Oxcal Bayesian cKDE results. The model generation for #####
#####this method is only available on the Oxcal desktop or browser client: #####
#####https://c14.arch.ox.ac.uk/oxcal/OxCal.html# See Bronk Ramsey 2017 for more #####
#####details regarding workflow. #####

#now add the bayesian ckde results exported from oxcal
OXCAL = read.csv("Master Paper/csv/New_arabia_KDE_model_export.csv")

#wrangle the relevant data and rename the columns.
OXCAL_post= OXCAL %>% 
  filter(op == "KDE_Model"&type == "posterior")%>%
  mutate(time = round(value),
         value = NULL)

##this uses base R functions to linearly interpolate between probabilities.##
##This will not impact the curve as the data resolution is already sub-decadal.##
##Also sum to unity for normalisation.##
OXCAL_post_approx = as.data.frame(approx(OXCAL_post$time, OXCAL_post$probability, xout = seq(min(OXCAL_post$time, na.rm = TRUE), max(OXCAL_post$time, na.rm = TRUE), 1))) %>% 
  mutate(sum = sum(y),
         time = x,
         probability = y/sum,
         x=NULL,
         y=NULL,
         sum=NULL)
#test to see if the ckde sums to 1
sum(OXCAL_post_approx$probability)


#####NOW IMPORT THE OTHER TEST RESULTS, STARTING WITH THE GAUSSIAN MIXTURES#####
#####The script for these can be found in elsewhere in the supplementary material#####
#read the Gaussian Mixture test results
GMM = read_rds("Master Paper/GMM/outputs/K_final/arabia_dens_Kbest.rds")

#join the GM model & cKDE to the normalised and unnormalised spds.
joined = GMM %>% 
  left_join(extracted_nspd) %>% 
  left_join(extracted_spd, join_by(calBP == calBP, time == time)) %>%
  left_join(cKDE.means, join_by(calBP == calBP)) %>%
  left_join(OXCAL_post_approx, join_by(time == time)) %>%
  rename(nspd_dens = PrDens.x, spd_dens = PrDens.y) %>% 
  filter(time >= -3700 & time <=-800) %>% 
  mutate(Rate_50 = Rate_50*100,
         Rate_2.5 = Rate_2.5*100,
         Rate_97.5 = Rate_97.5*100)
  

####save the resulting joined dataframe to file####
write_rds(joined,"Master Paper/rds/All_Curves.rds")





##### GMM GROWTH RATES
####Table 2: This uses the GMM 50% quantiles to create approximations of average####
####growth for periods, age, and the different millenium. ####
#periods (by the southeastern chronology)
GMM.Growth_Rates.Periods = joined %>% 
  mutate(period = ifelse((time>= -3200&time< -2500), "Hafit",
                  ifelse((time>= -2500&time< -2000), "Umm an-Nar",
                  ifelse((time>= -2000&time< -1300), "Wadi Suq", NA_real_)))) %>% 
  group_by(period) %>% 
  summarise(mean = mean(c(Rate_50, Rate_2.5, Rate_97.5))) %>% 
  filter(period != "NA") 
####save  dataframe to file####
write_csv(GMM.Growth_Rates.Periods,"Master Paper/csv/Results/Growth_Rates_by_Periods.csv")


#archaeological ages, stone, iron, bronze etc
GMM.Growth_Rates.Ages = joined %>% 
  mutate(period = ifelse((time>= -3200&time< -2000), "Early Bronze",
                  ifelse((time>= -2000&time< -1600), "Middle Bronze",
                  ifelse((time>= -1600&time< -1300), "Late Bronze", NA_real_)))) %>% 
  group_by(period) %>% 
  summarise(mean = mean(c(Rate_50, Rate_2.5, Rate_97.5))) %>% 
  filter(period != "NA") 
####save  dataframe to file####
write_csv(GMM.Growth_Rates.Ages,"Master Paper/csv/Results/Growth_Rates_by_Age.csv")

#by millennium
GMM.Growth_Rates.Millennium = joined %>% 
  mutate(period = ifelse((time>= -3200&time< -3000), "4thBC",
                  ifelse((time>= -3000&time< -2000), "3rdBC",
                  ifelse((time>= -2000&time< -1300), "2ndBC", NA_real_)))) %>% 
  group_by(period) %>% 
  summarise(mean = mean(c(Rate_50, Rate_2.5, Rate_97.5))) %>% 
  filter(period != "NA")  
####save  dataframe to file####
write_csv(GMM.Growth_Rates.Millennium,"Master Paper/csv/Results/Growth_Rates_by_Millennium.csv")

#Before and after Change point (SPD DOUBLE EXPONENTIAL MODEL CHANGEPOINT MEAN, SEE POSTERIOR IN FIGURE 2E!!)
GMM.Growth_Rates.Model_rates = joined %>% 
  mutate(period = ifelse((time>= -3200&time< (mean(c(-3593,-3922)*+1)+1950)), "GMM_r1",
                         ifelse((time>= (mean(c(-3593,-3922)*+1)+1950)&time< -1300), "GMM_r2", NA_real_))) %>% 
  group_by(period) %>% 
  summarise(mean = mean(c(Rate_50, Rate_2.5, Rate_97.5))) %>% 
  filter(period != "NA")  
####save  dataframe to file####
write_csv(GMM.Growth_Rates.Model_rates,"Master Paper/csv/Results/Growth_Rates_for_model_comparison.csv")
####


#Calculate the Pearson correlation coefficient between the unnormalised and#
#normalised SPDS, BayDem SPDs, and the midpoint of the Gaussian Mixture model#
#obtained calibrating radiocarbon samples#
pears.cor.table = cor(joined[c(3:5,10:11)], method = "pearson")
#now round to three decimal places
pears.cor.table = round(pears.cor.table, 3)

####Table 1--- save the correlation matrix as a .csv file, which will be####
#####wrangled in EXCEL for the product seen in the final table 1.####
write_csv(as.data.frame(pears.cor.table), "Master Paper/csv/Results/All_Models_Correlation.csv")


#Calculate the Pearson correlation coefficient between the unnormalised and normalised SPDS obtained calibrating radiocarbon samples
normalisation.cor = cor.test(spd$grid$PrDens, nspd$grid$PrDens, method = "pearson")
#the two curves are highly correlated. rho=0.8, p-value <0.00001


#For the later plot, calculate the Pearson correlation coefficient between the
#two KDE curves
KDE.cor = cor.test(joined$probability, joined$cKDE.mean, method = "pearson")
#the two curves are highly correlated. rho=0.857, p-value <0.00001

###before plotting, load up the SPD modelling posterior predictive check for tidyverse plotting!
Post.Pred.Check = read_rds("Master Paper/rds/Double_Exp_post.pred.check.rds")
###and the MCMC posteriors
MCMC_posteriors = read_rds("Master Paper/rds/MCMC_Double_Exp_posteriors.rds")


#calculating the mean growth rates 1 and 2 from the MCMC posterior distributions for comparison with GMM Growth rates
r1_mean = mean(c(MCMC_posteriors$r1)*100)
r2_mean = mean(c(MCMC_posteriors$r2)*100)

#Extract the 90% confidence envelope for r1, with median.
mcmc_r1_conf = quantile(MCMC_posteriors$r1, probs=c(.1,.5,.9))
#Extract the 90% confidence envelope for r2, with median.
mcmc_r2_conf = quantile(MCMC_posteriors$r2, probs=c(.1,.5,.9))
#Extract the 90% confidence envelope for chp, with median.
mcmc_chp_conf = quantile(MCMC_posteriors$chp, probs=c(.1,.5,.9))

####plotting!#####

######FIGURE2!!!!###
#sum to unity
SPDs_stu = joined[9:11] %>% 
  filter(calBP<=5200&calBP>=3200) %>% 
  mutate(sum = sum(nspd_dens),
         re_Norm_SPD = nspd_dens/sum,
         sum=NULL,
         sum = sum(spd_dens),
         re_Unnorm_SPD = spd_dens/sum,
         sum=NULL)
##### Figure 2a -- SPD curves.####
Figure_2a = ggplot(SPDs_stu, aes(x = calBP)) +
  geom_line(aes(y = re_Norm_SPD, colour = "Normalised SPD"), size=1)+
  geom_line(aes(y = re_Unnorm_SPD, colour = "Unnormalised SPD"), size=1)+
  scale_colour_manual(values=c("Unnormalised SPD"="#50786B","Normalised SPD"="#F4C431"))+
  #this is to show the pearsons correlation between normalised and unnormalised SPDs
  annotate("text", x=5200,y=0.000775, label = paste0("rho = ",round(normalisation.cor$estimate, 3)), color = "grey25", 
           angle = 0, hjust = 0, vjust = 0.5, size = 3.25, fontface = "bold", family = "sans")+
  annotate("text", x=5200,y=0.0007, label = "p-value = <0.05", color = "grey25", 
           angle = 0, hjust = 0, vjust = 0.5, size = 3.25, fontface = "bold", family = "sans")+
  scale_y_continuous(limits = c(0,0.0012), breaks=seq(0,0.0012, .0004))+
  scale_x_reverse(limits = c(5200,3200), breaks=seq(5200,3200, -200), sec.axis = sec_axis(~ . *-1+1950, name = "Years (calBC)", breaks = seq(-3300,-1300, 200))) +
  labs(x='Years (calBP)', y= 'Summed Probability')+
  guides(color = guide_legend(title = NULL), fill = guide_legend(title = NULL))+
  theme_bw()+
  theme(plot.margin = unit(c(0.1,0.25,0.1,0.25), "cm"), legend.key.size = unit(.5, 'cm'),legend.position = c(.1415, .925),legend.justification = c("right", "top"),legend.box.just = "left",legend.margin = margin(1, 1, 1, 1),legend.key.width = unit(.65,"cm"), axis.text.y = element_text(angle = 90,hjust=0.5))

#
#join the two dataframes
Post.Pred_joined = full_join(Post.Pred.Check, SPDs_stu)
#subset only those being used in correlation
Post.Pred.Cor = Post.Pred_joined[c(1:5,8)]
###This provides the correlation results for all proxies compared with 
cor.res = as.numeric(cor(Post.Pred.Cor, method = "pearson")[1,-1])
#calculate the min and max values to give the range of the correlation matrix.
cor.res.min = round(min(cor.res), 3)
cor.res.max = round(max(cor.res), 3)

####a bit of wrangling to make sure the ribbon behaves properly for red and green colouration
Post.Pred.Check_joined = Post.Pred_joined %>% 
  mutate(above = ifelse((re_Norm_SPD>=`90%`), re_Norm_SPD, NA_real_), below = ifelse((re_Norm_SPD<=`10%`), re_Norm_SPD, NA_real_))
##### Figure 2b -- SPD Model####
Figure_2b = ggplot(Post.Pred.Check_joined, aes(x = calBP)) +
  geom_ribbon(aes(ymin=`10%`, ymax=`90%`, fill = "Double Exponential Model Envelope"), alpha = 0.5)+
  geom_ribbon(aes(ymin=`90%`, ymax=above, fill = ">90%"), alpha = 0.5)+
  geom_ribbon(aes(ymin=`10%`, ymax=below, fill = "<10%"), alpha = 0.5)+
  geom_line(aes(y = re_Norm_SPD, colour = "SPD"), size=1, linetype = "solid")+
  scale_colour_manual(values=c("SPD"="black"))+
  scale_fill_manual(values=c("Double Exponential Model Envelope"="grey75", ">90%" = "blue", "<10%" = "red"))+
 #this is to show the pearsons correlation between spd and model
  annotate("text", x=5200,y=0.000525, label = paste0("rho = ~",cor.res.min,"-",cor.res.max), color = "grey25", 
           angle = 0, hjust = 0, vjust = 0.5, size = 3.25, fontface = "bold", family = "sans")+
  annotate("text", x=5200,y=0.00045, label = "p-value = <0.05", color = "grey25", 
           angle = 0, hjust = 0, vjust = 0.5, size = 3.25, fontface = "bold", family = "sans")+
  scale_y_continuous(limits = c(0,0.0012), breaks=seq(0,0.0012, .0004))+
  scale_x_reverse(limits = c(5200,3200), breaks=seq(5200,3200, -200), sec.axis = sec_axis(~ . *-1+1950, name = "Years (calBC)", breaks = seq(-3300,-1300, 200))) +
  labs(x='Years (calBP)', y= 'Summed Probability')+
  guides(fill = guide_legend(title = NULL),color = guide_legend(title = NULL))+
  theme_bw()+
  theme(plot.margin = unit(c(0.1,0.25,0.1,0.25), "cm"), legend.key.size = unit(.225, 'cm'),legend.position = c(.225, .550),legend.justification = c("right", "bottom"),legend.box.just = "left",legend.margin = margin(1, 1, 1, 1),legend.key.width = unit(.55,"cm"), axis.text.y = element_text(angle = 90,hjust=0.5))

###### wrangle for the next three sub-figures####
MCMC_posteriors=MCMC_posteriors %>% 
  mutate(r1_quant = ifelse((r1>=quantile(r1, probs=.15)&r1<=quantile(r1, probs=.95)), r1,NA_real_),
         r2_quant = ifelse((r2>=quantile(r2, probs=.15)&r2<=quantile(r2, probs=.95)), r2,NA_real_),
         chp_quant = ifelse((chp>=quantile(chp, probs=.15)&chp<=quantile(chp, probs=.95)), chp,NA_real_))
##### Figure 2c -- MCMC posterior r1####
Figure_2c = ggplot(MCMC_posteriors)+
  geom_density(aes(r1),colour = "grey25", fill = "lightblue")+
  labs(x='Growth Rate 1', y= 'Frequency')+
  scale_y_continuous(limits = c(0,3000), breaks=seq(0,3000, 1000), guide = guide_axis(angle = 90))+
  scale_x_continuous(limits = c(0.0003,0.0016), breaks=seq(0.0003,0.0016, 0.0003))+
  annotate("text", x=0.00095,y=2900, label = paste0("90%HDPI: ",round(mcmc_r1_conf[1],5),"~",round(mcmc_r1_conf[3],5)), color = "grey25", 
           angle = 0, hjust = 0.5, vjust = 0.5, size = 3.25, fontface = "bold", family = "sans")+
  theme_bw()+
  theme(plot.margin = unit(c(0.1,0.25,0.25,0.25), "cm"))

##### Figure 2d -- MCMC posterior r2####
Figure_2d = ggplot(MCMC_posteriors)+
  geom_density(aes(r2),colour = "grey25", fill= "lightblue")+
  labs(x='Growth Rate 2', y= 'Frequency')+
  scale_y_continuous(limits = c(0,1500), breaks=seq(0,1500, 500), guide = guide_axis(angle = 90))+
  scale_x_continuous(limits = c(-0.00175,0.000275), breaks=seq(-0.0017,0.00027, 0.0004))+
  annotate("text", x=-0.0007375,y=1450, label = paste0("90%HDPI: ",round(mcmc_r2_conf[1],5),"~",round(mcmc_r2_conf[3],5)), color = "grey25", 
           angle = 0, hjust = 0.5, vjust = 0.5, size = 3.25, fontface = "bold", family = "sans")+
  theme_bw()+
  theme(plot.margin = unit(c(0.1,0.25,0.25,0.25), "cm"))

##### Figure 2e -- MCMC posterior chp####
Figure_2e = ggplot(MCMC_posteriors)+
  geom_density(aes(chp), colour = "grey25", fill= "lightblue")+
  labs(x='Changepoint (calBP)', y= 'Density')+
  scale_y_continuous(limits = c(0,0.0055), breaks=seq(0,0.006, 0.002), guide = guide_axis(angle = 90))+
  scale_x_continuous(limits = c(3400,4100), breaks=seq(3400,4200, 200))+
  annotate("text", x=3750,y=0.00525, label = paste0("90%HDPI: ",round(mcmc_chp_conf[1],0),"~",round(mcmc_chp_conf[3],0)), color = "grey25", 
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
       path = "Master Paper/paper pdf/", width = 3508, height = 2480,
       units = 'px', dpi = 300)
       




#####FIGURE 3######
##### Figure 3a -- Gaussian Mixture Model####
Fig3a = ggplot(joined, aes(x = calBP)) +
  geom_line(aes(y = GaussMix_50, colour = "GMM 50%"), size=1, linetype = "dashed")+
  geom_ribbon(aes(ymin=GaussMix_2.5, ymax=GaussMix_97.5, fill = "GMM 2.5%-97.5% Quantiles"), alpha = 0.5)+
  scale_colour_manual(values=c("GMM 50%"="#A91C72"))+
  scale_y_continuous(limits = c(0,0.001), breaks=seq(0,0.001, .0004))+
  scale_x_reverse(limits = c(5200,3200), breaks=seq(5200,3200, -200), sec.axis = sec_axis(~ . *-1+1950, name = "Years (calBC)", breaks = seq(-3300,-1300, 200))) +
  labs(x='Years (calBP)', y= 'Probability')+
  guides(color = guide_legend(title = NULL), fill = guide_legend(title = NULL))+
  theme_bw()+
  theme(plot.margin = unit(c(0.1,0.25,0.1,0.25), "cm"), legend.key.size = unit(.5, 'cm'),legend.position = c(.195, .85),legend.justification = c("right", "top"),legend.box.just = "left",legend.margin = margin(1, 1, 1, 1),legend.key.width = unit(.65,"cm"), axis.text.y = element_text(angle = 90,hjust=0.5))


##
####a bit of wrangling to make sure the ribbon behaves properly for red and green colouration
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
  theme(plot.margin = unit(c(0.1,0.25,0.1,0.25), "cm"), legend.key.size = unit(.35, 'cm'),legend.position = c(.75, .425),legend.justification = c("right", "top"),legend.box.just = "left",legend.margin = margin(1, 1, 1, 1),legend.key.width = unit(.65,"cm"), axis.text.y = element_text(angle = 90,hjust=0.5))

####connect figure 3a & 3b together into one.
Figure_3 = ggarrange(Fig3a,Fig3b,
                         labels = c("A", "B"),
                         ncol = 1, nrow = 2)
#plot figure 3a and 3b
Figure_3
#save figure 3  to file.
ggsave("Figure_3.pdf",  plot = last_plot(), device = "pdf",
       path = "Master Paper/paper pdf/", width = 3508, height = 1750,
       units = 'px', dpi = 300)


##### Figure 4 -- all methods plotted at once.####
ggplot(joined, aes(x = calBP)) +
  geom_ribbon(aes(ymin=GaussMix_2.5, ymax=GaussMix_97.5, fill = "GMM 2.5%-97.5% Quantiles"), alpha = 0.5)+
  geom_line(aes(y = GaussMix_50, colour = "GMM 50%"), size=1, linetype = "dashed")+
  geom_line(aes(y = spd_dens, colour = "Unnormalised SPD"), size=1)+
  geom_line(aes(y = nspd_dens, colour = "Normalised SPD"), size=1)+
  scale_colour_manual(values=c("GMM 50%"="#A91C72","Unnormalised SPD"="#50786B","Normalised SPD"="#F4C431"))+
  ##this is for the extra regional periods
  #Southeastern
  annotate("text", x=5175,y=0.000825, label = "Southeast", color = "navyblue", 
           angle = 0, hjust = 0.0, size = 4, fontface = "bold", family = "sans")+
  annotate("segment", x = 3200, y = 0.000785, xend = 5150, yend = 0.000785, color = "navyblue", linetype = "solid",linewidth = .5)+
  annotate("segment", x = 5150, y = 0.000785, xend = 5200, yend = 0.000785, color = "navyblue", linetype = "dotted",linewidth = .5)+
  #Hafit start
  annotate("text", x=4900,y=0.00075, label = "Hafit", color = "navyblue", 
           angle = 0, hjust = 0.5, size = 3, fontface = "bold", family = "sans")+
  annotate("segment", x = 5150, xend = 5150, y = 0.000785, yend = 0.000755,  color = "navyblue",linetype="solid",linewidth = .5) +
  #Hafit end
  annotate("segment", x = 4650, xend = 4650, y = 0.000785, yend = 0.000755, color = "navyblue",linetype="solid",linewidth = .5) +
  #Umm an-Nar start
  annotate("text", x=4300,y=0.00075, label = "Umm an-Nar", color = "navyblue", 
           angle = 0, hjust = 0.5, size = 3, fontface = "bold", family = "sans")+
  annotate("segment", x = 4650, xend = 4650, y = 0.000785, yend = 0.000755,  color = "navyblue",linetype="solid",linewidth = .5) +
  #Umm an-Nar end
  annotate("segment", x = 3950, xend = 3950, y = 0.000785, yend = 0.000755,  color = "navyblue",linetype="solid",linewidth = .5) +
  #Wadi Suq start
  annotate("text", x=3775,y=0.00075, label = "Wadi Suq", color = "navyblue", 
           angle = 0, hjust = 0.5, size = 3, fontface = "bold", family = "sans")+
  annotate("segment", x = 4650, xend = 4650, y = 0.000785, yend = 0.000755,  color = "navyblue",linetype="solid",linewidth = .5) +
  #Wadi Suq end
  annotate("segment", x = 3600, xend = 3600, y = 0.000785, yend = 0.000755,  color = "navyblue",linetype="solid",linewidth = .5) +
  #Late Bronze start
  annotate("text", x=3400,y=0.00075, label = "Late Bronze", color = "navyblue", 
           angle = 0, hjust = 0.5, size = 3, fontface = "bold", family = "sans")+
  annotate("segment", x = 3600, xend = 3600, y = 0.000785, yend = 0.000755,  color = "navyblue",linetype="solid",linewidth = .5) +
  #Late Bronze start
  annotate("segment", x = 3200, xend = 3200, y = 0.000785, yend = 0.000755,  color = "navyblue",linetype="solid",linewidth = .5) +
  ##
  #Southwestern
  annotate("text", x=5175,y=0.001125, label = "Southwest", color = "navyblue", 
           angle = 0, hjust = 0.0, size = 4, fontface = "bold", family = "sans")+
  annotate("segment", x = 5150, y = 0.001085, xend = 5200, yend = 0.001085, color = "navyblue", linetype = "dotted",linewidth = .5)+
  annotate("segment", x = 3200, y = 0.001085, xend = 5150, yend = 0.001085, color = "navyblue", linetype = "solid",linewidth = .5)+
  #Neolithic end
  annotate("text", x=5025,y=0.00105, label = "Neolithic", color = "navyblue", 
           angle = 0, hjust = 0.5, size = 3, fontface = "bold", family = "sans")+
  annotate("segment", x = 4850, xend = 4850, y = 0.001085, yend = 0.001055,  color = "navyblue",linetype="solid",linewidth = .5) +
  #Bronze age end
  annotate("text", x=4025,y=0.00105, label = "Bronze Age", color = "navyblue", 
           angle = 0, hjust = 0.5, size = 3, fontface = "bold", family = "sans")+
  annotate("segment", x = 3200, xend = 3200, y = 0.001085, yend = 0.001055,  color = "navyblue",linetype="solid",linewidth = .5) +
  ##
  #Northwest
  annotate("text", x=5175,y=0.001425, label = "Northwest", color = "navyblue", 
           angle = 0, hjust = 0.0, size = 4, fontface = "bold", family = "sans")+
  annotate("segment", x = 5150, y = 0.001385, xend = 5200, yend = 0.001385, color = "navyblue", linetype = "dotted",linewidth = .5)+
  annotate("segment", x = 3200, y = 0.001385, xend = 5150, yend = 0.001385, color = "navyblue", linetype = "solid",linewidth = .5)+
  #Late EB I end
  annotate("text", x=5025,y=0.00135, label = "Late EB I", color = "navyblue", 
           angle = 0, hjust = 0.5, size = 3, fontface = "bold", family = "sans")+
  annotate("segment", x = 4850, xend = 4850, y = 0.001385, yend = 0.001355,  color = "navyblue",linetype="solid",linewidth = .5) +
  #EB II end
  annotate("text", x=4750,y=0.00135, label = "EB II", color = "navyblue", 
           angle = 0, hjust = 0.5, size = 3, fontface = "bold", family = "sans")+
  annotate("segment", x = 4650, xend = 4650, y = 0.001385, yend = 0.001355,  color = "navyblue",linetype="solid",linewidth = .5) +
  #EB III end
  annotate("text", x=4450,y=0.00135, label = "EB III", color = "navyblue", 
           angle = 0, hjust = 0.5, size = 3, fontface = "bold", family = "sans")+
  annotate("segment", x = 4250, xend = 4250, y = 0.001385, yend = 0.001355,  color = "navyblue",linetype="solid",linewidth = .5) +
  #EB IV end
  annotate("text", x=4100,y=0.00135, label = "EB IV", color = "navyblue", 
           angle = 0, hjust = 0.5, size = 3, fontface = "bold", family = "sans")+
  annotate("segment", x = 3950, xend = 3950, y = 0.001385, yend = 0.001355,  color = "navyblue",linetype="solid",linewidth = .5) +
  #MBA end
  annotate("text", x=3750,y=0.00135, label = "Middle Bronze", color = "navyblue", 
           angle = 0, hjust = 0.5, size = 3, fontface = "bold", family = "sans")+
  annotate("segment", x = 3550, xend = 3550, y = 0.001385, yend = 0.001355,  color = "navyblue",linetype="solid",linewidth = .5) +
  #LBA end
  annotate("text", x=3375,y=0.00135, label = "Late Bronze", color = "navyblue", 
           angle = 0, hjust = 0.5, size = 3, fontface = "bold", family = "sans")+
  annotate("segment", x = 3200, xend = 3200, y = 0.001385, yend = 0.001355,  color = "navyblue",linetype="solid",linewidth = .5) +
  ##
  #Bahrain
  annotate("text", x=5175,y=0.001725, label = "Bahrain", color = "navyblue", 
           angle = 0, hjust = 0.0, size = 4, fontface = "bold", family = "sans")+
  annotate("segment", x = 5150, y = 0.001685, xend = 5200, yend = 0.001685, color = "navyblue", linetype = "dotted",linewidth = .5)+
  annotate("segment", x = 3250, y = 0.001685, xend = 5150, yend = 0.001685, color = "navyblue", linetype = "solid",linewidth = .5)+
  annotate("segment", x = 3250, y = 0.001685, xend = 3200, yend = 0.001685, color = "navyblue", linetype = "dotted",linewidth = .5)+
 #Early bronze start
  annotate("text", x=4700,y=0.00165, label = "Late Bronze", color = "navyblue", 
           angle = 0, hjust = 0.5, size = 3, fontface = "bold", family = "sans")+
  annotate("segment", x = 4200, xend = 4200, y = 0.001685, yend = 0.001655,  color = "navyblue",linetype="solid",linewidth = .5) +
  #Early Dilmun start
  annotate("text", x=3875,y=0.00165, label = "Late Dilmun", color = "navyblue", 
           angle = 0, hjust = 0.5, size = 3, fontface = "bold", family = "sans")+
  annotate("segment", x = 3550, xend = 3550, y = 0.001685, yend = 0.001655,  color = "navyblue",linetype="solid",linewidth = .5) +
  #Middle Dilmun start
  annotate("text", x=3375,y=0.00165, label = "Middle Dilmun", color = "navyblue", 
           angle = 0, hjust = 0.5, size = 3, fontface = "bold", family = "sans")+
  #end of captioning
  scale_y_continuous(limits = c(0,0.0018), breaks=seq(0,0.0008, .0004))+
  scale_x_reverse(limits = c(5200,3200), breaks=seq(5200,3200, -200), sec.axis = sec_axis(~ . *-1+1950, name = "Years (calBC)", breaks = seq(-3300,-1300, 200))) +
  labs(x='Years (calBP)', y= 'Probability')+
  guides(color = guide_legend(title = NULL), fill = guide_legend(title = NULL))+
  theme(plot.margin = unit(c(0,0,0,0), "cm"))+
  theme_bw()+
  theme(plot.margin = unit(c(0.1,0.25,0.1,0.25), "cm"),legend.key.size = unit(.35, 'cm'), legend.position = c(.19, .35),legend.justification = c("right", "top"),legend.box.just = "left",legend.margin = margin(1, 1, 1, 1),legend.key.width = unit(.6,"cm"), axis.text.y = element_text(angle = 90,hjust=0.5), axis.title.y = element_text(hjust=0.25))

ggsave("Figure_4.pdf",  plot = last_plot(), device = "pdf",
       path = "Master Paper/paper pdf/", width = 4350, height = 2150,
       units = 'px', dpi = 300)





##### Figure S1 -- KDE curves####
ggplot(joined, aes(x = calBP)) +
  geom_line(aes(y = cKDE.mean, colour = "Composite KDE"), size=1)+
  geom_line(aes(y = probability, colour = "Oxcal KDE"), size=1)+
  scale_colour_manual(values=c("Composite KDE"="#723BF5","Oxcal KDE"="#F4542C"))+
  #this is to show the pearsons correlation between the two KDE curves
  annotate("text", x=5200,y=0.000525, label = paste0("rho = ",round(KDE.cor$estimate, 3)), color = "grey25", 
           angle = 0, hjust = 0, vjust = 0.5, size = 3.25, fontface = "bold", family = "sans")+
  annotate("text",  x=5200,y=0.00045, label = "p-value = <0.05", color = "grey25", 
           angle = 0, hjust = 0, vjust = 0.5, size = 3.25, fontface = "bold", family = "sans")+
  scale_y_continuous(limits = c(0,0.0008), breaks=seq(0,0.0008, .0004))+
  scale_x_reverse(limits = c(5200,3200), breaks=seq(5200,3200, -200), sec.axis = sec_axis(~ . *-1+1950, name = "Years (calBC)", breaks = seq(-3300,-1300, 200))) +
  labs(x='Years (calBP)', y= 'Probability')+
  guides(color = guide_legend(title = NULL), fill = guide_legend(title = NULL))+
  theme(plot.margin = unit(c(0,0,0,0), "cm"))+
  theme_bw()+
  theme(plot.margin = unit(c(0.1,0.25,0.1,0.25), "cm"), legend.position = c(.1225, .925),legend.justification = c("right", "top"),legend.box.just = "left",legend.margin = margin(1, 1, 1, 1),legend.key.width = unit(.65,"cm"), axis.text.y = element_text(angle = 90,hjust=0.5))

#save figure S1 plot to file.
ggsave("Figure_S1.pdf",  plot = last_plot(), device = "pdf",
       path = "Master Paper/paper pdf/", width = 3508, height = 875,
       units = 'px', dpi = 300)



####a bit of wrangling to make sure the ribbon behaves properly for red and green colouration
Post.Pred.Check_joined = Post.Pred_joined %>% 
  mutate(above = ifelse((re_Norm_SPD>=`90%`), re_Norm_SPD, NA_real_), below = ifelse((re_Norm_SPD<=`10%`), re_Norm_SPD, NA_real_))
##### Figure S2a  -- SPD Model####
Figure_S2a = ggplot(Post.Pred.Check_joined, aes(x = calBP)) +
  geom_ribbon(aes(ymin=`10%`, ymax=`90%`, fill = "Double Exponential Model Envelope"), alpha = 0.5)+
  geom_ribbon(aes(ymin=`90%`, ymax=above, fill = ">90%"), alpha = 0.5)+
  geom_ribbon(aes(ymin=`10%`, ymax=below, fill = "<10%"), alpha = 0.5)+
  geom_line(aes(y = re_Norm_SPD, colour = "SPD"), size=1, linetype = "solid")+
  scale_colour_manual(values=c("SPD"="black"))+
  scale_fill_manual(values=c("Double Exponential Model Envelope"="grey75", ">90%" = "blue", "<10%" = "red"))+
  #this is to show the pearsons correlation between spd and model
  annotate("text", x=5200,y=0.000625, label = paste0("rho = ~",cor.res.min,"-",cor.res.max), color = "grey25", 
           angle = 0, hjust = 0, vjust = 0.5, size = 3.25, fontface = "bold", family = "sans")+
  annotate("text", x=5200,y=0.000555, label = "p-value = <0.05", color = "grey25", 
           angle = 0, hjust = 0, vjust = 0.5, size = 3.25, fontface = "bold", family = "sans")+
  scale_y_continuous(limits = c(0,0.0012), breaks=seq(0,0.0012, .0004))+
  scale_x_reverse(limits = c(5200,3200), breaks=seq(5200,3200, -200), sec.axis = sec_axis(~ . *-1+1950, name = "Years (calBC)", breaks = seq(-3300,-1300, 200))) +
  labs(x='Years (calBP)', y= 'Summed Probability')+
  guides(fill = guide_legend(title = NULL),color = guide_legend(title = NULL))+
  theme_bw()+
  theme(plot.margin = unit(c(0.1,0.25,0.1,0.25), "cm"), legend.key.size = unit(.225, 'cm'),legend.position = c(.21, .565),legend.justification = c("right", "bottom"),legend.box.just = "left",legend.margin = margin(1, 1, 1, 1),legend.key.width = unit(.4,"cm"), axis.text.y = element_text(angle = 90,hjust=0.5))
####a bit of wrangling to make sure the ribbon behaves properly for red and green colouration
Log_Post.Pred.Check = read_rds("Master Paper/rds/Logistic_post_pred_check.rds")
###and the MCMC posteriors
Log_MCMC_posteriors = read_rds("Master Paper/rds/MCMC_Logistic_posteriors.rds")
#
#join the two dataframes
Log_Post.Pred_joined = full_join(Log_Post.Pred.Check, SPDs_stu)
#subset only those being used in correlation
Log_Post.Pred_cor = Log_Post.Pred_joined[c(1:5,8)]
###This provides the correlation results for all proxies compared with 
log.cor.res = as.numeric(cor(Log_Post.Pred_cor, method = "pearson")[1,-1])
#calculate the min and max values to give the range of the correlation matrix.
log_cor.res.min = round(min(log.cor.res), 3)
log_cor.res.max = round(max(log.cor.res), 3)
#
Log_Post.Pred.Check_joined = Log_Post.Pred_joined %>% 
  mutate(above = ifelse((re_Norm_SPD>=`90%`), re_Norm_SPD, NA_real_), below = ifelse((re_Norm_SPD<=`10%`), re_Norm_SPD, NA_real_))
##### Figure s2b -- SPD Model####
Figure_S2b = ggplot(Log_Post.Pred.Check_joined, aes(x = calBP)) +
  geom_ribbon(aes(ymin=`10%`, ymax=`90%`, fill = "Double Exponential Model Envelope"), alpha = 0.5)+
  geom_ribbon(aes(ymin=`90%`, ymax=above, fill = ">90%"), alpha = 0.5)+
  geom_ribbon(aes(ymin=`10%`, ymax=below, fill = "<10%"), alpha = 0.5)+
  geom_line(aes(y = re_Norm_SPD, colour = "SPD"), size=1, linetype = "solid")+
  scale_colour_manual(values=c("SPD"="black"))+
  scale_fill_manual(values=c("Double Exponential Model Envelope"="grey75", ">90%" = "blue", "<10%" = "red"))+
  #this is to show the pearsons correlation between spd and model
  annotate("text", x=5200,y=0.000925, label = paste0("rho = ~",log_cor.res.min,"-",log_cor.res.max), color = "grey25", 
           angle = 0, hjust = 0, vjust = 0.5, size = 3.25, fontface = "bold", family = "sans")+
  annotate("text", x=5200,y=0.000825, label = "p-value = <0.05", color = "grey25", 
           angle = 0, hjust = 0, vjust = 0.5, size = 3.25, fontface = "bold", family = "sans")+
  scale_y_continuous(limits = c(0,0.0018), breaks=seq(0,0.0018, .0004))+
  scale_x_reverse(limits = c(5200,3200), breaks=seq(5200,3200, -200), sec.axis = sec_axis(~ . *-1+1950, name = "Years (calBC)", breaks = seq(-3300,-1300, 200))) +
  labs(x='Years (calBP)', y= 'Summed Probability')+
  guides(fill = guide_legend(title = NULL),color = guide_legend(title = NULL))+
  theme_bw()+
  theme(plot.margin = unit(c(0.1,0.25,0.1,0.25), "cm"), legend.key.size = unit(.225, 'cm'),legend.position = c(.42, .565),legend.justification = c("right", "bottom"),legend.box.just = "left",legend.margin = margin(1, 1, 1, 1),legend.key.width = unit(.4,"cm"), axis.text.y = element_text(angle = 90,hjust=0.5))

####connect figure 2a & 2b together into one.
Figure_S2 = ggarrange(Figure_S2a,Figure_S2b,
                       labels = c("A", "B"),
                       ncol = 1, nrow = 2)
#plot figure 2a and 2b
Figure_S2

#save figure S2 plot to file.
ggsave("Figure_S2.pdf",  plot = last_plot(), device = "pdf",
       path = "Master Paper/paper pdf/", width = 3508, height = 1750,
       units = 'px', dpi = 300)



#wrangle for S3
#calculating the mean growth rates 1 and 2 from the MCMC posterior distributions for comparison with GMM Growth rates
r_mean = mean(c(Log_MCMC_posteriors$r)*100)
k_mean = mean(c(Log_MCMC_posteriors$k)*100)

#Extract the 90% confidence envelope for r1, with median.
mcmc_r_conf = quantile(Log_MCMC_posteriors$r, probs=c(.1,.5,.9))
#Extract the 90% confidence envelope for r2, with median.
mcmc_k_conf = quantile(Log_MCMC_posteriors$k, probs=c(.1,.5,.9))
##### Figure S3a -- MCMC posterior r####
Figure_S3a = ggplot(Log_MCMC_posteriors)+
  geom_density(aes(r),colour = "grey25", fill= "lightblue")+
  labs(x='Growth Rate', y= 'Frequency')+
  scale_y_continuous(limits = c(0,3000), breaks=seq(0,3000, 500), guide = guide_axis(angle = 90))+
  scale_x_continuous(limits = c(-0.0002,0.001), breaks=seq(-0.0002,0.0014, 0.0004))+
  annotate("text", x=0.0004,y=2950, label = paste0("90%HDPI: ",round(mcmc_r_conf[1],5),"~",round(mcmc_r_conf[3],5)), color = "grey25", 
           angle = 0, hjust = 0.5, vjust = 0.5, size = 3.25, fontface = "bold", family = "sans")+
  theme_bw()+
  theme(plot.margin = unit(c(0.1,0.25,0.25,0.25), "cm"))

##### Figure S3b -- MCMC posterior k####
Figure_S3b = ggplot(Log_MCMC_posteriors)+
  geom_density(aes(as.numeric(k)),colour = "grey25", fill= "lightblue")+
  labs(x='Carrying Capacity (k)', y= 'Frequency')+
  scale_y_continuous(limits = c(0,15), breaks=seq(0,15, 4), guide = guide_axis(angle = 90))+
  scale_x_continuous(limits = c(0,0.2), breaks=seq(0,0.2, 0.06))+
  annotate("text", x=0.1,y=14.75, label = paste0("90%HDPI: ",round(mcmc_k_conf[1],5),"~",round(mcmc_k_conf[3],5)), color = "grey25", 
           angle = 0, hjust = 0.5, vjust = 0.5, size = 3.25, fontface = "bold", family = "sans")+
  theme_bw()+
  theme(plot.margin = unit(c(0.1,0.25,0.25,0.25), "cm"))


####connect figure 3a & 3b together into one.
Figure_S3 = ggarrange(Figure_S3a,Figure_S3b,
                      labels = c("A", "B"),
                      ncol = 2, nrow = 1)
#plot figure 3a and 3b
Figure_S3

#save figure S3 plot to file.
ggsave("Figure_S3.pdf",  plot = last_plot(), device = "pdf",
       path = "Master Paper/paper pdf/", width = 3508, height = 1750,
       units = 'px', dpi = 300)


###plot the oxcal results for supplementary materials
#extract just the posterior values from the Oxcal Bayesian KDE ready for plotting.
OXCAL_plot1= OXCAL %>% 
  filter(type == "posterior")
#plot the posterior values.
ggplot(OXCAL_plot1)+
  geom_line(aes(x = value, y = probability))+
  facet_wrap(~op, scales = "free")

#extract just the likelihood values from the Oxcal Bayesian KDE ready for plotting.
OXCAL_plot2= OXCAL %>% 
  filter(type == "likelihood")
#plot the likelihood values.
ggplot(OXCAL_plot2)+
  geom_line(aes(x = value, y = probability))+
  facet_wrap(~op, scales = "free")

