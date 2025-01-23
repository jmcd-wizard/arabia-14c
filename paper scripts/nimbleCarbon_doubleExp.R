####Load libraries####
library(rcarbon)
library(nimbleCarbon)
library(truncnorm)
library(coda)
library(tidyverse)

#brute force stops scientific numbers in favour of full digits.
options(scipen=999)

####load data####
#load c14 dataset
arabia_data = read.csv("csv/arabia_dates_20.09.24.csv")%>% 
  filter(Country=="AE"|Country=="OM")
#prepare the 'bins' for thindate
#calibrate the dates
obs.caldates = calibrate(x=arabia_data$CRA,errors=arabia_data$Error,calCurves='intcal20',verbose=FALSE) #calibration
# Generate bins grouping dates within 100 yrs
#cal.bins = binPrep(sites=arabia_data$SiteID,ages=arabia_data$CRA,h=50) 
# Sample 1 date from each bin, selecting randomly the sample with the smallest error
#obs.caldates = obs.caldates[thinDates(ages=arabia_data$CRA,  errors=arabia_data$Error, bins=cal.bins, size=1, thresh=1,seed=123,method='splitsample')]

obs.spd = spd(obs.caldates,timeRange=c(5200,3200),verbose=FALSE)
plot(obs.spd)

#create an SPD from the observed c14 dates
#obs.spd.all = spd(obs.caldates,timeRange=c(8000,0),verbose=FALSE)
#obs.spd = spd(obs.caldates,timeRange=c(5150,3250),verbose=FALSE)
#par(mfrow=c(1,1))
#plot(obs.spd.all)
#par(mfrow=c(1,1))
#plot(obs.spd)

# Consider samples that have a probability mass over 0.5 within 5150 and 3250, This is the Bronze age period in southeastern Oman(n=118)
index = which.CalDates(obs.caldates,BP <=5200 & BP >=3200,p=0.5)

#subset these
obs.caldates = obs.caldates[index]

# Extract relevant CRA and CRAErrors
CRA=arabia_data$CRA[index]
Errors=arabia_data$Error[index]
SiteID = arabia_data$SiteID[index]
LabCode = arabia_data$LabID[index]
# Extract Median Calibrated Dates
medDates = medCal(obs.caldates)
# Save to R Image
obs.data = data.frame(LabCode=LabCode,CRA=CRA,Error=Errors,MedCalDate=medDates,SiteID=SiteID)

# Define Model ####
p = list(r1=0.0006,r2=-0.001,mu=3800)
modelPlot(model = dDoubleExponentialGrowth,a=5200,b=3200,params=p,alpha = 1)

modelPlot(model=dDoubleExponentialGrowth,a=5200,b=3200,params=list(r1=rnorm(500, sd=0.0004),r2=rnorm(500, sd=0.0004),mu=round(runif(500,3200,5200))), alpha = 1) 

#using as sd rate the average growth rate from Zahid et al 2015
mDE <- nimbleCode({
  for (i in 1:N){
    theta[i] ~ dDoubleExponentialGrowth(a=5200,b=3200,r1=r1,r2=r2,mu=changept);
    mu[i] <- interpLin(z=theta[i], x=calBP[], y=C14BP[]);
    sigmaCurve[i] <- interpLin(z=theta[i], x=calBP[], y=C14err[]);
    sd[i] <- (sigma[i]^2+sigmaCurve[i]^2)^(1/2);
    X[i] ~ dnorm(mean=mu[i],sd=sd[i]);
  }
  r1 ~ dnorm(0,sd=0.0004); 
  r2 ~ dnorm(0,sd=0.0004);
  chp ~ T(dnorm(3700,sd=200),3200,5200);
  changept <- round(chp);
})  


# Define Constants and Data ####
data(intcal20)
constants <- list(N=length(obs.caldates),calBP=intcal20$CalBP,C14BP=intcal20$C14Age,C14err=intcal20$C14Age.sigma)
data <- list(X=obs.data$CRA,sigma=obs.data$Error)

# Define Initialisation Function
inits =  list(r1=rnorm(1,sd=0.0004),r2 = rnorm(1,sd=0.0004), theta=as.numeric(obs.data$MedCalDate), chp=round(rtruncnorm(1,mean=4200,a=3200,b=5200)))

# Run MCMC 
mcmc.mDE.samples<- nimbleMCMC(code = mDE,constants = constants,data = data,niter = 10000, nchains = 3, thin=1, nburnin = 3000, summary = FALSE, monitors=c('r1','r2','chp','theta'),WAIC=TRUE,samplesAsCodaMCMC=TRUE,inits=inits,setSeed=c(1,2,3))

# Quick Summaries
gelman.diag(mcmc.mDE.samples$samples)$psrf[1:3,]


par(mfrow=c(3,3))
#chain 1
plot(as.numeric(mcmc.mDE.samples$samples$chain1[,'r1']),type='l',xlab='MCMC Iteration',ylab='r1',main='mDE r1 chain 1')
plot(as.numeric(mcmc.mDE.samples$samples$chain1[,'r2']),type='l',xlab='MCMC Iteration',ylab='r2',main='mDE r2 chain 1')
plot(as.numeric(mcmc.mDE.samples$samples$chain1[,'chp']),type='l',xlab='MCMC Iteration',ylab='chp',main='mDE chp chain 1')
#chain 2
plot(as.numeric(mcmc.mDE.samples$samples$chain2[,'r1']),type='l',xlab='MCMC Iteration',ylab='r1',main='mDE r1 chain 2')
plot(as.numeric(mcmc.mDE.samples$samples$chain2[,'r2']),type='l',xlab='MCMC Iteration',ylab='r2',main='mDE r2 chain 2')
plot(as.numeric(mcmc.mDE.samples$samples$chain2[,'chp']),type='l',xlab='MCMC Iteration',ylab='chp',main='mDE chp chain 2')
#chain 3
plot(as.numeric(mcmc.mDE.samples$samples$chain3[,'r1']),type='l',xlab='MCMC Iteration',ylab='r1',main='mDE r1 chain 3')
plot(as.numeric(mcmc.mDE.samples$samples$chain3[,'r2']),type='l',xlab='MCMC Iteration',ylab='r2',main='mDE r2 chain 3')
plot(as.numeric(mcmc.mDE.samples$samples$chain3[,'chp']),type='l',xlab='MCMC Iteration',ylab='chp',main='mDE chp chain 3')


rhat.mDE=gelman.diag(mcmc.mDE.samples$samples)
head(rhat.mDE$psrf)

ess.mDE=effectiveSize(mcmc.mDE.samples$samples)
ess.mDE[1:2]


#extract parameters for post pred check below
params.mDE = list(r1 = c(mcmc.mDE.samples$samples$chain1[,'r1'],mcmc.mDE.samples$samples$chain2[,'r1'],mcmc.mDE.samples$samples$chain3[,'r1']),
                  r2 = c(mcmc.mDE.samples$samples$chain1[,'r2'],mcmc.mDE.samples$samples$chain2[,'r2'],mcmc.mDE.samples$samples$chain3[,'r2']),
                  mu = round(c(mcmc.mDE.samples$samples$chain1[,'chp'],mcmc.mDE.samples$samples$chain2[,'chp'],mcmc.mDE.samples$samples$chain3[,'chp'])))

#this calculates the posterior predictive check
pp.check.mDE.calsample=postPredSPD(obs.data$CRA,obs.data$Error,calCurve = 'intcal20',model = dDoubleExponentialGrowth,a = 5200,b=3200,params=params.mDE,nsim = 500,ncores = 5,verbose=FALSE,method='calsample', spdnormalised = TRUE, datenormalised = TRUE)

#export the exact SPD curve created for comparison with other proxies
Post.Pred.Check.SPD = as.data.frame(pp.check.mDE.calsample$obs)
#export the simulated values and put them into quantiles.
Post.Pred.Check = apply(pp.check.mDE.calsample$simmatrix, 1, quantile, c(0,.05,.5,.95,1))
#turn into a tidy-friendly dataframe for export to the main arabia analysis
Post.Pred.Check = as.data.frame(t(Post.Pred.Check))
#now add dates for joining below
Post.Pred.Check = Post.Pred.Check %>% 
  mutate(calBP = 5200:3200)
#join the two dataframes
Post.Pred.Check = full_join(Post.Pred.Check, Post.Pred.Check.SPD)

#now write the post pred results to file!
write_rds(Post.Pred.Check,"rds/Double_Exp_post.pred.check.rds")


#this calculates the posterior predictive check but for an unnormalised date
pp.check.mDE.calsample.unnorm=postPredSPD(obs.data$CRA,obs.data$Error,calCurve = 'intcal20',model = dDoubleExponentialGrowth,a = 5200,b=3200,params=params.mDE,nsim = 500,ncores = 5,verbose=FALSE,method='calsample', spdnormalised = TRUE, datenormalised = FALSE)

#export the exact SPD curve created for comparison with other proxies
Post.Pred.Check.unnorm.SPD = as.data.frame(pp.check.mDE.calsample.unnorm$obs)
#export the simulated values and put them into quantiles.
Post.Pred.Check.unnorm = apply(pp.check.mDE.calsample.unnorm$simmatrix, 1, quantile, c(0,.05,.5,.95,1))
#turn into a tidy-friendly dataframe for export to the main arabia analysis
Post.Pred.Check.unnorm = as.data.frame(t(Post.Pred.Check.unnorm))
#now add dates for joining below
Post.Pred.Check.unnorm = Post.Pred.Check.unnorm %>% 
  mutate(calBP = 5200:3200)
#join the two dataframes
Post.Pred.Check.unnorm = full_join(Post.Pred.Check.unnorm, Post.Pred.Check.unnorm.SPD)

#now write the post pred results to file!
write_rds(Post.Pred.Check.unnorm,"rds/Double_Exp_post.pred.check.unnorm.rds")





##this plots the collection between posterior predictive check samples (using the
##rcarbon 'calsample' date random sampling) and SPDs
par(mfrow=c(1,1))
postHPDplot(postPredCor(pp.check.mDE.calsample),xlab="Pearson's Correlation coefficient",ylab='Density',main='mDE calsample goodness-of-fit')
#this plots the posterior predictive check samples and SPDS
par(mar = c(4.1,4,1,1))
plot(pp.check.mDE.calsample, xlab = "Years (calBP)", ylab = "Summed Probability", interval = .9)
grid(lty = 1, lwd = .5)


#Marginal posterior distributions#
par(mfrow=c(1,3))
postHPDplot(mcmc.mDE.samples$samples$chain2[,'r1'],rnd=5,xlab='Growth Rate 1',ylab='Frequency',prob = 0.95) #Double Exponential Model Growth Rate 1'
title(main='A', adj = 0)
postHPDplot(mcmc.mDE.samples$samples$chain2[,'r2'],rnd=5,xlab='Growth Rate 2',ylab='Frequency',prob = 0.95) #Double Exponential Model Growth Rate 2'
title(main='B', adj = 0)
postHPDplot(mcmc.mDE.samples$samples$chain2[,'chp'],rnd=5,xlab='Changepoint',ylab='Density',prob = 0.95)#Double Exponential Model Changepoint
title(main='C', adj = 0)

#extract the MCMC posteriors - I have arbiratilry chosen chain 2 as there is little need to plot more than one chain...
MCMC_posteriors = as.data.frame(mcmc.mDE.samples$samples$chain2[,'r1']) %>% 
  mutate(r1 = var1,
         var1=NULL,
         r2 = mcmc.mDE.samples$samples$chain2[,'r2'],
         chp = mcmc.mDE.samples$samples$chain2[,'chp'])
#now write the exported MCMC posteriors to file!
write_rds(MCMC_posteriors,"rds/MCMC_Double_Exp_posteriors.rds")


#show the prior in comparison to the posterior distribution
par(mfrow=c(1,1))
postHPDplot(mcmc.mDE.samples$samples$chain1[,'r1'],rnd=5,xlab='r1',ylab='Density',prob = 0.95)
abline(v=0.0006,lty=2)
axis(3,at=0.0006,label='True value of r1')


#extract mean (and also median for goodness-of-fits) growth rates and changepoints from each chain
estimated.r1.c1=mean(mcmc.mDE.samples$samples$chain1[,'r1'])
estimated.r1.c2=mean(mcmc.mDE.samples$samples$chain2[,'r1'])
estimated.r1.c3=mean(mcmc.mDE.samples$samples$chain3[,'r1'])

estimated.r2.c1=mean(mcmc.mDE.samples$samples$chain1[,'r2'])
estimated.r2.c2=mean(mcmc.mDE.samples$samples$chain2[,'r2'])
estimated.r2.c3=mean(mcmc.mDE.samples$samples$chain3[,'r2'])

estimated.mu.c1 = mean(mcmc.mDE.samples$samples$chain1[,'chp'])
estimated.mu.c2 = mean(mcmc.mDE.samples$samples$chain2[,'chp'])
estimated.mu.c3 = mean(mcmc.mDE.samples$samples$chain3[,'chp'])


estimated.gof.cal.median = median(postPredCor(pp.check.mDE.calsample))
estimated.gof.cal.mean = mean(postPredCor(pp.check.mDE.calsample))
