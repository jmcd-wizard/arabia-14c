library(nimble)
library(nimbleCarbon)
library(rcarbon)
library(coda)
library(tidyverse)


#load c14 dataset
arabia_data = read.csv("csv/arabia_dates_20.09.24.csv")%>% 
  filter(Country=="AE"|Country=="OM")

#calibrate the dates
caldates=calibrate(x=arabia_data$CRA,errors=arabia_data$Error,calCurves='intcal20',verbose=FALSE) #calibration

#subset the caldates to the start and end date of our study period. p= means the cumulative prob density of the date is above 0.5.
caldates = subset(caldates, BP <=5200 & BP >=3200, p=0.5)


#create an SPD from the observed c14 dates
obs.spd = spd(caldates,timeRange=c(5200,3200),verbose=FALSE)
#test plot it
plot(obs.spd)

#extract the observed BP value from the above subset
obs.CRA = caldates$metadata$CRA
#and SD/Error
obs.Errors = caldates$metadata$Error

#now define the start and end for the growth models:
#start:
a = 5200
#end
b = 3200

#extract the input data as lists
constants <- list(N=length(obs.CRA),calBP=intcal20$CalBP,C14BP=intcal20$C14Age,C14err=intcal20$C14Age.sigma,start=a,end=b)
data <- list(X=obs.CRA,sigma=obs.Errors)

#define the length of the list of CRA objects (ie. the number of dates)
N = length(obs.CRA)

###Logistic Growth Model: mL.a (crema)
mL.a <- nimbleCode({
  for (i in 1:N){
    # Growth Model Likelihood
    theta[i] ~ dLogisticGrowth(a=start,b=end,k=k,r=r);
    # Calibration
    mu[i] <- interpLin(z=theta[i], x=calBP[], y=C14BP[]);
    sigmaCurve[i] <- interpLin(z=theta[i], x=calBP[], y=C14err[]);
    sd[i] <- (sigma[i]^2+sigmaCurve[i]^2)^(1/2);
    X[i] ~ dnorm(mean=mu[i],sd=sd[i]);
  }
  # Prior
  r ~ dexp(1/0.004); # Prior
  k ~ T(dexp(1/0.05),0.001,0.2)
})  



#initialisation functions for the two models so that chains have different starting parameter values for the growth models. 
m.dates = medCal(caldates)
if(any(m.dates>3200|m.dates<3200)){m.dates[m.dates>5200]=5200;m.dates[m.dates<3200]=3200}
inits.function.mL.a = function() list(r=rexp(1,1/0.0004),k=runif(1,0.0001,0.2),theta=m.dates)

#define & run mL.a/LOG-warning! this takes hours!
mcmc.samples.mL.a<- nimbleMCMC(code = mL.a,constants = constants,data = data,niter = 10000, nchains = 2, thin=1, nburnin = 3000, progressBar = FALSE, monitors=c('r','k','theta'), inits=inits.function.mL.a, samplesAsCodaMCMC=TRUE,setSeed=c(123,456),WAIC=TRUE)


par(mfrow=c(2,2))
plot(as.numeric(mcmc.samples.mL.a$samples$chain1[,'r']),type='l',xlab='MCMC Iteration',ylab='r',main='mL.a r chain 1')
plot(as.numeric(mcmc.samples.mL.a$samples$chain2[,'r']),type='l',xlab='MCMC Iteration',ylab='r',main='mL.a r chain 2')
plot(as.numeric(mcmc.samples.mL.a$samples$chain1[,'k']),type='l',xlab='MCMC Iteration',ylab='k',main='mL.a k chain 1')
plot(as.numeric(mcmc.samples.mL.a$samples$chain2[,'k']),type='l',xlab='MCMC Iteration',ylab='k',main='mL.a k chain 2')


mL.a.rhat=gelman.diag(mcmc.samples.mL.a$samples)
mL.a.ess=effectiveSize(mcmc.samples.mL.a$samples)
head(mL.a.rhat$psrf)

mL.a.ess[1:2]


par(mfrow=c(1,2))
postHPDplot(mcmc.samples.mL.a$samples$chain2[,'r'],rnd=5,xlab='r',ylab='Density',prob = 0.95,main='Logistic Model: r')
postHPDplot(mcmc.samples.mL.a$samples$chain2[,'k'],rnd=5,xlab='k',ylab='Density',prob = 0.95,main='Logistic Model: k')

#extract the MCMC posteriors - I have arbiratilry chosen chain 2 as there is little need to plot more than one chain...
MCMC_log_posteriors = as.data.frame(mcmc.samples.mL.a$samples$chain2[,'r']) %>% 
  mutate(r = var1,
         var1=NULL,
         k = mcmc.samples.mL.a$samples$chain2[,'k'])
#now write the exported MCMC posteriors to file!
write_rds(MCMC_log_posteriors,"rds/MCMC_Logistic_posteriors.rds")


params.mL.a = list(r=c(mcmc.samples.mL.a$samples$chain1[,'r'],mcmc.samples.mL.a$samples$chain2[,'r']),k=c(mcmc.samples.mL.a$samples$chain1[,'k'],mcmc.samples.mL.a$samples$chain2[,'k']))

par(mfrow=c(1,2))
set.seed(123)
modelPlot(dLogisticGrowth,a=5200,b=3200,params=params.mL.a,nsample=100,alpha = 0.4,ylim=c(0,0.001),main='mL.a: Logistic',type='envelope')



#plot mL.a to see how it fits visually with the observed spds
set.seed(123)
pp.check.mL.a=postPredSPD(obs.CRA,errors = obs.Errors,calCurve = 'intcal20',model = dLogisticGrowth,a = 5200,b=3200,params=list(r=mcmc.samples.mL.a$samples$chain2[,'r'],k=mcmc.samples.mL.a$samples$chain2[,'k']),method='calsample',nsim = 100,ncores = 2,verbose=FALSE)
par(mfrow=c(1,1))
plot(pp.check.mL.a, main = "mL.a")

postHPDplot(postPredCor(pp.check.mL.a),xlab="Pearson's Correlation coefficient",ylab='Density',main='Logistic Model goodness-of-fit',xlim=c(0,1))

#export the exact SPD curve created for comparison with other proxies
Post.Pred.Check.SPD = as.data.frame(pp.check.mL.a$obs)
#export the simulated values and put them into quantiles.
Post.Pred.Check = apply(pp.check.mL.a$simmatrix, 1, quantile, c(0,.05,.5,.95,1))
#turn into a tidy-friendly dataframe for export to the main arabia analysis
Post.Pred.Check = as.data.frame(t(Post.Pred.Check))
#now add dates for joining below
Post.Pred.Check = Post.Pred.Check %>% 
  mutate(calBP = 5200:3200)
#join the two dataframes
Post.Pred.Check.Log = full_join(Post.Pred.Check, Post.Pred.Check.SPD)


#now write the post pred results to file!
write_rds(Post.Pred.Check.Log,"rds/Logistic_post_pred_check.rds")



plot(obs.spd,spdnormalised = TRUE)
highest.cor.index = which.max(postPredCor(pp.check.mL.a))
lines(5200:3200,pp.check.mL.a$simmatrix[,highest.cor.index],lty=2)
legend('topleft',legend=c('observed SPD','Posterior Predictive SPD with the highest correlation'),col=c('lightgrey','black'),lwd=c(4,1),lty=c(1,2),bty='n')


#now compare the two models
compare.models(mcmc.mDE.samples,mcmc.samples.mL.a)
