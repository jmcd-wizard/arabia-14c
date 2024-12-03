rm(list = ls())

#library(devtools)
#install_github("eehh-stanford/baydem")
library(baydem)
library(magrittr)
library(rstan)
library(tidyverse)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

#Load the radiocarbon dates
arabia_dates <- read.csv("Master Paper/csv/arabia_dates_20.09.24.csv",
                         stringsAsFactors=FALSE, check.names=FALSE)
#wrangle the period subset with 1k years edge to avoid negative edge effects.
arabia = as.data.frame(arabia_dates)%>% 
  filter(CRA <=5650 &CRA >= 2750) %>%
  #mutate(CRA = (CRA*-1)+1950)%>% 
  write_csv("Master Paper/csv/arabia_dates_20.09.24_filtered.csv")

# Add columns to arabia with the uncalibrated years BP and corresponding error
# using the baydem naming conventions.
arabia$trc_m <- arabia[,"CRA"]
arabia$sig_trc_m <- arabia[,"Error"]

# Calculate the fraction modern and associated uncertainty for all data
arabia$phi_m <- exp(-arabia$trc_m/8033)
arabia$sig_m <- arabia$sig_trc_m * arabia$phi_m / 8033

data_dir <- "Master Paper/GMM/outputs/K_final"

hp <-
  list(
    # Parameter for the dirichlet draw of the mixture probabilities
    alpha_d = 1,
    # The gamma distribution shape parameter for sigma
    alpha_s = 10,
    # The gamma distribution rate parameter for sigma, yielding a mode of 100
    alpha_r = (10 - 1) / 100,
    # Spacing for the measurement matrix (years)
    dtau = 1
  )
saveRDS(hp, file.path("Master Paper/GMM/outputs/K_final","arabia_hp.rds"))

# Do Bayesian inference with 2 to 14 mixtures:
Kvect <- 2:14
num_models <- length(Kvect)
# Define random number seeds for (1) initializing the Bayesian inference and
# (2) the sampling in Stan. The following seeds were created via the following
# set of commands:

base_seed <- 826775 # from random.org, between 1 and 1,000,000 (test = 72826775)
set.seed(base_seed)
seed_mat <- matrix(sample.int(1000000,2*num_models),ncol=2)
init_seed_vect <- seed_mat[,1]
stan_seed_vect <- seed_mat[,2]

# To reduce memory usage, use a for loop here to do Bayesian inference rather
# than using the standard pipeline functions. To allow interupted runs to be
# continued, only do the inference if a save file does not exist.

# Use the intcal20 calibration curve
calib_curve <- "intcal20"
calib_df <- load_calib_curve(calib_curve)


# Use only the observations from the dataset which was wrangled above
rc_meas <- arabia

# Define the base density model (without K)
density_model0 <- list(type="trunc_gauss_mix",
                       tau_min=-3700,
                       tau_max=-800)

# Save a list containing the variables used for the inference in the following
# for loop
arabia_inference_inputs <- list(Kvect=Kvect,
                               density_model0=density_model0,
                               rc_meas=rc_meas,
                               hp=hp,
                               calib_df=calib_df,
                               init_seed_vect=init_seed_vect,
                               stan_seed_vect=stan_seed_vect)
saveRDS(arabia_inference_inputs, file.path("Master Paper/GMM/outputs/K_final","arabia_inference_inputs.rds"))
for (m_K in 1:length(Kvect)) {
  K <- Kvect[m_K]
  save_file <- file.path("Master Paper/GMM/outputs/K_final", paste0("arabia_K",K,".rds"))
  if(!file.exists(save_file)) {
    density_model <- density_model0
    density_model$K <- K
    t0 <- Sys.time() # start time
    bayesian_solution <- sample_theta(rc_meas,
                                      density_model,
                                      hp,
                                      calib_df,
                                      init_seed=init_seed_vect[m_K],
                                      stan_seed=stan_seed_vect[m_K])
    t1 <- Sys.time() # end time
    run_time_sec <- as.numeric(difftime(t1,t0,units="secs"))
    # Calculate PSIS-LOO CV
    log_lik_mat <- rstan::extract(bayesian_solution$fit,"logh")[[1]]
    loo_analysis <- loo::loo(log_lik_mat)
    loo_value <- loo_analysis["estimates"][[1]]["elpd_loo","Estimate"]
    # Save the results to file
    saveRDS(list(bayesian_solution=bayesian_solution,
                 run_time_sec=run_time_sec,
                 loo_value=loo_value),
            save_file)
  }
}


# Create a vector of loo values
loo_vect <- rep(NA,length(Kvect))
for (m_K in 1:length(Kvect)) {
  K <- Kvect[m_K]
  save_file <- file.path("Master Paper/GMM/outputs/K_final", paste0("arabia_K",K,".rds"))
  if(!file.exists(save_file)) {
    stop(paste0("Missing save file for K=",K))
  }
  bayesian_soln <- readRDS(save_file)
  loo_vect[m_K] <- bayesian_soln$loo
}

# Identify the best model and write a summary of the loo as a .yaml file
m_K_best <- which.max(loo_vect)
K_best <- Kvect[m_K_best]
loo_summary <- list(m_K_best=m_K_best,K_best=K_best,loo_vect=loo_vect)
yaml::write_yaml(loo_summary, file.path("Master Paper/GMM/outputs/K_final","loo_summary.yaml"))

# Do the Baysian summary for the best model
density_model <- density_model0
density_model$K <- K_best
save_file <- file.path("Master Paper/GMM/outputs/K_final", paste0("arabia_K",K_best,".rds"))
#save_file <- file.path("End-to-end/outputs/K_test", paste0("arabia_K","X",".rds")) #for this, KX is the value with the highest loo value. this is if you want to exactly specify the Kx
bayesian_soln <- readRDS(save_file)
bayesian_soln <- bayesian_soln$bayesian_solution
bayesian_summ <- summarize_bayesian_inference(bayesian_soln,rc_meas,density_model,calib_df,hp$dtau)

##extract values from the bayesian summary completed on the best fitting Gaussian Mixture
#density values
extracted_Qdens = as.data.frame(bayesian_summ$Qdens)
extracted_Qdens = as.data.frame(t(extracted_Qdens))
#"growth" rates
extracted_Qrate = as.data.frame(bayesian_summ$Qrate)
extracted_Qrate = as.data.frame(t(extracted_Qrate))


#extra wrangling before exporting
bayesian_df_KX = as.data.frame(bayesian_summ$tau) %>% 
  mutate(time = bayesian_summ$tau,
         BayDEM_SPD = bayesian_summ$f_spdf,
         GaussMix_50 = extracted_Qdens$V2,
         GaussMix_2.5 = extracted_Qdens$V1,
         GaussMix_97.5 = extracted_Qdens$V3,
         Rate_2.5 = extracted_Qrate$V1,
         Rate_97.5 = extracted_Qrate$V3,
         Rate_50 = extracted_Qrate$V2,
         "bayesian_summ$tau" = NULL)



# a quick plot to check it has all worked fine
ggplot(bayesian_df_KX, aes(x = time)) +
  geom_line(aes(y = GaussMix_50, colour = "Bayesian"), linetype = 2)+
  geom_line(aes(y = BayDEM_SPD, colour = "BayDem SPD"))+
  geom_ribbon(aes(ymin=GaussMix_2.5, ymax=GaussMix_97.5, fill = "band"), alpha = 0.3)+
  scale_y_continuous(limits = c(0,.00125), breaks=seq(0,.00125, .0005))+
  scale_x_continuous(limits = c(-3300,-1200), breaks=seq(-3200,-1200, 200), sec.axis = sec_axis(~ . *-1+1950, name = "BP", breaks = seq(1000,8000, 500))) +
  theme_bw()


##write the data frame to file for use in 'final_analysis'
#write_rds(bayesian_df_KX,"Master Paper/GMM/outputs/K_final/arabia_dens_Kbest.rds")





