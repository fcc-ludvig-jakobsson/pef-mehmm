rm(list = ls())

library(tidyverse)
library(arm)
library(reshape2)
library(parallel)
library(matrixStats)

source("./MEHMM_SAEM.R")

log_pef <- FALSE
# The inverse of the parameter transformation function g
inv_g <- function(phi_i) return(c(exp(phi_i[1]), invlogit(phi_i[2]), exp(phi_i[3]), invlogit(phi_i[4])/2, invlogit(phi_i[5])/2))


## Simulate population, select true population parameters
mu_0_mean <- log(250)
d_mean <- 1
sigma2_mean <- log(300)
p_mean <- c(logit(0.1), logit(0.1))
mu_var <- 0.1
d_var <- 0.25
sigma2_var <- 0.05
p_var <- c(0.4,0.3)

mu_true <- c(mu_0_mean, d_mean, sigma2_mean, p_mean)
omega_true <- c(mu_var, d_var, sigma2_var, p_var)

# Select dose types per parameter
dose_types <- c('categorical', 'categorical', 'categorical', 'categorical', 'categorical')
beta_true <- list()
beta_true[[1]] <- 0.25
beta_true[[2]] <- 0
beta_true[[3]] <- 0
beta_true[[4]] <- -0.5
beta_true[[5]] <- 0

N <- 200
T <- 200
data <- matrix(, nrow=N, ncol=T)
phi_true <- matrix(, nrow=N, ncol=length(mu_true))
true_states <- matrix(, nrow=N, ncol=T+1)

x <- sample(0:1,N,replace = TRUE) # Binary doses
#x <- sample(c(0,1,3,10),N,replace=TRUE) # Continuous doses
for (i in 1:N) {
  sim <- simulate_two_state_model(T, x[i], mu_true, beta_true, omega_true, dose_types, inv_g, log_pef=log_pef)
  data[i,] <- sim$y
  phi_true[i,] <- sim$phi_i
  true_states[i,] <- sim$states
}


## Run estimation

K <- 700 # SAEM iterations
min_explore <- 150 # Minimum iterations of explore phase
convergence_window <- 50
convergence_tol <- 10 # Convergence indicator parameters for phase switch
burn_in <- 5
M1 <- 2
M2 <- 2
M <- M1+M2 # MCMC iteration within SAEM
M1_se <- 50
M2_se <- 150 # MCMC iterations after SAEM

# Initial fixed effects
init_mu <- c(5.5,0,6,-1,-1)
init_beta <- list(0,0,0,0,0)

target_rate <- 0.4

est_params <- c(1,2,3,4,5) # Which parameters to estimate
init_omega <- c(0,0,0,0,0)
init_omega[est_params] <- 1 # Initial random effect variance

start_time <- Sys.time()
# Run SAEM
saem_result <- SAEM(data, x, K, burn_in, M1, M2, init_mu, init_beta, init_omega, convergence_window, convergence_tol, min_explore, target_rate, dose_types, est_params, inv_g, log_pef=log_pef)

mu_hat <- saem_result$mu
beta_hat <- saem_result$beta
omega_hat <- saem_result$omega


## Extract final estimates
final_phi <- saem_result$phi[,K+1,]
final_mu <- mu_hat[K+1,]
final_beta <- list(NA,NA,NA,NA,NA)
final_beta[[1]] <- beta_hat[[1]][K+1]
final_beta[[2]] <- beta_hat[[2]][K+1]
final_beta[[3]] <- beta_hat[[3]][K+1]
final_beta[[4]] <- beta_hat[[4]][K+1]
final_beta[[5]] <- beta_hat[[5]][K+1]
final_omega <- omega_hat[K+1,]
scaling_factor <- saem_result$scaling_factor
phase_switch <- saem_result$phase_switch
convergence_indicator <- saem_result$convergence_indicator

print(final_mu)
print(final_beta)
print(final_omega)

# Draw posterior samples given theta_hat (parallel)
cl <- makeCluster(detectCores()-1)
clusterExport(cl, c("K",  "data", "x", "M1_se", "M2_se", "mu_true", "beta_true", "omega_true", "final_mu", "beta_hat", "final_beta", "final_omega", "final_phi", "scaling_factor", "dose_types", "est_params", "log_pef", "inv_g", "MCMC", "forward_algorithm", "IMP_log_likelihood", "invlogit"))
MCMC_samples <- parLapply(cl, 1:N, function(i) MCMC(K+1,data[i,],x[i],M1_se,M2_se,final_mu,beta_hat,final_omega,final_phi[i,],scaling_factor,dose_types,est_params,inv_g,log_pef=log_pef))
stopCluster(cl)

# log_likelihood, MAP, Viterbi, standard errors
log_likelihood <- sum(sapply(1:N, function(i) IMP_log_likelihood(MCMC_samples[[i]])))

MAP_samples <- do.call(rbind,lapply(lapply(1:N, function(i) MAP(MCMC_samples[[i]])$MAP_sample), unlist))

viterbi_paths <- do.call(rbind,lapply(1:N, function(i) viterbi(data[i,], inv_g(MAP_samples[i,]))))

residuals <- t(sapply(1:N, function(i) calculate_residuals(data[i,], exp(MAP_samples[i,1]), invlogit(MAP_samples[i,2]), viterbi_paths[i,])))

confusion_matrix <- Reduce('+', lapply(1:N, function(i) calc_confusion_matrix(T, true_states[i,], viterbi_paths[i,]))) / N

info_matrix <- Reduce('+',lapply(1:N, function(i) info_matrix_estimation(x[i], final_mu, final_beta, final_omega, MCMC_samples[[i]], dose_types, est_params)))
se <- sqrt(diag(solve(info_matrix)))

elapsed_time <- Sys.time() - start_time
