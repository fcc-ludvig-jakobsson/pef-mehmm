library(tidyverse)
library(arm)
library(reshape2)
library(parallel)
library(matrixStats)
library(zoo)
library(minpack.lm)

simulate_two_state_model <- function(T, x_i, mu, beta, omega, dose_types, inv_g, log_pef=FALSE) {
  # Compute means based on assumed dose-response
  loc <- numeric(5)
  for (c in 1:5) {
    if (dose_types[c] == 'none') {
      loc[c] = mu[c]
    }
    else if ( dose_types[c] == 'const') {
      loc[c] <- mu[c] + beta[[c]]*as.numeric(x_i > 0)
    }
    else if (dose_types[c] == 'categorical') {
      beta_bin <- append(beta[[c]], 0, after=0)
      
      loc[c] <- mu[c] + beta_bin[x_i+1]
    } 
    else if (dose_types[c] == 'linear') {
      loc[c] <- mu[c] + beta[[c]]*x_i
    }
    else if (dose_types[c] == 'emax') {
      emax <- beta[[c]][1]
      ed50 <- exp(beta[[c]][2])
      loc[c] <- mu[c] + emax*x_i / (ed50+x_i)
    }
  }
  

  # Draw individual model parameters
  phi_i <- rnorm(length(loc), loc, sqrt(omega))
  psi_i <- inv_g(phi_i)
  means <- c(psi_i[1], psi_i[1]*psi_i[2])

  if (log_pef==TRUE) means = log(means)
  
  # Build HMM time series
  y <- numeric(T)
  states = numeric(T+1)
  states[1] <- 0
  for (t in 1:T) {
    if (runif(1) < psi_i[4+states[t]]) {
      states[t+1] = ifelse(states[t]==0, 1, 0)
    } else {
      states[t+1] <- states[t]
    }
    y[t] = rnorm(1, means[states[t+1]+1], sqrt(psi_i[3]))
  }
  return(list(phi_i=phi_i, y=y, states=states))
}

forward_algorithm <- function(data_i, psi_i, log_pef=FALSE) {
  
  log_sum_exp <- function(v) {
    a <- v[1]
    b <- v[2]
    max <- max(a,b)
    return(max + log(exp(a-max)+exp(b-max)))
  }
  
  T <- length(data_i)
  mu_i <- psi_i[1]
  mu_i[2] <- mu_i*psi_i[2]
  if (log_pef==TRUE) mu_i <- log(mu_i)
  sigma2_i <- psi_i[3]
  p_i <- psi_i[4:5]
  P <- matrix(c(1-p_i[1],p_i[1],p_i[2],1-p_i[2]), nrow=2, byrow=T)
  log_P <- log(P)
  
  alpha <- matrix(,nrow = T+1, ncol=2)
  log_alpha <- matrix(-Inf, nrow = T+1, ncol = 2)
  
  alpha[1,] <- c(1,0)
  log_alpha[1,] <- c(0,-Inf)
  log_Z <- 0
  for (t in 1:T) {
    for (k in 1:2) {
      log_alpha[t+1,k] <- dnorm(data_i[t], mu_i[k], sqrt(sigma2_i),log=T) + log_sum_exp(log_alpha[t,]+log_P[,k])
    }
    log_normalizing_constant <- log_sum_exp(log_alpha[t+1,])
    log_Z <- log_Z + log_normalizing_constant # Update log_likelihood
    log_alpha[t+1,] <- log_alpha[t+1,] - log_normalizing_constant # Normalized log-scale
    alpha[t+1,] <- exp(log_alpha[t+1,]) / sum(exp(log_alpha[t+1,])) # convert to normalized prob-scale
  }
  return(list(alpha=alpha[2:(T+1),],log_alpha=log_alpha[2:(T+1),],log_Z=log_Z))
}

MCMC <- function(k, data_i, x_i, M1, M2, mu_k, beta, omega_k, init_phi_i, scaling_factor, dose_types, est_params, inv_g, log_pef=FALSE) {
  n_params <- length(mu_k)
  M <- M1+M2
  y_max <- ifelse(log_pef==FALSE, max(data_i), max(exp(data_i)))
  
  phi_i <- matrix(, nrow=M+1, ncol=n_params)
  phi_i[1,] <- init_phi_i
  
  log_Z <- numeric(M+1)
  component_log_pdf <- matrix(, nrow=M+1, ncol=length(est_params))
  log_likelihood <- numeric(M+1)
  
  log_Z_old <- forward_algorithm(data_i, inv_g(phi_i[1,]),log_pef=log_pef)$log_Z
  log_Z[1] <- log_Z_old
  
  loc <- numeric(5)
  loc[-est_params] <- mu_k[-est_params]
  for (i in 1:length(est_params)) {
    c <- est_params[i]
    if (dose_types[i] == 'none') {
      loc[c] = mu_k[c]
    }
    else if ( dose_types[i] == 'const') {
      loc[c] <- mu_k[c] + beta[[c]][k]*as.numeric(x_i > 0)
    }
    else if (dose_types[i] == 'categorical') {
      beta_bin <- append(beta[[c]][k,], 0, after=0)
      loc[c] <- mu_k[c] + beta_bin[x_i+1]
    } 
    else if (dose_types[i] == 'linear') {
      loc[c] <- mu_k[c] + beta[[c]][k]*x_i
    }
    else if (dose_types[i] == 'emax') {
      emax <- beta[[c]][k,1]
      ed50 <- exp(beta[[c]][k,2])
      loc[c] <- mu_k[c] + emax*x_i / (ed50+x_i)
    }
  }
  
  scale <- sqrt(omega_k)
  
  component_log_pdf[1,] <- dnorm(phi_i[1,], mean=loc, sd=scale,log=T)[est_params]
  log_likelihood[1] <- log_Z[1] + sum(component_log_pdf[1,])
  num_accepts <- numeric(n_params)
  if (M1 > 0) {
    for (m in 1:M1) {
      phi_i[m+1,] <- rnorm(n_params, mean=loc, sd=scale)
      phi_i[m+1,1] <- min(phi_i[m+1,1], log(y_max))
      log_Z_new <- forward_algorithm(data_i, inv_g(phi_i[m+1,]),log_pef=log_pef)$log_Z
      component_log_pdf[m+1,] <- dnorm(phi_i[m+1,], mean=loc, sd=scale,log=T)[est_params]
      log_alpha <- min(0, log_Z_new - log_Z_old + sum(component_log_pdf[m+1,]) - sum(component_log_pdf[m,]))
      if (runif(1) > exp(log_alpha)) {
        phi_i[m+1,] <- phi_i[m,] # Reject proposal
        component_log_pdf[m+1,] <- component_log_pdf[m,]
      } else {
        log_Z_old <- log_Z_new # Accept proposal
      }
      log_Z[m+1] <- log_Z_old
      log_likelihood[m+1] <- log_Z[m+1] + sum(component_log_pdf[m+1,])
    }
  }
  for (m in (M1+1):M) {
    phi_i[m+1,] <- phi_i[m,]
    component_log_pdf[m+1,] <- component_log_pdf[m,]
    phi_i_proposal <- rnorm(n_params, phi_i[m,], sqrt(scaling_factor)*scale)
    phi_i_proposal[1] <- min(phi_i_proposal[1], log(y_max))
    order <- sample(1:length(est_params),length(est_params))
    for (i in order) {
      c <- est_params[i]
      phi_i[m+1,c] <- phi_i_proposal[c]
      log_Z_new <- forward_algorithm(data_i, inv_g(phi_i[m+1,]),log_pef=log_pef)$log_Z
      component_log_pdf[m+1,i] <- dnorm(phi_i[m+1,c], mean=loc[c], sd=scale[c],log=T)
      log_alpha <- min(0, log_Z_new - log_Z_old + component_log_pdf[m+1,i] - component_log_pdf[m,i])
      if (runif(1) > exp(log_alpha)) {
        phi_i[m+1,c] <- phi_i[m,c] # Reject proposal
        component_log_pdf[m+1,i] <- component_log_pdf[m,i]
      } else {
        log_Z_old <- log_Z_new # Accept proposal
        num_accepts[c] <- num_accepts[c] + 1
      }
      log_Z[m+1] <- log_Z_old
    }
    log_likelihood[m+1] <- log_Z[m+1] + sum(component_log_pdf[m+1,])
  }
  return(list(phi_i = phi_i, num_accepts = num_accepts, log_Z=log_Z, log_likelihood = log_likelihood))
}

SAEM <- function(data, x, K, burn_in, M1, M2, init_mu, init_beta, init_omega, convergence_window, convergence_tol, min_explore, target_rate, dose_types, est_params, inv_g, data_type='array', log_pef=FALSE) {
  n_params <- length(init_mu)

  phase <- 'explore'
  gamma <- numeric(K)+1
  acceptance_rate <- matrix(, nrow=K, ncol=n_params)
  scaling_factor <- numeric(n_params) + 1
  convergence_indicator <- numeric(K)
  
  mu <- matrix(rep(init_mu,K+1), nrow=K+1, ncol=n_params, byrow=T)
  beta <- list()
  init_omega[-est_params] <- 0
  omega <- matrix(rep(init_omega,K+1), nrow=K+1, ncol=n_params, byrow=T)
  
  for (i in 1:length(est_params)) {
    c <- est_params[i]
    if (dose_types[i] == 'none') {
      # No dose parameters
    }
    else if (dose_types[i] == 'const') {
      beta[[c]] <- rep(init_beta[[c]],K+1)
    }
    else if (dose_types[i] == 'categorical') {
      beta[[c]] <- matrix(rep(init_beta[[c]],K+1), nrow=K+1, byrow=TRUE)
    }
    else if (dose_types[i] == 'linear') {
      beta[[c]] <- rep(init_beta[[c]],K+1)
    }
    else if (dose_types[i] == 'emax') {
      beta[[c]] <- matrix(rep(init_beta[[c]],K+1), nrow=K+1, byrow=TRUE)
      beta[[c]][,2] <- log(median(x))
    }
  }
  
  N <- length(x)
  phi <- array(dim=c(N,K+1,n_params))
  phi[,1,] <- matrix(rep(init_mu,N), ncol=n_params, byrow=T)
  
  cl <- makeCluster(detectCores()-1)
  clusterExport(cl, c("data", "x", "M1", "M2", "dose_types", "est_params", "log_pef", "inv_g", "MCMC", "forward_algorithm", "IMP_log_likelihood"))
  clusterEvalQ(cl, library(arm))
  
  mh_individual_update <- function(i,k,mu,beta,omega,phi,scaling_factor) {
    data_i <- if (data_type=='table') data[data$ID == i,]$PEF_interp else data[i,]
    result <- MCMC(k, data_i, x[i], M1, M2, mu[k,], beta, omega[k,], phi[i,k,], scaling_factor, dose_types, est_params, inv_g, log_pef=log_pef)
    return(list(phi = result$phi[M1+M2+1,], num_accepts = result$num_accepts, log_likelihood = result$log_likelihood[M1+M2+1]))
  }
  
  for (k in 1:K) {
    results <- parLapply(cl, 1:N, mh_individual_update, k, mu, beta, omega, phi, scaling_factor)
    phi[1:N,k+1,] <- do.call(rbind, lapply(results, function(res) res$phi))
    acceptance_rate[k,] <- colSums(do.call(rbind, lapply(results, function(res) res$num_accepts))) / (N*M2)
    convergence_indicator[k] <- Reduce('+',lapply(results, function(res) res$log_likelihood))
    
    if (k > burn_in) {
      if (phase=='explore' & k > min_explore) {
        avg_convergence <- rollmean(convergence_indicator,convergence_window,na.pad=TRUE,align='right')
        if ( (max(avg_convergence[(k-convergence_window):k])-min(avg_convergence[(k-convergence_window):k]) < convergence_tol) 
             | (k > K/2) ) {
          phase <- 'smoothing'
          phase_switch <- k
          gamma[k:K] <- 1e-2^(1:(K+1-k)/(K+1-k))
        }
      }
      for (i in 1:length(est_params)) {
        c <- est_params[i]
        phi_c <- phi[,k+1,c]
        if (dose_types[i] == 'none') {
          mu[k+1,c] <- mu[k,c] + gamma[k]*(mean(phi_c) - mu[k,c])
          omega[k+1,c] <- omega[k,c] + gamma[k]*(var(phi_c) - omega[k,c])
        }
        else if (dose_types[i] == 'const') {
          mod <- lm(phi_c ~ factor(as.numeric(x>0)))
          coefs <- coefficients(mod)
          mu[k+1,c] <- mu[k,c] + gamma[k]*(coefs[1] - mu[k,c])
          beta[[c]][k+1] <- beta[[c]][k] + gamma[k]*(coefs[2] - beta[[c]][k])
          omega[k+1,c] <- omega[k,c] + gamma[k]*(sum(residuals(mod)^2)/(N-2) - omega[k,c])
        }
        else if (dose_types[i] == 'categorical') {
          mod <- lm(phi_c ~ factor(x))
          p <- nlevels(factor(x))
          coefs <- coefficients(mod)
          mu[k+1,c] <- mu[k,c] + gamma[k]*(coefs[1] - mu[k,c])
          beta[[c]][k+1,] <- beta[[c]][k,] + gamma[k]*(coefs[-1] - beta[[c]][k,])
          omega[k+1,c] <- omega[k,c] + gamma[k]*(sum(residuals(mod)^2)/(N-p) - omega[k,c])
        }
        else if (dose_types[i] == 'linear') {
          mod <- lm(phi_c ~ x)
          coefs <- coefficients(mod)
          mu[k+1,c] <- mu[k,c] + gamma[k]*(coefs[1] - mu[k,c])
          beta[[c]][k+1] <- beta[[c]][k] + gamma[k]*(coefs[2] - beta[[c]][k])
          omega[k+1,c] <- omega[k,c] + gamma[k]*(sum(residuals(mod)^2)/(N-2) - omega[k,c])
        }
        else if (dose_types[i] == 'emax') {
          fcn <- function(par, y, x) {
            ed50 <- exp(par[3])
            pred <- par[1] + par[2]*x / (ed50+x)
            return(y-pred)
          }
          jac <- function(par, y, x) {
            d_e0 <- -1
            d_emax <- -x / (exp(par[3])+x)
            d_ed50 <- par[2]*x / (exp(par[3])+x)^2 * exp(par[3])
            return(cbind(d_e0, d_emax, d_ed50))
          }
          init_par <- c(mu[k,c], beta[[c]][k,1], beta[[c]][k,2])
          res <- nls.lm(par = init_par, lower=c(-Inf, -Inf, log(1e-4)), upper = c(Inf, Inf, log(5*max(x))), fn = fcn, jac = jac, y = phi_c, x = x)
          #print(paste0("Termination: ", res$info))
          residuals <- res$fvec
          
          mu[k+1,c] <- mu[k,c] + gamma[k]*(res$par[1] - mu[k,c])
          beta[[c]][k+1,] <- beta[[c]][k,] + gamma[k]*(res$par[2:3] - beta[[c]][k,])
          omega[k+1,c] <- omega[k,c] + gamma[k]*(sum(residuals^2)/(N-3) - omega[k,c])
        }
      }
     
      # Update proposal variance scaling every tenth iteration
      if (k %% 10 == 0) {
        scaling_factor[est_params] <- scaling_factor[est_params] * exp(colMeans(acceptance_rate[(k-9):k,est_params,drop=F]) - target_rate)
      }
    
    }
  }
  stopCluster(cl)
  return(list(mu=mu, beta=beta, omega=omega, phi=phi, acceptance_rate = acceptance_rate, scaling_factor=scaling_factor, phase_switch=phase_switch, convergence_indicator=convergence_indicator))
}

MAP <- function(MCMC_samples_i) {
  phi_i_samples <- MCMC_samples_i$phi_i

  MAP_sample <- colMeans(phi_i_samples)
  
  return(list(MAP_sample = MAP_sample))
  
}

viterbi <- function(data_i, psi_i, log_pef=FALSE) {
  T <- length(data_i)
  mu0 <- psi_i[1]
  d <- psi_i[2]
  sigma <- sqrt(psi_i[3])
  p01 <- psi_i[4]
  p10 <- psi_i[5]
  
  # Derived parameters
  mu1 <- mu0 * d
  p00 <- 1 - p01
  p11 <- 1 - p10
  
  if (log_pef==TRUE) {
    mu0 <- log(mu0)
    mu1 <- log(mu1)
  }
  
  # Log probabilities for numerical stability
  log_emission <- matrix(0, nrow = 2, ncol = T)
  log_emission[1, ] <- dnorm(data_i, mean = mu0, sd = sigma, log = TRUE)
  log_emission[2, ] <- dnorm(data_i, mean = mu1, sd = sigma, log = TRUE)
  
  log_trans <- matrix(c(log(p00), log(p01), 
                        log(p10), log(p11)), nrow = 2, byrow = TRUE)
  
  # Initialization
  delta <- matrix(-Inf, nrow = 2, ncol = T)
  backpointer <- matrix(0, nrow = 2, ncol = T)
  
  delta[,1] <- log(0.5) + log_emission[,1]
  # Recursion
  for (t in 2:T) {
    for (j in 1:2) {
      probs <- delta[, t-1] + log_trans[, j]
      backpointer[j, t] <- which.max(probs)
      delta[j, t] <- max(probs) + log_emission[j, t]
    }
  }
  # Backtracking
  states <- numeric(T)
  states[T] <- which.max(delta[, T]) - 1  # Convert to 0/1
  for (t in (T-1):1) {
    states[t] <- backpointer[states[t+1] + 1, t + 1] - 1
  }
  
  return(states)
}

calc_confusion_matrix <- function(T, true_states, viterbi_path) {
  state_df <- data.frame(true_states = true_states[2:(T+1)], viterbi_path=viterbi_path)
  
  TPs <- nrow(state_df %>% filter(true_states == 1, viterbi_path == 1))
  FNs <- nrow(state_df %>% filter(true_states == 1, viterbi_path == 0))
  FPs <- nrow(state_df %>% filter(true_states == 0, viterbi_path == 1))
  TNs <- nrow(state_df %>% filter(true_states == 0, viterbi_path == 0))
  confusion_matrix <- round(matrix(c(TPs,FNs,FPs,TNs)/T, nrow=2, ncol=2), 2)
  
  return(confusion_matrix)
}

info_matrix_estimation <- function(x_i,mu,beta,omega,MCMC_samples_i,dose_types,est_params,est_beta=1:5) {
  analytical_grad <- function(phi_i) {
    grad <- numeric()
    for (i in 1:length(est_params)) {
      c <- est_params[i]
      phi_ic <- phi_i[c]
      mu_c <- mu[c]
      beta_c <- beta[[c]]
      omega_c <- omega[c]

      if (dose_types[i] == 'none') {
        r <- phi_ic - mu_c
        d_mu <- r / omega_c
        d_omega <- r^2 / (2*omega_c^2) - 1/(2*omega_c)
        grad <- append(grad, c(d_mu, d_omega, use.names=FALSE))
      }
      else if (dose_types[i] == 'const') {
        r <- phi_ic - mu_c - beta_c*as.numeric(x_i > 0)
        d_mu <- r / omega_c
        d_beta <- r * as.numeric(x_i > 0) / omega_c
        d_omega <- r^2 / (2*omega_c^2) - 1/(2*omega_c)
        grad <- append(grad, c(d_mu, d_beta, d_omega, use.names=FALSE))
      }
      else if (dose_types[i] == 'categorical') {
        v_i <- rep(0,length(beta_c))
        if (x_i != 0) v_i[x_i] <- 1
        r <- phi_ic - mu_c - as.numeric(beta_c%*%v_i)
        d_mu <- r / omega_c
        d_beta <- as.vector(v_i %o% (r / omega_c))
        d_omega <- r^2 / (2*omega_c^2) - 1/(2*omega_c)
        grad <- append(grad, c(d_mu, d_beta, d_omega, use.names=FALSE))
      }
      else if (dose_types[i] == 'linear') {
        r <- phi_ic - mu_c - beta_c*x_i
        d_mu <- r / omega_c
        d_beta <- r*x_i / omega_c
        d_omega <- r^2 / (2*omega_c^2) - 1/(2*omega_c)
        grad <- append(grad, c(d_mu, d_beta, d_omega, use.names=FALSE))
      }
      else if (dose_types[i] == 'emax') {
        dose_effects <- beta_c[1]*x_i / (exp(beta_c[2])+x_i)
        r <- phi_ic - mu_c - dose_effects
        D <- exp(beta_c[2])+x_i
        d_mu <- r / omega_c
        d_emax <- x_i / D * r / omega_c
        d_ed50 <- - exp(beta_c[2]) * beta_c[1]*x_i / D^2 * r / omega_c
        d_omega <- r^2 / (2*omega_c^2) - 1/(2*omega_c)
        grad <- append(grad, c(d_mu, d_emax, d_ed50, d_omega, use.names=FALSE))
      }
    }
    return(grad)
  }
  analytical_hessian <- function(phi_i) {
    sub_hessians <- list()
    for (i in 1:length(est_params)) {
      c <- est_params[i]
      phi_ic <- phi_i[c]
      mu_c <- mu[c]
      beta_c <- beta[[c]]
      omega_c <- omega[c]
      
      if (dose_types[i] == 'none') {
        r <- phi_ic - mu_c
        d_mu_d_mu <- - 1 / omega_c
        d_omega_d_omega <- - r^2/(omega_c^3) + 1/(2*omega_c^2)
        d_mu_d_omega <- - r / omega_c^2
        sub_hessians[[i]] <- matrix(c(d_mu_d_mu, d_mu_d_omega, 0, d_omega_d_omega), byrow=T, nrow=2)
      }
      if (dose_types[i] == 'const') {
        r <- phi_ic - mu_c - beta_c*as.numeric(x_i > 0)
        d_mu_d_mu <- - 1 / omega_c
        d_beta_d_beta <- - as.numeric(x_i > 0)^2 / omega_c
        d_omega_d_omega <- - r^2/(omega_c^3) + 1/(2*omega_c^2)
        d_mu_d_beta <- - as.numeric(x_i > 0) / omega_c
        d_mu_d_omega <- - r / omega_c^2
        d_beta_d_omega <- - r * as.numeric(x_i > 0) / omega_c^2
        sub_hessians[[i]] <- matrix(c(d_mu_d_mu, d_mu_d_beta, d_mu_d_omega, 0, d_beta_d_beta, d_beta_d_omega, 0, 0, d_omega_d_omega), byrow=T, nrow=3)
      }
      if (dose_types[i] == 'categorical') {
        v_i <- rep(0,length(beta_c))
        if (x_i != 0) v_i[x_i] <- 1
        r <- phi_ic - mu_c - as.numeric(beta_c%*%v_i)
        d_mu_d_mu <- - 1 / omega_c
        d_beta_d_beta <- -as.vector(outer(v_i, omega_c, '/'))
        d_omega_d_omega <- - r^2/(omega_c^3) + 1/(2*omega_c^2)
        d_mu_d_beta <- -as.vector(outer(v_i,omega_c, '/'))
        d_mu_d_omega <- - r / omega_c^2
        d_beta_d_omega <- - as.vector(outer(v_i, r / omega_c^2))
        sub_hessians[[i]] <- diag(c(d_mu_d_mu, d_beta_d_beta, d_omega_d_omega))
        sub_hessians[[i]][1, 2:(1+length(v_i))] <- d_mu_d_beta
        sub_hessians[[i]][1,ncol(sub_hessians[[i]])] <- d_mu_d_omega
        sub_hessians[[i]][2:(1+length(v_i)),ncol(sub_hessians[[i]])] <- d_beta_d_omega
      }
      if (dose_types[i] == 'linear') {
        r <- phi_ic - mu_c - beta_c*x_i
        d_mu_d_mu <- - 1 / omega_c
        d_beta_d_beta <- - x_i^2 / omega_c
        d_omega_d_omega <- - r^2/(omega_c^3) + 1/(2*omega_c^2)
        d_mu_d_beta <- - x_i / omega_c
        d_mu_d_omega <- - r / omega_c^2
        d_beta_d_omega <- - r*x_i / omega_c^2
        sub_hessians[[i]] <- matrix(c(d_mu_d_mu, d_mu_d_beta, d_mu_d_omega, 0, d_beta_d_beta, d_beta_d_omega, 0, 0, d_omega_d_omega), byrow=T, nrow=3)
      }
      if (dose_types[i] == 'emax') {
        dose_effects <- beta_c[1]*x_i / (exp(beta_c[2])+x_i)
        r <- phi_ic - mu_c - dose_effects
        D <- exp(beta_c[2])+x_i
        d_mu <- r / omega_c
        
        d_mu_d_mu <- - 1 / omega_c
        d_emax_d_emax <- - x_i^2 / D^2 / omega_c
        d_ed50_d_ed50 <- - exp(beta_c[2]) * beta_c[1]*x_i / D^2 / omega_c * (exp(beta_c[2])*beta_c[1]*x_i / D^2 - r + 2 * r / D)
        d_omega_d_omega <- - r^2/(omega_c^3) + 1/(2*omega_c^2)
        d_mu_d_emax <- - x_i / D / omega_c
        d_mu_d_ed50 <- exp(beta_c[2]) * beta_c[1]*x_i / D^2 / omega_c
        d_mu_d_omega <- - r / (omega_c^2)
        d_emax_d_ed50 <- - exp(beta_c[2])*x_i / D^2 / omega_c * (r - dose_effects)
        d_emax_d_omega <- - x_i / D * r / (omega_c^2)
        d_ed50_d_omega <- exp(beta_c[2]) * beta_c[1]*x_i / D^2 * r / (omega_c^2)
        
        sub_hessians[[i]] <- matrix(c(d_mu_d_mu, d_mu_d_emax, d_mu_d_ed50, d_mu_d_omega, 0, d_emax_d_emax, d_emax_d_ed50, d_emax_d_omega, 0, 0, d_ed50_d_ed50, d_ed50_d_omega, 0, 0, 0, d_omega_d_omega), byrow=T, nrow=4)
      }
    }
    hessian <- as.matrix(bdiag(sub_hessians))
    hessian[lower.tri(hessian)] <- t(hessian)[lower.tri(hessian)]
    return(hessian)
  }
  
  phi_i <- MCMC_samples_i$phi_i
  M <- dim(phi_i)[1]
  
  gradients <- do.call(rbind, lapply(1:M, function(m) analytical_grad(phi_i[m,])))
  hessians <- lapply(1:M, function(m) analytical_hessian(phi_i[m,]))
    
  hessian_mean <- Reduce('+', hessians) / M
  gradient_covariance <- cov(gradients)

  info_matrix_i <- - hessian_mean - gradient_covariance

  return(info_matrix_i)
}

IMP_log_likelihood <- function(data_i, x_i, M, mu, beta, omega, dose_types, est_params, inv_g, log_pef=FALSE) {
  
  log_sum_exp <- function(v) {
    max <- max(v)
    return(max + log(sum(exp(v-max))))
  }
  
  loc <- numeric(5)
  loc[-est_params] <- mu[-est_params]
  for (i in 1:length(est_params)) {
    c <- est_params[i]
    if (dose_types[i] == 'none') {
      loc[c] = mu[c]
    }
    else if ( dose_types[i] == 'const') {
      loc[c] <- mu[c] + beta[[c]]*as.numeric(x_i > 0)
    }
    else if (dose_types[i] == 'categorical') {
      beta_bin <- append(beta[[c]], 0, after=0)
      loc[c] <- mu[c] + beta_bin[x_i+1]
    } 
    else if (dose_types[i] == 'linear') {
      loc[c] <- mu[c] + beta[[c]]*x_i
    }
    else if (dose_types[i] == 'emax') {
      emax <- beta[[c]][1]
      ed50 <- exp(beta[[c]][2])
      loc[c] <- mu[c] + emax*x_i / (ed50+x_i)
    }
  }
  
  scale <- sqrt(omega)
  
  phi_samples <- matrix(rnorm(5*M, loc, scale),nrow=M,byrow=T)
  
  log_Z <- sapply(1:M, function(m) forward_algorithm(data_i,inv_g(phi_samples[m,]))$log_Z, log_pef=log_pef)
  
  log_likelihood <- log_sum_exp(log_Z) - log(M)
  
  return(log_likelihood)
}

calculate_residuals <- function(data_i, m_i, d_i, viterbi_path_i) {
  mean_path <- m_i*(1 - (1-d_i)*viterbi_path_i)
  residuals <- data_i - mean_path
}
