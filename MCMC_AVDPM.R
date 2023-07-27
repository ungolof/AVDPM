library(mvtnorm)
library(invgamma)
library(truncnorm)
library(LaplacesDemon)
library(dplyr)
library(MASS)
library(TruncatedDistributions)


##===================== - Final dataset for inference - =============================================
### Comment: The only covariate here is the age of the couple member (1) for male and (2) for female
# - version of the implementation with the baseline alpha

# - Additional cleaning of data
data <- data[data$agediff!=0,]
data <- data[data$TObs>0,]

data <- data[data$ageM+data$TObs>50,]
data <- data[data$ageF+data$TObs>50,]


sample_size <- nrow(data)

data_shuffle <- data[sample(1:nrow(data)),]
data_final <- data_shuffle[1:round(nrow(data_shuffle)*0.75,0),]

X <- as.data.frame(cbind(data_final$ageM, data_final$ageF) - data_final$trunc)[3,]
#X <-  as.data.frame(cbind(data$ageM, data$ageF))
colnames(X) <- c("male", "female")
Z <- cbind(data_final$agediff2, data_final$maleelder)
colnames(Z) <- c("agediff", "maleelder")
Z <- as.data.frame(Z)
Z$lAD <- log(abs(Z$agediff))

t_i <- data_final$TObs
a_i <- data_final$trunc
d_ci <- as.data.frame(cbind(data_final$dmale, data_final$dfemale))
colnames(d_ci) <- c("male", "female")

x_bar <- 70
X <- X - x_bar

sample_size <- nrow(data_final)

#================================ - Parameters storage - =============================================

n_iter <- 50000 #number of iterations, set by the researcher

K = 25 # Dunson (2010) claims a number of pieces between 20 and 50
M <- 2

# - beta regression parameters
alpha <- matrix(0, n_iter, 2)
beta <- matrix(0, n_iter, 2)
colnames(beta) <- colnames(alpha) <- colnames(X) <- c("male", "female")

# - gamma frailty parameters
gamma_star <- matrix(0, n_iter, M*K)
colnames(gamma_star) <- rep(c(sprintf("gamma_s_%d", c(1:M))), K)

z_ad_v <- rep(0, n_iter) # - zeta for the common variance of the age difference distribution
z_star_ad <- matrix(0, n_iter, K) # - zeta* for the age difference - conjugate normal
z_star_mo <- matrix(0, n_iter, K) # - zeta* for the elder male (Bernoulli distributed, with Beta(3,1) baseline)

# - pi mixture weight
pi_mw <- matrix(NA, n_iter, K)
colnames(pi_mw) <- c(sprintf("pi_%d", c(1:K)))

# - mu_gamma (parameter of the base distribution)
#mu_gamma <- matrix(NA, n_iter, M)
#colnames(mu_gamma) <- c(sprintf("mu_gamma_%d", c(1:M)))

# - Covariance_gamma
Sigma_gamma <- matrix(NA, n_iter, (M^2))

# - mu_zeta (for the normal distribution of the age-difference)
mu_zeta <- rep(0, n_iter)
sigma2_zeta <- rep(1, n_iter)

# - psi from stick-breaking procedure
psi_sbp <- matrix(NA, n_iter, K)
colnames(psi_sbp) <- c(sprintf("psi_%d", c(1:K)))

# phi concentration parameter of the Dirichlet Process
phi <- matrix(NA, n_iter, 1)

# - Number of classes
N_k_tb <- matrix(NA, n_iter, K)
N_c_tb <- rep(0, n_iter)


############################### - Set prior distributions - #####################################

# - beta
mu_beta_prior <- rep(0, M)

Sigma_beta_prior <- rep(9, M) #rep(1, 2) # assuming independently distributed beta_cl

# - gamma
## - mu_gamma - set to zero!
#lambda1 <- rep(0, M) # - indicated in the paper as m
#lambda2 <-  9 * matrix(diag(1, M), M, M)# - 0.0001 * matrix(diag(1, M), M, M) # - indicated in the paper as B

## - Covariance gamma
lambda3 <- M + 5
lambda4 <- 0.001 * diag(1, M) %*% matrix(c(1,0.5,0.5,1), M, M, byrow=TRUE) %*% diag(1, M)

# - phi (reperform lambda5=0.5, lambda6=10)
lambda5 <- 1
lambda6 <- 1




####################################### - Starting values - ################################################

# beta
beta[1,] <- runif(2, 0.03, 0.2)# c(0.1, 0.08)#runif(M, 10^(-10), 0.3)# c(0, 0.14, 0.8, -0.2, -0.07, -0.3, 0.01) #st_val
alpha[1,] <- runif(2, -20, -2) #c(-10, -10)#runif(M, 10^(-10), 0.3)# c(0, 0.14, 0.8, -0.2, -0.07, -0.3, 0.01) #st_val

# - gamma frailty parameters
gamma_star[1,] <- rnorm(M*K, 0, 0.01)# c(0.2, 0.3, 0.1) #st_val

z_star_ad[1,] <- 0
z_ad_v[1] <- 1 / rgamma(1, 1, 1)

z_star_mo[1,] <- runif(K, 0, 1)#matrix(0, n_iter, K) # - zeta* for the elder male (Bernoulli distributed, with Beta(3,1) baseline)

# - mu_gamma (parameter of the base distribution)
# mu_gamma[1,] <- rnorm(M, 0, 1)# c(0.2, 0.3, 0.1) #st_val

# - Covariance_gamma
Sigma_gamma[1,] <- c(1,0.5,0.5,1)

# - psi from stick-breaking procedure
# - Initial value for the stick-breaking weights V (random)
phi[1,1] <- 1

psi_sbp[1,K] <- 1
for(k in 1:(K-1)){
  psi_sbp[1,k] <- rbeta(1, shape1 = 1, shape2 = phi[1,1], ncp = 0)
}

# - pi mixture weight
pi_mw[1, 1] <- psi_sbp[1]
for(k in 2:K){
  pi_mw[1,k] <- psi_sbp[1,k] * prod(1-psi_sbp[1,1:(k-1)])
}

## - Adaptive Metropolis-H. Vihola (2012)
Sigma_ga <- matrix(rep(diag(1, M), K), M, M*K) # see if possible to do rep... diag(1, ncol(post_par_MH)) * 0.0001
Sigma_ch_ga <- matrix(rep(t(chol(Sigma_ga[1:M,1:M])), K), M, M*K)# matrix(0, M, M*K) # <- t(chol(Sigma_ga)) # since chol return the upper triangular matrix
gamma_Vihola <- 0.6 # - should be between 0.5 and 1
opt_accept <- 0.234 # Efficient acceptance rate from Roberts and Rosenthal (2001)

sigma_beta_pr <- rep(1, M) * 0.001
beta_tg_accept <- 0.234


# - Preparatory
N_k <- rep(0, K)
prob_alloc_matrix <- matrix(NA, 1, K) # storage of the posterior weights from which to sample K

S_unit <- rep(0, sample_size)
S_unit <- sample(c(1:K), sample_size, prob=rep(1/K, K), replace=TRUE)

S_unit_fc <- function(unit){
  
  # - use log sum exp trick to deal with very small numbers which R inappropriately set to zero
  pf_vector <- rep(0, K)
  ga_k_pos <- c(1:M) # position vector for gamma
  
  pf_vector <- log(pi_mw[i-1,]) - exp(alpha[i-1,1] + beta[i-1,1] * X$male[unit]   + gamma_s1_lp) * (exp(beta[i-1,1] * (t_i[unit] + a_i[unit])) - exp(beta[i-1,1] * a_i[unit])) / beta[i-1,1] + d_ci$male[unit]   * gamma_s1_lp - 
                                  exp(alpha[i-1,2] + beta[i-1,2] * X$female[unit] + gamma_s2_lp) * (exp(beta[i-1,2] * (t_i[unit] + a_i[unit])) - exp(beta[i-1,2] * a_i[unit])) / beta[i-1,2] + d_ci$female[unit] * gamma_s2_lp + 
    dnorm(Z$lAD[unit], mean=gamma_s3_lp, sd=sqrt(z_ad_v[i-1]), log = TRUE) + Z$maleelder[unit] * log(z_star_mo[i-1,]) + (1 - Z$maleelder[unit]) * log(1 - z_star_mo[i-1,])
  
  max_pf <- max(pf_vector)
  
  denom_calc_lse <- max_pf + log(sum(exp(pf_vector-max_pf)))
  
  prob_alloc <- exp(pf_vector - denom_calc_lse) # - in vectorized form
  
  # - Sample cluster allocation
  S_unit_output <- sample(c(1:K), 1, prob=prob_alloc)
  
  return(S_unit_output)
}


i <- 2
for(i in i:n_iter){
  #  start_time <- Sys.time()
  
  # - Step 1: Sample new mixture component
  
  gamma_s1_lp <- gamma_star[i-1,seq(1, M*K-1, M)]
  gamma_s2_lp <- gamma_star[i-1,seq(2, M*K, M)]
  gamma_s3_lp <- z_star_ad[i-1,]
  
  S_unit <- sapply(c(1:sample_size), FUN = S_unit_fc)
  
  for(unit in 1:sample_size){
    S_unit[unit] <- S_unit_fc(unit)
  }
  
  #    profvis({
  
  # - Step 2.1: Sample stick-breaking weights
  for(k in 1:(K-1)){
    psi_sbp[i,k] <- rbeta(1, shape1 = 1 + sum(ifelse(S_unit==k,1,0)), shape2 = phi[i-1,1] + sum(ifelse(S_unit > k,1,0)), ncp = 0)
  }
  psi_sbp[,K] <- 1
  # - Step 2.2: Update mixture distribution pi
  ## - Computation of pi
  pi_mw[i, 1] <- psi_sbp[i,1]
  for(k in 2:K){
    pi_mw[i,k] <- psi_sbp[i,k] * prod(1-psi_sbp[i,1:(k-1)])
  }
  
  # - Step 2.9: Sample alpha_c
  alpha[i,1] <- log(rgamma(1, shape=sum(d_ci$male  ), rate=1 + sum(exp(beta[i-1,1] * X$male   + gamma_star[i-1,(1 + M * (S_unit - 1))]) * (exp((t_i + a_i) * beta[i-1,1]) - exp(a_i * beta[i-1,1])) / beta[i-1,1])))
  alpha[i,2] <- log(rgamma(1, shape=sum(d_ci$female), rate=1 + sum(exp(beta[i-1,2] * X$female + gamma_star[i-1,(2 + M * (S_unit - 1))]) * (exp((t_i + a_i) * beta[i-1,2]) - exp(a_i * beta[i-1,2])) / beta[i-1,2])))
  
  # - Step 3: Sample Beta_c
  # beta*
  
  beta_star <- rtruncnorm(1, a=(10^(-5)), mean = beta[i-1,1], sd = sigma_beta_pr[1], b=5)# rnorm(n = 1, beta[i-1,1], sigma_beta_pr[1])  #
  
  # Acceptance step
  #start_time <- Sys.time()
  # - Get energy functions + proposals
  log_ef <- log(dtruncnorm(beta[i-1,1], mean=beta_star, sd=sigma_beta_pr[1], a=(10^(-5)), b=5)) - log(dtruncnorm(beta_star, mean=beta[i-1,1], sd=sigma_beta_pr[1], a=(10^(-5)), b=5)) #  - dnorm(beta_star, mean=beta[i-1,1], sd=sigma_beta_pr[1], log=TRUE) + dnorm(beta[i-1,1], mean=beta_star, sd=sigma_beta_pr[1], log=TRUE)# log(dtruncnorm(beta[i-1,1], mean=beta_star, sd=sigma_beta_pr[1], a=(10^(-5)), b=5)) - log(dtruncnorm(beta_star, mean=beta[i-1,1], sd=sigma_beta_pr[1], a=(10^(-5)), b=5)) # 
  
  log_ef <- log_ef + sum( - exp(alpha[i,1] + beta_star   * X$male + gamma_star[i-1,(1 + M * (S_unit - 1))]) * (exp(beta_star   * (t_i + a_i)) - exp(beta_star   * a_i)) / beta_star   + d_ci$male * (X$male + t_i + a_i) * (beta_star - beta[i-1,1]) + 
                            exp(alpha[i,1] + beta[i-1,1] * X$male + gamma_star[i-1,(1 + M * (S_unit - 1))]) * (exp(beta[i-1,1] * (t_i + a_i)) - exp(beta[i-1,1] * a_i)) / beta[i-1,1])
  
  
  accept_beta <- min(1, exp(log_ef))
  beta_unif <- runif(1, 0, 1)
  beta[i,1] <- ifelse(beta_unif > accept_beta, beta[i-1,1], beta_star)
  # - Adjust sigma_beta
  sigma_beta_pr[1] <- sigma_beta_pr[1] * (1 + (i^(-gamma_Vihola)) * (accept_beta - beta_tg_accept))^0.5
  
  # - Female
  beta_star <- rtruncnorm(1, a=(10^(-5)), mean = beta[i-1,2], sd = sigma_beta_pr[2], b=5)# rnorm(n = 1, beta[i-1,2], sigma_beta_pr[2])  #
  # Acceptance step
  #start_time <- Sys.time()
  # - Get energy functions + proposals
  log_ef <- log(dtruncnorm(beta[i-1,2], mean=beta_star, sd=sigma_beta_pr[2], a=(10^(-5)), b=5)) - log(dtruncnorm(beta_star, mean=beta[i-1,2], sd=sigma_beta_pr[2], a=(10^(-5)), b=5)) #  - dnorm(beta_star, mean=beta[i-1,2], sd=sigma_beta_pr[2], log=TRUE) + dnorm(beta[i-1,2], mean=beta_star, sd=sigma_beta_pr[2], log=TRUE)# 
  
  log_ef <- log_ef + sum( - exp(alpha[i,2] + beta_star   * X$female + gamma_star[i-1,(2 + M * (S_unit - 1))]) * (exp(beta_star   * (t_i + a_i)) - exp(beta_star   * a_i)) / beta_star   + d_ci$female * (X$female + t_i + a_i) * (beta_star - beta[i-1,2]) + 
                            exp(alpha[i,2] + beta[i-1,2] * X$female + gamma_star[i-1,(2 + M * (S_unit - 1))]) * (exp(beta[i-1,2] * (t_i + a_i)) - exp(beta[i-1,2] * a_i)) / beta[i-1,2])
  
  accept_beta <- min(1, exp(log_ef))
  beta_unif <- runif(1, 0, 1)
  beta[i,2] <- ifelse(beta_unif > accept_beta, beta[i-1,2], beta_star)
  # - Adjust sigma_beta
  sigma_beta_pr[2] <- sigma_beta_pr[2] * (1 + (i^(-gamma_Vihola)) * (accept_beta - beta_tg_accept))^0.5
  
  # - Step 4: Sample gamma*
  ga_k_pos <- c(1:M)
  
  for(k in 1:K){
    log_sum_Si <- 0
    N_k[k] <- sum(ifelse(S_unit==k,1,0)) # number of people with the same k - or within the same cluster
    
    # - Sample gamma* using the Adaptive Metropolis Hastings of Vihola (2012)
    r_i <- mvrnorm(n = 1, rep(0, M), diag(1, M))
    gamma_proposed <- matrix(gamma_star[i-1,ga_k_pos],M,1) + Sigma_ch_ga[, ((k-1) * M + 1):(k * M)] %*% matrix(r_i, M,1)
    
    # - Get likelihood of the data conditional to the clustering variable S_i
    log_sum_Si <- log_sum_Si - dmvnorm(t(gamma_proposed), mean = gamma_star[i-1,ga_k_pos], sigma = Sigma_ga[, ((k-1) * M + 1):(k * M)], log = TRUE) +
                               dmvnorm(t(gamma_proposed), mean = rep(0,M), sigma = matrix(Sigma_gamma[i-1,], M,M), log = TRUE) + # this line corresponds to the prior
                               dmvnorm(gamma_star[i-1,ga_k_pos], mean = t(gamma_proposed), sigma = Sigma_ga[, ((k-1) * M + 1):(k * M)], log = TRUE) -
                               dmvnorm(gamma_star[i-1,ga_k_pos], mean = rep(0,M), sigma = matrix(Sigma_gamma[i-1,], M,M), log = TRUE) # this line corresponds to the prior
    
    log_sum_Si <- log_sum_Si + sum(ifelse(S_unit==k, (- exp(alpha[i,1] + beta[i,1] * X$male  ) * (exp(beta[i,1] * (t_i + a_i)) - exp(beta[i,1] * a_i)) / beta[i,1]) * (exp(gamma_proposed[1]) - exp(gamma_star[i-1,ga_k_pos[1]])) + 
                                                     (- exp(alpha[i,2] + beta[i,2] * X$female) * (exp(beta[i,2] * (t_i + a_i)) - exp(beta[i,2] * a_i)) / beta[i,2]) * (exp(gamma_proposed[2]) - exp(gamma_star[i-1,ga_k_pos[2]])) + d_ci$male * (gamma_proposed[1] - gamma_star[i-1,ga_k_pos[1]]) + d_ci$female * (gamma_proposed[2] - gamma_star[i-1,ga_k_pos[2]]), 0))
    
    # Acceptance step
    accept_i <- min(1, exp(log_sum_Si))
    
    u_unif <- runif(1, 0, 1)
    for(j in 1:M){
      gamma_star[i,ga_k_pos[j]] <- ifelse(u_unif > accept_i, gamma_star[i-1,ga_k_pos[j]], gamma_proposed[j])
    }
    
    # Adjsut Sigma_delta_chol
    Sigma_ga[, ((k-1) * M + 1):(k * M)] <- Sigma_ch_ga[, ((k-1) * M + 1):(k * M)] %*% (diag(1, M) + (i^(-gamma_Vihola)) * (accept_i - opt_accept) * (r_i %*% t(r_i)) / (sum(r_i^2))) %*% t(Sigma_ch_ga[, ((k-1) * M + 1):(k * M)])
    Sigma_ch_ga[, ((k-1) * M + 1):(k * M)] <- t(chol(Sigma_ga[, ((k-1) * M + 1):(k * M)]))
    
    ga_k_pos <- ga_k_pos + M
    
    # - Update zeta_star_ad (conjugate normal for the mean - add unknown base distribution, with random mean and variance. We have enough data to learn these two)
    v_zeta_ad_p <- 1/ (N_k[k] / z_ad_v[i-1] + 1 / sigma2_zeta[i-1])
    mean_zeta_ad_p <- (sum(ifelse(S_unit==k,1,0) * Z$lAD) / z_ad_v[i-1] + mu_zeta[i-1] / sigma2_zeta[i-1]) * v_zeta_ad_p  #  - beta1[i,1] * X1_cov[,1] - beta1[i,2] * X1_cov[,2] - beta1[i,3] * X1_cov[,3] - beta1[i,4] * X1_cov[,4] - beta1[i,5] * X1_cov[,5] - beta1[i,6] * X1_cov[,6] - beta1[i,7] * X1_cov[,7])
    
    z_star_ad[i,k] <- rnorm(1, mean=mean_zeta_ad_p, sd=sqrt(v_zeta_ad_p))
    
    # - Update zeta_star_mo
    z_star_mo[i,k] <- rbeta(1, shape1 = (3 + sum(Z$maleelder * ifelse(S_unit==k,1,0))), shape2 = (1 + sum((1 - Z$maleelder) * ifelse(S_unit==k,1,0))))
    
  }
  N_k_tb[i,] <- N_k
  
  # - Step 4 and 5: sample mu_gamma and Sigma_gamma
  
  N_c <- sum(ifelse(N_k[1:K]>0,1,0)) # Number of clusters with at least one element
  
  # - Posterior variance of mu_0
  #  Lambda2 <- solve(solve(lambda2) + N_c * solve(matrix(Sigma_gamma[i-1,], M,M)))
  
  # - Posterior mean of mu_0 and preparatory calculation for the posterior parameter of the Wishart matrix
  #  mu_gamma_sum_over_k <- rep(0, M)
  #  ga_k_pos <- c(1:M)
  
  #  for(k_count in 1:K){ ## - check if it performs matrix operations (unlikely)
  #    mu_gamma_sum_over_k[1] <- mu_gamma_sum_over_k[1] + ifelse(N_k[k_count] > 0, gamma_star[i,ga_k_pos[1]], 0)
  #    mu_gamma_sum_over_k[2] <- mu_gamma_sum_over_k[2] + ifelse(N_k[k_count] > 0, gamma_star[i,ga_k_pos[2]], 0)
  #    ga_k_pos <- ga_k_pos + M
  #  }
  
  #  Lambda1 <- Lambda2 %*% (solve(lambda2) %*% lambda1 + solve(matrix(Sigma_gamma[i-1,], M,M)) %*% mu_gamma_sum_over_k)
  #  mu_gamma[i,] <- rmvnorm(1, mean = Lambda1, sigma = Lambda2) # - update of mu_0
  
  # - Update of Sigma_gamma
  means_cross_prod <- matrix(0, M, M)
  ga_k_pos <- c(1:M)
  
  for(k in 1:K){
    if(N_k[k] > 0){
      means_cross_prod <- means_cross_prod + gamma_star[i,ga_k_pos] %*% t(gamma_star[i,ga_k_pos])
    }
    ga_k_pos <- ga_k_pos + M
  }
  
  Lambda3 <- N_c + lambda3
  Lambda4 <- lambda3 * lambda4 + means_cross_prod
  
  Sigma_gamma[i,] <- matrix(rinvwishart(nu = Lambda3, S = Lambda4), nrow=1, ncol=ncol(Sigma_gamma), byrow = TRUE)
  
  # - Update of the baseline distribution of zeta*star age_difference mean. This is conjugate
  ## - Posterior variance of mu_zeta
  Lambda8 <- 1/(1 + N_c / sigma2_zeta[i-1])
  
  mu_zeta_sum_over_k <- 0
  
  for(k_count in 1:K){ ## - check if it performs matrix operations (unlikely)
    mu_zeta_sum_over_k <- mu_zeta_sum_over_k + ifelse(N_k[k_count] > 0, z_star_ad[i,k], 0)
  }
  
  Lambda7 <- Lambda8 * (mu_zeta_sum_over_k / sigma2_zeta[i-1])
  mu_zeta[i] <- rnorm(1, mean = Lambda7, sd = sqrt(Lambda8)) # - update of mu_0
  
  ## - Update of the variance of the baseline distribution of zeta  
  means_cross_prod <- 0
  
  for(k in 1:K){
    if(N_k[k] > 0){
      means_cross_prod <- means_cross_prod + ((z_star_ad[i,k] - mu_zeta[i])^2)
    }
  }
  
  Lambda9 <- N_c/2 + 2
  Lambda10 <- 2 + 0.5 * means_cross_prod
  
  sigma2_zeta[i] <- 1 / rgamma(1, shape = Lambda9, rate = Lambda10)
  
  # - Update of sigma of the normal distribution of the age difference covariate (the inverse gamma)
  suff_stat <- 0.5 * sum((Z$lAD - z_star_ad[i,S_unit])^2)
  z_ad_v[i] <- 1 / rgamma(1, shape = (1 + 0.5 * sample_size), rate = (1 + suff_stat) )
  
  # - Step 9: Update concentration parameter phi
  N_c <- sum(ifelse(N_k>0,1,0)) # Number of clusters with at least one element
  N_c_tb[i] <- N_c
  
  zeta <- rbeta(1, shape1 = phi[i-1,1] + 1, shape2 = sample_size, ncp=0)
  pi_zeta <- (lambda5 + N_c - 1) / (lambda5 + N_c - 1 + sample_size * (lambda6 - log(zeta) ))
  
  mix_subs <- sample(c(1, 0), 1, prob = c(pi_zeta, 1 - pi_zeta))
  phi[i,1] <- rgamma(1, shape = (lambda5 + N_c - ifelse(mix_subs==0,0,1)), rate = lambda6 - log(zeta) )
  
  
  #     })
  #  end_time <- Sys.time()
  #  end_time - start_time
  print(paste(i, "% iteration"))
  
  if(i %in% seq(5000, n_iter, 5000)){
    save.image("~/Desktop/AVM/Canlif_MCMC.RData")
  }
  
}

plot(c(1:i), beta[1:i,1], type='l')
plot(c(1:i), beta[1:i,2], type='l')
