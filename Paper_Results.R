#============================= - Get final draws - ========================
## - 1) Discard first 40,000 iterations (burn-in) and thinning every 50
burn_in <- 80000
thinning <- 20

seq_bi_thn <- seq((burn_in+thinning), n_iter, thinning)

# - beta regression parameters
alpha_thn <- alpha[seq_bi_thn,]
beta_thn <- beta[seq_bi_thn,]

# - gamma frailty parameters
gamma_star_thn <- gamma_star[seq_bi_thn,]

# - zeta age difference
z_star_ad_thn <- z_star_ad[seq_bi_thn,]
z_ad_v_thn <- z_ad_v[seq_bi_thn]

# - zeta oldest male
z_star_mo_thn <- z_star_mo[seq_bi_thn,]

# - pi mixture weight
pi_mw_thn <- pi_mw[seq_bi_thn,]

# - mu_theta (parameter of the base distribution)
#mu_gamma_thn <- mu_gamma[seq_bi_thn,]

# - Covariance_theta
Sigma_gamma_thn <- Sigma_gamma[seq_bi_thn,]

# - psi from stick-breaking procedure
psi_sbp_thn <- psi_sbp[seq_bi_thn,]

# phi concentration parameter of the Dirichlet Process
phi_thn <- phi[seq_bi_thn,]

# - Get posterior mean of the parameters
alpha_post_m <- colMeans(alpha_thn)
beta_post_m <- colMeans(beta_thn)
gamma_star_post_m <- colMeans(gamma_star_thn)
z_star_ad_post_m <- colMeans(z_star_ad_thn)
z_star_mo_post_m <- colMeans(z_star_mo_thn)
z_ad_v_post_m <- mean(z_ad_v[seq_bi_thn])
pi_mw_post_m <- colMeans(pi_mw_thn)

gamma_star_post_m_M <- gamma_star_post_m[seq(1,M*K-1,M)]
gamma_star_post_m_F <- gamma_star_post_m[seq(2,M*K  ,M)]

colQuantiles(alpha_thn, probs=c(0.025, 0.975))
colQuantiles(beta_thn, probs=c(0.025, 0.975))

colMeans(gamma_star_thn)[c(1:8,11,12,17,18)]
colMeans(z_star_ad_thn)[c(1,2,3,4,6,9)]
colMeans(z_star_mo_thn)[c(1,2,3,4,6,9)]

round(colQuantiles(gamma_star_thn[,c(1:8,11,12,17,18)], probs=c(0.025, 0.975)), digits=2)
round(colQuantiles(z_star_ad_thn[,c(1,2,3,4,6,9)], probs=c(0.025, 0.975)), digits=2)
round(colQuantiles(z_star_mo_thn[,c(1,2,3,4,6,9)], probs=c(0.025, 0.975)), digits=2)

round(pi_mw_post_m[c(1,2,3,4,6,9)], digits=4)
round(colQuantiles(pi_mw_thn[,c(1,2,3,4,6,9)], probs=c(0.025, 0.975)), digits=4)

mean(phi_thn)
quantile(phi_thn, probs=c(0.025, 0.975))

library(matrixStats)

par(mfrow=c(1,1), mar = c(4.5,4,1,0.5))
hist(N_c_tb[seq_bi_thn])
hist(N_c_tb[seq_bi_thn], breaks=100)
hist(N_c_tb, breaks=100, main='', xlab='N. of occupied classes')

plot(c(1:(i-1)), alpha[c(1:(i-1)),1], type='l', ylim=c(-20,0))
plot(c(1:(i-1)), alpha[c(1:(i-1)),2], type='l', ylim=c(-20,0))

plot(c(1:(i-1)), beta[c(1:(i-1)),1], type='l', ylim=c(0.04,0.20))
plot(c(1:(i-1)), beta[c(1:(i-1)),2], type='l', ylim=c(0.04,0.30))

## - phi
plot(c(1:(i-1)), phi[c(1:(i-1)),1], type='l', ylim=c(0,10))

density_zA <- function(z_AD){
  f_zA_k <- matrix(NA, nrow=length(seq_bi_thn), ncol=K)
  f_za_m <- rep(0, length(seq_bi_thn))
  for(k in 1:K){
    f_zA_k[,k] <- pi_mw_thn[,k] * dnorm(z_AD, mean=z_star_ad_thn[,k], sd=z_ad_v_thn)
  }
  f_za_m <- rowSums(f_zA_k)
  return(mean(f_za_m))
}

density(Z$agediff)
plot(hist(Z$agediff, breaks=100), freq=FALSE, ylim=c(0,0.5), xlim=c(min(Z$agediff), 4))
density_zA <- sapply(log(seq(0.1, 20, 0.1)), FUN=density_zA)
lines(log(seq(0.1, 20, 0.1)), density_zA)
plot(log(seq(0.1, 20, 0.1)), density_zA, type='l')

length(Z$agediff[Z$agediff>2]) / sample_size

## - gamma* as drawn from the Dirichlet Process
plot(c(1:(i-1)), gamma_star[c(1:(i-1)),1], type='l')
plot(density(gamma_star[seq_bi_thn,1]))
plot(c(1:(i-1)), gamma_star[c(1:(i-1)),2], type='l')
plot(density(gamma_star[seq_bi_thn,2]))

plot(c(1:(i-1)), gamma_star[c(1:(i-1)),3], type='l')
plot(density(gamma_star[seq_bi_thn,3]))
plot(c(1:(i-1)), gamma_star[c(1:(i-1)),4], type='l')
plot(density(gamma_star[seq_bi_thn,4]))

plot(c(1:(i-1)), gamma_star[c(1:(i-1)),5], type='l')
plot(density(gamma_star[seq_bi_thn,5]))
plot(c(1:(i-1)), gamma_star[c(1:(i-1)),6], type='l')
plot(density(gamma_star[seq_bi_thn,6]))

plot(c(1:(i-1)), gamma_star[c(1:(i-1)),7], type='l')
plot(c(1:(i-1)), gamma_star[c(1:(i-1)),8], type='l')

plot(c(1:(i-1)), gamma_star[c(1:(i-1)),9], type='l')
plot(c(1:(i-1)), gamma_star[c(1:(i-1)),10], type='l')

plot(c(1:(i-1)), gamma_star[c(1:(i-1)),11], type='l')
plot(c(1:(i-1)), gamma_star[c(1:(i-1)),12], type='l')

plot(c(1:(i-1)), gamma_star[c(1:(i-1)),13], type='l')
plot(c(1:(i-1)), gamma_star[c(1:(i-1)),14], type='l')

plot(c(1:(i-1)), gamma_star[c(1:(i-1)),15], type='l')
plot(c(1:(i-1)), gamma_star[c(1:(i-1)),16], type='l')

plot(c(1:(i-1)), gamma_star[c(1:(i-1)),17], type='l')
plot(c(1:(i-1)), gamma_star[c(1:(i-1)),18], type='l')

plot(c(1:(i-1)), gamma_star[c(1:(i-1)),19], type='l')
plot(c(1:(i-1)), gamma_star[c(1:(i-1)),20], type='l')

plot(c(1:(i-1)), gamma_star[c(1:(i-1)),21], type='l')
plot(c(1:(i-1)), gamma_star[c(1:(i-1)),22], type='l')

plot(c(1:(i-1)), gamma_star[c(1:(i-1)),23], type='l')
plot(c(1:(i-1)), gamma_star[c(1:(i-1)),24], type='l')

plot(c(1:(i-1)), z_star_ad[c(1:(i-1)),1], type='l')
plot(c(1:(i-1)), z_star_ad[c(1:(i-1)),2], type='l')
plot(c(1:(i-1)), z_star_ad[c(1:(i-1)),3], type='l')
plot(c(1:(i-1)), z_star_ad[c(1:(i-1)),4], type='l')
plot(c(1:(i-1)), z_star_ad[c(1:(i-1)),5], type='l')
plot(c(1:(i-1)), z_star_ad[c(1:(i-1)),6], type='l')
plot(c(1:(i-1)), z_star_ad[c(1:(i-1)),7], type='l')
plot(c(1:(i-1)), z_star_ad[c(1:(i-1)),8], type='l')
plot(c(1:(i-1)), z_star_ad[c(1:(i-1)),9], type='l')
plot(c(1:(i-1)), z_star_ad[c(1:(i-1)),10], type='l')
plot(c(1:(i-1)), z_star_ad[c(1:(i-1)),11], type='l')
plot(c(1:(i-1)), z_star_ad[c(1:(i-1)),12], type='l')

plot(c(1:(i-1)), z_star_mo[c(1:(i-1)),1], type='l')
plot(c(1:(i-1)), z_star_mo[c(1:(i-1)),2], type='l')
plot(c(1:(i-1)), z_star_mo[c(1:(i-1)),3], type='l')
plot(c(1:(i-1)), z_star_mo[c(1:(i-1)),4], type='l')
plot(c(1:(i-1)), z_star_mo[c(1:(i-1)),5], type='l')
plot(c(1:(i-1)), z_star_mo[c(1:(i-1)),6], type='l')
plot(c(1:(i-1)), z_star_mo[c(1:(i-1)),7], type='l')
plot(c(1:(i-1)), z_star_mo[c(1:(i-1)),8], type='l')
plot(c(1:(i-1)), z_star_mo[c(1:(i-1)),9], type='l')
plot(c(1:(i-1)), z_star_mo[c(1:(i-1)),10], type='l')
plot(c(1:(i-1)), z_star_mo[c(1:(i-1)),11], type='l')
plot(c(1:(i-1)), z_star_mo[c(1:(i-1)),12], type='l')

log_abs_z_ad <- log(abs(Z$agediff))
hist(Z$agediff, breaks=1000)
plot(density(Z$agediff))


# - pi (mixture weights from the Dirichlet Process)
plot(c(1:(i-1)), pi_mw[c(1:(i-1)),1], type='l',ylim=c(0,1))
plot(c(1:(i-1)), pi_mw[c(1:(i-1)),2], type='l',ylim=c(0,1))
plot(c(1:(i-1)), pi_mw[c(1:(i-1)),3], type='l',ylim=c(0,1))
plot(c(1:(i-1)), pi_mw[c(1:(i-1)),4], type='l',ylim=c(0,1))
plot(c(1:(i-1)), pi_mw[c(1:(i-1)),5], type='l',ylim=c(0,1))
plot(c(1:(i-1)), pi_mw[c(1:(i-1)),6], type='l',ylim=c(0,1))
plot(c(1:(i-1)), pi_mw[c(1:(i-1)),7], type='l',ylim=c(0,1))
plot(c(1:(i-1)), pi_mw[c(1:(i-1)),8], type='l',ylim=c(0,1))
plot(c(1:(i-1)), pi_mw[c(1:(i-1)),9], type='l',ylim=c(0,1))
plot(c(1:(i-1)), pi_mw[c(1:(i-1)),10], type='l',ylim=c(0,1))
plot(c(1:(i-1)), pi_mw[c(1:(i-1)),11], type='l',ylim=c(0,1))

## - Number of clusters and cluster occupancy
plot(c(1:(i-1)), N_c_tb[c(1:(i-1))], type='l')
plot(c(1:(i-1)), N_k_tb[c(1:(i-1)),1], type='l', ylim=c(0, sample_size))
plot(c(1:(i-1)), N_k_tb[c(1:(i-1)),2], type='l', ylim=c(0, sample_size))
plot(c(1:(i-1)), N_k_tb[c(1:(i-1)),3], type='l', ylim=c(0, sample_size))
plot(c(1:(i-1)), N_k_tb[c(1:(i-1)),4], type='l', ylim=c(0, sample_size))
plot(c(1:(i-1)), N_k_tb[c(1:(i-1)),5], type='l', ylim=c(0, sample_size))
plot(c(1:(i-1)), N_k_tb[c(1:(i-1)),6], type='l', ylim=c(0, sample_size))


#============================= - log-hazard rates and comparison males/females - =======================

## - Males

log_hf_m <- function(t, ageM, mo, ad){
  
  piece1 <- log(pi_mw_post_m) - exp(beta_post_m[1] * (ageM - x_bar) + gamma_star_post_m_M) * (exp(beta_post_m[1] * t) - 1) / beta_post_m[1] + beta_post_m[1] * ((ageM - x_bar) + t) + gamma_star_post_m_M
  
  max_piece1 <- max(piece1)
  
  num_calc_lse <- max_piece1 + log(sum(exp(piece1-max_piece1)))
  
  piece2 <- log(pi_mw_post_m) - exp(beta_post_m[1] * (ageM - x_bar) + gamma_star_post_m_M) * (exp(beta_post_m[1] * t) - 1) / beta_post_m[1] 
  
  max_piece2 <- max(piece2)
  
  den_calc_lse <- max_piece2 + log(sum(exp(piece2-max_piece2)))
  
  return(num_calc_lse - den_calc_lse)
}

log_hf_f <- function(t, ageF, mo, ad){
  
  piece1 <- log(pi_mw_post_m) - exp(beta_post_m[2] * (ageF - x_bar) + gamma_star_post_m_F) * (exp(beta_post_m[2] * t) - 1) / beta_post_m[2] + beta_post_m[2] * ((ageF - x_bar) + t) + gamma_star_post_m_F
  
  max_piece1 <- max(piece1)
  
  num_calc_lse <- max_piece1 + log(sum(exp(piece1-max_piece1)))
  
  piece2 <- log(pi_mw_post_m) - exp(beta_post_m[2] * (ageF - x_bar) + gamma_star_post_m_F) * (exp(beta_post_m[2] * t) - 1) / beta_post_m[2] 
  
  max_piece2 <- max(piece2)
  
  den_calc_lse <- max_piece2 + log(sum(exp(piece2-max_piece2)))
  
  return(num_calc_lse - den_calc_lse)
}



data$ageDM <- data$ageM + data$TObs1
data$ageDF <- data$ageF + data$TObs2



#========== - log - hazard (by age and age difference) - ==============


l_hf <- function(t, ageM, ageF, j, mo, ad){
  
  # - use logsumexp trick
  pf_vector0 <- log(pi_mw_post_m) - ((exp(t * beta_post_m[1]) - 1) / beta_post_m[1]) * exp(beta_post_m[1] * (ageM - x_bar) + gamma_star_post_m[seq(1,49,2)]) - ((exp(t * beta_post_m[2]) - 1) / beta_post_m[2]) * exp(beta_post_m[2] * (ageF - x_bar) + gamma_star_post_m[seq(2,50,2)])  + (beta_post_m[j] * (ifelse(j==1,ageM,ageF) - x_bar + t) + gamma_star_post_m[seq(j,ifelse(j==1,49,50),2)]) + mo * log(z_star_mo_post_m) + log(1-z_star_mo_post_m) * (1-mo) + dnorm(log(ad), mean = z_star_ad_post_m, sd = sqrt(z_ad_v_post_m), log=TRUE)
  max_pf0 <- max(pf_vector0)
  
  num_calc_lse <- max_pf0 + log(sum(exp(pf_vector0-max_pf0)))
  
  pf_vector1 <- log(pi_mw_post_m) - ((exp(t * beta_post_m[1]) - 1) / beta_post_m[1]) * exp(beta_post_m[1] * (ageM - x_bar) + gamma_star_post_m[seq(1,49,2)]) - ((exp(t * beta_post_m[2]) - 1) / beta_post_m[2]) * exp(beta_post_m[2] * (ageF - x_bar) + gamma_star_post_m[seq(2,50,2)])                                                                                                             + mo * log(z_star_mo_post_m) + log(1-z_star_mo_post_m) * (1-mo) + dnorm(log(ad), mean = z_star_ad_post_m, sd = sqrt(z_ad_v_post_m), log=TRUE)
  max_pf1 <- max(pf_vector1)
  
  denom_calc_lse <- max_pf1 + log(sum(exp(pf_vector1-max_pf1)))
  
  return(num_calc_lse - denom_calc_lse)
  
}


### - Using the posterior mean
l_hf <- function(t, age, j, mo, ad){
  
  # - use logsumexp trick
  pf_vector0 <- log(pi_mw_post_m) - ((exp(t * beta_post_m[j]) - 1) / beta_post_m[j]) * exp(alpha_post_m[j] + beta_post_m[j] * (age - x_bar) + gamma_star_post_m[seq(j,ifelse(j==1,49,50),2)]) + (alpha_post_m[j] + beta_post_m[j] * (age - x_bar + t) + gamma_star_post_m[seq(j,ifelse(j==1,49,50),2)]) + mo * log(z_star_mo_post_m) + log(1-z_star_mo_post_m) * (1-mo) + dnorm(log(ad), mean = z_star_ad_post_m, sd = sqrt(z_ad_v_post_m), log=TRUE)
  max_pf0 <- max(pf_vector0)
  
  num_calc_lse <- max_pf0 + log(sum(exp(pf_vector0-max_pf0)))
  
  pf_vector1 <- log(pi_mw_post_m) - ((exp(t * beta_post_m[j]) - 1) / beta_post_m[j]) * exp(alpha_post_m[j] + beta_post_m[j] * (age - x_bar) + gamma_star_post_m[seq(j,ifelse(j==1,49,50),2)])                                                                                         + mo * log(z_star_mo_post_m) + log(1-z_star_mo_post_m) * (1-mo) + dnorm(log(ad), mean = z_star_ad_post_m, sd = sqrt(z_ad_v_post_m), log=TRUE)
  max_pf1 <- max(pf_vector1)
  
  denom_calc_lse <- max_pf1 + log(sum(exp(pf_vector1-max_pf1)))
  
  return(num_calc_lse - denom_calc_lse)
  
}


### - Averaged over the posterior draws
l_hf <- function(t, age, j, mo, ad){
  M <- length(seq_bi_thn)
  pf_vector0 <- pf_vector1 <- matrix(NA, M, K) #rep(0, length(seq_bi_thn))
  
  for(m in 1:M){
#    for(k in 1:K){
  # - use logsumexp trick
  pf_vector0[m,] <- log(pi_mw_thn[m,]) - ((exp(t * beta_thn[m,j]) - 1) / beta_thn[m,j]) * exp(alpha_thn[m,j] + beta_thn[m,j] * (age - x_bar) + gamma_star_thn[m,seq(j,ifelse(j==1,49,50),2)]) + (alpha_thn[m,j] + beta_thn[m,j] * (age - x_bar + t) + gamma_star_thn[m,seq(j,ifelse(j==1,49,50),2)]) + mo * log(z_star_mo_thn[m,]) + log(1-z_star_mo_thn[m,]) * (1-mo) + dnorm(log(ad), mean = z_star_ad_thn[m,], sd = sqrt(z_ad_v_thn[m]), log=TRUE)
  
  ## - denominator with logsumexp
  pf_vector1[m,] <- log(pi_mw_thn[m,]) - ((exp(t * beta_thn[m,j]) - 1) / beta_thn[m,j]) * exp(alpha_thn[m,j] + beta_thn[m,j] * (age - x_bar) + gamma_star_thn[m,seq(j,ifelse(j==1,49,50),2)])                                                                                                        + mo * log(z_star_mo_thn[m,]) + log(1-z_star_mo_thn[m,]) * (1-mo) + dnorm(log(ad), mean = z_star_ad_thn[m,], sd = sqrt(z_ad_v_thn[m]), log=TRUE)
    }  #}
  
  max_pf0 <- max(pf_vector0)
  num_calc_lse <- max_pf0 + log(sum(exp(pf_vector0-max_pf0)))
  
  max_pf1 <- max(pf_vector1)
  denom_calc_lse <- max_pf1 + log(sum(exp(pf_vector1-max_pf1)))
  
  return(num_calc_lse - denom_calc_lse)
  
}




# - male older (Yes) and ad=1
log_hf_males_mo1_ad2 <- rep(0, 40)
log_hf_females_mo1_ad2 <- rep(0, 40)
log_hf_males_mo1_ad5 <- rep(0, 40)
log_hf_females_mo1_ad5 <- rep(0, 40)
log_hf_males_mo1_ad10 <- rep(0, 40)
log_hf_females_mo1_ad10 <- rep(0, 40)

log_hf_males_mo0_ad2 <- rep(0, 40)
log_hf_females_mo0_ad2 <- rep(0, 40)
log_hf_males_mo0_ad5 <- rep(0, 40)
log_hf_females_mo0_ad5 <- rep(0, 40)
log_hf_males_mo0_ad10 <- rep(0, 40)
log_hf_females_mo0_ad10 <- rep(0, 40)


for(t in 1:40){
  log_hf_males_mo1_ad2[t] <- l_hf(t, 60, 1, 1, 2)
  log_hf_females_mo1_ad2[t] <- l_hf(t, 60, 2, 1, 2)
  log_hf_males_mo1_ad5[t] <- l_hf(t, 60, 1, 1, 5)
  log_hf_females_mo1_ad5[t] <- l_hf(t, 60, 2, 1, 5)
  log_hf_males_mo1_ad10[t] <- l_hf(t, 60, 1, 1, 10)
  log_hf_females_mo1_ad10[t] <- l_hf(t, 60, 2, 1, 10)
  
  log_hf_males_mo0_ad2[t] <- l_hf(t, 60, 1, 0, 2)
  log_hf_females_mo0_ad2[t] <- l_hf(t, 60, 2, 0, 2)
  log_hf_males_mo0_ad5[t] <- l_hf(t, 60, 1, 0, 5)
  log_hf_females_mo0_ad5[t] <- l_hf(t, 60, 2, 0, 5)
  log_hf_males_mo0_ad10[t] <- l_hf(t, 60, 1, 0, 10)
  log_hf_females_mo0_ad10[t] <- l_hf(t, 60, 2, 0, 10)
  
}


par(mfrow=c(1,2), mar = c(4.5,4,1,0.5))
plot(c(61:100), log_hf_males_mo1_ad2, type='l', ylim=c(-7,-1), xlab="Age", ylab='log-hazard', main="Male older", lwd=1.5)
lines(c(61:100), log_hf_females_mo1_ad2, lty=1, col='grey', lwd=1.5)
lines(c(61:100), log_hf_males_mo1_ad5, lty=2, lwd=1.5)
lines(c(61:100), log_hf_females_mo1_ad5, lty=2, col='grey', lwd=1.5)
lines(c(61:100), log_hf_males_mo1_ad10, lty=3, lwd=1.5)
lines(c(61:100), log_hf_females_mo1_ad10, lty=3, col='grey', lwd=1.5)
legend('topleft', legend=c("Males", "Females"), col=c('black', 'grey'), lwd=1.5, bty='n')

plot(c(61:100), log_hf_males_mo0_ad2, type='l', ylim=c(-7,-1), xlab="Age", ylab='', main="Female older", lwd=1.5)
lines(c(61:100), log_hf_females_mo0_ad2, lty=1, col='grey', lwd=1.5)
lines(c(61:100), log_hf_males_mo0_ad5, lty=2, lwd=1.5)
lines(c(61:100), log_hf_females_mo0_ad5, lty=2, col='grey', lwd=1.5)
lines(c(61:100), log_hf_males_mo0_ad10, lty=3, lwd=1.5)
lines(c(61:100), log_hf_females_mo0_ad10, lty=3, col='grey', lwd=1.5)
legend('bottomright', legend=c("Age diff.=2", "Age diff.=5", "Age diff.=10"), lwd=2, lty=c(1,2,3), bty='n')

log_hf_males_mo1_ad2 - log_hf_males_mo0_ad2
log_hf_males_mo1_ad5 - log_hf_males_mo0_ad5
log_hf_males_mo1_ad10 - log_hf_males_mo0_ad10

log_hf_females_mo1_ad2 - log_hf_females_mo0_ad2
log_hf_females_mo1_ad5 - log_hf_females_mo0_ad5
log_hf_females_mo1_ad10 - log_hf_females_mo0_ad10


# - Comparison with empirical rates
age_baseline <- floor(min(X + a_i) + x_bar) # 60
max_age <- max(ceiling(X + c(t1_i, t2_i) + a_i + x_bar)) #90


age_baseline <- 50 #round(min(X), digits=0)  + x_bar
max_age <- 99# round(max(X), digits=0) + x_bar # age_baseline + n_ages
n_ages <- max_age - age_baseline # 40

X_mo1 <- X[Z$maleelder==1,] + x_bar
t1_i_mo1 <- t1_i[Z$maleelder==1]
t2_i_mo1 <- t2_i[Z$maleelder==1]
a_i_mo1 <- data_final$trunc[Z$maleelder==1]
d_ci_mo1 <- as.data.frame(cbind(data_final$dmale[Z$maleelder==1], data_final$dfemale[Z$maleelder==1]))
colnames(d_ci_mo1) <- c("male", "female")

E_cxM_mo1 <- matrix(0, nrow(Z[Z$maleelder==1,]), n_ages) #rep(0, n_ages)
d_xM_mo1 <- matrix(0, nrow(Z[Z$maleelder==1,]), n_ages)
E_cxF_mo1 <- matrix(0, nrow(Z[Z$maleelder==1,]), n_ages)
d_xF_mo1 <- matrix(0, nrow(Z[Z$maleelder==1,]), n_ages)

X_mo0 <- X[Z$maleelder==0,] + x_bar
t1_i_mo0 <- t1_i[Z$maleelder==0]
t2_i_mo0 <- t2_i[Z$maleelder==0]
a_i_mo0 <- data_final$trunc[Z$maleelder==0]
d_ci_mo0 <- as.data.frame(cbind(data_final$dmale[Z$maleelder==0], data_final$dfemale[Z$maleelder==0]))
colnames(d_ci_mo0) <- c("male", "female")

E_cxM_mo0 <- matrix(0, nrow(Z[Z$maleelder==0,]), n_ages) #rep(0, n_ages)
d_xM_mo0 <- matrix(0, nrow(Z[Z$maleelder==0,]), n_ages)
E_cxF_mo0 <- matrix(0, nrow(Z[Z$maleelder==0,]), n_ages)
d_xF_mo0 <- matrix(0, nrow(Z[Z$maleelder==0,]), n_ages)


for(age in age_baseline:(max_age-1)){

#  E_cxM_mo1[age - age_baseline] <- 0
#  d_xM_mo1[age - age_baseline] <- 0
#  E_cxF_mo1[age - age_baseline] <- 0
#  d_xF_mo1[age - age_baseline] <- 0
  for(i in 1:nrow(Z[Z$maleelder==1,])){
    E_cxM_mo1[i,age + 1 - age_baseline] <- max(0, min(age+1, X_mo1$male[i] + a_i_mo1[i] + t1_i_mo1[i]) - max(X_mo1$male[i] + a_i_mo1[i], age))
    d_xM_mo1[i,age + 1 - age_baseline] <- ifelse(d_ci_mo1$male[i]==1 & ((X_mo1$male[i] + a_i_mo1[i] + t1_i_mo1[i]) >= age) & ((X_mo1$male[i] + a_i_mo1[i] + t1_i_mo1[i]) < age+1), 1, 0)
    E_cxF_mo1[i,age + 1 - age_baseline] <- max(0, min(age+1, X_mo1$female[i] + a_i_mo1[i] + t2_i_mo1[i]) - max(X_mo1$female[i] + a_i_mo1[i], age))
    d_xF_mo1[i,age + 1 - age_baseline] <- ifelse(d_ci_mo1$female[i]==1 & ((X_mo1$female[i] + a_i_mo1[i] + t2_i_mo1[i]) >= age) & ((X_mo1$female[i] + a_i_mo1[i] + t2_i_mo1[i]) < age+1), 1, 0)
  }
  
}

for(age in age_baseline:(max_age-1)){
#  E_cxM_mo0[age - age_baseline] <- 0
#  d_xM_mo0[age - age_baseline] <- 0
#  E_cxF_mo0[age - age_baseline] <- 0
#  d_xF_mo0[age - age_baseline] <- 0
  for(i in 1:nrow(Z[Z$maleelder==0,])){
    E_cxM_mo0[i,age - age_baseline] <- max(0, min(age+1, X_mo0$male[i] + a_i_mo0[i] + t1_i_mo0[i]) - max(X_mo0$male[i] + a_i_mo0[i], age))
    d_xM_mo0[i,age - age_baseline] <- ifelse(d_ci_mo0$male[i]==1 & ((X_mo0$male[i] + a_i_mo0[i] + t1_i_mo0[i]) >= age) & ((X_mo0$male[i] + a_i_mo0[i] + t1_i_mo0[i]) < age+1), 1, 0)
    E_cxF_mo0[i,age - age_baseline] <- max(0, min(age+1, X_mo0$female[i] + a_i_mo0[i] + t2_i_mo0[i]) - max(X_mo0$female[i] + a_i_mo0[i], age))
    d_xF_mo0[i,age - age_baseline] <- ifelse(d_ci_mo0$female[i]==1 & ((X_mo0$female[i] + a_i_mo0[i] + t2_i_mo0[i]) >= age) & ((X_mo0$female[i] + a_i_mo0[i] + t2_i_mo0[i]) < age+1), 1, 0)
  }
  
}

colSums(E_cxF_mo1)
colSums(d_xF_mo0)
sum(d_xF_mo0)

colSums(E_cxM_mo1)
colSums(d_xM_mo0)
sum(d_xM_mo0)


par(mfrow=c(1,1), mar = c(4.5,4,1,0.5))
plot(c(age_baseline:(max_age-1)), log(colSums(d_xM_mo1) / colSums(E_cxM_mo1)), pch=20, ylim=c(-6,0), xlab="Age", ylab="log-death rates", cex=1)
points(c(age_baseline:(max_age-1)), log(colSums(d_xF_mo1) / colSums(E_cxF_mo1)), pch=20, col='grey', cex=1)
points(c(age_baseline:(max_age-1)), log(colSums(d_xM_mo0) / colSums(E_cxM_mo0)), pch=4, cex=0.6)
points(c(age_baseline:(max_age-1)), log(colSums(d_xF_mo0) / colSums(E_cxF_mo0)), pch=4, col='grey', cex=0.6)
legend('bottomright', legend=c("Males", "Females"), pch=20, col=c('black', 'grey'), bty='n', cex=0.8)
legend('topleft', legend=c("Male older = 1", "Male older = 0"), pch=c(20,4), bty='n', cex=0.8)

## - plot by males-females
par(mfrow=c(1,1), mar = c(4.5,4,1,0.5))
plot(c(age_baseline:(max_age-1)), log((colSums(d_xM_mo1) + colSums(d_xM_mo0)) / (colSums(E_cxM_mo1) + colSums(E_cxM_mo0))), pch=20, ylim=c(-8,0), xlab="Age", ylab="log-hazard function", cex=1)
points(c(age_baseline:(max_age-1)), log((colSums(d_xF_mo1) + colSums(d_xF_mo0)) / (colSums(E_cxF_mo1) + colSums(E_cxF_mo0))), pch=20, col='grey', cex=1)
legend('bottomright', legend=c("Males", "Females"), pch=20, col=c('black', 'grey'), bty='n', cex=0.8)
legend('topleft', legend=c("AVDPM", "Base Gompertz (BG)", "Prop. Hazard (PH)"), lty=c(1,2,3), bty='n', lwd=2)
### - add Gompertz
lines(c(age_baseline:(max_age-1)), alpha_BG_m + beta_BG_m * (c(age_baseline:(max_age-1)) - x_bar), lwd=2, lty=2)
lines(c(age_baseline:(max_age-1)), alpha_BG_f + beta_BG_f * (c(age_baseline:(max_age-1)) - x_bar), lwd=2, col='grey', lty=2)
lines(c(age_baseline:(max_age-1)), alpha_G_m + beta_G_m * (c(age_baseline:(max_age-1)) - x_bar), lwd=2, lty=3)
lines(c(age_baseline:(max_age-1)), alpha_G_f + beta_G_f * (c(age_baseline:(max_age-1)) - x_bar), lwd=2, col='grey', lty=3)

ln_haz_AVDPM_m <- ln_haz_AVDPM_f <- rep(0, (max_age - age_baseline))
for(i in 1:(max_age - age_baseline)){
  ln_haz_AVDPM_m[i] <- l_hf(i, age_baseline, 1, 0, 1)
  ln_haz_AVDPM_f[i] <- l_hf(i, age_baseline, 2, 0, 1)
}

lines(c(age_baseline:(max_age-1)), ln_haz_AVDPM_m, lwd=2)
lines(c(age_baseline:(max_age-1)), ln_haz_AVDPM_f, col='grey', lwd=2)


# - Density plot (hence produce an 4x4 plot)
n_samples <- 10000
row_sample <- sample(seq_bi_thn, n_samples, replace=TRUE)
pi_sample <- pi_mw[row_sample,]
alphaM_sample <- alpha[row_sample,1]
alphaF_sample <- alpha[row_sample,2]
betaM_sample <- beta[row_sample,1]
betaF_sample <- beta[row_sample,2]

gammaM_sample <- gamma_star[row_sample, seq(1, 49, 2)]
gammaF_sample <- gamma_star[row_sample, seq(2, 50, 2)]

z_star_ad_sample <- z_star_ad[row_sample,]
z_star_mo_sample <- z_star_mo[row_sample,]
z_ad_v_sample <- z_ad_v[row_sample]

s_sample <- rep(0, n_samples)
gammaM_sample2 <- rep(0, n_samples) #matrix(NA, nrow=n_samples, col=K)
gammaF_sample2 <- rep(0, n_samples) #matrix(NA, nrow=n_samples, col=K)


t_density <- function(age, Z_ad, Z_mo){
  y_sample <- matrix(NA, nrow=n_samples, ncol=2)
  
  for(sample in 1:n_samples){
    
    prob_vector0 <- pi_sample[sample,] * dnorm(log(Z_ad), mean=z_star_ad_sample[sample,], sd=sqrt(z_ad_v_sample[sample])) * (z_star_mo_sample[sample,]^Z_mo) * ((1 - z_star_mo_sample[sample,])^(1 - Z_mo))
    prob_vector1 <- prob_vector0 / sum(prob_vector0)
    
    s_sample[sample] <- sample(c(1:K), 1, prob=pi_sample[sample,])
    
    gammaM_sample2[sample] <- gammaM_sample[sample, s_sample[sample]]
    gammaF_sample2[sample] <- gammaF_sample[sample, s_sample[sample]]
    
    y_sample[sample,1] <- log(1 + betaM_sample[sample] * rexp(1, rate=exp(alphaM_sample[sample] + gammaM_sample2[sample] + betaM_sample[sample] * (age-x_bar)))) / betaM_sample[sample]   # - rnorm(1, mean=theta1_bar[s_sample[sample]], sd=sqrt(sigma2_c_bar[1]))
    y_sample[sample,2] <- log(1 + betaF_sample[sample] * rexp(1, rate=exp(alphaF_sample[sample] + gammaF_sample2[sample] + betaF_sample[sample] * (age-x_bar)))) / betaF_sample[sample]   # rnorm(1, mean=theta2_bar[s_sample[sample]], sd=sqrt(sigma2_c_bar[2]))
  }
  return(y_sample)
}

t_ad2_mo1 <- t_density(60, 2, 1)
t_ad5_mo1 <- t_density(60, 5, 1)
t_ad10_mo1 <- t_density(60, 10, 1)
t_ad2_mo0 <- t_density(60, 2, 0)
t_ad5_mo0 <- t_density(60, 5, 0)
t_ad10_mo0 <- t_density(60, 10, 0)

colMeans(t_ad2_mo1)
colMeans(t_ad5_mo1)
colMeans(t_ad10_mo1)
colMeans(t_ad2_mo0)
colMeans(t_ad5_mo0)
colMeans(t_ad10_mo0)

plot(density(t_ad2_mo1[,1]), xlab="t", ylab='f(t)', main="Male older=1", ylim=c(0,0.05), lwd=0.5)  
lines(density(t_ad5_mo1[,1]), lty=2, lwd=0.5)
lines(density(t_ad10_mo1[,1]), lty=3, lwd=0.5)
lines(density(t_ad2_mo1[,2]), lty=1, col='grey', lwd=0.5)
lines(density(t_ad5_mo1[,2]), lty=2, col='grey', lwd=0.5)
lines(density(t_ad10_mo1[,2]), lty=3, col='grey', lwd=0.5)

plot(density(t_ad2_mo0[,1]), xlab="t", ylab='f(t)', main="Male older=0", ylim=c(0,0.05), lwd=0.5)  
lines(density(t_ad5_mo0[,1]), lty=2, lwd=0.5)
lines(density(t_ad10_mo0[,1]), lty=3, lwd=0.5)
lines(density(t_ad2_mo0[,2]), lty=1, col='grey', lwd=0.5)
lines(density(t_ad5_mo0[,2]), lty=2, col='grey', lwd=0.5)
lines(density(t_ad10_mo0[,2]), lty=3, col='grey', lwd=0.5)

#======================== - Effect of the covariates on gamma - =========================

Cov_effect_male <- matrix(NA, 2, 10)
Cov_effect_female <- matrix(NA, 2, 10)

# - MO = 1
for(log_AD in -6:3){
  prob_vector0 <- pi_mw_post_m * dnorm(log_AD, mean=z_star_ad_post_m, sd=sqrt(z_ad_v_post_m)) * z_star_mo_post_m
  prob_vector1 <- prob_vector0 / sum(prob_vector0)
  
  Cov_effect_male[1,log_AD+7] <- sum(prob_vector1 * gamma_star_post_m_M)
  Cov_effect_female[1,log_AD+7] <- sum(prob_vector1 * gamma_star_post_m_F)
}

for(log_AD in -6:3){
  prob_vector0 <- pi_mw_post_m * dnorm(log_AD, mean=z_star_ad_post_m, sd=sqrt(z_ad_v_post_m)) * (1 - z_star_mo_post_m)
  prob_vector1 <- prob_vector0 / sum(prob_vector0)
  
  Cov_effect_male[2,log_AD+7] <- sum(prob_vector1 * gamma_star_post_m_M)
  Cov_effect_female[2,log_AD+7] <- sum(prob_vector1 * gamma_star_post_m_F)
  
}

par(mfrow=c(1,1), mar = c(4.5,4,1,0.5))
plot(exp(-6:3), Cov_effect_male[1,], type='l', ylim=c(-2.5,3), xlab='Age difference', ylab='Average gamma', lwd=2)
lines(exp(-6:3), Cov_effect_male[2,], lty=2, lwd=2)
lines(exp(-6:3), Cov_effect_female[1,], lty=1, col='grey', lwd=2)
lines(exp(-6:3), Cov_effect_female[2,], lty=2, col='grey', lwd=2)
legend('topright', legend=c("Males", "Females"), col=c('black', 'grey'), lwd=2, bty='n', cex=0.75)
legend('bottomright', legend=c("Male Older = 1", "Male Older = 0"), lty=c(1,2), lwd=2, bty='n', cex=0.75)


xtable(Cov_effect[1:4,], digits = rep(3, 21))

# - Effect of the age difference
AD_effect <- matrix(NA, 2, 10)

for(log_AD in -6:3){
  prob_vector0 <- pi_mw_post_m * dnorm(log_AD, mean=z_star_ad_post_m, sd=sqrt(z_ad_v_post_m))
  prob_vector1 <- prob_vector0 / sum(prob_vector0)
  
  AD_effect[1,log_AD+7] <- sum(prob_vector1 * gamma_star_post_m_M)
  AD_effect[2,log_AD+7] <- sum(prob_vector1 * gamma_star_post_m_F)
}

plot(exp(-6:3), AD_effect[1,], type='l', ylim=c(-2,2.5), xlab='Age difference', ylab='Average gamma', lwd=2)
lines(exp(-6:3), AD_effect[2,], col='grey', lwd=2)

# - Effect of the male older
MO_effect <- matrix(NA, 2, 2)

for(mo in 0:1){
  prob_vector0 <- pi_mw_post_m * (z_star_mo_post_m^(mo)) * ((1 - z_star_mo_post_m)^(1 - mo))
  prob_vector1 <- prob_vector0 / sum(prob_vector0)
  
  MO_effect[1,mo+1] <- sum(prob_vector1 * gamma_star_post_m_M)
  MO_effect[2,mo+1] <- sum(prob_vector1 * gamma_star_post_m_F)
}

colnames(MO_effect) <- c("MO=0", "MO=1")
rownames(MO_effect) <- c("Male", "Females")

MO_effect[,1] - MO_effect[,2]

#======================== - Analyis of dependent time to event - =========================

# - Comment: eliminating the age-dependence (i.e. if we random draw it) then we note an increase in the correlation between
# - the time to events in males and females

n_samples <- 10000


# -  This function uses the posterior mean
t_dependent <- function(n_samples, age_base, Z_ad, Z_mo){#function(age_base, Z_ad, Z_mo){
  s_sample <- rep(0, n_samples)
  t_dep_sample <- matrix(NA, n_samples, 2)
  
  # - Step 2: choose the conditional probability distribution
  prob_vector0 <- pi_mw_post_m * dnorm(log(Z_ad), mean=z_star_ad_post_m, sd=sqrt(z_ad_v_post_m)) * (z_star_mo_post_m^Z_mo) * ((1 - z_star_mo_post_m)^(1 - Z_mo))
  prob_vector1 <- prob_vector0 / sum(prob_vector0)
  
  for(sample in 1:n_samples){
    s_sample[sample] <- sample(c(1:K), 1, prob=prob_vector1)
    age <- runif(1, age_base, 100)
    #    t_dep_sample[sample,1] <- log(1 + beta_post_m[1] * rexp(1, rate=exp(alpha_post_m[1] + gamma_star_post_m_M[s_sample[sample]] + beta_post_m[1] * (age +  Z_ad *      Z_mo  )  ))) / beta_post_m[1]   # - rnorm(1, mean=theta1_bar[s_sample[sample]], sd=sqrt(sigma2_c_bar[1]))
    #    t_dep_sample[sample,2] <- log(1 + beta_post_m[2] * rexp(1, rate=exp(alpha_post_m[2] + gamma_star_post_m_F[s_sample[sample]] + beta_post_m[2] * (age - (Z_ad * (1 - Z_mo)))  ))) / beta_post_m[2]   # rnorm(1, mean=theta2_bar[s_sample[sample]], sd=sqrt(sigma2_c_bar[2]))
    t_dep_sample[sample,1] <- log(1 + beta_post_m[1] * rexp(1, rate=exp(alpha_post_m[1] + gamma_star_post_m_M[s_sample[sample]] + beta_post_m[1] * (age - x_bar +  Z_ad *      Z_mo  )  ))) / beta_post_m[1]   # - rnorm(1, mean=theta1_bar[s_sample[sample]], sd=sqrt(sigma2_c_bar[1]))
    t_dep_sample[sample,2] <- log(1 + beta_post_m[2] * rexp(1, rate=exp(alpha_post_m[2] + gamma_star_post_m_F[s_sample[sample]] + beta_post_m[2] * (age - x_bar - (Z_ad * (1 - Z_mo)))  ))) / beta_post_m[2]   # rnorm(1, mean=theta2_bar[s_sample[sample]], sd=sqrt(sigma2_c_bar[2]))
  }
  return(t_dep_sample)
}

t_dependent <- function(n_samples, age_base, Z_ad, Z_mo){#function(age_base, Z_ad, Z_mo){
  s_sample <- rep(0, n_samples)
  t_dep_sample <- matrix(NA, n_samples, 2)
  
  # - Step 2: choose the conditional probability distribution
  prob_vector0 <- pi_mw_post_m * dnorm(log(Z_ad), mean=z_star_ad_post_m, sd=sqrt(z_ad_v_post_m)) * (z_star_mo_post_m^Z_mo) * ((1 - z_star_mo_post_m)^(1 - Z_mo))
  prob_vector1 <- prob_vector0 / sum(prob_vector0)
  
  for(sample in 1:n_samples){
    s_sample[sample] <- sample(c(1:K), 1, prob=prob_vector1)
    age <- runif(1, age_base-20, age_base+20)
    #    t_dep_sample[sample,1] <- log(1 + beta_post_m[1] * rexp(1, rate=exp(alpha_post_m[1] + gamma_star_post_m_M[s_sample[sample]] + beta_post_m[1] * (age +  Z_ad *      Z_mo  )  ))) / beta_post_m[1]   # - rnorm(1, mean=theta1_bar[s_sample[sample]], sd=sqrt(sigma2_c_bar[1]))
    #    t_dep_sample[sample,2] <- log(1 + beta_post_m[2] * rexp(1, rate=exp(alpha_post_m[2] + gamma_star_post_m_F[s_sample[sample]] + beta_post_m[2] * (age - (Z_ad * (1 - Z_mo)))  ))) / beta_post_m[2]   # rnorm(1, mean=theta2_bar[s_sample[sample]], sd=sqrt(sigma2_c_bar[2]))
    t_dep_sample[sample,1] <- log(1 + beta_post_m[1] * rexp(1, rate=exp(alpha_post_m[1] + gamma_star_post_m_M[s_sample[sample]] + beta_post_m[1] * (age - x_bar +  Z_ad *      Z_mo  )  ))) / beta_post_m[1]   # - rnorm(1, mean=theta1_bar[s_sample[sample]], sd=sqrt(sigma2_c_bar[1]))
    t_dep_sample[sample,2] <- log(1 + beta_post_m[2] * rexp(1, rate=exp(alpha_post_m[2] + gamma_star_post_m_F[s_sample[sample]] + beta_post_m[2] * (age - x_bar - (Z_ad * (1 - Z_mo)))  ))) / beta_post_m[2]   # rnorm(1, mean=theta2_bar[s_sample[sample]], sd=sqrt(sigma2_c_bar[2]))
  }
  return(t_dep_sample)
}


# -  This function uses 1 sample from the posterior distribution
t_dependent <- function(n_samples, age_base, Z_ad, Z_mo){#function(age_base, Z_ad, Z_mo){
  s_sample <- rep(0, n_samples)
  t_dep_sample <- matrix(NA, n_samples, 2)
  
  sample_thn <- sample(c(1:nrow(alpha_thn)), 1, replace=FALSE)
  # - Step 2: choose the conditional probability distribution
  prob_vector0 <- pi_mw_thn[sample_thn,] * dnorm(log(Z_ad), mean=z_star_ad_thn[sample_thn,], sd=sqrt(z_ad_v_thn[sample_thn])) * (z_star_mo_thn[sample_thn,]^Z_mo) * ((1 - z_star_mo_thn[sample_thn,])^(1 - Z_mo))
  prob_vector1 <- prob_vector0 / sum(prob_vector0)
  
  for(sample in 1:n_samples){
    s_sample[sample] <- sample(c(1:K), 1, prob=prob_vector1)
    age <- runif(1, age_base, 100)
    #    t_dep_sample[sample,1] <- log(1 + beta_post_m[1] * rexp(1, rate=exp(alpha_post_m[1] + gamma_star_post_m_M[s_sample[sample]] + beta_post_m[1] * (age +  Z_ad *      Z_mo  )  ))) / beta_post_m[1]   # - rnorm(1, mean=theta1_bar[s_sample[sample]], sd=sqrt(sigma2_c_bar[1]))
    #    t_dep_sample[sample,2] <- log(1 + beta_post_m[2] * rexp(1, rate=exp(alpha_post_m[2] + gamma_star_post_m_F[s_sample[sample]] + beta_post_m[2] * (age - (Z_ad * (1 - Z_mo)))  ))) / beta_post_m[2]   # rnorm(1, mean=theta2_bar[s_sample[sample]], sd=sqrt(sigma2_c_bar[2]))
    t_dep_sample[sample,1] <- log(1 + beta_thn[sample_thn,1] * rexp(1, rate=exp(alpha_thn[sample_thn,1] + gamma_star_thn[sample_thn,s_sample[sample]*2-1] + beta_thn[sample_thn,1] * (age - x_bar +  Z_ad *      Z_mo  )  ))) / beta_thn[sample_thn,1]   # - rnorm(1, mean=theta1_bar[s_sample[sample]], sd=sqrt(sigma2_c_bar[1]))
    t_dep_sample[sample,2] <- log(1 + beta_thn[sample_thn,2] * rexp(1, rate=exp(alpha_thn[sample_thn,2] + gamma_star_thn[sample_thn,s_sample[sample]*2  ] + beta_thn[sample_thn,2] * (age - x_bar - (Z_ad * (1 - Z_mo)))  ))) / beta_thn[sample_thn,2]   # rnorm(1, mean=theta2_bar[s_sample[sample]], sd=sqrt(sigma2_c_bar[2]))
  }
  return(t_dep_sample)
}

# -  This function the average t sampled using all the samples
t_dependent <- function(n_samples, age_base, Z_ad, Z_mo){#function(age_base, Z_ad, Z_mo){
  s_sample <- rep(0, n_samples)
  t_dep_sample <- matrix(NA, n_samples, 2)
  
  sample_thn <- sample(c(1:nrow(alpha_thn)), 1, replace=FALSE)
  # - Step 2: choose the conditional probability distribution
  prob_vector0 <- pi_mw_thn[sample_thn,] * dnorm(log(Z_ad), mean=z_star_ad_thn[sample_thn,], sd=sqrt(z_ad_v_thn[sample_thn])) * (z_star_mo_thn[sample_thn,]^Z_mo) * ((1 - z_star_mo_thn[sample_thn,])^(1 - Z_mo))
  prob_vector1 <- prob_vector0 / sum(prob_vector0)
  
  for(sample in 1:n_samples){
    s_sample[sample] <- sample(c(1:K), 1, prob=prob_vector1)
    age <- runif(1, age_base-5, age_base+5)
    #    t_dep_sample[sample,1] <- log(1 + beta_post_m[1] * rexp(1, rate=exp(alpha_post_m[1] + gamma_star_post_m_M[s_sample[sample]] + beta_post_m[1] * (age +  Z_ad *      Z_mo  )  ))) / beta_post_m[1]   # - rnorm(1, mean=theta1_bar[s_sample[sample]], sd=sqrt(sigma2_c_bar[1]))
    #    t_dep_sample[sample,2] <- log(1 + beta_post_m[2] * rexp(1, rate=exp(alpha_post_m[2] + gamma_star_post_m_F[s_sample[sample]] + beta_post_m[2] * (age - (Z_ad * (1 - Z_mo)))  ))) / beta_post_m[2]   # rnorm(1, mean=theta2_bar[s_sample[sample]], sd=sqrt(sigma2_c_bar[2]))
    t_dep_sample[sample,1] <- mean(log(1 + beta_thn[,1] * rexp(1, rate=exp(alpha_thn[,1] + gamma_star_thn[,s_sample[sample]*2-1] + beta_thn[,1] * (age - x_bar +  Z_ad *      Z_mo  )  ))) / beta_thn[,1])   # - rnorm(1, mean=theta1_bar[s_sample[sample]], sd=sqrt(sigma2_c_bar[1]))
    t_dep_sample[sample,2] <- mean(log(1 + beta_thn[,2] * rexp(1, rate=exp(alpha_thn[,2] + gamma_star_thn[,s_sample[sample]*2  ] + beta_thn[,2] * (age - x_bar - (Z_ad * (1 - Z_mo)))  ))) / beta_thn[,2])   # rnorm(1, mean=theta2_bar[s_sample[sample]], sd=sqrt(sigma2_c_bar[2]))
  }
  return(t_dep_sample)
}




# - table
sequence <-  seq(1,20,1) # c(seq(0.1,1,0.1), seq(1.5, 20, 0.5))# c(0.5, 20, 0.5)
spearman_tb60 <- matrix(NA, 2, length(sequence))
kendall_tb60 <- matrix(NA, 2, length(sequence))

spearman_tb70 <- matrix(NA, 2, length(sequence))
kendall_tb70 <- matrix(NA, 2, length(sequence))

spearman_tb80 <- matrix(NA, 2, length(sequence))
kendall_tb80 <- matrix(NA, 2, length(sequence))

t_dependent <- function(n_samples, age_base, range, Z_ad, Z_mo){#function(age_base, Z_ad, Z_mo){
  s_sample <- rep(0, n_samples)
  t_dep_sample <- matrix(NA, n_samples, 2)
  
  # - Step 2: choose the conditional probability distribution
  
  for(sample in 1:n_samples){
    m <- sample(c(1:length(seq_bi_thn)), 1)
    
    prob_vector0 <- pi_mw_thn[m,] * dnorm(log(Z_ad), mean=z_star_ad_thn[m,], sd=sqrt(z_ad_v_thn[m])) * (z_star_mo_thn[m,]^Z_mo) * ((1 - z_star_mo_thn[m,])^(1 - Z_mo))
    prob_vector1 <- prob_vector0 / sum(prob_vector0)
    s_sample[sample] <- sample(c(1:K), 1, prob=prob_vector1)
    
    age <- runif(1, age_base-range, age_base+range)
    t_dep_sample[sample,1] <- log(1 + beta_thn[m,1] * rexp(1, rate=exp(alpha_thn[m,1] + gamma_star_thn[m,2 * s_sample[sample]-1] + beta_thn[m,1] * (age - x_bar +  Z_ad *      Z_mo  )  ))) / beta_thn[m,1]   # - rnorm(1, mean=theta1_bar[s_sample[sample]], sd=sqrt(sigma2_c_bar[1]))
    t_dep_sample[sample,2] <- log(1 + beta_thn[m,2] * rexp(1, rate=exp(alpha_thn[m,2] + gamma_star_thn[m,2 * s_sample[sample]  ] + beta_thn[m,2] * (age - x_bar - (Z_ad * (1 - Z_mo)))  ))) / beta_thn[m,2]   # rnorm(1, mean=theta2_bar[s_sample[sample]], sd=sqrt(sigma2_c_bar[2]))
  }
  return(t_dep_sample)
}



ab <- 60
n_samples <- 20000
col<- 1
for(ad in sequence){
  row <- 1
  for(mo in c(1,0)){
#    t_dep60 <- t_dependent(n_samples, 60, ad, mo)
#    t_dep70 <- t_dependent(n_samples, 70, ad, mo)
    t_dep60 <- t_dependent(n_samples, 60, 20, ad, mo)
    t_dep70 <- t_dependent(n_samples, 60, 5, ad, mo)
    #    t_dep80 <- t_dependent(n_samples, 80, 5, ad, mo)
    
    spearman_tb60[row,col] <- cor.test(x=t_dep60[,1], y=t_dep60[,2], method = 'spearman')$estimate
    kendall_tb60[row,col] <-  cor.test(x=t_dep60[,1], y=t_dep60[,2], method = 'kendall')$estimate
    
    spearman_tb70[row,col] <- cor.test(x=t_dep70[,1], y=t_dep70[,2], method = 'spearman')$estimate
    kendall_tb70[row,col] <-  cor.test(x=t_dep70[,1], y=t_dep70[,2], method = 'kendall')$estimate
  
#    spearman_tb80[row,col] <- cor.test(x=t_dep80[,1], y=t_dep80[,2], method = 'spearman')$estimate
#    kendall_tb80[row,col] <-  cor.test(x=t_dep80[,1], y=t_dep80[,2], method = 'kendall')$estimate

    row <- row + 1
  }
  col <- col + 1
}


par(mfrow=c(1,2), mar = c(4.5,4,1,0.5))
plot(sequence, spearman_tb60[1,], type='l', ylim=c(0,1), xlab="Age difference", ylab='Spearman rho', lwd=2)
lines(sequence, spearman_tb60[2,], col='grey', lwd=2)
lines(sequence, spearman_tb70[1,], lwd=2, lty=2)
lines(sequence, spearman_tb70[2,], col='grey', lwd=2, lty=2)
#lines(sequence, spearman_tb80[1,], lwd=2, lty=3)
#lines(sequence, spearman_tb80[2,], col='grey', lwd=2, lty=3)
legend('topleft', legend=c("MO=1", "MO=0"), col=c('black', 'grey'), lwd=2, bty='n', cex=0.75)
legend('topright', legend=c("40-80", "55-65"), lty=c(1,2), lwd=2, bty='n', cex=0.75)

plot(sequence, kendall_tb60[1,], type='l', ylim=c(0,1), xlab="Age difference", ylab='kendall tau', lwd=2)
lines(sequence, kendall_tb60[2,], col='grey', lwd=2)
lines(sequence, kendall_tb70[1,], lwd=2, lty=2)
lines(sequence, kendall_tb70[2,], col='grey', lwd=2, lty=2)
#lines(sequence, kendall_tb80[1,], lwd=2, lty=3)
#lines(sequence, kendall_tb80[2,], col='grey', lwd=2, lty=3)
legend('topleft', legend=c("MO=1", "MO=0"), col=c('black', 'grey'), lwd=2, bty='n', cex=0.75)
legend('topright', legend=c("40-80", "55-65"), lty=c(1,2), lwd=2, bty='n', cex=0.75)

# - Below there's the older part about the dependence analysis (using AD=2,5,10)

ab <- 60
t_dep_ad2_mo1 <- t_dependent(1000, ab, 2, 1)
t_dep_ad5_mo1 <- t_dependent(1000, ab, 5, 1)
t_dep_ad10_mo1 <- t_dependent(1000, ab, 10, 1)
t_dep_ad2_mo0 <- t_dependent(1000, ab, 2, 0)
t_dep_ad5_mo0 <- t_dependent(1000, ab, 5, 0)
t_dep_ad10_mo0 <- t_dependent(1000, ab, 10, 0)

ab <- 80
t_dep_ad2_mo1 <- t_dependent(1000, ab, 2, 1)
t_dep_ad5_mo1 <- t_dependent(1000, ab, 5, 1)
t_dep_ad10_mo1 <- t_dependent(1000, ab, 10, 1)
t_dep_ad2_mo0 <- t_dependent(1000, ab, 2, 0)
t_dep_ad5_mo0 <- t_dependent(1000, ab, 5, 0)
t_dep_ad10_mo0 <- t_dependent(1000, ab, 10, 0)

par(mfrow=c(1,2), mar = c(4.5,4,1,0.5))
plot(log(t_dep_ad2_mo1[,1]),log(t_dep_ad2_mo1[,2]), pch=20, cex=0.1, xlab='log(T) males', ylab='log(T) females', main='Male Older', xlim=c(-10,4), ylim=c(-8,4))
points(log(t_dep_ad5_mo1[,1]),log(t_dep_ad5_mo1[,2]), col='gray', pch=20, cex=0.1)
points(log(t_dep_ad10_mo1[,1]),log(t_dep_ad10_mo1[,2]), col='azure2', pch=20, cex=0.1)
legend('bottomleft', legend=c("AD=2", "AD=5", 'AD=10'), col=c('black', 'grey', 'azure2'), pch=16, bty='n', cex=0.75)

plot(log(t_dep_ad2_mo0[,1]),log(t_dep_ad2_mo0[,2]), pch=20, cex=0.1, xlab='log(T) males', ylab='log(T) females', main='Female Older', xlim=c(-10,4), ylim=c(-8,4))
points(log(t_dep_ad5_mo0[,1]),log(t_dep_ad5_mo0[,2]), col='gray', pch=20, cex=0.1)
points(log(t_dep_ad10_mo0[,1]),log(t_dep_ad10_mo0[,2]), col='azure2', pch=20, cex=0.1)
legend('bottomleft', legend=c("AD=2", "AD=5", 'AD=10'), col=c('black', 'grey', 'azure2'), pch=16, bty='n', cex=0.75)

par(mfrow=c(1,2), mar = c(4.5,4,1,0.5))
plot(t_dep_ad2_mo1[,1],t_dep_ad2_mo1[,2], pch=20, cex=0.1, xlab='T males', ylab='T females', main='Male Older')#, xlim=c(-10,4), ylim=c(-8,4))
points(t_dep_ad5_mo1[,1],t_dep_ad5_mo1[,2], col='gray', pch=20, cex=0.1)
points(t_dep_ad10_mo1[,1],t_dep_ad10_mo1[,2], col='azure2', pch=20, cex=0.1)
legend('topright', legend=c("AD=2", "AD=5", 'AD=10'), col=c('black', 'grey', 'azure2'), pch=16, bty='n', cex=0.75)

plot(t_dep_ad2_mo0[,1],t_dep_ad2_mo0[,2], pch=20, cex=0.1, xlab='log(T) males', ylab='T females', main='Female Older')#, xlim=c(-10,4), ylim=c(-8,4))
points(t_dep_ad5_mo0[,1],t_dep_ad5_mo0[,2], col='gray', pch=20, cex=0.1)
points(t_dep_ad10_mo0[,1],t_dep_ad10_mo0[,2], col='azure2', pch=20, cex=0.1)
legend('topright', legend=c("AD=2", "AD=5", 'AD=10'), col=c('black', 'grey', 'azure2'), pch=16, bty='n', cex=0.75)


plot(density(t_dep_ad2_mo1[,1]), xlab='T males', main='Male Older', ylim=c(0,1))#, xlim=c(-10,4), ylim=c(-8,4))
lines(density(t_dep_ad5_mo1[,1]), xlab='T males', main='Male Older')#, xlim=c(-10,4), ylim=c(-8,4))
lines(density(t_dep_ad10_mo1[,1]), xlab='T males', main='Male Older')#, xlim=c(-10,4), ylim=c(-8,4))
points(t_dep_ad5_mo1[,1],t_dep_ad5_mo1[,2], col='gray', pch=20, cex=0.1)
points(t_dep_ad10_mo1[,1],t_dep_ad10_mo1[,2], col='azure2', pch=20, cex=0.1)
legend('topright', legend=c("AD=2", "AD=5", 'AD=10'), col=c('black', 'grey', 'azure2'), pch=16, bty='n', cex=0.75)




cor(log(t_dep_ad2_mo1[,1]), log(t_dep_ad2_mo1[,2]))
cor(log(t_dep_ad5_mo1[,1]), log(t_dep_ad5_mo1[,2]))
cor(log(t_dep_ad10_mo1[,1]), log(t_dep_ad10_mo1[,2]))

cor(t_dep_ad2_mo0[,1], y=t_dep_ad2_mo0[,2])
cor(t_dep_ad5_mo0[,1], y=t_dep_ad5_mo0[,2])
cor(t_dep_ad10_mo0[,1], y=t_dep_ad10_mo0[,2])

rhosphea <- cor.test(x=t_dep_ad2_mo1[,1], y=t_dep_ad2_mo1[,2], method = 'spearman')$estimate
cor.test(x=t_dep_ad5_mo1[,1], y=t_dep_ad5_mo1[,2], method = 'spearman')
cor.test(x=t_dep_ad10_mo1[,1], y=t_dep_ad10_mo1[,2], method = 'spearman')

cor.test(x=t_dep_ad2_mo0[,1], y=t_dep_ad2_mo0[,2], method = 'spearman')
cor.test(x=t_dep_ad5_mo0[,1], y=t_dep_ad5_mo0[,2], method = 'spearman')
cor.test(x=t_dep_ad10_mo0[,1], y=t_dep_ad10_mo0[,2], method = 'spearman')

cor.test(x=t_dep_ad2_mo1[,1], y=t_dep_ad2_mo1[,2], method = 'kendall')
cor.test(x=t_dep_ad5_mo1[,1], y=t_dep_ad5_mo1[,2], method = 'kendall')
cor.test(x=t_dep_ad10_mo1[,1], y=t_dep_ad10_mo1[,2], method = 'kendall')

cor.test(x=t_dep_ad2_mo0[,1], y=t_dep_ad2_mo0[,2], method = 'kendall')
cor.test(x=t_dep_ad5_mo0[,1], y=t_dep_ad5_mo0[,2], method = 'kendall')
cor.test(x=t_dep_ad10_mo0[,1], y=t_dep_ad10_mo0[,2], method = 'kendall')

cor.test(x=log(t_dep_ad2_mo1[,1]), y=log(t_dep_ad2_mo1[,2]), method = 'kendall')
cor.test(x=log(t_dep_ad5_mo1[,1]), y=log(t_dep_ad5_mo1[,2]), method = 'kendall')
cor.test(x=log(t_dep_ad10_mo1[,1]), y=log(t_dep_ad10_mo1[,2]), method = 'kendall')

cor.test(x=log(t_dep_ad2_mo0[,1]), y=log(t_dep_ad2_mo0[,2]), method = 'kendall')
cor.test(x=log(t_dep_ad5_mo0[,1]), y=log(t_dep_ad5_mo0[,2]), method = 'kendall')
cor.test(x=log(t_dep_ad10_mo0[,1]), y=log(t_dep_ad10_mo0[,2]), method = 'kendall')


Correlation_table_60 <- matrix(NA, 2, 10)
#rho_table_60 <- matrix(NA, 2, 10)
#tau_table_60 <- matrix(NA, 2, 10)
age_base1 <- 60
for(Z_mo in 0:1){
  for(Z_ad in 1:10){
    Correlation_table_60[Z_mo+1,Z_ad] <- cor(log(t_dependent(age_base1, Z_ad, Z_mo)))[1,2]
  }
}

Correlation_table_80 <- matrix(NA, 2, 10)
age_base1 <- 80
for(Z_mo in 0:1){
  for(Z_ad in 1:10){
    Correlation_table_80[Z_mo+1,Z_ad] <- cor(log(t_dependent(age_base1, Z_ad, Z_mo)))[1,2]
  }
}

plot(c(1:10), Correlation_table_60[1,], type='l', lwd=1.5, ylab = "Correlation log(T)", xlab='Age Difference', ylim=c(0,0.5))
lines(c(1:10), Correlation_table_60[2,], type='l', col='grey', lwd=1.5)
lines(c(1:10), Correlation_table_80[1,], type='l', col='black', lwd=1.5, lty=2)
lines(c(1:10), Correlation_table_80[2,], type='l', col='grey', lwd=1.5, lty=2)

# - Comment: since the analysis of correlated time to event to not provide any meaningful result, it is the case
# - to analyse the correlated gamma, given the value of zeta. We have still to figure out whether we should consider its
# - predictive distribution in terms of the MCMC output (given the poor mixing), or just the posterior mean

n_samples <- 10000

gamma_star_post_m_M <- gamma_star_post_m[seq(1,M*K-1,M)]
gamma_star_post_m_F <- gamma_star_post_m[seq(2,M*K  ,M)]

gamma_dependent <- function(Z_ad, Z_mo){
  s_sample <- rep(0, n_samples)
  gamma_dep_sample <- matrix(NA, n_samples, 2)
  #  age_base_cor <- age_base - x_bar
  
  # - Step 2: choose the conditional probability distribution
  prob_vector0 <- pi_mw_post_m * dnorm(log(Z_ad), mean=z_star_ad_post_m, sd=sqrt(z_ad_v_post_m)) * (z_star_mo_post_m^Z_mo) * ((1 - z_star_mo_post_m)^(1 - Z_mo))
  prob_vector1 <- prob_vector0 / sum(prob_vector0)
  
  for(sample in 1:n_samples){
    s_sample[sample] <- sample(c(1:K), 1, prob=prob_vector1)
    gamma_dep_sample[sample,1] <- gamma_star_post_m_M[s_sample[sample]]
    gamma_dep_sample[sample,2] <- gamma_star_post_m_F[s_sample[sample]]
  }
  return(gamma_dep_sample)
}

ga_dep_ad2_mo1 <- gamma_dependent(2, 1)
ga_dep_ad5_mo1 <- gamma_dependent(5, 1)
ga_dep_ad10_mo1 <- gamma_dependent(10, 1)
ga_dep_ad2_mo0 <- gamma_dependent(2, 0)
ga_dep_ad5_mo0 <- gamma_dependent(5, 0)
ga_dep_ad10_mo0 <- gamma_dependent(10, 0)

plot(ga_dep_ad2_mo1[,1], ga_dep_ad2_mo1[,2], pch=16, cex=0.1)
points(ga_dep_ad5_mo1[,1], ga_dep_ad5_mo1[,2], col='gray', pch=16, cex=0.1)
points(ga_dep_ad10_mo1[,1], ga_dep_ad10_mo1[,2], col='azure2', pch=16, cex=0.1)


Correlation_table_gamma <- matrix(NA, 2, 10)
rho_table_60 <- matrix(NA, 2, 10)
tau_table_60 <- matrix(NA, 2, 10)
age_base <- 60
for(Z_mo in 0:1){
  for(Z_ad in 1:10){
    Correlation_table_gamma[Z_mo+1,Z_ad] <- cor(gamma_dependent(Z_ad, Z_mo))[1,2]
  }
}

plot(c(1:10), Correlation_table_gamma[1,], type='l', lwd=1.5, ylab = "Correlation log(T)", xlab='Age Difference', ylim=c(-1,1))
lines(c(1:10), Correlation_table_gamma[2,], type='l', col='grey', lwd=1.5)


#===================== - Group analysis & Analysis of clusters by Bayes' rule - ==========================

N_c_tb_thn <- mean(N_c_tb[seq_bi_thn])
N_k_tb_thn <- colMeans(N_k_tb[seq_bi_thn,])
round(100 * N_k_tb_thn / sample_size,2)

sum(N_k_tb_thn[c(1,2,3,4,6,9)]) / sample_size# sum(N_k_tb_thn[c(2,3,4,6,8)]) / sample_size

par(mfrow=c(1,1), mar = c(4.5,4,1,0.5))
barplot(N_k_tb_thn, main="", horiz=FALSE, xlab = "Mixture component (k)", ylab="Post. mean n_k", cex.names=0.7, names.arg=c(1:K))
abline(h=round(0.025 * sample_size,0))


post_prob_fc <- function(unit){
  
  # - use log sum exp trick to deal with very small numbers which R inappropriately set to zero
  pf_vector <- rep(0, K)
  th_k_pos <- c(1:M) # position vector for theta
  
  pf_vector <- log(pi_mw_post_m) - exp(alpha_post_m[1] + beta_post_m[1] * X$male[unit]   + gamma_star_post_m[seq(1, M*K-1, M)]) * (exp(beta_post_m[1] * (t1_i[unit] + a_i[unit])) - exp(beta_post_m[1] * a_i[unit])) / beta_post_m[1] + d_ci$male[unit]   * gamma_star_post_m[seq(1, M*K-1, M)] - 
                                   exp(alpha_post_m[2] + beta_post_m[2] * X$female[unit] + gamma_star_post_m[seq(2, M*K  , M)]) * (exp(beta_post_m[2] * (t2_i[unit] + a_i[unit])) - exp(beta_post_m[2] * a_i[unit])) / beta_post_m[2] + d_ci$female[unit] * gamma_star_post_m[seq(2, M*K  , M)] + 
           dnorm(Z$agediff[unit], mean=z_star_ad_post_m, sd=sqrt(z_ad_v_post_m), log = TRUE) + Z$maleelder[unit] * log(z_star_mo_post_m) + (1 - Z$maleelder[unit]) * log(1 - z_star_mo_post_m)
  
  max_bayes <- max(pf_vector)
  bayes_alloc <- which(pf_vector==max_bayes)
  
  return(bayes_alloc)
}


post_prob_fc <- function(unit){
  
  # - use log sum exp trick to deal with very small numbers which R inappropriately set to zero
  M <- length(seq_bi_thn)
  pf_vector <- matrix(0, M, K) # rep(0, K)

  for(m in 1:M){
  pf_vector[m,] <- log(pi_mw_thn[m,]) - exp(alpha_thn[m,1] + beta_thn[m,1] * X$male[unit]   + gamma_star_thn[m,seq(1, 2*K-1, 2)]) * (exp(beta_thn[m,1] * (t1_i[unit] + a_i[unit])) - exp(beta_thn[m,1] * a_i[unit])) / beta_thn[m,1] + d_ci$male[unit]   * gamma_star_thn[m,seq(1, 2*K-1, 2)] - 
                                        exp(alpha_thn[m,2] + beta_thn[m,2] * X$female[unit] + gamma_star_thn[m,seq(2, 2*K  , 2)]) * (exp(beta_thn[m,2] * (t2_i[unit] + a_i[unit])) - exp(beta_thn[m,2] * a_i[unit])) / beta_thn[m,2] + d_ci$female[unit] * gamma_star_thn[m,seq(2, 2*K  , 2)] + 
    dnorm(Z$agediff[unit], mean=z_star_ad_thn[m,], sd=sqrt(z_ad_v_thn[m]), log = TRUE) + Z$maleelder[unit] * log(z_star_mo_thn[m,]) + (1 - Z$maleelder[unit]) * log(1 - z_star_mo_thn[m,])
  }
  pf_vector_m <- colMeans(pf_vector)
  max_bayes <- max(pf_vector_m)
  bayes_alloc <- which(pf_vector_m==max_bayes)
  
  return(bayes_alloc)
}


#Bayes_allocation <- sapply(c(1:sample_size), FUN = post_prob_fc)

Bayes_allocation <- rep(0, sample_size)
for(unit in 1:sample_size){
  Bayes_allocation[unit] <- post_prob_fc(unit)
  print(unit)
}


Group_analysis <- cbind(X, Z, d_ci, Bayes_allocation)
colnames(Group_analysis) <- c("AgeM", "AgeF", "AD", "MO", "DeathM", "DeathF", "BayesGroup")

class <- c(100 * nrow(Group_analysis[Group_analysis$BayesGroup==1,]) / sample_size,
           100 * nrow(Group_analysis[Group_analysis$BayesGroup==2,]) / sample_size,
           100 * nrow(Group_analysis[Group_analysis$BayesGroup==3,]) / sample_size,
           100 * nrow(Group_analysis[Group_analysis$BayesGroup==4,]) / sample_size,
           100 * nrow(Group_analysis[Group_analysis$BayesGroup==5,]) / sample_size,
           100 * nrow(Group_analysis[Group_analysis$BayesGroup==6,]) / sample_size,
           100 * nrow(Group_analysis[Group_analysis$BayesGroup==7,]) / sample_size,
           100 * nrow(Group_analysis[Group_analysis$BayesGroup==8,]) / sample_size,
           100 * nrow(Group_analysis[Group_analysis$BayesGroup==9,]) / sample_size,
           100 * nrow(Group_analysis[Group_analysis$BayesGroup==10,]) / sample_size,
           100 * nrow(Group_analysis[Group_analysis$BayesGroup==11,]) / sample_size,
           100 * nrow(Group_analysis[Group_analysis$BayesGroup==12,]) / sample_size,
           100 * nrow(Group_analysis[Group_analysis$BayesGroup==13,]) / sample_size,
           100 * nrow(Group_analysis[Group_analysis$BayesGroup==14,]) / sample_size,
           100 * nrow(Group_analysis[Group_analysis$BayesGroup==15,]) / sample_size,
           100 * nrow(Group_analysis[Group_analysis$BayesGroup==16,]) / sample_size,
           100 * nrow(Group_analysis[Group_analysis$BayesGroup==17,]) / sample_size,
           100 * nrow(Group_analysis[Group_analysis$BayesGroup==18,]) / sample_size,
           100 * nrow(Group_analysis[Group_analysis$BayesGroup==19,]) / sample_size,
           100 * nrow(Group_analysis[Group_analysis$BayesGroup==20,]) / sample_size,
           100 * nrow(Group_analysis[Group_analysis$BayesGroup==21,]) / sample_size,
           100 * nrow(Group_analysis[Group_analysis$BayesGroup==22,]) / sample_size,
           100 * nrow(Group_analysis[Group_analysis$BayesGroup==23,]) / sample_size,
           100 * nrow(Group_analysis[Group_analysis$BayesGroup==24,]) / sample_size,
           100 * nrow(Group_analysis[Group_analysis$BayesGroup==25,]) / sample_size)


sum(class[c(1,2,4,9)])

## - Group selection Analysis
group_select <- 1
Gr_select <- Group_analysis[Group_analysis$BayesGroup==group_select,]

100*nrow(Gr_select) / nrow(Group_analysis)

mean(Gr_select$AgeM) + x_bar
mean(Gr_select$AgeF) + x_bar
mean(Gr_select$AD)
mean(log(abs(Gr_select$AgeM-Gr_select$AgeF)))
mean(Gr_select$MO)

c(gamma_star_post_m_M[group_select], gamma_star_post_m_F[group_select])
c(z_star_ad_post_m[group_select], z_star_mo_post_m[group_select])

plot(gamma_star_thn[,seq(1, M*K-1,M)[group_select]], gamma_star_thn[,seq(2, M*K,M)[group_select]], pch=20, cex=0.1)
cor(gamma_star_thn[,seq(1, M*K-1,M)[group_select]], gamma_star_thn[,seq(2, M*K,M)[group_select]]) # - Careful, because this is a matter of posterior parameter correlation

group_select <- 2
Gr_select <- Group_analysis[Group_analysis$BayesGroup==group_select,]

100*nrow(Gr_select) / nrow(Group_analysis)

mean(Gr_select$AgeM) + x_bar
mean(Gr_select$AgeF) + x_bar
mean(Gr_select$AD)
mean(log(abs(Gr_select$AgeM-Gr_select$AgeF)))
mean(Gr_select$MO)

c(gamma_star_post_m_M[group_select], gamma_star_post_m_F[group_select])
c(z_star_ad_post_m[group_select], z_star_mo_post_m[group_select])


group_select <- 4
Gr_select <- Group_analysis[Group_analysis$BayesGroup==group_select,]

100*nrow(Gr_select) / nrow(Group_analysis)

mean(Gr_select$AgeM) + x_bar
mean(Gr_select$AgeF) + x_bar
mean(Gr_select$AD)
mean(log(abs(Gr_select$AgeM-Gr_select$AgeF)))
mean(Gr_select$MO)

c(gamma_star_post_m_M[group_select], gamma_star_post_m_F[group_select])
c(z_star_ad_post_m[group_select], z_star_mo_post_m[group_select])


group_select <- 9
Gr_select <- Group_analysis[Group_analysis$BayesGroup==group_select,]

100*nrow(Gr_select) / nrow(Group_analysis)

mean(Gr_select$AgeM) + x_bar
mean(Gr_select$AgeF) + x_bar
mean(Gr_select$AD)
mean(log(abs(Gr_select$AgeM-Gr_select$AgeF)))
mean(Gr_select$MO)

c(gamma_star_post_m_M[group_select], gamma_star_post_m_F[group_select])
c(z_star_ad_post_m[group_select], z_star_mo_post_m[group_select])


# - In the training sample
mean(Group_analysis$AgeM)
mean(Group_analysis$AgeF)
mean(Group_analysis$AD)
mean(Group_analysis$MO)

#=============================== - Model comparison - =========================

# - Model fitting
log_L_m <- function(vdParameters){
  alpha <- vdParameters[1]
  beta <- vdParameters[2]
  
  logL <- -exp(alpha + beta * X$male) * (exp(beta * (t1_i + a_i)) - exp(beta * a_i)) / beta + d_ci$male * (alpha + beta * (X$male + t1_i + a_i)) 
  
  return(-sum(logL))
}

par_init <- c(-10, 0.1)
par_est_male_s <- nlm(log_L_m, p=par_init, typsize=par_init, hessian=T, iterlim=10000)
sqrt(diag(solve(par_est_male_s$hessian)))

### - female
log_L_f <- function(vdParameters){
  alpha <- vdParameters[1]
  beta <- vdParameters[2]
  
  logL <- -exp(alpha + beta * X$female) * (exp(beta * (t2_i + a_i)) - exp(beta * a_i)) / beta + d_ci$female * (alpha + beta * (X$female + t2_i + a_i)) 
  
  return(-sum(logL))
}

par_init <- c(-10, 0.1)
par_est_female_s <- nlm(log_L_f, p=par_init, typsize=par_init, hessian=T, iterlim=10000)

sqrt(diag(solve(par_est_female_s$hessian)))

log_L_m_cov <- function(vdParameters){
  alpha <- vdParameters[1]
  beta <- vdParameters[2]
  gamma_AD <- vdParameters[3]
  gamma_MO <- vdParameters[4]
  
  logL <- -exp(alpha + beta * X$male + gamma_AD * Z$agediff + gamma_MO * Z$maleelder) * (exp(beta * (t1_i + a_i)) - exp(beta * a_i)) / beta + d_ci$male * (alpha + beta * (X$male + t1_i + a_i) + gamma_AD * Z$agediff + gamma_MO * Z$maleelder)
  
  return(-sum(logL))
}

par_init <- c(-10, 0.1, 0.1, 0.1)
par_est_male_s_cov <- nlm(log_L_m_cov, p=par_init, typsize=par_init, hessian=T, iterlim=10000)
par_est_male_s_cov$estimate / (0.5 * sqrt(diag(solve(par_est_male_s_cov$hessian))))

log_L_f_cov <- function(vdParameters){
  alpha <- vdParameters[1]
  beta <- vdParameters[2]
  gamma_AD <- vdParameters[3]
  gamma_MO <- vdParameters[4]
  
  logL <- -exp(alpha + beta * X$female + gamma_AD * Z$agediff + gamma_MO * Z$maleelder) * (exp(beta * (t2_i + a_i)) - exp(beta * a_i)) / beta + d_ci$female * (alpha + beta * (X$female + t2_i + a_i) + gamma_AD * Z$agediff + gamma_MO * Z$maleelder)
  
  return(-sum(logL))
}

par_init <- c(-10, 0.1, 0.1, 0.1)
par_est_female_s_cov <- nlm(log_L_f_cov, p=par_init, typsize=par_init, hessian=T, iterlim=10000)
par_est_female_s_cov$estimate / (sqrt(diag(solve(par_est_female_s_cov$hessian))))



# - in sample performance

## - AIC base Gompertz
alpha_BG_m <- par_est_male_s$estimate[1]
beta_BG_m <- par_est_male_s$estimate[2]

alpha_BG_f <- par_est_female_s$estimate[1]
beta_BG_f <- par_est_female_s$estimate[2]

log_L_BG_is <- sum( - exp(alpha_BG_m + beta_BG_m * X$male  ) * (exp(beta_BG_m * (t1_i + a_i)) - exp(beta_BG_m * a_i)) / beta_BG_m + d_ci$male   * (alpha_BG_m + beta_BG_m * (X$male   + a_i + t1_i)) -
                      exp(alpha_BG_f + beta_BG_f * X$female) * (exp(beta_BG_f * (t2_i + a_i)) - exp(beta_BG_f * a_i)) / beta_BG_f + d_ci$female * (alpha_BG_f + beta_BG_f * (X$female + a_i + t2_i)))

AIC_BG_is <- - 2 * log_L_BG_is + 2 * 4


## - AIC Gompertz independent
par_init <- c(-10, 0.1, 0.1, 0.1)
par_est_male_s_cov <- nlm(log_L_m_cov, p=par_init, typsize=par_init, hessian=T, iterlim=10000)
par_est_female_s_cov <- nlm(log_L_f_cov, p=par_init, typsize=par_init, hessian=T, iterlim=10000)

alpha_G_m <- par_est_male_s_cov$estimate[1]
beta_G_m <- par_est_male_s_cov$estimate[2]
gamma_AD_G_m <- par_est_male_s_cov$estimate[3]
gamma_MO_G_m <- par_est_male_s_cov$estimate[4]

alpha_G_f <- par_est_female_s_cov$estimate[1]
beta_G_f <- par_est_female_s_cov$estimate[2]
gamma_AD_G_f <- par_est_female_s_cov$estimate[3]
gamma_MO_G_f <- par_est_female_s_cov$estimate[4]

log_L_G_is <- sum( - exp(alpha_G_m + beta_G_m * X$male   + gamma_AD_G_m * Z$agediff + gamma_MO_G_m * Z$maleelder) * (exp(beta_G_m * (t1_i + a_i)) - exp(beta_G_m * a_i)) / beta_G_m + d_ci$male   * (alpha_G_m + beta_G_m * (X$male   + a_i + t1_i) + gamma_AD_G_m * Z$agediff + gamma_MO_G_m * Z$maleelder) -
                     exp(alpha_G_f + beta_G_f * X$female + gamma_AD_G_f * Z$agediff + gamma_MO_G_f * Z$maleelder) * (exp(beta_G_f * (t2_i + a_i)) - exp(beta_G_f * a_i)) / beta_G_f + d_ci$female * (alpha_G_f + beta_G_f * (X$female + a_i + t2_i) + gamma_AD_G_f * Z$agediff + gamma_MO_G_f * Z$maleelder))

AIC_G_is <- - 2 * log_L_G_is + 2 * 8

### - WAIC computation
WAIC_p1_i <- function(unit){
  
  M_theta <- length(seq_bi_thn)
  
  sum_over_m <- 0
  
  for(m in 1:M_theta){
    pf_vector_num <- log(pi_mw_thn[m,]) - exp(alpha_thn[m,1] + beta_thn[m,1] * X$male[unit]   + gamma_star_thn[m,seq(1,M*K-1,M)]) * (exp(beta_thn[m,1] * (t1_i[unit] + a_i[unit])) - exp(beta_thn[m,1] * a_i[unit])) / beta_thn[m,1] + d_ci$male[unit]   * (alpha_thn[m,1] + beta_thn[m,1] * (t1_i[unit] + X$male[unit]   + a_i[unit]) + gamma_star_thn[m,seq(1,M*K-1,M)]) - 
                                          exp(alpha_thn[m,2] + beta_thn[m,2] * X$female[unit] + gamma_star_thn[m,seq(2,M*K  ,M)]) * (exp(beta_thn[m,2] * (t2_i[unit] + a_i[unit])) - exp(beta_thn[m,2] * a_i[unit])) / beta_thn[m,2] + d_ci$female[unit] * (alpha_thn[m,2] + beta_thn[m,2] * (t2_i[unit] + X$female[unit] + a_i[unit]) + gamma_star_thn[m,seq(2,M*K  ,M)]) + 
      dnorm(Z$agediff[unit], mean=z_star_ad_thn[m,], sd=sqrt(z_ad_v_thn[m]), log = TRUE) + Z$maleelder[unit] * log(z_star_mo_thn[m,]) + (1 - Z$maleelder[unit]) * log(1 - z_star_mo_thn[m,])
    
    max_pf_n <- max(pf_vector_num)
    
    pf_vector_den <- log(pi_mw_thn[m,]) + dnorm(Z$agediff[unit], mean=z_star_ad_thn[m,], sd=sqrt(z_ad_v_thn[m]), log = TRUE) + Z$maleelder[unit] * log(z_star_mo_thn[m,]) + (1 - Z$maleelder[unit]) * log(1 - z_star_mo_thn[m,])
    
    max_pf_d <- max(pf_vector_den)
    
    f_MF <- exp(max_pf_n + log(sum(exp(pf_vector_num - max_pf_n))) - max_pf_d - log(sum(exp(pf_vector_den - max_pf_d))))
    
    sum_over_m <- sum_over_m + f_MF
    
  }
  
  return(-log(M_theta) + log(sum_over_m))
  
}

WAIC_p2_i <- function(unit){
  
  M_theta <- length(seq_bi_thn)
  
  sum_over_m <- 0
  
  for(m in 1:M_theta){
    pf_vector_num <- log(pi_mw_thn[m,]) - exp(alpha_thn[m,1] + beta_thn[m,1] * X$male[unit]   + gamma_star_thn[m,seq(1,M*K-1,M)]) * (exp(beta_thn[m,1] * (t1_i[unit] + a_i[unit])) - exp(beta_thn[m,1] * a_i[unit])) / beta_thn[m,1] + d_ci$male[unit]   * (alpha_thn[m,1] + beta_thn[m,1] * (t1_i[unit] + X$male[unit]   + a_i[unit]) + gamma_star_thn[m,seq(1,M*K-1,M)]) - 
                                          exp(alpha_thn[m,2] + beta_thn[m,2] * X$female[unit] + gamma_star_thn[m,seq(2,M*K  ,M)]) * (exp(beta_thn[m,2] * (t2_i[unit] + a_i[unit])) - exp(beta_thn[m,2] * a_i[unit])) / beta_thn[m,2] + d_ci$female[unit] * (alpha_thn[m,2] + beta_thn[m,2] * (t2_i[unit] + X$female[unit] + a_i[unit]) + gamma_star_thn[m,seq(2,M*K  ,M)]) + 
      dnorm(Z$agediff[unit], mean=z_star_ad_thn[m,], sd=sqrt(z_ad_v_thn[m]), log = TRUE) + Z$maleelder[unit] * log(z_star_mo_thn[m,]) + (1 - Z$maleelder[unit]) * log(1 - z_star_mo_thn[m,])
    
    max_pf_n <- max(pf_vector_num)
    
    pf_vector_den <- log(pi_mw_thn[m,]) + dnorm(Z$agediff[unit], mean=z_star_ad_thn[m,], sd=sqrt(z_ad_v_thn[m]), log = TRUE) + Z$maleelder[unit] * log(z_star_mo_thn[m,]) + (1 - Z$maleelder[unit]) * log(1 - z_star_mo_thn[m,])
    
    max_pf_d <- max(pf_vector_den)
    
    log_f_MF <- max_pf_n + log(sum(exp(pf_vector_num - max_pf_n))) - max_pf_d - log(sum(exp(pf_vector_den - max_pf_d)))
    
    sum_over_m <- sum_over_m + log_f_MF
    
  }
  
  return(sum_over_m/M_theta)
  
}

WAIC_p1_vector <- sapply(c(1:sample_size), WAIC_p1_i)
WAIC_p2_vector <- sapply(c(1:sample_size), WAIC_p2_i)

p_WAIC <- 2 * (sum(WAIC_p1_vector) - sum(WAIC_p2_vector))
WAIC <- 2*(- sum(WAIC_p1_vector) + p_WAIC)


# - out of sample

## - Held out data set
data_test <- data_shuffle[(round(nrow(data_shuffle)*0.75,0)+1):nrow(data_shuffle),]

X <- as.data.frame(cbind(data_test$ageM, data_test$ageF)) - x_bar
#X <-  as.data.frame(cbind(data$ageM, data$ageF))
colnames(X) <- c("male", "female")
Z <- cbind(log(as.numeric(data_test$agediff)), data_test$maleelder)
colnames(Z) <- c("agediff", "maleelder")
Z <- as.data.frame(Z)
#Z$lAD <- log(abs(Z$agediff))

t1_i <- data_test$TObs1
t2_i <- data_test$TObs2
a_i <- data_test$trunc
d_ci <- as.data.frame(cbind(data_test$dmale, data_test$dfemale))
colnames(d_ci) <- c("male", "female")


sample_size_test <- nrow(data_test)

## - Base Gompertz
log_L_BG_os <- sum( - exp(alpha_BG_m + beta_BG_m * X$male  ) * (exp(beta_BG_m * (t1_i + a_i)) - exp(beta_BG_m * a_i)) / beta_BG_m + d_ci$male   * (alpha_BG_m + beta_BG_m * (X$male   + a_i + t1_i)) -
                      exp(alpha_BG_f + beta_BG_f * X$female) * (exp(beta_BG_f * (t2_i + a_i)) - exp(beta_BG_f * a_i)) / beta_BG_f + d_ci$female * (alpha_BG_f + beta_BG_f * (X$female + a_i + t2_i)))

log_L_BG_m_os <- sum( - exp(alpha_BG_m + beta_BG_m * X$male  ) * (exp(beta_BG_m * (t1_i + a_i)) - exp(beta_BG_m * a_i)) / beta_BG_m + d_ci$male   * (alpha_BG_m + beta_BG_m * (X$male   + a_i + t1_i)))

log_L_BG_f_os <- sum( - exp(alpha_BG_f + beta_BG_f * X$female) * (exp(beta_BG_f * (t2_i + a_i)) - exp(beta_BG_f * a_i)) / beta_BG_f + d_ci$female * (alpha_BG_f + beta_BG_f * (X$female + a_i + t2_i)))

AIC_BG_os <- - 2 * log_L_BG_os + 2 * 4
AIC_BG_m_os <- - 2 * log_L_BG_m_os + 2 * 2
AIC_BG_f_os <- - 2 * log_L_BG_f_os + 2 * 2

## - Proportional hazard
log_L_G_os <- sum( - exp(alpha_G_m + beta_G_m * X$male   + gamma_AD_G_m * Z$agediff + gamma_MO_G_m * Z$maleelder) * (exp(beta_G_m * (t1_i + a_i)) - exp(beta_G_m * a_i)) / beta_G_m + d_ci$male   * (alpha_G_m + beta_G_m * (X$male   + a_i + t1_i) + gamma_AD_G_m * Z$agediff + gamma_MO_G_m * Z$maleelder) -
                     exp(alpha_G_f + beta_G_f * X$female + gamma_AD_G_f * Z$agediff + gamma_MO_G_f * Z$maleelder) * (exp(beta_G_f * (t2_i + a_i)) - exp(beta_G_f * a_i)) / beta_G_f + d_ci$female * (alpha_G_f + beta_G_f * (X$female + a_i + t2_i) + gamma_AD_G_f * Z$agediff + gamma_MO_G_f * Z$maleelder))

log_L_G_m_os <- sum( - exp(alpha_G_m + beta_G_m * X$male   + gamma_AD_G_m * Z$agediff + gamma_MO_G_m * Z$maleelder) * (exp(beta_G_m * (t1_i + a_i)) - exp(beta_G_m * a_i)) / beta_G_m + d_ci$male   * (alpha_G_m + beta_G_m * (X$male   + a_i + t1_i) + gamma_AD_G_m * Z$agediff + gamma_MO_G_m * Z$maleelder))

log_L_G_f_os <- sum( - exp(alpha_G_f + beta_G_f * X$female + gamma_AD_G_f * Z$agediff + gamma_MO_G_f * Z$maleelder) * (exp(beta_G_f * (t2_i + a_i)) - exp(beta_G_f * a_i)) / beta_G_f + d_ci$female * (alpha_G_f + beta_G_f * (X$female + a_i + t2_i) + gamma_AD_G_f * Z$agediff + gamma_MO_G_f * Z$maleelder))

AIC_G_os <- - 2 * log_L_G_os + 2 * 8
AIC_G_m_os <- - 2 * log_L_G_m_os + 2 * 4
AIC_G_f_os <- - 2 * log_L_G_f_os + 2 * 4

## - AVM
WAIC_p1_vector_os <- sapply(c(1:sample_size_test), WAIC_p1_i)
WAIC_p2_vector_os <- sapply(c(1:sample_size_test), WAIC_p2_i)

p_WAIC_os <- 2 * (sum(WAIC_p1_vector_os) - sum(WAIC_p2_vector_os))
WAIC_os <- 2*(- sum(WAIC_p1_vector_os) + p_WAIC_os)

#================================== - Actuarial analysis - =========================

ann_last_indG <- function(iota, x_male_ann, ad, mo){
  
  integrand <- function(t){
    #      exp(-iota * t) * (1 - (1 - exp(-((exp(betaM_ann_iG * t) - 1)/betaM_ann_iG) * exp(alphaM_ann_iG + betaM_ann_iG * x_male_ann + deltaM_ad_ann_iG * log(ad) + deltaM_mo_ann_iG * mo))) * (1 - exp(-((exp(betaF_ann_iG * t) - 1)/betaF_ann_iG) * exp(alphaF_ann_iG + betaF_ann_iG * (x_male_ann - ad * mo + ad * (1 - mo)) + deltaF_ad_ann_iG * log(ad) + deltaF_mo_ann_iG * mo))))
    exp(-iota * t) * (1 - (1 - exp(-((exp(betaM_ann_iG * t) - 1)/betaM_ann_iG) * exp(alphaM_ann_iG + betaM_ann_iG * (x_male_ann - x_bar - (1 - mo) * ad)))) * (1 - exp(-((exp(betaF_ann_iG * t) - 1)/betaF_ann_iG) * exp(alphaF_ann_iG + betaF_ann_iG * (x_male_ann - x_bar - ad * mo)))))
  }
  
  #    annuity = integrate(integrand, lower=0, upper=(120-min(x_male_ann, (x_male_ann - ad * mo + ad * (1 - mo)))))$value
  annuity = integrate(integrand, lower=0, upper=100)$value
  
  return(annuity)
  
}

ann_joint_indG <- function(iota, x_male_ann, ad, mo){
  
  integrand <- function(t){
    #      exp(-iota * t - ((exp(betaM_ann_iG * t) - 1)/betaM_ann_iG) * exp(alphaM_ann_iG + betaM_ann_iG * x_male_ann + deltaM_ad_ann_iG * log(ad) + deltaM_mo_ann_iG * mo) - ((exp(betaF_ann_iG * t) - 1)/betaF_ann_iG) * exp(alphaF_ann_iG + betaF_ann_iG * (x_male_ann - ad * mo + ad * (1 - mo)) + deltaF_ad_ann_iG * log(ad) + deltaF_mo_ann_iG * mo)) 
    exp(-iota * t - ((exp(betaM_ann_iG * t) - 1)/betaM_ann_iG) * exp(alphaM_ann_iG + betaM_ann_iG * (x_male_ann - x_bar - (1 - mo) * ad)) - ((exp(betaF_ann_iG * t) - 1)/betaF_ann_iG) * exp(alphaF_ann_iG + betaF_ann_iG * (x_male_ann - x_bar - ad * mo))) 
  }
  
  annuity = integrate(integrand, lower=0, upper=100)$value
  
  return(annuity)
}



# - Annuity factor for the last survivor
ann_last_indPH <- function(iota, x_male_ann, ad, mo){
  
  integrand <- function(t){
    #      exp(-iota * t) * (1 - (1 - exp(-((exp(betaM_ann_iPH * t) - 1)/betaM_ann_iPH) * exp(alphaM_ann_iPH + betaM_ann_iPH * x_male_ann + deltaM_ad_ann_iPH * log(ad) + deltaM_mo_ann_iPH * mo))) * (1 - exp(-((exp(betaF_ann_iPH * t) - 1)/betaF_ann_iPH) * exp(alphaF_ann_iPH + betaF_ann_iPH * (x_male_ann - ad * mo + ad * (1 - mo)) + deltaF_ad_ann_iPH * log(ad) + deltaF_mo_ann_iPH * mo))))
    exp(-iota * t) * (1 - (1 - exp(-((exp(betaM_ann_iPH * t) - 1)/betaM_ann_iPH) * exp(alphaM_ann_iPH + betaM_ann_iPH * (x_male_ann - x_bar - (1 - mo) * ad) + deltaM_ad_ann_iPH * log(ad) + deltaM_mo_ann_iPH * mo))) * (1 - exp(-((exp(betaF_ann_iPH * t) - 1)/betaF_ann_iPH) * exp(alphaF_ann_iPH + betaF_ann_iPH * (x_male_ann - x_bar - ad * mo) + deltaF_ad_ann_iPH * log(ad) + deltaF_mo_ann_iPH * mo))))
  }
  
  #    annuity = integrate(integrand, lower=0, upper=(120-min(x_male_ann, (x_male_ann - ad * mo + ad * (1 - mo)))))$value
  annuity = integrate(integrand, lower=0, upper=100)$value
  
  return(annuity)
  
}

ann_joint_indPH <- function(iota, x_male_ann, ad, mo){
  
  integrand <- function(t){
    #      exp(-iota * t - ((exp(betaM_ann_iPH * t) - 1)/betaM_ann_iPH) * exp(alphaM_ann_iPH + betaM_ann_iPH * x_male_ann + deltaM_ad_ann_iPH * log(ad) + deltaM_mo_ann_iPH * mo) - ((exp(betaF_ann_iPH * t) - 1)/betaF_ann_iPH) * exp(alphaF_ann_iPH + betaF_ann_iPH * (x_male_ann - ad * mo + ad * (1 - mo)) + deltaF_ad_ann_iPH * log(ad) + deltaF_mo_ann_iPH * mo)) 
    exp(-iota * t - ((exp(betaM_ann_iPH * t) - 1)/betaM_ann_iPH) * exp(alphaM_ann_iPH + betaM_ann_iPH * (x_male_ann - x_bar - (1 - mo) * ad) + deltaM_ad_ann_iPH * log(ad) + deltaM_mo_ann_iPH * mo) - ((exp(betaF_ann_iPH * t) - 1)/betaF_ann_iPH) * exp(alphaF_ann_iPH + betaF_ann_iPH * (x_male_ann - x_bar - ad * mo) + deltaF_ad_ann_iPH * log(ad) + deltaF_mo_ann_iPH * mo)) 
  }
  
  annuity = integrate(integrand, lower=0, upper=100)$value
  
  return(annuity)
}




## - Function of interest rates, ad and mo (+1/0)
ann_last <- function(iota, x_male_ann, ad, mo){
  den_annuity = sum(pi_mw_ann * dnorm(log(ad), mean=zeta_AD_ann, sd=sqrt(var_AD_ann)) * (zeta_MO_ann^(mo)) * ((1 - zeta_MO_ann)^(1 - mo)))
  annuity <- rep(0, K)
  
  for(k in 1:K){
    integrand <- function(t){
      #      exp(-iota * t) * (1 - (1 - exp(-((exp(betaM_ann * t) - 1)/betaM_ann) * exp(alphaM_ann + betaM_ann * x_male_ann + gammaM_ann[k]))) * (1 - exp(-((exp(betaF_ann * t) - 1)/betaF_ann) * exp(alphaF_ann + betaF_ann * (x_male_ann - ad * mo + ad * (1 - mo)) + gammaF_ann[k]))))
      exp(-iota * t) * (1 - (1 - exp(-((exp(betaM_ann * t) - 1)/betaM_ann) * exp(alphaM_ann + betaM_ann * (x_male_ann - x_bar - (1 - mo) * ad) + gammaM_ann[k]))) * (1 - exp(-((exp(betaF_ann * t) - 1)/betaF_ann) * exp(alphaF_ann + betaF_ann * (x_male_ann - x_bar - ad * mo) + gammaF_ann[k]))))
    }
    
    #    annuity[k] = integrate(integrand, lower=0, upper=(120-min(x_male_ann, (x_male_ann - ad * mo + ad * (1 - mo)))))$value * dnorm(log(ad), mean=zeta_AD_ann[k], sd=sqrt(var_AD_ann)) * (zeta_MO_ann[k]^(mo)) * ((1 - zeta_MO_ann[k])^(1 - mo))
    annuity[k] = integrate(integrand, lower=0, upper=100)$value * dnorm(log(ad), mean=zeta_AD_ann[k], sd=sqrt(var_AD_ann)) * (zeta_MO_ann[k]^(mo)) * ((1 - zeta_MO_ann[k])^(1 - mo))
  }
  annuity_value <- sum(pi_mw_ann * annuity) / den_annuity
  return(annuity_value)
  
}

ann_last <- function(iota, x_male_ann, ad, mo){
  
  M <- length(seq_bi_thn) # - eventually modify to just take a smaller sample from the posterior
  annuity <- matrix(0, M, K)
  annuity_value <- rep(0, M)
  
  for(m in 1:M){
    den_annuity = sum(pi_mw_thn[m,] * dnorm(log(ad), mean=z_star_ad_thn[m,], sd=sqrt(z_ad_v_thn[m])) * (z_star_mo_thn[m,]^(mo)) * ((1 - z_star_mo_thn[m,])^(1 - mo)))
    
    for(k in 1:K){
    integrand <- function(t){
      #      exp(-iota * t) * (1 - (1 - exp(-((exp(betaM_ann * t) - 1)/betaM_ann) * exp(alphaM_ann + betaM_ann * x_male_ann + gammaM_ann[k]))) * (1 - exp(-((exp(betaF_ann * t) - 1)/betaF_ann) * exp(alphaF_ann + betaF_ann * (x_male_ann - ad * mo + ad * (1 - mo)) + gammaF_ann[k]))))
      exp(-iota * t) * (1 - (1 - exp(-((exp(beta_thn[m,1] * t) - 1)/beta_thn[m,1]) * exp(alpha_thn[m,1] + beta_thn[m,1] * (x_male_ann - x_bar - (1 - mo) * ad) + gamma_star_thn[m,k*2-1]))) * (1 - exp(-((exp(beta_thn[m,2] * t) - 1)/beta_thn[m,2]) * exp(alpha_thn[m,2] + beta_thn[m,2] * (x_male_ann - x_bar - ad * mo) + gamma_star_thn[m,k*2]))))
    }
    
    #    annuity[k] = integrate(integrand, lower=0, upper=(120-min(x_male_ann, (x_male_ann - ad * mo + ad * (1 - mo)))))$value * dnorm(log(ad), mean=zeta_AD_ann[k], sd=sqrt(var_AD_ann)) * (zeta_MO_ann[k]^(mo)) * ((1 - zeta_MO_ann[k])^(1 - mo))
    annuity[m,k] = integrate(integrand, lower=0, upper=100)$value * dnorm(log(ad), mean=z_star_ad_thn[m,k], sd=sqrt(z_ad_v_thn[m])) * (z_star_mo_thn[m,k]^(mo)) * ((1 - z_star_mo_thn[m,k])^(1 - mo))
  }
  annuity_value[m] <- sum(pi_mw_thn[m,] * annuity[m,]) / den_annuity
  }
  return(mean(annuity_value))
  
}



ann_joint <- function(iota, x_male_ann, ad, mo){
  den_annuity = sum(pi_mw_ann * dnorm(log(ad), mean=zeta_AD_ann, sd=sqrt(var_AD_ann)) * (zeta_MO_ann^(mo)) * ((1 - zeta_MO_ann)^(1 - mo)))
  annuity <- rep(0, K)
  
  for(k in 1:K){
    integrand <- function(t){
      #    exp(-iota * t - ((exp(betaM_ann * t) - 1)/betaM_ann) * exp(alphaM_ann + betaM_ann * x_male_ann + gammaM_ann[k]) - ((exp(betaF_ann * t) - 1)/betaF_ann) * exp(alphaF_ann + betaF_ann * (x_male_ann - ad * mo + ad * (1 - mo)) + gammaF_ann[k])) 
      exp(-iota * t - ((exp(betaM_ann * t) - 1)/betaM_ann) * exp(alphaM_ann + betaM_ann * (x_male_ann - x_bar - (1 - mo) * ad) + gammaM_ann[k]) - ((exp(betaF_ann * t) - 1)/betaF_ann) * exp(alphaF_ann + betaF_ann * (x_male_ann - x_bar - ad * mo) + gammaF_ann[k])) 
    }
    
    annuity[k] = integrate(integrand, lower=0, upper=100)$value * dnorm(log(ad), mean=zeta_AD_ann[k], sd=sqrt(var_AD_ann)) * (zeta_MO_ann[k]^(mo)) * ((1 - zeta_MO_ann[k])^(1 - mo))
  }
  annuity_value <- sum(pi_mw_ann * annuity) / den_annuity
  return(annuity_value)
}


ann_joint <- function(iota, x_male_ann, ad, mo){
  M <- length(seq_bi_thn) # - eventually modify to just take a smaller sample from the posterior
  annuity <- matrix(0, M, K)
  annuity_value <- rep(0, M)
  
  for(m in 1:M){
    den_annuity = sum(pi_mw_thn[m,] * dnorm(log(ad), mean=z_star_ad_thn[m,], sd=sqrt(z_ad_v_thn[m])) * (z_star_mo_thn[m,]^(mo)) * ((1 - z_star_mo_thn[m,])^(1 - mo)))

  for(k in 1:K){
    integrand <- function(t){
      #    exp(-iota * t - ((exp(betaM_ann * t) - 1)/betaM_ann) * exp(alphaM_ann + betaM_ann * x_male_ann + gammaM_ann[k]) - ((exp(betaF_ann * t) - 1)/betaF_ann) * exp(alphaF_ann + betaF_ann * (x_male_ann - ad * mo + ad * (1 - mo)) + gammaF_ann[k])) 
      exp(-iota * t - ((exp(beta_thn[m,1] * t) - 1)/beta_thn[m,1]) * exp(alpha_thn[m,1] + beta_thn[m,1] * (x_male_ann - x_bar - (1 - mo) * ad) + gamma_star_thn[m,k*2-1]) - ((exp(beta_thn[m,2] * t) - 1)/beta_thn[m,2]) * exp(alpha_thn[m,2] + beta_thn[m,2] * (x_male_ann - x_bar - ad * mo) + gamma_star_thn[m,k*2])) 
    }
    
    annuity[m,k] = integrate(integrand, lower=0, upper=100)$value * dnorm(log(ad), mean=z_star_ad_thn[m,k], sd=sqrt(z_ad_v_thn[m])) * (z_star_mo_thn[m,k]^(mo)) * ((1 - z_star_mo_thn[m,k])^(1 - mo))
  }
  annuity_value[m] <- sum(pi_mw_thn[m,] * annuity[m,]) / den_annuity
  }
  return(mean(annuity_value))
}


## - Computations

# - At maximum likelihood

alphaM_ann_iG <- par_est_male_s$estimate[1]
alphaF_ann_iG <- par_est_female_s$estimate[1]

betaM_ann_iG <- par_est_male_s$estimate[2]
betaF_ann_iG <- par_est_female_s$estimate[2]


alphaM_ann_iPH <- par_est_male_s_cov$estimate[1]
alphaF_ann_iPH <- par_est_female_s_cov$estimate[1]

betaM_ann_iPH <- par_est_male_s_cov$estimate[2]
betaF_ann_iPH <- par_est_female_s_cov$estimate[2]

deltaM_ad_ann_iPH <- par_est_male_s_cov$estimate[3]
deltaF_ad_ann_iPH <- par_est_female_s_cov$estimate[3]

deltaM_mo_ann_iPH <- par_est_male_s_cov$estimate[4]
deltaF_mo_ann_iPH <- par_est_female_s_cov$estimate[4]


# - At the posterior mean of the parameters
#alphaM_ann <- alpha_post_m[1]
#alphaF_ann <- alpha_post_m[2]

#betaM_ann <- beta_post_m[1]
#betaF_ann <- beta_post_m[2]

#gammaM_ann <- gamma_star_post_m_M
#gammaF_ann <- gamma_star_post_m_F

#pi_mw_ann <- pi_mw_post_m

#zeta_AD_ann <- z_star_ad_post_m
#zeta_MO_ann <- z_star_mo_post_m

#var_AD_ann <- z_ad_v_post_m

age_ann <- c(60, 70)
iota_ann <- c(0.01, 0.05)
sequence <- seq(1,20,1) # c(seq(0.1,1,0.1), seq(2, 20, 1))# seq(0.5,20,0.5)

# - produce annuity tables as continuum with respect to AD (understand the impact of dependence on joint and last survivor annuities)
ann_J_tb_G_i1_x60 <- matrix(NA, 2, length(sequence))
ann_J_tb_PH_i1_x60 <- matrix(NA, 2, length(sequence))
ann_J_tb_AV_i1_x60 <- matrix(NA, 2, length(sequence))

ann_L_tb_G_i1_x60 <- matrix(NA, 2, length(sequence))
ann_L_tb_PH_i1_x60 <- matrix(NA, 2, length(sequence))
ann_L_tb_AV_i1_x60 <- matrix(NA, 2, length(sequence))


ann_J_tb_G_i5_x60 <- matrix(NA, 2, length(sequence))
ann_J_tb_PH_i5_x60 <- matrix(NA, 2, length(sequence))
ann_J_tb_AV_i5_x60 <- matrix(NA, 2, length(sequence))

ann_L_tb_G_i5_x60 <- matrix(NA, 2, length(sequence))
ann_L_tb_PH_i5_x60 <- matrix(NA, 2, length(sequence))
ann_L_tb_AV_i5_x60 <- matrix(NA, 2, length(sequence))


ann_J_tb_G_i1_x70 <- matrix(NA, 2, length(sequence))
ann_J_tb_PH_i1_x70 <- matrix(NA, 2, length(sequence))
ann_J_tb_AV_i1_x70 <- matrix(NA, 2, length(sequence))

ann_L_tb_G_i1_x70 <- matrix(NA, 2, length(sequence))
ann_L_tb_PH_i1_x70 <- matrix(NA, 2, length(sequence))
ann_L_tb_AV_i1_x70 <- matrix(NA, 2, length(sequence))


ann_J_tb_G_i5_x70 <- matrix(NA, 2, length(sequence))
ann_J_tb_PH_i5_x70 <- matrix(NA, 2, length(sequence))
ann_J_tb_AV_i5_x70 <- matrix(NA, 2, length(sequence))

ann_L_tb_G_i5_x70 <- matrix(NA, 2, length(sequence))
ann_L_tb_PH_i5_x70 <- matrix(NA, 2, length(sequence))
ann_L_tb_AV_i5_x70 <- matrix(NA, 2, length(sequence))

# - 1
iota <- 0.01
age <- 60

col <- 1
for(ad in sequence){
  row <- 1
  for(mo in c(1,0)){
    ann_J_tb_G_i1_x60[row,col] <- ann_joint_indG(iota, age, ad, mo)
    ann_J_tb_PH_i1_x60[row,col] <- ann_joint_indPH(iota, age, ad, mo)
    ann_J_tb_AV_i1_x60[row,col] <- ann_joint(iota, age, ad, mo)
    ann_L_tb_G_i1_x60[row,col] <- ann_last_indG(iota, age, ad, mo)
    ann_L_tb_PH_i1_x60[row,col] <- ann_last_indPH(iota, age, ad, mo)
    ann_L_tb_AV_i1_x60[row,col] <- ann_last(iota, age, ad, mo)
    row <- row + 1
  }
  col <- col + 1
}

# - 2
iota <- 0.05

col <- 1
for(ad in sequence){
  row <- 1
  for(mo in c(1,0)){
    ann_J_tb_G_i5_x60[row,col] <- ann_joint_indG(iota, age, ad, mo)
    ann_J_tb_PH_i5_x60[row,col] <- ann_joint_indPH(iota, age, ad, mo)
    ann_J_tb_AV_i5_x60[row,col] <- ann_joint(iota, age, ad, mo)
    ann_L_tb_G_i5_x60[row,col] <- ann_last_indG(iota, age, ad, mo)
    ann_L_tb_PH_i5_x60[row,col] <- ann_last_indPH(iota, age, ad, mo)
    ann_L_tb_AV_i5_x60[row,col] <- ann_last(iota, age, ad, mo)
    row <- row + 1
  }
  col <- col + 1
}

# - 3
iota <- 0.01
age <- 70

col <- 1
for(ad in sequence){
  row <- 1
  for(mo in c(1,0)){
    ann_J_tb_G_i1_x70[row,col] <- ann_joint_indG(iota, age, ad, mo)
    ann_J_tb_PH_i1_x70[row,col] <- ann_joint_indPH(iota, age, ad, mo)
    ann_J_tb_AV_i1_x70[row,col] <- ann_joint(iota, age, ad, mo)
    ann_L_tb_G_i1_x70[row,col] <- ann_last_indG(iota, age, ad, mo)
    ann_L_tb_PH_i1_x70[row,col] <- ann_last_indPH(iota, age, ad, mo)
    ann_L_tb_AV_i1_x70[row,col] <- ann_last(iota, age, ad, mo)
    row <- row + 1
  }
  col <- col + 1
}

# - 4
iota <- 0.05

col <- 1
for(ad in sequence){
  row <- 1
  for(mo in c(1,0)){
    ann_J_tb_G_i5_x70[row,col] <- ann_joint_indG(iota, age, ad, mo)
    ann_J_tb_PH_i5_x70[row,col] <- ann_joint_indPH(iota, age, ad, mo)
    ann_J_tb_AV_i5_x70[row,col] <- ann_joint(iota, age, ad, mo)
    ann_L_tb_G_i5_x70[row,col] <- ann_last_indG(iota, age, ad, mo)
    ann_L_tb_PH_i5_x70[row,col] <- ann_last_indPH(iota, age, ad, mo)
    ann_L_tb_AV_i5_x70[row,col] <- ann_last(iota, age, ad, mo)
    row <- row + 1
  }
  col <- col + 1
}

# - Percentage difference plots
par(mfrow=c(2,2), mar = c(4.5,4,1,0.5))
plot(sequence, 100 * (ann_J_tb_AV_i1_x60[1,] / ann_J_tb_PH_i1_x60[1,] - 1), type='l', xlab="Age difference", ylim=c(-1,20), ylab='Joint L. % Diff (Age 60)', lwd=2, main="iota=0.01")
abline(h=0, lwd=0.5)
lines(sequence, 100 * (ann_J_tb_AV_i1_x60[2,] / ann_J_tb_PH_i1_x60[2,] - 1), col='grey', lwd=2)
lines(sequence, 100 * (ann_J_tb_AV_i1_x60[1,] / ann_J_tb_G_i1_x60[1,] - 1), lty=3, lwd=1.5)
lines(sequence, 100 * (ann_J_tb_AV_i1_x60[2,] / ann_J_tb_G_i1_x60[2,] - 1), col='grey', lty=3, lwd=1.5)
legend('topright', legend=c("MO=1", "MO=0", "Covariates", "No covariates"), col=c('black', 'grey', "black", "black"), lwd=c(2,2,1.5,1.5), lty=c(1,1,1,3), bty='n', cex=0.75)

plot(sequence, 100 * (ann_J_tb_AV_i5_x60[1,] / ann_J_tb_PH_i5_x60[1,] - 1), type='l', xlab="Age difference", ylim=c(-1,20), ylab='', lwd=2, main="iota=0.05")
abline(h=0, lwd=0.5)
lines(sequence, 100 * (ann_J_tb_AV_i5_x60[2,] / ann_J_tb_PH_i5_x60[2,] - 1), col='grey', lwd=2)
lines(sequence, 100 * (ann_J_tb_AV_i5_x60[1,] / ann_J_tb_G_i5_x60[1,] - 1), lty=3, lwd=1.5)
lines(sequence, 100 * (ann_J_tb_AV_i5_x60[2,] / ann_J_tb_G_i5_x60[2,] - 1), col='grey', lty=3, lwd=1.5)
legend('topright', legend=c("MO=1", "MO=0", "Covariates", "No covariates"), col=c('black', 'grey', "black", "black"), lwd=c(2,2,1.5,1.5), lty=c(1,1,1,3), bty='n', cex=0.75)

plot(sequence, 100 * (ann_J_tb_AV_i1_x70[1,] / ann_J_tb_PH_i1_x70[1,] - 1), type='l', xlab="Age difference", ylim=c(-1,20), ylab='Joint L. % Diff (Age 70)', lwd=2)
abline(h=0, lwd=0.5)
lines(sequence, 100 * (ann_J_tb_AV_i1_x70[2,] / ann_J_tb_PH_i1_x70[2,] - 1), col='grey', lwd=2)
lines(sequence, 100 * (ann_J_tb_AV_i1_x70[1,] / ann_J_tb_G_i1_x70[1,] - 1), lty=3, lwd=1.5)
lines(sequence, 100 * (ann_J_tb_AV_i1_x70[2,] / ann_J_tb_G_i1_x70[2,] - 1), col='grey', lty=3, lwd=1.5)
legend('topright', legend=c("MO=1", "MO=0", "Covariates", "No covariates"), col=c('black', 'grey', "black", "black"), lwd=c(2,2,1.5,1.5), lty=c(1,1,1,3), bty='n', cex=0.75)

plot(sequence, 100 * (ann_J_tb_AV_i5_x70[1,] / ann_J_tb_PH_i5_x70[1,] - 1), type='l', xlab="Age difference", ylim=c(-1,20), ylab='', lwd=2)
abline(h=0, lwd=0.5)
lines(sequence, 100 * (ann_J_tb_AV_i5_x70[2,] / ann_J_tb_PH_i5_x70[2,] - 1), col='grey', lwd=2)
lines(sequence, 100 * (ann_J_tb_AV_i5_x70[1,] / ann_J_tb_G_i5_x70[1,] - 1), lty=3, lwd=1.5)
lines(sequence, 100 * (ann_J_tb_AV_i5_x70[2,] / ann_J_tb_G_i5_x70[2,] - 1), col='grey', lty=3, lwd=1.5)
legend('topright', legend=c("MO=1", "MO=0", "Covariates", "No covariates"), col=c('black', 'grey', "black", "black"), lwd=c(2,2,1.5,1.5), lty=c(1,1,1,3), bty='n', cex=0.75)

## - Last survivor annuity
par(mfrow=c(2,2), mar = c(4.5,4,1,0.5))
plot(sequence, 100 * (ann_L_tb_AV_i1_x60[1,] / ann_L_tb_PH_i1_x60[1,] - 1), type='l', xlab="Age difference", ylim=c(-3,20), ylab='Last S. % Diff (Age 60)', lwd=2, main="iota=0.01")
abline(h=0, lwd=0.5)
lines(sequence, 100 * (ann_L_tb_AV_i1_x60[2,] / ann_L_tb_PH_i1_x60[2,] - 1), col='grey', lwd=2)
lines(sequence, 100 * (ann_L_tb_AV_i1_x60[1,] / ann_L_tb_G_i1_x60[1,] - 1), lty=3, lwd=1.5)
lines(sequence, 100 * (ann_L_tb_AV_i1_x60[2,] / ann_L_tb_G_i1_x60[2,] - 1), col='grey', lty=3, lwd=1.5)
legend('topright', legend=c("MO=1", "MO=0", "Covariates", "No covariates"), col=c('black', 'grey', "black", "black"), lwd=c(2,2,1.5,1.5), lty=c(1,1,1,3), bty='n', cex=0.75)

plot(sequence, 100 * (ann_L_tb_AV_i5_x60[1,] / ann_L_tb_PH_i5_x60[1,] - 1), type='l', xlab="Age difference", ylim=c(-3,20), ylab='', lwd=2, main="iota=0.05")
abline(h=0, lwd=0.5)
lines(sequence, 100 * (ann_L_tb_AV_i5_x60[2,] / ann_L_tb_PH_i5_x60[2,] - 1), col='grey', lwd=2)
lines(sequence, 100 * (ann_L_tb_AV_i5_x60[1,] / ann_L_tb_G_i5_x60[1,] - 1), lty=3, lwd=1.5)
lines(sequence, 100 * (ann_L_tb_AV_i5_x60[2,] / ann_L_tb_G_i5_x60[2,] - 1), col='grey', lty=3, lwd=1.5)
legend('topright', legend=c("MO=1", "MO=0", "Covariates", "No covariates"), col=c('black', 'grey', "black", "black"), lwd=c(2,2,1.5,1.5), lty=c(1,1,1,3), bty='n', cex=0.75)

plot(sequence, 100 * (ann_L_tb_AV_i1_x70[1,] / ann_L_tb_PH_i1_x70[1,] - 1), type='l', xlab="Age difference", ylim=c(-3,20), ylab='Last S. % Diff (Age 70)', lwd=2)
abline(h=0, lwd=0.5)
lines(sequence, 100 * (ann_L_tb_AV_i1_x70[2,] / ann_L_tb_PH_i1_x70[2,] - 1), col='grey', lwd=2)
lines(sequence, 100 * (ann_L_tb_AV_i1_x70[1,] / ann_L_tb_G_i1_x70[1,] - 1), lty=3, lwd=1.5)
lines(sequence, 100 * (ann_L_tb_AV_i1_x70[2,] / ann_L_tb_G_i1_x70[2,] - 1), col='grey', lty=3, lwd=1.5)
legend('topright', legend=c("MO=1", "MO=0", "Covariates", "No covariates"), col=c('black', 'grey', "black", "black"), lwd=c(2,2,1.5,1.5), lty=c(1,1,1,3), bty='n', cex=0.75)

plot(sequence, 100 * (ann_L_tb_AV_i5_x70[1,] / ann_L_tb_PH_i5_x70[1,] - 1), type='l', xlab="Age difference", ylim=c(-3,20), ylab='', lwd=2)
abline(h=0, lwd=0.5)
lines(sequence, 100 * (ann_L_tb_AV_i5_x70[2,] / ann_L_tb_PH_i5_x70[2,] - 1), col='grey', lwd=2)
lines(sequence, 100 * (ann_L_tb_AV_i5_x70[1,] / ann_L_tb_G_i5_x70[1,] - 1), lty=3, lwd=1.5)
lines(sequence, 100 * (ann_L_tb_AV_i5_x70[2,] / ann_L_tb_G_i5_x70[2,] - 1), col='grey', lty=3, lwd=1.5)
legend('topright', legend=c("MO=1", "MO=0", "Covariates", "No covariates"), col=c('black', 'grey', "black", "black"), lwd=c(2,2,1.5,1.5), lty=c(1,1,1,3), bty='n', cex=0.75)




# - Annuity value plots
par(mfrow=c(2,2), mar = c(4.5,4,1,0.5))
plot(sequence, ann_J_tb_AV_i1_x60[1,], type='l', ylim=c(15,25), xlab="Age difference", ylab='Ann. f. (Age 60)', lwd=2, main="iota=0.01")
lines(sequence, ann_J_tb_AV_i1_x60[2,], col='grey', lwd=2)
lines(sequence, ann_J_tb_PH_i1_x60[1,], lwd=2, lty=2)
lines(sequence, ann_J_tb_PH_i1_x60[2,], col='grey', lwd=2, lty=2)
lines(sequence, ann_J_tb_G_i1_x60[1,], lwd=2, lty=3)
lines(sequence, ann_J_tb_G_i1_x60[2,], col='grey', lwd=2, lty=3)
legend('bottomright', legend=c("MO=1", "MO=0"), col=c('black', 'grey'), lwd=2, bty='n', cex=0.75)
legend('topleft', legend=c("AVDPM", "PH", "BG"), lty=c(1,2,3), lwd=2, bty='n', cex=0.75)

plot(sequence, ann_J_tb_AV_i5_x60[1,], type='l', ylim=c(11,15), xlab="Age difference", ylab='Ann. f. (Age 60)', lwd=2, main="iota=0.05")
lines(sequence, ann_J_tb_AV_i5_x60[2,], col='grey', lwd=2)
lines(sequence, ann_J_tb_PH_i5_x60[1,], lwd=2, lty=2)
lines(sequence, ann_J_tb_PH_i5_x60[2,], col='grey', lwd=2, lty=2)
lines(sequence, ann_J_tb_G_i5_x60[1,], lwd=2, lty=3)
lines(sequence, ann_J_tb_G_i5_x60[2,], col='grey', lwd=2, lty=3)
legend('bottomright', legend=c("MO=1", "MO=0"), col=c('black', 'grey'), lwd=2, bty='n', cex=0.75)
legend('topleft', legend=c("AVDPM", "PH", 'BG'), lty=c(1,2,3), lwd=2, bty='n', cex=0.75)

plot(sequence, ann_J_tb_AV_i1_x70[1,], type='l', ylim=c(10,18), xlab="Age difference", ylab='Ann. f. (Age 70)', lwd=2)
lines(sequence, ann_J_tb_AV_i1_x70[2,], col='grey', lwd=2)
lines(sequence, ann_J_tb_PH_i1_x70[1,], lwd=2, lty=2)
lines(sequence, ann_J_tb_PH_i1_x70[2,], col='grey', lwd=2, lty=2)
lines(sequence, ann_J_tb_G_i1_x70[1,], lwd=2, lty=3)
lines(sequence, ann_J_tb_G_i1_x70[2,], col='grey', lwd=2, lty=3)
legend('bottomright', legend=c("MO=1", "MO=0"), col=c('black', 'grey'), lwd=2, bty='n', cex=0.75)
legend('topleft', legend=c("AVDPM", "PH", "BG"), lty=c(1,2,3), lwd=2, bty='n', cex=0.75)

plot(sequence,ann_J_tb_AV_i5_x70[1,], type='l', ylim=c(8,12), xlab="Age difference", ylab='Ann. f. (Age 70)', lwd=2)
lines(sequence,ann_J_tb_AV_i5_x70[2,], col='grey', lwd=2)
lines(sequence,ann_J_tb_PH_i5_x70[1,], lwd=2, lty=2)
lines(sequence,ann_J_tb_PH_i5_x70[2,], col='grey', lwd=2, lty=2)
lines(sequence, ann_J_tb_G_i5_x70[1,], lwd=2, lty=3)
lines(sequence, ann_J_tb_G_i5_x70[2,], col='grey', lwd=2, lty=3)
legend('bottomright', legend=c("MO=1", "MO=0"), col=c('black', 'grey'), lwd=2, bty='n', cex=0.75)
legend('topleft', legend=c("AVDPM", "PH", "BG"), lty=c(1,2,3), lwd=2, bty='n', cex=0.75)

par(mfrow=c(2,2), mar = c(4.5,4,1,0.5))
plot(sequence,ann_L_tb_AV_i1_x60[1,], type='l', ylim=c(25,40), xlab="Age difference", ylab='Ann. f. (Age 60)', lwd=2, main="iota=0.01")
lines(sequence,ann_L_tb_AV_i1_x60[2,], col='grey', lwd=2)
lines(sequence,ann_L_tb_PH_i1_x60[1,], lwd=2, lty=2)
lines(sequence,ann_L_tb_PH_i1_x60[2,], col='grey', lwd=2, lty=2)
lines(sequence, ann_L_tb_G_i1_x60[1,], lwd=2, lty=3)
lines(sequence, ann_L_tb_G_i1_x60[2,], col='grey', lwd=2, lty=3)
legend('bottomright', legend=c("MO=1", "MO=0"), col=c('black', 'grey'), lwd=2, bty='n', cex=0.75)
legend('topleft', legend=c("AVDPM", "PH", "BG"), lty=c(1,2,3), lwd=2, bty='n', cex=0.75)

plot(sequence,ann_L_tb_AV_i5_x60[1,], type='l', ylim=c(15,18), xlab="Age difference", ylab='Ann. f. (Age 60)', lwd=2, main="iota=0.05")
lines(sequence,ann_L_tb_AV_i5_x60[2,], col='grey', lwd=2)
lines(sequence,ann_L_tb_PH_i5_x60[1,], lwd=2, lty=2)
lines(sequence,ann_L_tb_PH_i5_x60[2,], col='grey', lwd=2, lty=2)
lines(sequence, ann_L_tb_G_i5_x60[1,], lwd=2, lty=3)
lines(sequence, ann_L_tb_G_i5_x60[2,], col='grey', lwd=2, lty=3)
legend('bottomright', legend=c("MO=1", "MO=0"), col=c('black', 'grey'), lwd=2, bty='n', cex=0.75)
legend('topleft', legend=c("AVDPM", "PH", "BG"), lty=c(1,2,3), lwd=2, bty='n', cex=0.75)

plot(sequence,ann_L_tb_AV_i1_x70[1,], type='l', ylim=c(17,32), xlab="Age difference", ylab='Ann. f. (Age 70)', lwd=2)
lines(sequence,ann_L_tb_AV_i1_x70[2,], col='grey', lwd=2)
lines(sequence,ann_L_tb_PH_i1_x70[1,], lwd=2, lty=2)
lines(sequence,ann_L_tb_PH_i1_x70[2,], col='grey', lwd=2, lty=2)
lines(sequence, ann_L_tb_G_i1_x70[1,], lwd=2, lty=3)
lines(sequence, ann_L_tb_G_i1_x70[2,], col='grey', lwd=2, lty=3)
legend('bottomright', legend=c("MO=1", "MO=0"), col=c('black', 'grey'), lwd=2, bty='n', cex=0.75)
legend('topleft', legend=c("AVDPM", "PH", "BG"), lty=c(1,2,3), lwd=2, bty='n', cex=0.75)

plot(sequence,ann_L_tb_AV_i5_x70[1,], type='l', ylim=c(12.5,17), xlab="Age difference", ylab='Ann. f. (Age 70)', lwd=2)
lines(sequence,ann_L_tb_AV_i5_x70[2,], col='grey', lwd=2)
lines(sequence,ann_L_tb_PH_i5_x70[1,], lwd=2, lty=2)
lines(sequence,ann_L_tb_PH_i5_x70[2,], col='grey', lwd=2, lty=2)
lines(sequence, ann_L_tb_G_i5_x70[1,], lwd=2, lty=3)
lines(sequence, ann_L_tb_G_i5_x70[2,], col='grey', lwd=2, lty=3)
legend('bottomright', legend=c("MO=1", "MO=0"), col=c('black', 'grey'), lwd=2, bty='n', cex=0.75)
legend('topleft', legend=c("AVDPM", "PH", "BG"), lty=c(1,2,3), lwd=2, bty='n', cex=0.75)



ad_mo_tb <- matrix(NA, 6, 2)
ad_mo_tb[1,] <- c(2, 0)
ad_mo_tb[2,] <- c(2, 1)
ad_mo_tb[3,] <- c(5, 0)
ad_mo_tb[4,] <- c(5, 1)
ad_mo_tb[5,] <- c(10, 0)
ad_mo_tb[6,] <- c(10, 1)

age_count <- 1
admo_count <- 1
iota_count <- 1

annuity_table <- matrix(NA, 12, 8)

# - age = 60
## - iota = 0.01
for(admo in 1:6){
  annuity_table[admo,1] <- ann_joint(0.01, 60, ad_mo_tb[admo,1], ad_mo_tb[admo,2])
  annuity_table[admo,2] <- ann_joint_indG(0.01, 60, ad_mo_tb[admo,1], ad_mo_tb[admo,2])
  annuity_table[admo,3] <- ann_joint(0.05, 60, ad_mo_tb[admo,1], ad_mo_tb[admo,2])
  annuity_table[admo,4] <- ann_joint_indG(0.05, 60, ad_mo_tb[admo,1], ad_mo_tb[admo,2])
  annuity_table[admo,5] <- ann_joint(0.01, 70, ad_mo_tb[admo,1], ad_mo_tb[admo,2])
  annuity_table[admo,6] <- ann_joint_indG(0.01, 70, ad_mo_tb[admo,1], ad_mo_tb[admo,2])
  annuity_table[admo,7] <- ann_joint(0.05, 70, ad_mo_tb[admo,1], ad_mo_tb[admo,2])
  annuity_table[admo,8] <- ann_joint_indG(0.05, 70, ad_mo_tb[admo,1], ad_mo_tb[admo,2])
  
  annuity_table[admo + 6,1] <- ann_last(0.01, 60, ad_mo_tb[admo,1], ad_mo_tb[admo,2])
  annuity_table[admo + 6,2] <- ann_last_indG(0.01, 60, ad_mo_tb[admo,1], ad_mo_tb[admo,2])
  annuity_table[admo + 6,3] <- ann_last(0.05, 60, ad_mo_tb[admo,1], ad_mo_tb[admo,2])
  annuity_table[admo + 6,4] <- ann_last_indG(0.05, 60, ad_mo_tb[admo,1], ad_mo_tb[admo,2])
  annuity_table[admo + 6,5] <- ann_last(0.01, 70, ad_mo_tb[admo,1], ad_mo_tb[admo,2])
  annuity_table[admo + 6,6] <- ann_last_indG(0.01, 70, ad_mo_tb[admo,1], ad_mo_tb[admo,2])
  annuity_table[admo + 6,7] <- ann_last(0.05, 70, ad_mo_tb[admo,1], ad_mo_tb[admo,2])
  annuity_table[admo + 6,8] <- ann_last_indG(0.05, 70, ad_mo_tb[admo,1], ad_mo_tb[admo,2])
  
}

xtable(annuity_table, digits = c(0,2,2,2,2,2,2,2,2))

ann_joint(0.01, 65, 2, 1)
ann_joint_indG(0.01, 65, 2, 1)


ann_last(0.01, 65, 2, 1)
ann_last_indG(0.01, 65, 2, 1)


# - Trace plot mixing
alpha2 <- alpha
beta2 <- beta
gamma_star2 <- gamma_star

N_k_tb2 <- N_k_tb
N_c_tb2 <- N_c_tb
pi_mw2 <- pi_mw
z_ad_v2 <- z_ad_v
z_star_ad2 <- z_star_ad
z_star_mo2 <- z_star_mo

Sigma_gamma2 <- Sigma_gamma

# - mu_zeta (for the normal distribution of the age-difference)
mu_zeta2 <- mu_zeta
sigma2_zeta2 <- sigma2_zeta

# phi concentration parameter of the Dirichlet Process
phi2 <- phi

#------------------ Prior sensitivity --------------------

alpha3 <- alpha
beta3 <- beta
gamma_star3 <- gamma_star

N_k_tb3 <- N_k_tb
N_c_tb3 <- N_c_tb
pi_mw3 <- pi_mw
z_ad_v3 <- z_ad_v
z_star_ad3 <- z_star_ad
z_star_mo3 <- z_star_mo

Sigma_gamma3 <- Sigma_gamma

# - mu_zeta (for the normal distribution of the age-difference)
mu_zeta3 <- mu_zeta
sigma2_zeta3 <- sigma2_zeta

# phi concentration parameter of the Dirichlet Process
phi3 <- phi


alpha4 <- alpha
beta4 <- beta
gamma_star4 <- gamma_star

N_k_tb4 <- N_k_tb
N_c_tb4 <- N_c_tb
pi_mw4 <- pi_mw
z_ad_v4 <- z_ad_v
z_star_ad4 <- z_star_ad
z_star_mo4 <- z_star_mo

Sigma_gamma4 <- Sigma_gamma

# - mu_zeta (for the normal distribution of the age-difference)
mu_zeta4 <- mu_zeta
sigma2_zeta4 <- sigma2_zeta

# phi concentration parameter of the Dirichlet Process
phi4 <- phi


#---------------------------------------------------------



par(mfrow=c(1,1), mar = c(4.5,4,1,0.5))

plot(c(1:n_iter), rowSums(N_k_tb * log(ifelse(N_k_tb>0,N_k_tb,1))), type='l', ylim=c(0,95000), xlab="Iteration", ylab='Entropy')
lines(c(1:n_iter), rowSums(N_k_tb2 * log(ifelse(N_k_tb2>0, N_k_tb2, 1))), col='grey')

lines(c(1:n_iter), rowSums(N_k_tb3 * log(ifelse(N_k_tb3>0, N_k_tb3, 1))), col='yellow')
lines(c(1:n_iter), rowSums(N_k_tb4 * log(ifelse(N_k_tb4>0, N_k_tb4, 1))), col='green')


plot(c(1:(i-1)), alpha[c(1:(i-1)),1], type='l', ylim=c(-10,0))#)#, ylim=c(-10,0))
lines(c(1:(i-1)), alpha2[c(1:(i-1)),1], col='grey')
lines(c(1:(i-1)), alpha3[c(1:(i-1)),1], col='yellow')
lines(c(1:(i-1)), alpha4[c(1:(i-1)),1], col='green')

plot(density(alpha[seq(80020, 100000,20),1]), xlim=c(-10,0))
lines(density(alpha2[seq(80020, 100000,20),1]), col='grey')
lines(density(alpha3[seq(80020, 100000,20),1]), col='yellow')
lines(density(alpha4[seq(80020, 100000,20),1]), col='green')



#plot(c(1:(i-1)), 100*(alpha[c(1:(i-1)),1] - alpha2[c(1:(i-1)),1])/alpha[c(1:(i-1)),1], type='l')

plot(c(1:(i-1)), alpha[c(1:(i-1)),2], type='l', ylim=c(-20,0))#, ylim=c(-10,0))
lines(c(1:(i-1)), alpha2[c(1:(i-1)),2], col='grey')
lines(c(1:(i-1)), alpha3[c(1:(i-1)),2], col='yellow')
lines(c(1:(i-1)), alpha4[c(1:(i-1)),2], col='green')

plot(density(alpha[seq(80020, 100000,20),2]), xlim=c(-20,0))
lines(density(alpha2[seq(80020, 100000,20),2]), col='grey')
lines(density(alpha3[seq(80020, 100000,20),2]), col='yellow')
lines(density(alpha4[seq(80020, 100000,20),2]), col='green')


plot(c(1:(i-1)), beta[c(1:(i-1)),1], type='l', ylim=c(0,0.20))
lines(c(1:(i-1)), beta2[c(1:(i-1)),1], col='grey')
lines(c(1:(i-1)), beta3[c(1:(i-1)),1], col='yellow')
lines(c(1:(i-1)), beta4[c(1:(i-1)),1], col='green')

plot(density(beta[seq(80020, 100000,20),1]), xlim=c(0,0.20))
lines(density(beta2[seq(80020, 100000,20),1]), col='grey')
lines(density(beta3[seq(80020, 100000,20),1]), col='yellow')
lines(density(beta4[seq(80020, 100000,20),1]), col='green')


plot(c(1:(i-1)), beta[c(1:(i-1)),2], type='l', ylim=c(0,0.20))
lines(c(1:(i-1)), beta2[c(1:(i-1)),2], col='grey')
lines(c(1:(i-1)), beta3[c(1:(i-1)),2], col='yellow')
lines(c(1:(i-1)), beta4[c(1:(i-1)),2], col='green')

plot(density(beta[seq(80020, 100000,20),2]), xlim=c(0,0.20))
lines(density(beta2[seq(80020, 100000,20),2]), col='grey')
lines(density(beta3[seq(80020, 100000,20),2]), col='yellow')
lines(density(beta4[seq(80020, 100000,20),2]), col='green')

## - mu_zeta
plot(c(1:(i-1)), mu_zeta[c(1:(i-1))], type='l', ylim=c(-20,20))
lines(c(1:(i-1)), mu_zeta2[c(1:(i-1))], col='grey')
lines(c(1:(i-1)), mu_zeta3[c(1:(i-1))], col='yellow')
lines(c(1:(i-1)), mu_zeta4[c(1:(i-1))], col='green')

plot(density(mu_zeta[seq(80020, 100000,20)]), xlim=c(-20,100))
lines(density(mu_zeta2[seq(80020, 100000,20)]), col='grey')
lines(density(mu_zeta3[seq(80020, 100000,20)]), col='yellow')
lines(density(mu_zeta4[seq(80020, 100000,20)]), col='green')

## - sigma2_zeta
plot(c(1:(i-1)), sigma2_zeta[c(1:(i-1))], type='l', ylim=c(0,200))
lines(c(1:(i-1)), sigma2_zeta2[c(1:(i-1))], col='grey')
lines(c(1:(i-1)), sigma2_zeta3[c(1:(i-1))], col='yellow')
lines(c(1:(i-1)), sigma2_zeta4[c(1:(i-1))], col='green')

plot(density(sigma2_zeta[seq(80020, 100000,20)]), xlim=c(00,500))
lines(density(sigma2_zeta2[seq(80020, 100000,20)]), col='grey')
lines(density(sigma2_zeta3[seq(80020, 100000,20)]), col='yellow')
lines(density(sigma2_zeta4[seq(80020, 100000,20)]), col='green')


## - z_ad_v (variance Z^A)
plot(c(1:(i-1)), z_ad_v[c(1:(i-1))], type='l', ylim=c(0,1))
lines(c(1:(i-1)), z_ad_v2[c(1:(i-1))], col='grey')
lines(c(1:(i-1)), z_ad_v3[c(1:(i-1))], col='yellow')
lines(c(1:(i-1)), z_ad_v4[c(1:(i-1))], col='green')
plot(density(z_ad_v[seq(80020, 100000,20)]), xlim=c(00,0.6))
lines(density(z_ad_v2[seq(80020, 100000,20)]), col='grey')
lines(density(z_ad_v3[seq(80020, 100000,20)]), col='yellow')
lines(density(z_ad_v4[seq(80020, 100000,20)]), col='green')

## - phi
plot(c(1:(i-1)), phi[c(1:(i-1)),1], type='l', ylim=c(0,10))
lines(c(1:(i-1)), phi2[c(1:(i-1)),1], col='grey')
lines(c(1:(i-1)), phi3[c(1:(i-1)),1], col='yellow')
lines(c(1:(i-1)), phi4[c(1:(i-1)),1], col='green')

plot(density(phi[seq(80020, 100000,20),1]), xlim=c(00,10))
lines(density(phi2[seq(80020, 100000,20),1]), col='grey')
lines(density(phi3[seq(80020, 100000,20),1]), col='yellow')
lines(density(phi4[seq(80020, 100000,20),1]), col='green')


plot(c(1:i), N_c_tb[c(1:i)], type='l')
lines(c(1:i), rowSums(ifelse(N_k_tb2>0,1,0)), col='grey')
lines(c(1:i), rowSums(ifelse(N_k_tb3>0,1,0)), col='yellow')
lines(c(1:i), rowSums(ifelse(N_k_tb4>0,1,0)), col='green')



## - Sigma_gamma
plot(c(1:(i-1)), Sigma_gamma[c(1:(i-1)),1], type='l', ylim=c(0,30))
lines(c(1:(i-1)), Sigma_gamma2[c(1:(i-1)),1], col='grey')
lines(c(1:(i-1)), Sigma_gamma3[c(1:(i-1)),1], col='yellow')
lines(c(1:(i-1)), Sigma_gamma4[c(1:(i-1)),1], col='green')

plot(density(Sigma_gamma[seq(80020, 100000,20),1]), xlim=c(00,20))
lines(density(Sigma_gamma2[seq(80020, 100000,20),1]), col='grey')
lines(density(Sigma_gamma3[seq(80020, 100000,20),1]), col='yellow')
lines(density(Sigma_gamma4[seq(80020, 100000,20),1]), col='green')


## - Sigma_gamma
plot(c(1:(i-1)), Sigma_gamma[c(1:(i-1)),2], type='l', ylim=c(0,30))
lines(c(1:(i-1)), Sigma_gamma2[c(1:(i-1)),2], col='grey')
lines(c(1:(i-1)), Sigma_gamma3[c(1:(i-1)),2], col='yellow')
lines(c(1:(i-1)), Sigma_gamma4[c(1:(i-1)),2], col='green')

plot(density(Sigma_gamma[seq(80020, 100000,20),2]), xlim=c(00,20))
lines(density(Sigma_gamma2[seq(80020, 100000,20),2]), col='grey')
lines(density(Sigma_gamma3[seq(80020, 100000,20),2]), col='yellow')
lines(density(Sigma_gamma4[seq(80020, 100000,20),2]), col='green')

## - Sigma_gamma
plot(c(1:(i-1)), Sigma_gamma[c(1:(i-1)),4], type='l', ylim=c(0,30))
lines(c(1:(i-1)), Sigma_gamma2[c(1:(i-1)),4], col='grey')
lines(c(1:(i-1)), Sigma_gamma3[c(1:(i-1)),4], col='yellow')
lines(c(1:(i-1)), Sigma_gamma4[c(1:(i-1)),4], col='green')
plot(density(Sigma_gamma[seq(80020, 100000,20),4]), xlim=c(00,20))
lines(density(Sigma_gamma2[seq(80020, 100000,20),4]), col='grey')
lines(density(Sigma_gamma3[seq(80020, 100000,20),4]), col='yellow')
lines(density(Sigma_gamma4[seq(80020, 100000,20),4]), col='green')

## - correlation
plot(c(1:(i-1)), Sigma_gamma[c(1:(i-1)),2] / sqrt(Sigma_gamma[c(1:(i-1)),1]*Sigma_gamma[c(1:(i-1)),4]), type='l', ylim=c(-1,1))
lines(c(1:(i-1)), Sigma_gamma2[c(1:(i-1)),2] / sqrt(Sigma_gamma2[c(1:(i-1)),1]*Sigma_gamma2[c(1:(i-1)),4]), type='l', col='grey')
lines(c(1:(i-1)), Sigma_gamma3[c(1:(i-1)),2] / sqrt(Sigma_gamma3[c(1:(i-1)),1]*Sigma_gamma3[c(1:(i-1)),4]), type='l', col='yellow')
lines(c(1:(i-1)), Sigma_gamma4[c(1:(i-1)),2] / sqrt(Sigma_gamma4[c(1:(i-1)),1]*Sigma_gamma4[c(1:(i-1)),4]), type='l', col='green')

plot(density(Sigma_gamma[c(1:(i-1)),2] / sqrt(Sigma_gamma[c(1:(i-1)),1]*Sigma_gamma[c(1:(i-1)),4])), xlim=c(-1,1))
lines(density(Sigma_gamma2[c(1:(i-1)),2] / sqrt(Sigma_gamma2[c(1:(i-1)),1]*Sigma_gamma2[c(1:(i-1)),4])), col='grey')
lines(density(Sigma_gamma3[c(1:(i-1)),2] / sqrt(Sigma_gamma3[c(1:(i-1)),1]*Sigma_gamma3[c(1:(i-1)),4])), col='yellow')
lines(density(Sigma_gamma4[c(1:(i-1)),2] / sqrt(Sigma_gamma4[c(1:(i-1)),1]*Sigma_gamma4[c(1:(i-1)),4])), col='green')


mean(Sigma_gamma[seq(80020, 100000,20),2] / sqrt(Sigma_gamma[seq(80020, 100000,20),1]*Sigma_gamma[seq(80020, 100000,20),4]))
mean(Sigma_gamma2[seq(80020, 100000,20),2] / sqrt(Sigma_gamma2[seq(80020, 100000,20),1]*Sigma_gamma2[seq(80020, 100000,20),4]))
mean(Sigma_gamma3[seq(80020, 100000,20),2] / sqrt(Sigma_gamma3[seq(80020, 100000,20),1]*Sigma_gamma3[seq(80020, 100000,20),4]))
mean(Sigma_gamma4[seq(80020, 100000,20),2] / sqrt(Sigma_gamma4[seq(80020, 100000,20),1]*Sigma_gamma4[seq(80020, 100000,20),4]))


## - gamma* as drawn from the Dirichlet Process

### - Mixture component 1
plot(c(1:(i-1)), gamma_star[c(1:(i-1)),1] + alpha[c(1:(i-1)),1], type='l')
lines(c(1:(i-1)), gamma_star2[c(1:(i-1)),1] + alpha2[c(1:(i-1)),1], col='grey')

plot(c(1:(i-1)), gamma_star[c(1:(i-1)),2] + alpha[c(1:(i-1)),2], type='l')
lines(c(1:(i-1)), gamma_star2[c(1:(i-1)),2] + alpha2[c(1:(i-1)),2], col='grey')

plot(c(1:(i-1)), z_star_mo[c(1:(i-1)),1], type='l')
lines(c(1:(i-1)), z_star_mo2[c(1:(i-1)),1], col='grey')

plot(c(1:(i-1)), z_star_ad[c(1:(i-1)),1], type='l')
lines(c(1:(i-1)), z_star_ad2[c(1:(i-1)),1], col='grey')

plot(c(1:(i-1)), pi_mw[c(1:(i-1)),1], type='l')
lines(c(1:(i-1)), pi_mw2[c(1:(i-1)),1], col='grey')


### - Mixture component 2
plot(c(1:(i-1)), gamma_star[c(1:(i-1)),3] + alpha[c(1:(i-1)),1], type='l')
lines(c(1:(i-1)), gamma_star2[c(1:(i-1)),3] + alpha2[c(1:(i-1)),1], col='grey')

plot(c(1:(i-1)), gamma_star[c(1:(i-1)),4] + alpha[c(1:(i-1)),2], type='l')
lines(c(1:(i-1)), gamma_star2[c(1:(i-1)),4] + alpha2[c(1:(i-1)),2], col='grey')

plot(c(1:(i-1)), z_star_mo[c(1:(i-1)),2], type='l')
lines(c(1:(i-1)), z_star_mo2[c(1:(i-1)),2], col='grey')

plot(c(1:(i-1)), z_star_ad[c(1:(i-1)),2], type='l')
lines(c(1:(i-1)), z_star_ad2[c(1:(i-1)),2], col='grey')

plot(c(1:(i-1)), pi_mw[c(1:(i-1)),2], type='l')
lines(c(1:(i-1)), pi_mw2[c(1:(i-1)),2], col='grey')


### - Mixture component 3
plot(c(1:(i-1)), gamma_star[c(1:(i-1)),5] + alpha[c(1:(i-1)),1], type='l')
lines(c(1:(i-1)), gamma_star2[c(1:(i-1)),5] + alpha2[c(1:(i-1)),1], col='grey')

plot(c(1:(i-1)), gamma_star[c(1:(i-1)),6] + alpha[c(1:(i-1)),2], type='l')
lines(c(1:(i-1)), gamma_star2[c(1:(i-1)),6] + alpha2[c(1:(i-1)),2], col='grey')

plot(c(1:(i-1)), z_star_mo[c(1:(i-1)),3], type='l')
lines(c(1:(i-1)), z_star_mo2[c(1:(i-1)),3], col='grey')

plot(c(1:(i-1)), z_star_ad[c(1:(i-1)),3], type='l')
lines(c(1:(i-1)), z_star_ad2[c(1:(i-1)),3], col='grey')

plot(c(1:(i-1)), pi_mw[c(1:(i-1)),3], type='l')
lines(c(1:(i-1)), pi_mw2[c(1:(i-1)),3], col='grey')

### - Mixture component 4
plot(c(1:(i-1)), gamma_star[c(1:(i-1)),7] + alpha[c(1:(i-1)),1], type='l')
lines(c(1:(i-1)), gamma_star2[c(1:(i-1)),7] + alpha2[c(1:(i-1)),1], col='grey')

plot(c(1:(i-1)), gamma_star[c(1:(i-1)),8] + alpha[c(1:(i-1)),2], type='l')
lines(c(1:(i-1)), gamma_star2[c(1:(i-1)),8] + alpha2[c(1:(i-1)),2], col='grey')

plot(c(1:(i-1)), z_star_mo[c(1:(i-1)),4], type='l')
lines(c(1:(i-1)), z_star_mo2[c(1:(i-1)),4], col='grey')

plot(c(1:(i-1)), z_star_ad[c(1:(i-1)),4], type='l')
lines(c(1:(i-1)), z_star_ad2[c(1:(i-1)),4], col='grey')

plot(c(1:(i-1)), pi_mw[c(1:(i-1)),4], type='l')
lines(c(1:(i-1)), pi_mw2[c(1:(i-1)),4], col='grey')

### - Mixture component 6
plot(c(1:(i-1)), gamma_star[c(1:(i-1)),11] + alpha[c(1:(i-1)),1], type='l')
lines(c(1:(i-1)), gamma_star2[c(1:(i-1)),11] + alpha2[c(1:(i-1)),1], col='grey')

plot(c(1:(i-1)), gamma_star[c(1:(i-1)),12] + alpha[c(1:(i-1)),2], type='l')
lines(c(1:(i-1)), gamma_star2[c(1:(i-1)),12] + alpha2[c(1:(i-1)),2], col='grey')

plot(c(1:(i-1)), z_star_mo[c(1:(i-1)),6], type='l')
lines(c(1:(i-1)), z_star_mo2[c(1:(i-1)),6], col='grey')

plot(c(1:(i-1)), z_star_ad[c(1:(i-1)),6], type='l')
lines(c(1:(i-1)), z_star_ad2[c(1:(i-1)),6], col='grey')

plot(c(1:(i-1)), pi_mw[c(1:(i-1)),6], type='l')
lines(c(1:(i-1)), pi_mw2[c(1:(i-1)),6], col='grey')

### - Mixture component 8
plot(c(1:(i-1)), gamma_star[c(1:(i-1)),15] + alpha[c(1:(i-1)),1], type='l')
lines(c(1:(i-1)), gamma_star2[c(1:(i-1)),15] + alpha2[c(1:(i-1)),1], col='grey')

plot(c(1:(i-1)), gamma_star[c(1:(i-1)),16] + alpha[c(1:(i-1)),2], type='l')
lines(c(1:(i-1)), gamma_star2[c(1:(i-1)),16] + alpha2[c(1:(i-1)),2], col='grey')

plot(c(1:(i-1)), z_star_mo[c(1:(i-1)),8], type='l')
lines(c(1:(i-1)), z_star_mo2[c(1:(i-1)),8], col='grey')

plot(c(1:(i-1)), z_star_ad[c(1:(i-1)),8], type='l')
lines(c(1:(i-1)), z_star_ad2[c(1:(i-1)),8], col='grey')

plot(c(1:(i-1)), pi_mw[c(1:(i-1)),8], type='l')
lines(c(1:(i-1)), pi_mw2[c(1:(i-1)),8], col='grey')



log_abs_z_ad <- log(abs(Z$agediff))
hist(Z$agediff, breaks=1000)
plot(density(Z$agediff))


# - pi (mixture weights from the Dirichlet Process)
plot(c(1:(i-1)), pi_mw[c(1:(i-1)),1], type='l',ylim=c(0,1))
lines(c(1:(i-1)), pi_mw2[c(1:(i-1)),1], col='grey')

plot(c(1:(i-1)), pi_mw[c(1:(i-1)),2], type='l',ylim=c(0,1))
lines(c(1:(i-1)), pi_mw2[c(1:(i-1)),2], col='grey')

plot(c(1:(i-1)), pi_mw[c(1:(i-1)),3], type='l',ylim=c(0,1))
lines(c(1:(i-1)), pi_mw2[c(1:(i-1)),3], col='grey')

plot(c(1:(i-1)), pi_mw[c(1:(i-1)),4], type='l',ylim=c(0,1))
lines(c(1:(i-1)), pi_mw2[c(1:(i-1)),4], col='grey')

## - Number of clusters and cluster occupancy
plot(c(1:(i-1)), N_c_tb[c(1:(i-1))], type='l')
plot(c(1:(i-1)), N_k_tb[c(1:(i-1)),1], type='l', ylim=c(0, sample_size))
plot(c(1:(i-1)), N_k_tb[c(1:(i-1)),2], type='l', ylim=c(0, sample_size))
plot(c(1:(i-1)), N_k_tb[c(1:(i-1)),3], type='l', ylim=c(0, sample_size))
plot(c(1:(i-1)), N_k_tb[c(1:(i-1)),4], type='l', ylim=c(0, sample_size))
plot(c(1:(i-1)), N_k_tb[c(1:(i-1)),5], type='l', ylim=c(0, sample_size))
plot(c(1:(i-1)), N_k_tb[c(1:(i-1)),6], type='l', ylim=c(0, sample_size))



# - Predictive density
t_range <- seq(0, 50, 0.5)

post_pred_male <- post_pred_female <- rep(0, length(t_range))
f_t_k <- S_t_k <- matrix(0, length(t_range), K)
f_t <- rep(0, length(t_range))
f_t_k2 <- S_t_k2 <- matrix(0, length(t_range), K)
f_t2 <- rep(0, length(t_range))
S_t <- S_t2 <- rep(0, length(t_range))

age <- 80

# - beta regression parameters
alpha_thn2 <- alpha2[seq_bi_thn,]
beta_thn2 <- beta2[seq_bi_thn,]

# - gamma frailty parameters
gamma_star_thn2 <- gamma_star2[seq_bi_thn,]

# - zeta age difference
z_star_ad_thn2 <- z_star_ad2[seq_bi_thn,]
z_ad_v_thn2 <- z_ad_v2[seq_bi_thn]

# - zeta oldest male
z_star_mo_thn2 <- z_star_mo2[seq_bi_thn,]

# - pi mixture weight
pi_mw_thn2 <- pi_mw2[seq_bi_thn,]

# - Covariance_theta
Sigma_gamma_thn2 <- Sigma_gamma2[seq_bi_thn,]

# - psi from stick-breaking procedure
psi_sbp_thn2 <- psi_sbp2[seq_bi_thn,]

# phi concentration parameter of the Dirichlet Process
phi_thn2 <- phi2[seq_bi_thn,]



for(t in 1:length(t_range)){
  for(k in 1:K){
    f_t_k[t,k] <- mean(pi_mw_thn[,k] * exp(- exp(alpha_thn[,1] + beta_thn[,1] * (age - x_bar) + gamma_star_thn[,seq(1,49,2)[k]]) * (exp(beta_thn[,1] * t_range[t]) - 1)/beta_thn[,1] + alpha_thn[,1] + beta_thn[,1] * (age - x_bar + t_range[t]) + gamma_star_thn[,seq(1,49,2)[k]]))
    f_t_k2[t,k] <- mean(pi_mw_thn2[,k] * exp(- exp(alpha_thn2[,1] + beta_thn2[,1] * (age - x_bar) + gamma_star_thn2[,seq(1,49,2)[k]]) * (exp(beta_thn2[,1] * t_range[t]) - 1)/beta_thn2[,1] + alpha_thn2[,1] + beta_thn2[,1] * (age - x_bar + t_range[t]) + gamma_star_thn2[,seq(1,49,2)[k]]))
    S_t_k[t,k] <- mean(pi_mw_thn[,k] * exp(- exp(alpha_thn[,1] + beta_thn[,1] * (age - x_bar) + gamma_star_thn[,seq(1,49,2)[k]]) * (exp(beta_thn[,1] * t_range[t]) - 1)/beta_thn[,1]))
    S_t_k2[t,k] <- mean(pi_mw_thn2[,k] * exp(- exp(alpha_thn2[,1] + beta_thn2[,1] * (age - x_bar) + gamma_star_thn2[,seq(1,49,2)[k]]) * (exp(beta_thn2[,1] * t_range[t]) - 1)/beta_thn2[,1]))
    
  }
  f_t[t] <- sum(f_t_k[t,])
  f_t2[t] <- sum(f_t_k2[t,])
  S_t[t] <- sum(S_t_k[t,])
  S_t2[t] <- sum(S_t_k2[t,])
  
}

plot(t_range, f_t, type='l')
lines(t_range, f_t2, col='grey')
plot(t_range, S_t, type='l')
lines(t_range, S_t2, col='grey')

# - Use posterior mean (just as double check)
pi_mw_pm <- colMeans(pi_mw_thn)
alpha_pm <- colMeans(alpha_thn)
beta_pm <- colMeans(beta_thn)
gamma_s_pm <- colMeans(gamma_star_thn)


for(t in 1:length(t_range)){
  
  f_t_k[t,] <- pi_mw_pm * exp(- exp(alpha_pm[1] + beta_pm[1] * age + gamma_s_pm[seq(1,49,2)]) * (exp(beta_pm[1] * t_range[t]) - 1)/beta_pm[1] + alpha_pm[1] + beta_pm[1] * (age + t_range[t]) + gamma_s_pm[seq(1,49,2)])
  f_t[t] <- sum(f_t_k[t,])
}


# - Poisson Deviance residuals
data_test <- data_shuffle[(round(nrow(data_shuffle)*0.75,0)+1):nrow(data_shuffle),]

X <- as.data.frame(cbind(data_test$ageM, data_test$ageF)) - x_bar
#X <-  as.data.frame(cbind(data$ageM, data$ageF))
colnames(X) <- c("male", "female")
Z <- cbind(log(as.numeric(data_test$agediff)), data_test$maleelder)
colnames(Z) <- c("agediff", "maleelder")
Z <- as.data.frame(Z)

sample_size_test <- nrow(data_test)

t1_i <- data_test$TObs1
t2_i <- data_test$TObs2
a_i <- data_test$trunc
d_ci <- as.data.frame(cbind(data_test$dmale, data_test$dfemale))
colnames(d_ci) <- c("male", "female")

age_baseline <- floor(min(X + a_i) + x_bar) # 60
max_age <- max(ceiling(X + c(t1_i, t2_i) + a_i + x_bar)) #90

n_ages <- max_age - age_baseline
  
E_cxM <- matrix(0, sample_size, n_ages) #rep(0, n_ages)
d_xM <- matrix(0, sample_size, n_ages)
st_ageM <- matrix(0, sample_size, n_ages)
E_cxF <- matrix(0, sample_size, n_ages)
d_xF <- matrix(0, sample_size, n_ages)
st_ageF <- matrix(0, sample_size, n_ages)

for(age in age_baseline:(max_age-1)){
#  E_cxM[age + 1 - age_baseline] <- 0
# d_xM[age + 1 - age_baseline] <- 0
#  E_cxF[age + 1 - age_baseline] <- 0
#  d_xF[age + 1 - age_baseline] <- 0
  for(i in 1:sample_size){
    E_cxM[i,age + 1 - age_baseline] <- max(0, min(age + 1, X$male[i] + x_bar + a_i[i] + t1_i[i]) - max(X$male[i] + x_bar + a_i[i], age))
    d_xM[i,age + 1 - age_baseline] <- ifelse(d_ci$male[i]==1 & ((X$male[i] + x_bar + a_i[i] + t1_i[i]) >= age) & ((X$male[i] + x_bar + a_i[i] + t1_i[i]) < age+1), 1, 0)
    st_ageM[i,age + 1 - age_baseline] <- ifelse((E_cxM[i,age + 1 - age_baseline] > 0) & ((X$male[i] + x_bar + a_i[i]) < age), age, X$male[i] + x_bar + a_i[i])
    E_cxF[i,age + 1 - age_baseline] <- max(0, min(age + 1, X$female[i] + x_bar + a_i[i] + t2_i[i]) - max(X$female[i] + x_bar + a_i[i], age))
    d_xF[i,age + 1 - age_baseline] <- ifelse(d_ci$female[i]==1 & ((X$female[i] + x_bar + a_i[i] + t2_i[i]) >= age) & ((X$female[i] + x_bar + a_i[i] + t2_i[i]) < age+1), 1, 0)
    st_ageF[i,age + 1 - age_baseline] <- ifelse((E_cxF[i,age + 1 - age_baseline] > 0) & ((X$female[i] + x_bar + a_i[i]) < age), age, X$female[i] + x_bar + a_i[i])
  }
  
}

deaths_m <- colSums(d_xM)

## - Analysis of the integrated hazard function
### - Gompertz
Int_hf_Gom_m <- (exp(beta_BG_m * E_cxM) - 1) / beta_BG_m

for(row in 1:nrow(E_cxM)){
  for(col in 1:ncol(E_cxM)){
  Int_hf_Gom_m[row,col] <- Int_hf_Gom_m[row,col] * exp(alpha_BG_m + beta_BG_m * (st_ageM[row,col]-x_bar))
}}

### - AVDPM
integrand_m <- function(t){
  
 (pi_mw_post_m[1] * exp(- ((exp(t * beta_post_m[1]) - 1) / beta_post_m[1]) * exp(alpha_post_m[1] + beta_post_m[1] * (st_ageM[i,j] - x_bar) + gamma_star_post_m[1]) + (alpha_post_m[1] + beta_post_m[1] * (st_ageM[i,j] - x_bar + t) + gamma_star_post_m[1]) + Z$maleelder[i] * log(z_star_mo_post_m[1]) + log(1-z_star_mo_post_m[1]) * (1-Z$maleelder[i]) + dnorm(Z$agediff[i], mean = z_star_ad_post_m[1], sd = sqrt(z_ad_v_post_m), log=TRUE)) + 
  pi_mw_post_m[2] * exp(- ((exp(t * beta_post_m[1]) - 1) / beta_post_m[1]) * exp(alpha_post_m[1] + beta_post_m[1] * (st_ageM[i,j] - x_bar) + gamma_star_post_m[2]) + (alpha_post_m[1] + beta_post_m[1] * (st_ageM[i,j] - x_bar + t) + gamma_star_post_m[2]) + Z$maleelder[i] * log(z_star_mo_post_m[2]) + log(1-z_star_mo_post_m[2]) * (1-Z$maleelder[i]) + dnorm(Z$agediff[i], mean = z_star_ad_post_m[2], sd = sqrt(z_ad_v_post_m), log=TRUE)) + 
  pi_mw_post_m[3] * exp(- ((exp(t * beta_post_m[1]) - 1) / beta_post_m[1]) * exp(alpha_post_m[1] + beta_post_m[1] * (st_ageM[i,j] - x_bar) + gamma_star_post_m[3]) + (alpha_post_m[1] + beta_post_m[1] * (st_ageM[i,j] - x_bar + t) + gamma_star_post_m[3]) + Z$maleelder[i] * log(z_star_mo_post_m[3]) + log(1-z_star_mo_post_m[3]) * (1-Z$maleelder[i]) + dnorm(Z$agediff[i], mean = z_star_ad_post_m[3], sd = sqrt(z_ad_v_post_m), log=TRUE)) + 
  pi_mw_post_m[4] * exp(- ((exp(t * beta_post_m[1]) - 1) / beta_post_m[1]) * exp(alpha_post_m[1] + beta_post_m[1] * (st_ageM[i,j] - x_bar) + gamma_star_post_m[4]) + (alpha_post_m[1] + beta_post_m[1] * (st_ageM[i,j] - x_bar + t) + gamma_star_post_m[4]) + Z$maleelder[i] * log(z_star_mo_post_m[4]) + log(1-z_star_mo_post_m[4]) * (1-Z$maleelder[i]) + dnorm(Z$agediff[i], mean = z_star_ad_post_m[4], sd = sqrt(z_ad_v_post_m), log=TRUE)) + 
  pi_mw_post_m[5] * exp(- ((exp(t * beta_post_m[1]) - 1) / beta_post_m[1]) * exp(alpha_post_m[1] + beta_post_m[1] * (st_ageM[i,j] - x_bar) + gamma_star_post_m[5]) + (alpha_post_m[1] + beta_post_m[1] * (st_ageM[i,j] - x_bar + t) + gamma_star_post_m[5]) + Z$maleelder[i] * log(z_star_mo_post_m[5]) + log(1-z_star_mo_post_m[5]) * (1-Z$maleelder[i]) + dnorm(Z$agediff[i], mean = z_star_ad_post_m[5], sd = sqrt(z_ad_v_post_m), log=TRUE)) + 
  pi_mw_post_m[6] * exp(- ((exp(t * beta_post_m[1]) - 1) / beta_post_m[1]) * exp(alpha_post_m[1] + beta_post_m[1] * (st_ageM[i,j] - x_bar) + gamma_star_post_m[6]) + (alpha_post_m[1] + beta_post_m[1] * (st_ageM[i,j] - x_bar + t) + gamma_star_post_m[6]) + Z$maleelder[i] * log(z_star_mo_post_m[6]) + log(1-z_star_mo_post_m[6]) * (1-Z$maleelder[i]) + dnorm(Z$agediff[i], mean = z_star_ad_post_m[6], sd = sqrt(z_ad_v_post_m), log=TRUE)) + 
  pi_mw_post_m[7] * exp(- ((exp(t * beta_post_m[1]) - 1) / beta_post_m[1]) * exp(alpha_post_m[1] + beta_post_m[1] * (st_ageM[i,j] - x_bar) + gamma_star_post_m[7]) + (alpha_post_m[1] + beta_post_m[1] * (st_ageM[i,j] - x_bar + t) + gamma_star_post_m[7]) + Z$maleelder[i] * log(z_star_mo_post_m[7]) + log(1-z_star_mo_post_m[7]) * (1-Z$maleelder[i]) + dnorm(Z$agediff[i], mean = z_star_ad_post_m[7], sd = sqrt(z_ad_v_post_m), log=TRUE)) + 
  pi_mw_post_m[8] * exp(- ((exp(t * beta_post_m[1]) - 1) / beta_post_m[1]) * exp(alpha_post_m[1] + beta_post_m[1] * (st_ageM[i,j] - x_bar) + gamma_star_post_m[8]) + (alpha_post_m[1] + beta_post_m[1] * (st_ageM[i,j] - x_bar + t) + gamma_star_post_m[8]) + Z$maleelder[i] * log(z_star_mo_post_m[8]) + log(1-z_star_mo_post_m[8]) * (1-Z$maleelder[i]) + dnorm(Z$agediff[i], mean = z_star_ad_post_m[8], sd = sqrt(z_ad_v_post_m), log=TRUE)) + 
  pi_mw_post_m[9] * exp(- ((exp(t * beta_post_m[1]) - 1) / beta_post_m[1]) * exp(alpha_post_m[1] + beta_post_m[1] * (st_ageM[i,j] - x_bar) + gamma_star_post_m[9]) + (alpha_post_m[1] + beta_post_m[1] * (st_ageM[i,j] - x_bar + t) + gamma_star_post_m[9]) + Z$maleelder[i] * log(z_star_mo_post_m[9]) + log(1-z_star_mo_post_m[9]) * (1-Z$maleelder[i]) + dnorm(Z$agediff[i], mean = z_star_ad_post_m[9], sd = sqrt(z_ad_v_post_m), log=TRUE)) + 
  pi_mw_post_m[10] * exp(- ((exp(t * beta_post_m[1]) - 1) / beta_post_m[1]) * exp(alpha_post_m[1] + beta_post_m[1] * (st_ageM[i,j] - x_bar) + gamma_star_post_m[10]) + (alpha_post_m[1] + beta_post_m[1] * (st_ageM[i,j] - x_bar + t) + gamma_star_post_m[10]) + Z$maleelder[i] * log(z_star_mo_post_m[10]) + log(1-z_star_mo_post_m[10]) * (1-Z$maleelder[i]) + dnorm(Z$agediff[i], mean = z_star_ad_post_m[10], sd = sqrt(z_ad_v_post_m), log=TRUE)) + 
  pi_mw_post_m[11] * exp(- ((exp(t * beta_post_m[1]) - 1) / beta_post_m[1]) * exp(alpha_post_m[1] + beta_post_m[1] * (st_ageM[i,j] - x_bar) + gamma_star_post_m[11]) + (alpha_post_m[1] + beta_post_m[1] * (st_ageM[i,j] - x_bar + t) + gamma_star_post_m[11]) + Z$maleelder[i] * log(z_star_mo_post_m[11]) + log(1-z_star_mo_post_m[11]) * (1-Z$maleelder[i]) + dnorm(Z$agediff[i], mean = z_star_ad_post_m[11], sd = sqrt(z_ad_v_post_m), log=TRUE)) + 
  pi_mw_post_m[12] * exp(- ((exp(t * beta_post_m[1]) - 1) / beta_post_m[1]) * exp(alpha_post_m[1] + beta_post_m[1] * (st_ageM[i,j] - x_bar) + gamma_star_post_m[12]) + (alpha_post_m[1] + beta_post_m[1] * (st_ageM[i,j] - x_bar + t) + gamma_star_post_m[12]) + Z$maleelder[i] * log(z_star_mo_post_m[12]) + log(1-z_star_mo_post_m[12]) * (1-Z$maleelder[i]) + dnorm(Z$agediff[i], mean = z_star_ad_post_m[12], sd = sqrt(z_ad_v_post_m), log=TRUE)) + 
  pi_mw_post_m[13] * exp(- ((exp(t * beta_post_m[1]) - 1) / beta_post_m[1]) * exp(alpha_post_m[1] + beta_post_m[1] * (st_ageM[i,j] - x_bar) + gamma_star_post_m[13]) + (alpha_post_m[1] + beta_post_m[1] * (st_ageM[i,j] - x_bar + t) + gamma_star_post_m[13]) + Z$maleelder[i] * log(z_star_mo_post_m[13]) + log(1-z_star_mo_post_m[13]) * (1-Z$maleelder[i]) + dnorm(Z$agediff[i], mean = z_star_ad_post_m[13], sd = sqrt(z_ad_v_post_m), log=TRUE)) + 
  pi_mw_post_m[14] * exp(- ((exp(t * beta_post_m[1]) - 1) / beta_post_m[1]) * exp(alpha_post_m[1] + beta_post_m[1] * (st_ageM[i,j] - x_bar) + gamma_star_post_m[14]) + (alpha_post_m[1] + beta_post_m[1] * (st_ageM[i,j] - x_bar + t) + gamma_star_post_m[14]) + Z$maleelder[i] * log(z_star_mo_post_m[14]) + log(1-z_star_mo_post_m[14]) * (1-Z$maleelder[i]) + dnorm(Z$agediff[i], mean = z_star_ad_post_m[14], sd = sqrt(z_ad_v_post_m), log=TRUE)) + 
  pi_mw_post_m[15] * exp(- ((exp(t * beta_post_m[1]) - 1) / beta_post_m[1]) * exp(alpha_post_m[1] + beta_post_m[1] * (st_ageM[i,j] - x_bar) + gamma_star_post_m[15]) + (alpha_post_m[1] + beta_post_m[1] * (st_ageM[i,j] - x_bar + t) + gamma_star_post_m[15]) + Z$maleelder[i] * log(z_star_mo_post_m[15]) + log(1-z_star_mo_post_m[15]) * (1-Z$maleelder[i]) + dnorm(Z$agediff[i], mean = z_star_ad_post_m[15], sd = sqrt(z_ad_v_post_m), log=TRUE)) + 
  pi_mw_post_m[16] * exp(- ((exp(t * beta_post_m[1]) - 1) / beta_post_m[1]) * exp(alpha_post_m[1] + beta_post_m[1] * (st_ageM[i,j] - x_bar) + gamma_star_post_m[16]) + (alpha_post_m[1] + beta_post_m[1] * (st_ageM[i,j] - x_bar + t) + gamma_star_post_m[16]) + Z$maleelder[i] * log(z_star_mo_post_m[16]) + log(1-z_star_mo_post_m[16]) * (1-Z$maleelder[i]) + dnorm(Z$agediff[i], mean = z_star_ad_post_m[16], sd = sqrt(z_ad_v_post_m), log=TRUE)) + 
  pi_mw_post_m[17] * exp(- ((exp(t * beta_post_m[1]) - 1) / beta_post_m[1]) * exp(alpha_post_m[1] + beta_post_m[1] * (st_ageM[i,j] - x_bar) + gamma_star_post_m[17]) + (alpha_post_m[1] + beta_post_m[1] * (st_ageM[i,j] - x_bar + t) + gamma_star_post_m[17]) + Z$maleelder[i] * log(z_star_mo_post_m[17]) + log(1-z_star_mo_post_m[17]) * (1-Z$maleelder[i]) + dnorm(Z$agediff[i], mean = z_star_ad_post_m[17], sd = sqrt(z_ad_v_post_m), log=TRUE)) + 
  pi_mw_post_m[18] * exp(- ((exp(t * beta_post_m[1]) - 1) / beta_post_m[1]) * exp(alpha_post_m[1] + beta_post_m[1] * (st_ageM[i,j] - x_bar) + gamma_star_post_m[18]) + (alpha_post_m[1] + beta_post_m[1] * (st_ageM[i,j] - x_bar + t) + gamma_star_post_m[18]) + Z$maleelder[i] * log(z_star_mo_post_m[18]) + log(1-z_star_mo_post_m[18]) * (1-Z$maleelder[i]) + dnorm(Z$agediff[i], mean = z_star_ad_post_m[18], sd = sqrt(z_ad_v_post_m), log=TRUE)) + 
  pi_mw_post_m[19] * exp(- ((exp(t * beta_post_m[1]) - 1) / beta_post_m[1]) * exp(alpha_post_m[1] + beta_post_m[1] * (st_ageM[i,j] - x_bar) + gamma_star_post_m[19]) + (alpha_post_m[1] + beta_post_m[1] * (st_ageM[i,j] - x_bar + t) + gamma_star_post_m[19]) + Z$maleelder[i] * log(z_star_mo_post_m[19]) + log(1-z_star_mo_post_m[19]) * (1-Z$maleelder[i]) + dnorm(Z$agediff[i], mean = z_star_ad_post_m[19], sd = sqrt(z_ad_v_post_m), log=TRUE)) + 
  pi_mw_post_m[20] * exp(- ((exp(t * beta_post_m[1]) - 1) / beta_post_m[1]) * exp(alpha_post_m[1] + beta_post_m[1] * (st_ageM[i,j] - x_bar) + gamma_star_post_m[20]) + (alpha_post_m[1] + beta_post_m[1] * (st_ageM[i,j] - x_bar + t) + gamma_star_post_m[20]) + Z$maleelder[i] * log(z_star_mo_post_m[20]) + log(1-z_star_mo_post_m[20]) * (1-Z$maleelder[i]) + dnorm(Z$agediff[i], mean = z_star_ad_post_m[20], sd = sqrt(z_ad_v_post_m), log=TRUE)) + 
  pi_mw_post_m[21] * exp(- ((exp(t * beta_post_m[1]) - 1) / beta_post_m[1]) * exp(alpha_post_m[1] + beta_post_m[1] * (st_ageM[i,j] - x_bar) + gamma_star_post_m[21]) + (alpha_post_m[1] + beta_post_m[1] * (st_ageM[i,j] - x_bar + t) + gamma_star_post_m[21]) + Z$maleelder[i] * log(z_star_mo_post_m[21]) + log(1-z_star_mo_post_m[21]) * (1-Z$maleelder[i]) + dnorm(Z$agediff[i], mean = z_star_ad_post_m[21], sd = sqrt(z_ad_v_post_m), log=TRUE)) + 
  pi_mw_post_m[22] * exp(- ((exp(t * beta_post_m[1]) - 1) / beta_post_m[1]) * exp(alpha_post_m[1] + beta_post_m[1] * (st_ageM[i,j] - x_bar) + gamma_star_post_m[22]) + (alpha_post_m[1] + beta_post_m[1] * (st_ageM[i,j] - x_bar + t) + gamma_star_post_m[22]) + Z$maleelder[i] * log(z_star_mo_post_m[22]) + log(1-z_star_mo_post_m[22]) * (1-Z$maleelder[i]) + dnorm(Z$agediff[i], mean = z_star_ad_post_m[22], sd = sqrt(z_ad_v_post_m), log=TRUE)) + 
  pi_mw_post_m[23] * exp(- ((exp(t * beta_post_m[1]) - 1) / beta_post_m[1]) * exp(alpha_post_m[1] + beta_post_m[1] * (st_ageM[i,j] - x_bar) + gamma_star_post_m[23]) + (alpha_post_m[1] + beta_post_m[1] * (st_ageM[i,j] - x_bar + t) + gamma_star_post_m[23]) + Z$maleelder[i] * log(z_star_mo_post_m[23]) + log(1-z_star_mo_post_m[23]) * (1-Z$maleelder[i]) + dnorm(Z$agediff[i], mean = z_star_ad_post_m[23], sd = sqrt(z_ad_v_post_m), log=TRUE)) + 
  pi_mw_post_m[24] * exp(- ((exp(t * beta_post_m[1]) - 1) / beta_post_m[1]) * exp(alpha_post_m[1] + beta_post_m[1] * (st_ageM[i,j] - x_bar) + gamma_star_post_m[24]) + (alpha_post_m[1] + beta_post_m[1] * (st_ageM[i,j] - x_bar + t) + gamma_star_post_m[24]) + Z$maleelder[i] * log(z_star_mo_post_m[24]) + log(1-z_star_mo_post_m[24]) * (1-Z$maleelder[i]) + dnorm(Z$agediff[i], mean = z_star_ad_post_m[24], sd = sqrt(z_ad_v_post_m), log=TRUE)) + 
  pi_mw_post_m[25] * exp(- ((exp(t * beta_post_m[1]) - 1) / beta_post_m[1]) * exp(alpha_post_m[1] + beta_post_m[1] * (st_ageM[i,j] - x_bar) + gamma_star_post_m[25]) + (alpha_post_m[1] + beta_post_m[1] * (st_ageM[i,j] - x_bar + t) + gamma_star_post_m[25]) + Z$maleelder[i] * log(z_star_mo_post_m[25]) + log(1-z_star_mo_post_m[25]) * (1-Z$maleelder[i]) + dnorm(Z$agediff[i], mean = z_star_ad_post_m[25], sd = sqrt(z_ad_v_post_m), log=TRUE)) ) / 
 (pi_mw_post_m[1] * exp(- ((exp(t * beta_post_m[1]) - 1) / beta_post_m[1]) * exp(alpha_post_m[1] + beta_post_m[1] * (st_ageM[i,j] - x_bar) + gamma_star_post_m[1]) + Z$maleelder[i] * log(z_star_mo_post_m[1]) + log(1-z_star_mo_post_m[1]) * (1-Z$maleelder[i]) + dnorm(Z$agediff[i], mean = z_star_ad_post_m[1], sd = sqrt(z_ad_v_post_m), log=TRUE)) + 
  pi_mw_post_m[2] * exp(- ((exp(t * beta_post_m[1]) - 1) / beta_post_m[1]) * exp(alpha_post_m[1] + beta_post_m[1] * (st_ageM[i,j] - x_bar) + gamma_star_post_m[2]) + Z$maleelder[i] * log(z_star_mo_post_m[2]) + log(1-z_star_mo_post_m[2]) * (1-Z$maleelder[i]) + dnorm(Z$agediff[i], mean = z_star_ad_post_m[2], sd = sqrt(z_ad_v_post_m), log=TRUE)) + 
  pi_mw_post_m[3] * exp(- ((exp(t * beta_post_m[1]) - 1) / beta_post_m[1]) * exp(alpha_post_m[1] + beta_post_m[1] * (st_ageM[i,j] - x_bar) + gamma_star_post_m[3]) + Z$maleelder[i] * log(z_star_mo_post_m[3]) + log(1-z_star_mo_post_m[3]) * (1-Z$maleelder[i]) + dnorm(Z$agediff[i], mean = z_star_ad_post_m[3], sd = sqrt(z_ad_v_post_m), log=TRUE)) + 
  pi_mw_post_m[4] * exp(- ((exp(t * beta_post_m[1]) - 1) / beta_post_m[1]) * exp(alpha_post_m[1] + beta_post_m[1] * (st_ageM[i,j] - x_bar) + gamma_star_post_m[4]) + Z$maleelder[i] * log(z_star_mo_post_m[4]) + log(1-z_star_mo_post_m[4]) * (1-Z$maleelder[i]) + dnorm(Z$agediff[i], mean = z_star_ad_post_m[4], sd = sqrt(z_ad_v_post_m), log=TRUE)) + 
  pi_mw_post_m[5] * exp(- ((exp(t * beta_post_m[1]) - 1) / beta_post_m[1]) * exp(alpha_post_m[1] + beta_post_m[1] * (st_ageM[i,j] - x_bar) + gamma_star_post_m[5]) + Z$maleelder[i] * log(z_star_mo_post_m[5]) + log(1-z_star_mo_post_m[5]) * (1-Z$maleelder[i]) + dnorm(Z$agediff[i], mean = z_star_ad_post_m[5], sd = sqrt(z_ad_v_post_m), log=TRUE)) + 
  pi_mw_post_m[6] * exp(- ((exp(t * beta_post_m[1]) - 1) / beta_post_m[1]) * exp(alpha_post_m[1] + beta_post_m[1] * (st_ageM[i,j] - x_bar) + gamma_star_post_m[6]) + Z$maleelder[i] * log(z_star_mo_post_m[6]) + log(1-z_star_mo_post_m[6]) * (1-Z$maleelder[i]) + dnorm(Z$agediff[i], mean = z_star_ad_post_m[6], sd = sqrt(z_ad_v_post_m), log=TRUE)) + 
  pi_mw_post_m[7] * exp(- ((exp(t * beta_post_m[1]) - 1) / beta_post_m[1]) * exp(alpha_post_m[1] + beta_post_m[1] * (st_ageM[i,j] - x_bar) + gamma_star_post_m[7]) + Z$maleelder[i] * log(z_star_mo_post_m[7]) + log(1-z_star_mo_post_m[7]) * (1-Z$maleelder[i]) + dnorm(Z$agediff[i], mean = z_star_ad_post_m[7], sd = sqrt(z_ad_v_post_m), log=TRUE)) + 
  pi_mw_post_m[8] * exp(- ((exp(t * beta_post_m[1]) - 1) / beta_post_m[1]) * exp(alpha_post_m[1] + beta_post_m[1] * (st_ageM[i,j] - x_bar) + gamma_star_post_m[8]) + Z$maleelder[i] * log(z_star_mo_post_m[8]) + log(1-z_star_mo_post_m[8]) * (1-Z$maleelder[i]) + dnorm(Z$agediff[i], mean = z_star_ad_post_m[8], sd = sqrt(z_ad_v_post_m), log=TRUE)) + 
  pi_mw_post_m[9] * exp(- ((exp(t * beta_post_m[1]) - 1) / beta_post_m[1]) * exp(alpha_post_m[1] + beta_post_m[1] * (st_ageM[i,j] - x_bar) + gamma_star_post_m[9]) + Z$maleelder[i] * log(z_star_mo_post_m[9]) + log(1-z_star_mo_post_m[9]) * (1-Z$maleelder[i]) + dnorm(Z$agediff[i], mean = z_star_ad_post_m[9], sd = sqrt(z_ad_v_post_m), log=TRUE)) + 
  pi_mw_post_m[10] * exp(- ((exp(t * beta_post_m[1]) - 1) / beta_post_m[1]) * exp(alpha_post_m[1] + beta_post_m[1] * (st_ageM[i,j] - x_bar) + gamma_star_post_m[10]) + Z$maleelder[i] * log(z_star_mo_post_m[10]) + log(1-z_star_mo_post_m[10]) * (1-Z$maleelder[i]) + dnorm(Z$agediff[i], mean = z_star_ad_post_m[10], sd = sqrt(z_ad_v_post_m), log=TRUE)) + 
  pi_mw_post_m[11] * exp(- ((exp(t * beta_post_m[1]) - 1) / beta_post_m[1]) * exp(alpha_post_m[1] + beta_post_m[1] * (st_ageM[i,j] - x_bar) + gamma_star_post_m[11]) + Z$maleelder[i] * log(z_star_mo_post_m[11]) + log(1-z_star_mo_post_m[11]) * (1-Z$maleelder[i]) + dnorm(Z$agediff[i], mean = z_star_ad_post_m[11], sd = sqrt(z_ad_v_post_m), log=TRUE)) + 
  pi_mw_post_m[12] * exp(- ((exp(t * beta_post_m[1]) - 1) / beta_post_m[1]) * exp(alpha_post_m[1] + beta_post_m[1] * (st_ageM[i,j] - x_bar) + gamma_star_post_m[12]) + Z$maleelder[i] * log(z_star_mo_post_m[12]) + log(1-z_star_mo_post_m[12]) * (1-Z$maleelder[i]) + dnorm(Z$agediff[i], mean = z_star_ad_post_m[12], sd = sqrt(z_ad_v_post_m), log=TRUE)) + 
  pi_mw_post_m[13] * exp(- ((exp(t * beta_post_m[1]) - 1) / beta_post_m[1]) * exp(alpha_post_m[1] + beta_post_m[1] * (st_ageM[i,j] - x_bar) + gamma_star_post_m[13]) + Z$maleelder[i] * log(z_star_mo_post_m[13]) + log(1-z_star_mo_post_m[13]) * (1-Z$maleelder[i]) + dnorm(Z$agediff[i], mean = z_star_ad_post_m[13], sd = sqrt(z_ad_v_post_m), log=TRUE)) + 
  pi_mw_post_m[14] * exp(- ((exp(t * beta_post_m[1]) - 1) / beta_post_m[1]) * exp(alpha_post_m[1] + beta_post_m[1] * (st_ageM[i,j] - x_bar) + gamma_star_post_m[14]) + Z$maleelder[i] * log(z_star_mo_post_m[14]) + log(1-z_star_mo_post_m[14]) * (1-Z$maleelder[i]) + dnorm(Z$agediff[i], mean = z_star_ad_post_m[14], sd = sqrt(z_ad_v_post_m), log=TRUE)) + 
  pi_mw_post_m[15] * exp(- ((exp(t * beta_post_m[1]) - 1) / beta_post_m[1]) * exp(alpha_post_m[1] + beta_post_m[1] * (st_ageM[i,j] - x_bar) + gamma_star_post_m[15]) + Z$maleelder[i] * log(z_star_mo_post_m[15]) + log(1-z_star_mo_post_m[15]) * (1-Z$maleelder[i]) + dnorm(Z$agediff[i], mean = z_star_ad_post_m[15], sd = sqrt(z_ad_v_post_m), log=TRUE)) + 
  pi_mw_post_m[16] * exp(- ((exp(t * beta_post_m[1]) - 1) / beta_post_m[1]) * exp(alpha_post_m[1] + beta_post_m[1] * (st_ageM[i,j] - x_bar) + gamma_star_post_m[16]) + Z$maleelder[i] * log(z_star_mo_post_m[16]) + log(1-z_star_mo_post_m[16]) * (1-Z$maleelder[i]) + dnorm(Z$agediff[i], mean = z_star_ad_post_m[16], sd = sqrt(z_ad_v_post_m), log=TRUE)) + 
  pi_mw_post_m[17] * exp(- ((exp(t * beta_post_m[1]) - 1) / beta_post_m[1]) * exp(alpha_post_m[1] + beta_post_m[1] * (st_ageM[i,j] - x_bar) + gamma_star_post_m[17]) + Z$maleelder[i] * log(z_star_mo_post_m[17]) + log(1-z_star_mo_post_m[17]) * (1-Z$maleelder[i]) + dnorm(Z$agediff[i], mean = z_star_ad_post_m[17], sd = sqrt(z_ad_v_post_m), log=TRUE)) + 
  pi_mw_post_m[18] * exp(- ((exp(t * beta_post_m[1]) - 1) / beta_post_m[1]) * exp(alpha_post_m[1] + beta_post_m[1] * (st_ageM[i,j] - x_bar) + gamma_star_post_m[18]) + Z$maleelder[i] * log(z_star_mo_post_m[18]) + log(1-z_star_mo_post_m[18]) * (1-Z$maleelder[i]) + dnorm(Z$agediff[i], mean = z_star_ad_post_m[18], sd = sqrt(z_ad_v_post_m), log=TRUE)) + 
  pi_mw_post_m[19] * exp(- ((exp(t * beta_post_m[1]) - 1) / beta_post_m[1]) * exp(alpha_post_m[1] + beta_post_m[1] * (st_ageM[i,j] - x_bar) + gamma_star_post_m[19]) + Z$maleelder[i] * log(z_star_mo_post_m[19]) + log(1-z_star_mo_post_m[19]) * (1-Z$maleelder[i]) + dnorm(Z$agediff[i], mean = z_star_ad_post_m[19], sd = sqrt(z_ad_v_post_m), log=TRUE)) + 
  pi_mw_post_m[20] * exp(- ((exp(t * beta_post_m[1]) - 1) / beta_post_m[1]) * exp(alpha_post_m[1] + beta_post_m[1] * (st_ageM[i,j] - x_bar) + gamma_star_post_m[20]) + Z$maleelder[i] * log(z_star_mo_post_m[20]) + log(1-z_star_mo_post_m[20]) * (1-Z$maleelder[i]) + dnorm(Z$agediff[i], mean = z_star_ad_post_m[20], sd = sqrt(z_ad_v_post_m), log=TRUE)) + 
  pi_mw_post_m[21] * exp(- ((exp(t * beta_post_m[1]) - 1) / beta_post_m[1]) * exp(alpha_post_m[1] + beta_post_m[1] * (st_ageM[i,j] - x_bar) + gamma_star_post_m[21]) + Z$maleelder[i] * log(z_star_mo_post_m[21]) + log(1-z_star_mo_post_m[21]) * (1-Z$maleelder[i]) + dnorm(Z$agediff[i], mean = z_star_ad_post_m[21], sd = sqrt(z_ad_v_post_m), log=TRUE)) + 
  pi_mw_post_m[22] * exp(- ((exp(t * beta_post_m[1]) - 1) / beta_post_m[1]) * exp(alpha_post_m[1] + beta_post_m[1] * (st_ageM[i,j] - x_bar) + gamma_star_post_m[22]) + Z$maleelder[i] * log(z_star_mo_post_m[22]) + log(1-z_star_mo_post_m[22]) * (1-Z$maleelder[i]) + dnorm(Z$agediff[i], mean = z_star_ad_post_m[22], sd = sqrt(z_ad_v_post_m), log=TRUE)) + 
  pi_mw_post_m[23] * exp(- ((exp(t * beta_post_m[1]) - 1) / beta_post_m[1]) * exp(alpha_post_m[1] + beta_post_m[1] * (st_ageM[i,j] - x_bar) + gamma_star_post_m[23]) + Z$maleelder[i] * log(z_star_mo_post_m[23]) + log(1-z_star_mo_post_m[23]) * (1-Z$maleelder[i]) + dnorm(Z$agediff[i], mean = z_star_ad_post_m[23], sd = sqrt(z_ad_v_post_m), log=TRUE)) + 
  pi_mw_post_m[24] * exp(- ((exp(t * beta_post_m[1]) - 1) / beta_post_m[1]) * exp(alpha_post_m[1] + beta_post_m[1] * (st_ageM[i,j] - x_bar) + gamma_star_post_m[24]) + Z$maleelder[i] * log(z_star_mo_post_m[24]) + log(1-z_star_mo_post_m[24]) * (1-Z$maleelder[i]) + dnorm(Z$agediff[i], mean = z_star_ad_post_m[24], sd = sqrt(z_ad_v_post_m), log=TRUE)) + 
  pi_mw_post_m[25] * exp(- ((exp(t * beta_post_m[1]) - 1) / beta_post_m[1]) * exp(alpha_post_m[1] + beta_post_m[1] * (st_ageM[i,j] - x_bar) + gamma_star_post_m[25]) + Z$maleelder[i] * log(z_star_mo_post_m[25]) + log(1-z_star_mo_post_m[25]) * (1-Z$maleelder[i]) + dnorm(Z$agediff[i], mean = z_star_ad_post_m[25], sd = sqrt(z_ad_v_post_m), log=TRUE)))
  
}

# - Taking expectation over draws instead of posterior mean

integrand_f <- function(t){
  
  (pi_mw_post_m[1] * exp(- ((exp(t * beta_post_m[2]) - 1) / beta_post_m[2]) * exp(alpha_post_m[2] + beta_post_m[2] * (st_ageM[i,j] - x_bar) + gamma_star_post_m[1+1]) + (alpha_post_m[2] + beta_post_m[2] * (st_ageM[i,j] - x_bar + t) + gamma_star_post_m[1+1]) + Z$maleelder[i] * log(z_star_mo_post_m[1]) + log(1-z_star_mo_post_m[1]) * (1-Z$maleelder[i]) + dnorm(Z$agediff[i], mean = z_star_ad_post_m[1], sd = sqrt(z_ad_v_post_m), log=TRUE)) + 
   pi_mw_post_m[2] * exp(- ((exp(t * beta_post_m[2]) - 1) / beta_post_m[2]) * exp(alpha_post_m[2] + beta_post_m[2] * (st_ageM[i,j] - x_bar) + gamma_star_post_m[1+2]) + (alpha_post_m[2] + beta_post_m[2] * (st_ageM[i,j] - x_bar + t) + gamma_star_post_m[1+2]) + Z$maleelder[i] * log(z_star_mo_post_m[2]) + log(1-z_star_mo_post_m[2]) * (1-Z$maleelder[i]) + dnorm(Z$agediff[i], mean = z_star_ad_post_m[2], sd = sqrt(z_ad_v_post_m), log=TRUE)) + 
   pi_mw_post_m[3] * exp(- ((exp(t * beta_post_m[2]) - 1) / beta_post_m[2]) * exp(alpha_post_m[2] + beta_post_m[2] * (st_ageM[i,j] - x_bar) + gamma_star_post_m[1+3]) + (alpha_post_m[2] + beta_post_m[2] * (st_ageM[i,j] - x_bar + t) + gamma_star_post_m[1+3]) + Z$maleelder[i] * log(z_star_mo_post_m[3]) + log(1-z_star_mo_post_m[3]) * (1-Z$maleelder[i]) + dnorm(Z$agediff[i], mean = z_star_ad_post_m[3], sd = sqrt(z_ad_v_post_m), log=TRUE)) + 
   pi_mw_post_m[4] * exp(- ((exp(t * beta_post_m[2]) - 1) / beta_post_m[2]) * exp(alpha_post_m[2] + beta_post_m[2] * (st_ageM[i,j] - x_bar) + gamma_star_post_m[1+4]) + (alpha_post_m[2] + beta_post_m[2] * (st_ageM[i,j] - x_bar + t) + gamma_star_post_m[1+4]) + Z$maleelder[i] * log(z_star_mo_post_m[4]) + log(1-z_star_mo_post_m[4]) * (1-Z$maleelder[i]) + dnorm(Z$agediff[i], mean = z_star_ad_post_m[4], sd = sqrt(z_ad_v_post_m), log=TRUE)) + 
   pi_mw_post_m[5] * exp(- ((exp(t * beta_post_m[2]) - 1) / beta_post_m[2]) * exp(alpha_post_m[2] + beta_post_m[2] * (st_ageM[i,j] - x_bar) + gamma_star_post_m[1+5]) + (alpha_post_m[2] + beta_post_m[2] * (st_ageM[i,j] - x_bar + t) + gamma_star_post_m[1+5]) + Z$maleelder[i] * log(z_star_mo_post_m[5]) + log(1-z_star_mo_post_m[5]) * (1-Z$maleelder[i]) + dnorm(Z$agediff[i], mean = z_star_ad_post_m[5], sd = sqrt(z_ad_v_post_m), log=TRUE)) + 
   pi_mw_post_m[6] * exp(- ((exp(t * beta_post_m[2]) - 1) / beta_post_m[2]) * exp(alpha_post_m[2] + beta_post_m[2] * (st_ageM[i,j] - x_bar) + gamma_star_post_m[1+6]) + (alpha_post_m[2] + beta_post_m[2] * (st_ageM[i,j] - x_bar + t) + gamma_star_post_m[1+6]) + Z$maleelder[i] * log(z_star_mo_post_m[6]) + log(1-z_star_mo_post_m[6]) * (1-Z$maleelder[i]) + dnorm(Z$agediff[i], mean = z_star_ad_post_m[6], sd = sqrt(z_ad_v_post_m), log=TRUE)) + 
   pi_mw_post_m[7] * exp(- ((exp(t * beta_post_m[2]) - 1) / beta_post_m[2]) * exp(alpha_post_m[2] + beta_post_m[2] * (st_ageM[i,j] - x_bar) + gamma_star_post_m[1+7]) + (alpha_post_m[2] + beta_post_m[2] * (st_ageM[i,j] - x_bar + t) + gamma_star_post_m[1+7]) + Z$maleelder[i] * log(z_star_mo_post_m[7]) + log(1-z_star_mo_post_m[7]) * (1-Z$maleelder[i]) + dnorm(Z$agediff[i], mean = z_star_ad_post_m[7], sd = sqrt(z_ad_v_post_m), log=TRUE)) + 
   pi_mw_post_m[8] * exp(- ((exp(t * beta_post_m[2]) - 1) / beta_post_m[2]) * exp(alpha_post_m[2] + beta_post_m[2] * (st_ageM[i,j] - x_bar) + gamma_star_post_m[1+8]) + (alpha_post_m[2] + beta_post_m[2] * (st_ageM[i,j] - x_bar + t) + gamma_star_post_m[1+8]) + Z$maleelder[i] * log(z_star_mo_post_m[8]) + log(1-z_star_mo_post_m[8]) * (1-Z$maleelder[i]) + dnorm(Z$agediff[i], mean = z_star_ad_post_m[8], sd = sqrt(z_ad_v_post_m), log=TRUE)) + 
   pi_mw_post_m[9] * exp(- ((exp(t * beta_post_m[2]) - 1) / beta_post_m[2]) * exp(alpha_post_m[2] + beta_post_m[2] * (st_ageM[i,j] - x_bar) + gamma_star_post_m[1+9]) + (alpha_post_m[2] + beta_post_m[2] * (st_ageM[i,j] - x_bar + t) + gamma_star_post_m[1+9]) + Z$maleelder[i] * log(z_star_mo_post_m[9]) + log(1-z_star_mo_post_m[9]) * (1-Z$maleelder[i]) + dnorm(Z$agediff[i], mean = z_star_ad_post_m[9], sd = sqrt(z_ad_v_post_m), log=TRUE)) + 
   pi_mw_post_m[10] * exp(- ((exp(t * beta_post_m[2]) - 1) / beta_post_m[2]) * exp(alpha_post_m[2] + beta_post_m[2] * (st_ageM[i,j] - x_bar) + gamma_star_post_m[1+10]) + (alpha_post_m[2] + beta_post_m[2] * (st_ageM[i,j] - x_bar + t) + gamma_star_post_m[1+10]) + Z$maleelder[i] * log(z_star_mo_post_m[10]) + log(1-z_star_mo_post_m[10]) * (1-Z$maleelder[i]) + dnorm(Z$agediff[i], mean = z_star_ad_post_m[10], sd = sqrt(z_ad_v_post_m), log=TRUE)) + 
   pi_mw_post_m[11] * exp(- ((exp(t * beta_post_m[2]) - 1) / beta_post_m[2]) * exp(alpha_post_m[2] + beta_post_m[2] * (st_ageM[i,j] - x_bar) + gamma_star_post_m[1+11]) + (alpha_post_m[2] + beta_post_m[2] * (st_ageM[i,j] - x_bar + t) + gamma_star_post_m[1+11]) + Z$maleelder[i] * log(z_star_mo_post_m[11]) + log(1-z_star_mo_post_m[11]) * (1-Z$maleelder[i]) + dnorm(Z$agediff[i], mean = z_star_ad_post_m[11], sd = sqrt(z_ad_v_post_m), log=TRUE)) + 
   pi_mw_post_m[12] * exp(- ((exp(t * beta_post_m[2]) - 1) / beta_post_m[2]) * exp(alpha_post_m[2] + beta_post_m[2] * (st_ageM[i,j] - x_bar) + gamma_star_post_m[1+12]) + (alpha_post_m[2] + beta_post_m[2] * (st_ageM[i,j] - x_bar + t) + gamma_star_post_m[1+12]) + Z$maleelder[i] * log(z_star_mo_post_m[12]) + log(1-z_star_mo_post_m[12]) * (1-Z$maleelder[i]) + dnorm(Z$agediff[i], mean = z_star_ad_post_m[12], sd = sqrt(z_ad_v_post_m), log=TRUE)) + 
   pi_mw_post_m[13] * exp(- ((exp(t * beta_post_m[2]) - 1) / beta_post_m[2]) * exp(alpha_post_m[2] + beta_post_m[2] * (st_ageM[i,j] - x_bar) + gamma_star_post_m[1+13]) + (alpha_post_m[2] + beta_post_m[2] * (st_ageM[i,j] - x_bar + t) + gamma_star_post_m[1+13]) + Z$maleelder[i] * log(z_star_mo_post_m[13]) + log(1-z_star_mo_post_m[13]) * (1-Z$maleelder[i]) + dnorm(Z$agediff[i], mean = z_star_ad_post_m[13], sd = sqrt(z_ad_v_post_m), log=TRUE)) + 
   pi_mw_post_m[14] * exp(- ((exp(t * beta_post_m[2]) - 1) / beta_post_m[2]) * exp(alpha_post_m[2] + beta_post_m[2] * (st_ageM[i,j] - x_bar) + gamma_star_post_m[1+14]) + (alpha_post_m[2] + beta_post_m[2] * (st_ageM[i,j] - x_bar + t) + gamma_star_post_m[1+14]) + Z$maleelder[i] * log(z_star_mo_post_m[14]) + log(1-z_star_mo_post_m[14]) * (1-Z$maleelder[i]) + dnorm(Z$agediff[i], mean = z_star_ad_post_m[14], sd = sqrt(z_ad_v_post_m), log=TRUE)) + 
   pi_mw_post_m[15] * exp(- ((exp(t * beta_post_m[2]) - 1) / beta_post_m[2]) * exp(alpha_post_m[2] + beta_post_m[2] * (st_ageM[i,j] - x_bar) + gamma_star_post_m[1+15]) + (alpha_post_m[2] + beta_post_m[2] * (st_ageM[i,j] - x_bar + t) + gamma_star_post_m[1+15]) + Z$maleelder[i] * log(z_star_mo_post_m[15]) + log(1-z_star_mo_post_m[15]) * (1-Z$maleelder[i]) + dnorm(Z$agediff[i], mean = z_star_ad_post_m[15], sd = sqrt(z_ad_v_post_m), log=TRUE)) + 
   pi_mw_post_m[16] * exp(- ((exp(t * beta_post_m[2]) - 1) / beta_post_m[2]) * exp(alpha_post_m[2] + beta_post_m[2] * (st_ageM[i,j] - x_bar) + gamma_star_post_m[1+16]) + (alpha_post_m[2] + beta_post_m[2] * (st_ageM[i,j] - x_bar + t) + gamma_star_post_m[1+16]) + Z$maleelder[i] * log(z_star_mo_post_m[16]) + log(1-z_star_mo_post_m[16]) * (1-Z$maleelder[i]) + dnorm(Z$agediff[i], mean = z_star_ad_post_m[16], sd = sqrt(z_ad_v_post_m), log=TRUE)) + 
   pi_mw_post_m[17] * exp(- ((exp(t * beta_post_m[2]) - 1) / beta_post_m[2]) * exp(alpha_post_m[2] + beta_post_m[2] * (st_ageM[i,j] - x_bar) + gamma_star_post_m[1+17]) + (alpha_post_m[2] + beta_post_m[2] * (st_ageM[i,j] - x_bar + t) + gamma_star_post_m[1+17]) + Z$maleelder[i] * log(z_star_mo_post_m[17]) + log(1-z_star_mo_post_m[17]) * (1-Z$maleelder[i]) + dnorm(Z$agediff[i], mean = z_star_ad_post_m[17], sd = sqrt(z_ad_v_post_m), log=TRUE)) + 
   pi_mw_post_m[18] * exp(- ((exp(t * beta_post_m[2]) - 1) / beta_post_m[2]) * exp(alpha_post_m[2] + beta_post_m[2] * (st_ageM[i,j] - x_bar) + gamma_star_post_m[1+18]) + (alpha_post_m[2] + beta_post_m[2] * (st_ageM[i,j] - x_bar + t) + gamma_star_post_m[1+18]) + Z$maleelder[i] * log(z_star_mo_post_m[18]) + log(1-z_star_mo_post_m[18]) * (1-Z$maleelder[i]) + dnorm(Z$agediff[i], mean = z_star_ad_post_m[18], sd = sqrt(z_ad_v_post_m), log=TRUE)) + 
   pi_mw_post_m[19] * exp(- ((exp(t * beta_post_m[2]) - 1) / beta_post_m[2]) * exp(alpha_post_m[2] + beta_post_m[2] * (st_ageM[i,j] - x_bar) + gamma_star_post_m[1+19]) + (alpha_post_m[2] + beta_post_m[2] * (st_ageM[i,j] - x_bar + t) + gamma_star_post_m[1+19]) + Z$maleelder[i] * log(z_star_mo_post_m[19]) + log(1-z_star_mo_post_m[19]) * (1-Z$maleelder[i]) + dnorm(Z$agediff[i], mean = z_star_ad_post_m[19], sd = sqrt(z_ad_v_post_m), log=TRUE)) + 
   pi_mw_post_m[20] * exp(- ((exp(t * beta_post_m[2]) - 1) / beta_post_m[2]) * exp(alpha_post_m[2] + beta_post_m[2] * (st_ageM[i,j] - x_bar) + gamma_star_post_m[1+20]) + (alpha_post_m[2] + beta_post_m[2] * (st_ageM[i,j] - x_bar + t) + gamma_star_post_m[1+20]) + Z$maleelder[i] * log(z_star_mo_post_m[20]) + log(1-z_star_mo_post_m[20]) * (1-Z$maleelder[i]) + dnorm(Z$agediff[i], mean = z_star_ad_post_m[20], sd = sqrt(z_ad_v_post_m), log=TRUE)) + 
   pi_mw_post_m[21] * exp(- ((exp(t * beta_post_m[2]) - 1) / beta_post_m[2]) * exp(alpha_post_m[2] + beta_post_m[2] * (st_ageM[i,j] - x_bar) + gamma_star_post_m[1+21]) + (alpha_post_m[2] + beta_post_m[2] * (st_ageM[i,j] - x_bar + t) + gamma_star_post_m[1+21]) + Z$maleelder[i] * log(z_star_mo_post_m[21]) + log(1-z_star_mo_post_m[21]) * (1-Z$maleelder[i]) + dnorm(Z$agediff[i], mean = z_star_ad_post_m[21], sd = sqrt(z_ad_v_post_m), log=TRUE)) + 
   pi_mw_post_m[22] * exp(- ((exp(t * beta_post_m[2]) - 1) / beta_post_m[2]) * exp(alpha_post_m[2] + beta_post_m[2] * (st_ageM[i,j] - x_bar) + gamma_star_post_m[1+22]) + (alpha_post_m[2] + beta_post_m[2] * (st_ageM[i,j] - x_bar + t) + gamma_star_post_m[1+22]) + Z$maleelder[i] * log(z_star_mo_post_m[22]) + log(1-z_star_mo_post_m[22]) * (1-Z$maleelder[i]) + dnorm(Z$agediff[i], mean = z_star_ad_post_m[22], sd = sqrt(z_ad_v_post_m), log=TRUE)) + 
   pi_mw_post_m[23] * exp(- ((exp(t * beta_post_m[2]) - 1) / beta_post_m[2]) * exp(alpha_post_m[2] + beta_post_m[2] * (st_ageM[i,j] - x_bar) + gamma_star_post_m[1+23]) + (alpha_post_m[2] + beta_post_m[2] * (st_ageM[i,j] - x_bar + t) + gamma_star_post_m[1+23]) + Z$maleelder[i] * log(z_star_mo_post_m[23]) + log(1-z_star_mo_post_m[23]) * (1-Z$maleelder[i]) + dnorm(Z$agediff[i], mean = z_star_ad_post_m[23], sd = sqrt(z_ad_v_post_m), log=TRUE)) + 
   pi_mw_post_m[24] * exp(- ((exp(t * beta_post_m[2]) - 1) / beta_post_m[2]) * exp(alpha_post_m[2] + beta_post_m[2] * (st_ageM[i,j] - x_bar) + gamma_star_post_m[1+24]) + (alpha_post_m[2] + beta_post_m[2] * (st_ageM[i,j] - x_bar + t) + gamma_star_post_m[1+24]) + Z$maleelder[i] * log(z_star_mo_post_m[24]) + log(1-z_star_mo_post_m[24]) * (1-Z$maleelder[i]) + dnorm(Z$agediff[i], mean = z_star_ad_post_m[24], sd = sqrt(z_ad_v_post_m), log=TRUE)) + 
   pi_mw_post_m[25] * exp(- ((exp(t * beta_post_m[2]) - 1) / beta_post_m[2]) * exp(alpha_post_m[2] + beta_post_m[2] * (st_ageM[i,j] - x_bar) + gamma_star_post_m[1+25]) + (alpha_post_m[2] + beta_post_m[2] * (st_ageM[i,j] - x_bar + t) + gamma_star_post_m[1+25]) + Z$maleelder[i] * log(z_star_mo_post_m[25]) + log(1-z_star_mo_post_m[25]) * (1-Z$maleelder[i]) + dnorm(Z$agediff[i], mean = z_star_ad_post_m[25], sd = sqrt(z_ad_v_post_m), log=TRUE)) ) / 
    (pi_mw_post_m[1] * exp(- ((exp(t * beta_post_m[2]) - 1) / beta_post_m[2]) * exp(alpha_post_m[2] + beta_post_m[2] * (st_ageM[i,j] - x_bar) + gamma_star_post_m[1+1]) + Z$maleelder[i] * log(z_star_mo_post_m[1]) + log(1-z_star_mo_post_m[1]) * (1-Z$maleelder[i]) + dnorm(Z$agediff[i], mean = z_star_ad_post_m[1], sd = sqrt(z_ad_v_post_m), log=TRUE)) + 
       pi_mw_post_m[2] * exp(- ((exp(t * beta_post_m[2]) - 1) / beta_post_m[2]) * exp(alpha_post_m[2] + beta_post_m[2] * (st_ageM[i,j] - x_bar) + gamma_star_post_m[1+2]) + Z$maleelder[i] * log(z_star_mo_post_m[2]) + log(1-z_star_mo_post_m[2]) * (1-Z$maleelder[i]) + dnorm(Z$agediff[i], mean = z_star_ad_post_m[2], sd = sqrt(z_ad_v_post_m), log=TRUE)) + 
       pi_mw_post_m[3] * exp(- ((exp(t * beta_post_m[2]) - 1) / beta_post_m[2]) * exp(alpha_post_m[2] + beta_post_m[2] * (st_ageM[i,j] - x_bar) + gamma_star_post_m[1+3]) + Z$maleelder[i] * log(z_star_mo_post_m[3]) + log(1-z_star_mo_post_m[3]) * (1-Z$maleelder[i]) + dnorm(Z$agediff[i], mean = z_star_ad_post_m[3], sd = sqrt(z_ad_v_post_m), log=TRUE)) + 
       pi_mw_post_m[4] * exp(- ((exp(t * beta_post_m[2]) - 1) / beta_post_m[2]) * exp(alpha_post_m[2] + beta_post_m[2] * (st_ageM[i,j] - x_bar) + gamma_star_post_m[1+4]) + Z$maleelder[i] * log(z_star_mo_post_m[4]) + log(1-z_star_mo_post_m[4]) * (1-Z$maleelder[i]) + dnorm(Z$agediff[i], mean = z_star_ad_post_m[4], sd = sqrt(z_ad_v_post_m), log=TRUE)) + 
       pi_mw_post_m[5] * exp(- ((exp(t * beta_post_m[2]) - 1) / beta_post_m[2]) * exp(alpha_post_m[2] + beta_post_m[2] * (st_ageM[i,j] - x_bar) + gamma_star_post_m[1+5]) + Z$maleelder[i] * log(z_star_mo_post_m[5]) + log(1-z_star_mo_post_m[5]) * (1-Z$maleelder[i]) + dnorm(Z$agediff[i], mean = z_star_ad_post_m[5], sd = sqrt(z_ad_v_post_m), log=TRUE)) + 
       pi_mw_post_m[6] * exp(- ((exp(t * beta_post_m[2]) - 1) / beta_post_m[2]) * exp(alpha_post_m[2] + beta_post_m[2] * (st_ageM[i,j] - x_bar) + gamma_star_post_m[1+6]) + Z$maleelder[i] * log(z_star_mo_post_m[6]) + log(1-z_star_mo_post_m[6]) * (1-Z$maleelder[i]) + dnorm(Z$agediff[i], mean = z_star_ad_post_m[6], sd = sqrt(z_ad_v_post_m), log=TRUE)) + 
       pi_mw_post_m[7] * exp(- ((exp(t * beta_post_m[2]) - 1) / beta_post_m[2]) * exp(alpha_post_m[2] + beta_post_m[2] * (st_ageM[i,j] - x_bar) + gamma_star_post_m[1+7]) + Z$maleelder[i] * log(z_star_mo_post_m[7]) + log(1-z_star_mo_post_m[7]) * (1-Z$maleelder[i]) + dnorm(Z$agediff[i], mean = z_star_ad_post_m[7], sd = sqrt(z_ad_v_post_m), log=TRUE)) + 
       pi_mw_post_m[8] * exp(- ((exp(t * beta_post_m[2]) - 1) / beta_post_m[2]) * exp(alpha_post_m[2] + beta_post_m[2] * (st_ageM[i,j] - x_bar) + gamma_star_post_m[1+8]) + Z$maleelder[i] * log(z_star_mo_post_m[8]) + log(1-z_star_mo_post_m[8]) * (1-Z$maleelder[i]) + dnorm(Z$agediff[i], mean = z_star_ad_post_m[8], sd = sqrt(z_ad_v_post_m), log=TRUE)) + 
       pi_mw_post_m[9] * exp(- ((exp(t * beta_post_m[2]) - 1) / beta_post_m[2]) * exp(alpha_post_m[2] + beta_post_m[2] * (st_ageM[i,j] - x_bar) + gamma_star_post_m[1+9]) + Z$maleelder[i] * log(z_star_mo_post_m[9]) + log(1-z_star_mo_post_m[9]) * (1-Z$maleelder[i]) + dnorm(Z$agediff[i], mean = z_star_ad_post_m[9], sd = sqrt(z_ad_v_post_m), log=TRUE)) + 
       pi_mw_post_m[10] * exp(- ((exp(t * beta_post_m[2]) - 1) / beta_post_m[2]) * exp(alpha_post_m[2] + beta_post_m[2] * (st_ageM[i,j] - x_bar) + gamma_star_post_m[1+10]) + Z$maleelder[i] * log(z_star_mo_post_m[10]) + log(1-z_star_mo_post_m[10]) * (1-Z$maleelder[i]) + dnorm(Z$agediff[i], mean = z_star_ad_post_m[10], sd = sqrt(z_ad_v_post_m), log=TRUE)) + 
       pi_mw_post_m[11] * exp(- ((exp(t * beta_post_m[2]) - 1) / beta_post_m[2]) * exp(alpha_post_m[2] + beta_post_m[2] * (st_ageM[i,j] - x_bar) + gamma_star_post_m[1+11]) + Z$maleelder[i] * log(z_star_mo_post_m[11]) + log(1-z_star_mo_post_m[11]) * (1-Z$maleelder[i]) + dnorm(Z$agediff[i], mean = z_star_ad_post_m[11], sd = sqrt(z_ad_v_post_m), log=TRUE)) + 
       pi_mw_post_m[12] * exp(- ((exp(t * beta_post_m[2]) - 1) / beta_post_m[2]) * exp(alpha_post_m[2] + beta_post_m[2] * (st_ageM[i,j] - x_bar) + gamma_star_post_m[1+12]) + Z$maleelder[i] * log(z_star_mo_post_m[12]) + log(1-z_star_mo_post_m[12]) * (1-Z$maleelder[i]) + dnorm(Z$agediff[i], mean = z_star_ad_post_m[12], sd = sqrt(z_ad_v_post_m), log=TRUE)) + 
       pi_mw_post_m[13] * exp(- ((exp(t * beta_post_m[2]) - 1) / beta_post_m[2]) * exp(alpha_post_m[2] + beta_post_m[2] * (st_ageM[i,j] - x_bar) + gamma_star_post_m[1+13]) + Z$maleelder[i] * log(z_star_mo_post_m[13]) + log(1-z_star_mo_post_m[13]) * (1-Z$maleelder[i]) + dnorm(Z$agediff[i], mean = z_star_ad_post_m[13], sd = sqrt(z_ad_v_post_m), log=TRUE)) + 
       pi_mw_post_m[14] * exp(- ((exp(t * beta_post_m[2]) - 1) / beta_post_m[2]) * exp(alpha_post_m[2] + beta_post_m[2] * (st_ageM[i,j] - x_bar) + gamma_star_post_m[1+14]) + Z$maleelder[i] * log(z_star_mo_post_m[14]) + log(1-z_star_mo_post_m[14]) * (1-Z$maleelder[i]) + dnorm(Z$agediff[i], mean = z_star_ad_post_m[14], sd = sqrt(z_ad_v_post_m), log=TRUE)) + 
       pi_mw_post_m[15] * exp(- ((exp(t * beta_post_m[2]) - 1) / beta_post_m[2]) * exp(alpha_post_m[2] + beta_post_m[2] * (st_ageM[i,j] - x_bar) + gamma_star_post_m[1+15]) + Z$maleelder[i] * log(z_star_mo_post_m[15]) + log(1-z_star_mo_post_m[15]) * (1-Z$maleelder[i]) + dnorm(Z$agediff[i], mean = z_star_ad_post_m[15], sd = sqrt(z_ad_v_post_m), log=TRUE)) + 
       pi_mw_post_m[16] * exp(- ((exp(t * beta_post_m[2]) - 1) / beta_post_m[2]) * exp(alpha_post_m[2] + beta_post_m[2] * (st_ageM[i,j] - x_bar) + gamma_star_post_m[1+16]) + Z$maleelder[i] * log(z_star_mo_post_m[16]) + log(1-z_star_mo_post_m[16]) * (1-Z$maleelder[i]) + dnorm(Z$agediff[i], mean = z_star_ad_post_m[16], sd = sqrt(z_ad_v_post_m), log=TRUE)) + 
       pi_mw_post_m[17] * exp(- ((exp(t * beta_post_m[2]) - 1) / beta_post_m[2]) * exp(alpha_post_m[2] + beta_post_m[2] * (st_ageM[i,j] - x_bar) + gamma_star_post_m[1+17]) + Z$maleelder[i] * log(z_star_mo_post_m[17]) + log(1-z_star_mo_post_m[17]) * (1-Z$maleelder[i]) + dnorm(Z$agediff[i], mean = z_star_ad_post_m[17], sd = sqrt(z_ad_v_post_m), log=TRUE)) + 
       pi_mw_post_m[18] * exp(- ((exp(t * beta_post_m[2]) - 1) / beta_post_m[2]) * exp(alpha_post_m[2] + beta_post_m[2] * (st_ageM[i,j] - x_bar) + gamma_star_post_m[1+18]) + Z$maleelder[i] * log(z_star_mo_post_m[18]) + log(1-z_star_mo_post_m[18]) * (1-Z$maleelder[i]) + dnorm(Z$agediff[i], mean = z_star_ad_post_m[18], sd = sqrt(z_ad_v_post_m), log=TRUE)) + 
       pi_mw_post_m[19] * exp(- ((exp(t * beta_post_m[2]) - 1) / beta_post_m[2]) * exp(alpha_post_m[2] + beta_post_m[2] * (st_ageM[i,j] - x_bar) + gamma_star_post_m[1+19]) + Z$maleelder[i] * log(z_star_mo_post_m[19]) + log(1-z_star_mo_post_m[19]) * (1-Z$maleelder[i]) + dnorm(Z$agediff[i], mean = z_star_ad_post_m[19], sd = sqrt(z_ad_v_post_m), log=TRUE)) + 
       pi_mw_post_m[20] * exp(- ((exp(t * beta_post_m[2]) - 1) / beta_post_m[2]) * exp(alpha_post_m[2] + beta_post_m[2] * (st_ageM[i,j] - x_bar) + gamma_star_post_m[1+20]) + Z$maleelder[i] * log(z_star_mo_post_m[20]) + log(1-z_star_mo_post_m[20]) * (1-Z$maleelder[i]) + dnorm(Z$agediff[i], mean = z_star_ad_post_m[20], sd = sqrt(z_ad_v_post_m), log=TRUE)) + 
       pi_mw_post_m[21] * exp(- ((exp(t * beta_post_m[2]) - 1) / beta_post_m[2]) * exp(alpha_post_m[2] + beta_post_m[2] * (st_ageM[i,j] - x_bar) + gamma_star_post_m[1+21]) + Z$maleelder[i] * log(z_star_mo_post_m[21]) + log(1-z_star_mo_post_m[21]) * (1-Z$maleelder[i]) + dnorm(Z$agediff[i], mean = z_star_ad_post_m[21], sd = sqrt(z_ad_v_post_m), log=TRUE)) + 
       pi_mw_post_m[22] * exp(- ((exp(t * beta_post_m[2]) - 1) / beta_post_m[2]) * exp(alpha_post_m[2] + beta_post_m[2] * (st_ageM[i,j] - x_bar) + gamma_star_post_m[1+22]) + Z$maleelder[i] * log(z_star_mo_post_m[22]) + log(1-z_star_mo_post_m[22]) * (1-Z$maleelder[i]) + dnorm(Z$agediff[i], mean = z_star_ad_post_m[22], sd = sqrt(z_ad_v_post_m), log=TRUE)) + 
       pi_mw_post_m[23] * exp(- ((exp(t * beta_post_m[2]) - 1) / beta_post_m[2]) * exp(alpha_post_m[2] + beta_post_m[2] * (st_ageM[i,j] - x_bar) + gamma_star_post_m[1+23]) + Z$maleelder[i] * log(z_star_mo_post_m[23]) + log(1-z_star_mo_post_m[23]) * (1-Z$maleelder[i]) + dnorm(Z$agediff[i], mean = z_star_ad_post_m[23], sd = sqrt(z_ad_v_post_m), log=TRUE)) + 
       pi_mw_post_m[24] * exp(- ((exp(t * beta_post_m[2]) - 1) / beta_post_m[2]) * exp(alpha_post_m[2] + beta_post_m[2] * (st_ageM[i,j] - x_bar) + gamma_star_post_m[1+24]) + Z$maleelder[i] * log(z_star_mo_post_m[24]) + log(1-z_star_mo_post_m[24]) * (1-Z$maleelder[i]) + dnorm(Z$agediff[i], mean = z_star_ad_post_m[24], sd = sqrt(z_ad_v_post_m), log=TRUE)) + 
       pi_mw_post_m[25] * exp(- ((exp(t * beta_post_m[2]) - 1) / beta_post_m[2]) * exp(alpha_post_m[2] + beta_post_m[2] * (st_ageM[i,j] - x_bar) + gamma_star_post_m[1+25]) + Z$maleelder[i] * log(z_star_mo_post_m[25]) + log(1-z_star_mo_post_m[25]) * (1-Z$maleelder[i]) + dnorm(Z$agediff[i], mean = z_star_ad_post_m[25], sd = sqrt(z_ad_v_post_m), log=TRUE)))
  
}


## - males
Int_hf_AVDPM_m <- Int_hf_Gom_m

i <- 1
for(i in i:nrow(E_cxM)){
  for(j in 1:ncol(E_cxM)){
Int_hf_AVDPM_m[i,j] <- ifelse(E_cxM[i,j]==0, 0, integrate(integrand_m, lower=(st_ageM[i,j]- (age_baseline + j - 1)), upper=E_cxM[i,j])$value)
}
}

Dev_Gom_m <- 2 * (deaths_m * (log(deaths_m) - log(colSums(Int_hf_Gom_m))) - (deaths_m - colSums(Int_hf_Gom_m)))
Dev_AVDPM_m <- 2 * (deaths_m * (log(deaths_m) - log(colSums(Int_hf_AVDPM_m))) - (deaths_m - colSums(Int_hf_AVDPM_m)))

pdr_Gom_m <- sign(deaths_m - colSums(Int_hf_Gom_m)) * sqrt(Dev_Gom_m)
pdr_AVDPM_m <- sign(deaths_m - colSums(Int_hf_AVDPM_m)) * sqrt(Dev_AVDPM_m)

plot(c(age_baseline:(max_age-1)), pdr_Gom_m)
abline(h=2)
abline(h=-2)
points(c(age_baseline:(max_age-1)), pdr_AVDPM_m, pch=16)

chi2_Gom_m <- sum(na.omit(pdr_Gom_m)^2)
chi2_AVDPM_m <- sum(na.omit(pdr_AVDPM_m)^2)

pearson_Gom_m <- (deaths_m - colSums(Int_hf_Gom_m)) / sqrt(colSums(Int_hf_Gom_m))
pearson_AVDPM_m <- (deaths_m - colSums(Int_hf_AVDPM_m)) / sqrt(colSums(Int_hf_AVDPM_m))

qqnorm(pearson_Gom_m, pch = 1, frame = FALSE)
qqline(pearson_Gom_m, col = "steelblue", lwd = 2)
qqnorm(pearson_AVDPM_m, pch = 1, frame = FALSE)
qqline(pearson_AVDPM_m, col = "steelblue", lwd = 2)

## - females
deaths_f <- colSums(d_xF)

## - Analysis of the integrated hazard function
### - Gompertz
Int_hf_Gom_f <- (exp(beta_BG_f * E_cxF) - 1) / beta_BG_f

for(row in 1:nrow(E_cxF)){
  for(col in 1:ncol(E_cxF)){
    Int_hf_Gom_f[row,col] <- Int_hf_Gom_f[row,col] * exp(alpha_BG_f + beta_BG_f * (st_ageF[row,col]-x_bar))
  }}

Int_hf_AVDPM_f <- Int_hf_Gom_f

i <- 1
for(i in i:nrow(E_cxF)){
  for(j in 1:ncol(E_cxF)){
    Int_hf_AVDPM_f[i,j] <- ifelse(E_cxF[i,j]==0, 0, integrate(integrand_f, lower=(st_ageF[i,j]- (age_baseline + j - 1)), upper=E_cxF[i,j])$value)
  }
}


Dev_Gom_f <- 2 * (deaths_f * (log(deaths_f) - log(colSums(Int_hf_Gom_f))) - (deaths_f - colSums(Int_hf_Gom_f)))
Dev_AVDPM_f <- 2 * (deaths_f * (log(deaths_f) - log(colSums(Int_hf_AVDPM_f))) - (deaths_f - colSums(Int_hf_AVDPM_f)))

pdr_Gom_f <- sign(deaths_f - colSums(Int_hf_Gom_f)) * sqrt(Dev_Gom_f)
pdr_AVDPM_f <- sign(deaths_f - colSums(Int_hf_AVDPM_f)) * sqrt(Dev_AVDPM_f)

plot(c(age_baseline:(max_age - 1)), pdr_Gom_f, ylim=c(-4,4))
abline(h=2)
abline(h=-2)
points(c(age_baseline:(max_age - 1)), pdr_AVDPM_f, pch=16)

chi2_Gom_f <- sum(na.omit(pdr_Gom_f)^2)
chi2_AVDPM_f <- sum(na.omit(pdr_AVDPM_f)^2)

pearson_Gom_f <- (deaths_f - colSums(Int_hf_Gom_f)) / sqrt(colSums(Int_hf_Gom_f))
pearson_AVDPM_f <- (deaths_f - colSums(Int_hf_AVDPM_f)) / sqrt(colSums(Int_hf_AVDPM_f))

qqnorm(pearson_Gom_f, pch = 1, frame = FALSE)
qqline(pearson_Gom_f, col = "steelblue", lwd = 2)
qqnorm(pearson_AVDPM_f, pch = 1, frame = FALSE)
qqline(pearson_AVDPM_f, col = "steelblue", lwd = 2)


# - Calculate mean squared error of the hazard rates

## - In terms of integrated hazard functions, to compare with the number of deaths (advantage to weight by number of exposures)
RMSE_Gom_m <- mean((deaths_m - colSums(Int_hf_Gom_m))^2, na.rm=TRUE)
RMSE_AVDPM_m <- mean((deaths_m - colSums(Int_hf_AVDPM_m))^2, na.rm=TRUE)

RMSE_Gom_f <- mean((deaths_f - colSums(Int_hf_Gom_f))^2, na.rm=TRUE)
RMSE_AVDPM_f <- mean((deaths_f - colSums(Int_hf_AVDPM_f))^2, na.rm=TRUE)

R_RMSE_Gom_m_i <- R_RMSE_Gom_m <- rep(0,max_age - age_baseline)
R_RMSE_Gom_f_i <- R_RMSE_Gom_f <- rep(0,max_age - age_baseline)
R_RMSE_AVDPM_m_i <- R_RMSE_AVDPM_m <- rep(0,max_age - age_baseline)
R_RMSE_AVDPM_f_i <- R_RMSE_AVDPM_f <- rep(0,max_age - age_baseline)

RMSE_Gom_m_i <- (deaths_m - colSums(Int_hf_Gom_m))^2
RMSE_Gom_f_i <- (deaths_f - colSums(Int_hf_Gom_f))^2
RMSE_AVDPM_m_i <- (deaths_m - colSums(Int_hf_AVDPM_m))^2
RMSE_AVDPM_f_i <- (deaths_f - colSums(Int_hf_AVDPM_f))^2

for(i in 1:(max_age - age_baseline)){
  R_RMSE_Gom_m[i] <- mean(RMSE_Gom_m_i[1:i])
  R_RMSE_Gom_f[i] <- mean(RMSE_Gom_f_i[1:i])
  R_RMSE_AVDPM_m[i] <- mean(RMSE_AVDPM_m_i[1:i])
  R_RMSE_AVDPM_f[i] <- mean(RMSE_AVDPM_f_i[1:i])
}

plot(c(age_baseline:(max_age - 1)), R_RMSE_Gom_m, type='l', ylim=c(0, 50))
lines(c(age_baseline:(max_age - 1)), R_RMSE_AVDPM_m, col='grey')
plot(c(age_baseline:(max_age - 1)), R_RMSE_Gom_f, type='l', ylim=c(0, 50))
lines(c(age_baseline:(max_age - 1)), R_RMSE_AVDPM_f, col='grey')

# - Weighted (exposure) RMSE

WRMSE_Gom_m <- mean(colSums(E_cxM) * ((deaths_m - colSums(Int_hf_Gom_m))^2), na.rm=TRUE) / sum(E_cxM)
WRMSE_AVDPM_m <- mean(colSums(E_cxM) * ((deaths_m - colSums(Int_hf_AVDPM_m))^2), na.rm=TRUE) / sum(E_cxM)

WRMSE_Gom_f <- mean(colSums(E_cxF) * ((deaths_f - colSums(Int_hf_Gom_f))^2), na.rm=TRUE) / sum(E_cxF)
WRMSE_AVDPM_f <- mean(colSums(E_cxF) * ((deaths_f - colSums(Int_hf_AVDPM_f))^2), na.rm=TRUE) / sum(E_cxF)

R_WRMSE_Gom_m_i <- R_WRMSE_Gom_m <- rep(0,max_age - age_baseline)
R_WRMSE_Gom_f_i <- R_WRMSE_Gom_f <- rep(0,max_age - age_baseline)
R_WRMSE_AVDPM_m_i <- R_WRMSE_AVDPM_m <- rep(0,max_age - age_baseline)
R_WRMSE_AVDPM_f_i <- R_WRMSE_AVDPM_f <- rep(0,max_age - age_baseline)

WRMSE_Gom_m_i <- colSums(E_cxM) * ((deaths_m - colSums(Int_hf_Gom_m))^2) / sum(E_cxM)
WRMSE_Gom_f_i <- colSums(E_cxF) * ((deaths_f - colSums(Int_hf_Gom_f))^2) / sum(E_cxF)
WRMSE_AVDPM_m_i <- colSums(E_cxM) * ((deaths_m - colSums(Int_hf_AVDPM_m))^2) / sum(E_cxM)
WRMSE_AVDPM_f_i <- colSums(E_cxF) * ((deaths_f - colSums(Int_hf_AVDPM_f))^2) / sum(E_cxF)

for(i in 1:(max_age - age_baseline)){
  R_WRMSE_Gom_m[i] <- mean(WRMSE_Gom_m_i[1:i])
  R_WRMSE_Gom_f[i] <- mean(WRMSE_Gom_f_i[1:i])
  R_WRMSE_AVDPM_m[i] <- mean(WRMSE_AVDPM_m_i[1:i])
  R_WRMSE_AVDPM_f[i] <- mean(WRMSE_AVDPM_f_i[1:i])
}

plot(c(age_baseline:(max_age - 1)), R_WRMSE_Gom_m, type='l', ylim=c(0, 1))
lines(c(age_baseline:(max_age - 1)), R_WRMSE_AVDPM_m, col='grey')

plot(c(age_baseline:(max_age - 1)), R_WRMSE_Gom_f, type='l', ylim=c(0, 1))
lines(c(age_baseline:(max_age - 1)), R_WRMSE_AVDPM_f, col='grey')


#### - Overall assessment
deaths_mf <- colSums(d_xM) + colSums(d_xF)

## - Analysis of the integrated hazard function
### - Gompertz
Int_hf_Gom_mf <- Int_hf_Gom_m + Int_hf_Gom_f
Int_hf_AVDPM_mf <- Int_hf_AVDPM_m + Int_hf_AVDPM_f

Dev_Gom_mf <- 2 * (deaths_mf * (log(deaths_mf) - log(colSums(Int_hf_Gom_mf))) - (deaths_mf - colSums(Int_hf_Gom_mf)))
Dev_AVDPM_mf <- 2 * (deaths_mf * (log(deaths_mf) - log(colSums(Int_hf_AVDPM_mf))) - (deaths_mf - colSums(Int_hf_AVDPM_mf)))

pdr_Gom_mf <- sign(deaths_mf - colSums(Int_hf_Gom_mf)) * sqrt(Dev_Gom_mf)
pdr_AVDPM_mf <- sign(deaths_mf - colSums(Int_hf_AVDPM_mf)) * sqrt(Dev_AVDPM_mf)

plot(c(age_baseline:(max_age-1)), pdr_Gom_mf)
abline(h=2)
abline(h=-2)
points(c(age_baseline:(max_age-1)), pdr_AVDPM_mf, pch=16)

chi2_Gom_mf <- sum(na.omit(pdr_Gom_mf)^2)
chi2_AVDPM_mf <- sum(na.omit(pdr_AVDPM_mf)^2)

pearson_Gom_mf <- (deaths_mf - colSums(Int_hf_Gom_mf)) / sqrt(colSums(Int_hf_Gom_mf))
pearson_AVDPM_mf <- (deaths_mf - colSums(Int_hf_AVDPM_mf)) / sqrt(colSums(Int_hf_AVDPM_mf))

qqnorm(pearson_Gom_mf, pch = 1, frame = FALSE)
qqline(pearson_Gom_mf, col = "steelblue", lwd = 2)
qqnorm(pearson_AVDPM_mf, pch = 1, frame = FALSE)
qqline(pearson_AVDPM_mf, col = "steelblue", lwd = 2)


# - Calculate mean squared error of the hazard rates

## - In terms of integrated hazard functions, to compare with the number of deaths (advantage to weight by number of exposures)
RMSE_Gom_mf <- mean((deaths_mf - colSums(Int_hf_Gom_mf))^2, na.rm=TRUE)
RMSE_AVDPM_mf <- mean((deaths_mf - colSums(Int_hf_AVDPM_mf))^2, na.rm=TRUE)

R_RMSE_Gom_mf_i <- R_RMSE_Gom_mf <- rep(0,max_age - age_baseline)
R_RMSE_AVDPM_mf_i <- R_RMSE_AVDPM_mf <- rep(0,max_age - age_baseline)

RMSE_Gom_mf_i <- (deaths_mf - colSums(Int_hf_Gom_mf))^2
RMSE_AVDPM_mf_i <- (deaths_mf - colSums(Int_hf_AVDPM_mf))^2

for(i in 1:(max_age - age_baseline)){
  R_RMSE_Gom_mf[i] <- mean(RMSE_Gom_mf_i[1:i])
  R_RMSE_AVDPM_mf[i] <- mean(RMSE_AVDPM_mf_i[1:i])
}

plot(c(age_baseline:(max_age - 1)), R_RMSE_Gom_mf, type='l')#, ylim=c(0, 25))
lines(c(age_baseline:(max_age - 1)), R_RMSE_AVDPM_mf, col='grey')

# - Weighted (exposure) RMSE

WRMSE_Gom_mf <- mean(colSums(E_cxM + E_cxF) * ((deaths_mf - colSums(Int_hf_Gom_mf))^2), na.rm=TRUE) / sum(E_cxM + E_cxF)
WRMSE_AVDPM_mf <- mean(colSums(E_cxM + E_cxF) * ((deaths_mf - colSums(Int_hf_AVDPM_mf))^2), na.rm=TRUE) / sum(E_cxM + E_cxF)

R_WRMSE_Gom_mf_i <- R_WRMSE_Gom_mf <- rep(0,max_age - age_baseline)
R_WRMSE_AVDPM_mf_i <- R_WRMSE_AVDPM_mf <- rep(0,max_age - age_baseline)

WRMSE_Gom_mf_i <- colSums(E_cxM + E_cxF) * ((deaths_mf - colSums(Int_hf_Gom_mf))^2) / sum(E_cxM + E_cxF)
WRMSE_AVDPM_mf_i <- colSums(E_cxM + E_cxF) * ((deaths_mf - colSums(Int_hf_AVDPM_mf))^2) / sum(E_cxM + E_cxF)

for(i in 1:(max_age - age_baseline)){
  R_WRMSE_Gom_mf[i] <- mean(WRMSE_Gom_mf_i[1:i])
  R_WRMSE_AVDPM_mf[i] <- mean(WRMSE_AVDPM_mf_i[1:i])
}

plot(c(age_baseline:(max_age - 1)), R_WRMSE_Gom_mf, type='l')#, ylim=c(0, 1))
lines(c(age_baseline:(max_age - 1)), R_WRMSE_AVDPM_mf, col='grey')

# - Plot of number of deaths for each age
plot(c(age_baseline:(max_age - 1)), deaths_mf, type='l', ylim=c(0, 35))
lines(c(age_baseline:(max_age - 1)), colSums(Int_hf_Gom_mf), col='grey')
lines(c(age_baseline:(max_age - 1)), colSums(Int_hf_AVDPM_mf), col='red')


#### - 5 years bands age
deaths_m5 <- deaths_f5 <- deaths_mf5 <- Int_hf_Gom_m5 <- Int_hf_Gom_f5 <- Int_hf_AVDPM_m5 <- Int_hf_AVDPM_f5 <- Int_hf_Gom_mf5 <- Int_hf_AVDPM_mf5 <- rep(0,((max_age - age_baseline) / 5))

for(ag in 1:((max_age - age_baseline) / 5)){
  deaths_m5[ag] <- sum(deaths_m[((ag-1)*5 + 1) :(ag * 5)])
  deaths_f5[ag] <- sum(deaths_f[((ag-1)*5 + 1) :(ag * 5)])
  deaths_mf5[ag] <- sum(deaths_mf[((ag-1)*5 + 1) :(ag * 5)])
  Int_hf_Gom_m5[ag] <- sum(colSums(Int_hf_Gom_m)[((ag-1)*5 + 1) :(ag * 5)])
  Int_hf_Gom_f5[ag] <- sum(colSums(Int_hf_Gom_f)[((ag-1)*5 + 1) :(ag * 5)])
  Int_hf_Gom_mf5[ag] <- sum(colSums(Int_hf_Gom_mf)[((ag-1)*5 + 1) :(ag * 5)])
  Int_hf_AVDPM_m5[ag] <- sum(colSums(Int_hf_AVDPM_m)[((ag-1)*5 + 1) :(ag * 5)])
  Int_hf_AVDPM_f5[ag] <- sum(colSums(Int_hf_AVDPM_f)[((ag-1)*5 + 1) :(ag * 5)])
  Int_hf_AVDPM_mf5[ag] <- sum(colSums(Int_hf_AVDPM_mf)[((ag-1)*5 + 1) :(ag * 5)])
}

Dev_Gom_m5   <- 2 * (deaths_m5 * (log(deaths_m5) - log(Int_hf_Gom_m5  )) - (deaths_m5 - Int_hf_Gom_m5  ))
Dev_AVDPM_m5 <- 2 * (deaths_m5 * (log(deaths_m5) - log(Int_hf_AVDPM_m5)) - (deaths_m5 - Int_hf_AVDPM_m5))

pdr_Gom_m5   <- sign(deaths_m5 - Int_hf_Gom_m5  ) * sqrt(Dev_Gom_m5  )
pdr_AVDPM_m5 <- sign(deaths_m5 - Int_hf_AVDPM_m5) * sqrt(Dev_AVDPM_m5)

plot(1:6, pdr_Gom_m5)
abline(h=2)
abline(h=-2)
points(1:6, pdr_AVDPM_m5, pch=16)

chi2_Gom_m5 <- sum((pdr_Gom_m5)^2)
chi2_AVDPM_m5 <- sum((pdr_AVDPM_m5)^2)

pearson_Gom_m5 <- (deaths_m5 - Int_hf_Gom_m5) / sqrt(Int_hf_Gom_m5)
pearson_AVDPM_m5 <- (deaths_m5 - Int_hf_AVDPM_m5) / sqrt(Int_hf_AVDPM_m5)

qqnorm(pearson_Gom_m5, pch = 1, frame = FALSE)
qqline(pearson_Gom_m5, col = "steelblue", lwd = 2)
qqnorm(pearson_AVDPM_m5, pch = 1, frame = FALSE)
qqline(pearson_AVDPM_m5, col = "steelblue", lwd = 2)

## - females
Dev_Gom_f5   <- 2 * (deaths_f5 * (log(deaths_f5) - log(Int_hf_Gom_f5  )) - (deaths_f5 - Int_hf_Gom_f5  ))
Dev_AVDPM_f5 <- 2 * (deaths_f5 * (log(deaths_f5) - log(Int_hf_AVDPM_f5)) - (deaths_f5 - Int_hf_AVDPM_f5))

pdr_Gom_f5   <- sign(deaths_f5 - Int_hf_Gom_f5  ) * sqrt(Dev_Gom_f5  )
pdr_AVDPM_f5 <- sign(deaths_f5 - Int_hf_AVDPM_f5) * sqrt(Dev_AVDPM_f5)

plot(1:6, pdr_Gom_f5)
abline(h=2)
abline(h=-2)
points(1:6, pdr_AVDPM_f5, pch=16)

chi2_Gom_f5 <- sum(pdr_Gom_f5^2)
chi2_AVDPM_f5 <- sum(pdr_AVDPM_f5^2)

pearson_Gom_f5 <- (deaths_f5 - Int_hf_Gom_f5) / sqrt(Int_hf_Gom_f5)
pearson_AVDPM_f5 <- (deaths_f5 - Int_hf_AVDPM_f5) / sqrt(Int_hf_AVDPM_f5)

qqnorm(pearson_Gom_f5, pch = 1, frame = FALSE)
qqline(pearson_Gom_f5, col = "steelblue", lwd = 2)
qqnorm(pearson_AVDPM_f5, pch = 1, frame = FALSE)
qqline(pearson_AVDPM_f5, col = "steelblue", lwd = 2)


## - males + females
Dev_Gom_mf5   <- 2 * (deaths_mf5 * (log(deaths_mf5) - log(Int_hf_Gom_mf5  )) - (deaths_mf5 - Int_hf_Gom_mf5  ))
Dev_AVDPM_mf5 <- 2 * (deaths_mf5 * (log(deaths_mf5) - log(Int_hf_AVDPM_mf5)) - (deaths_mf5 - Int_hf_AVDPM_mf5))

pdr_Gom_mf5   <- sign(deaths_mf5 - Int_hf_Gom_mf5  ) * sqrt(Dev_Gom_mf5  )
pdr_AVDPM_mf5 <- sign(deaths_mf5 - Int_hf_AVDPM_mf5) * sqrt(Dev_AVDPM_mf5)

plot(1:6, pdr_Gom_mf5)
abline(h=2)
abline(h=-2)
points(1:6, pdr_AVDPM_mf5, pch=16)

chi2_Gom_mf5 <- sum(na.omit(pdr_Gom_mf5)^2)
chi2_AVDPM_mf5 <- sum(na.omit(pdr_AVDPM_mf5)^2)

pearson_Gom_mf5 <- (deaths_mf5 - Int_hf_Gom_mf5) / sqrt(Int_hf_Gom_mf5)
pearson_AVDPM_mf5 <- (deaths_mf5 - Int_hf_AVDPM_mf5) / sqrt(Int_hf_AVDPM_mf5)

qqnorm(pearson_Gom_mf5, pch = 1, frame = FALSE)
qqline(pearson_Gom_mf5, col = "steelblue", lwd = 2)
qqnorm(pearson_AVDPM_mf5, pch = 1, frame = FALSE)
qqline(pearson_AVDPM_mf5, col = "steelblue", lwd = 2)


# - Calculate mean squared error of the hazard rates

## - In terms of integrated hazard functions, to compare with the number of deaths (advantage to weight by number of exposures)
RMSE_Gom_mf5 <- mean((deaths_mf5 - Int_hf_Gom_mf5)^2, na.rm=TRUE)
RMSE_AVDPM_mf5 <- mean((deaths_mf5 - Int_hf_AVDPM_mf5)^2, na.rm=TRUE)

R_RMSE_Gom_mf5_i <- R_RMSE_Gom_mf5 <- rep(0, 6)
R_RMSE_AVDPM_mf5_i <- R_RMSE_AVDPM_mf5 <- rep(0, 6)

RMSE_Gom_mf5_i <- (deaths_mf5 - Int_hf_Gom_mf5)^2
RMSE_AVDPM_mf5_i <- (deaths_mf5 - Int_hf_AVDPM_mf5)^2

for(i in 1:6){
  R_RMSE_Gom_mf5[i] <- mean(RMSE_Gom_mf5_i[1:i])
  R_RMSE_AVDPM_mf5[i] <- mean(RMSE_AVDPM_mf5_i[1:i])
}

plot(1:6, R_RMSE_Gom_mf5, type='l')#, ylim=c(0, 25))
lines(1:6, R_RMSE_AVDPM_mf5, col='grey')

# - Weighted (exposure) RMSE

WRMSE_Gom_mf5 <- mean(colSums(E_cxM + E_cxF) * ((deaths_mf5 - colSums(Int_hf_Gom_mf5))^2), na.rm=TRUE) / sum(E_cxM + E_cxF)
WRMSE_AVDPM_mf5 <- mean(colSums(E_cxM + E_cxF) * ((deaths_mf5 - colSums(Int_hf_AVDPM_mf5))^2), na.rm=TRUE) / sum(E_cxM + E_cxF)

R_WRMSE_Gom_mf5_i <- R_WRMSE_Gom_mf5 <- rep(0,max_age - age_baseline)
R_WRMSE_AVDPM_mf5_i <- R_WRMSE_AVDPM_mf5 <- rep(0,max_age - age_baseline)

WRMSE_Gom_mf5_i <- colSums(E_cxM + E_cxF) * ((deaths_mf5 - colSums(Int_hf_Gom_mf5))^2) / sum(E_cxM + E_cxF)
WRMSE_AVDPM_mf5_i <- colSums(E_cxM + E_cxF) * ((deaths_mf5 - colSums(Int_hf_AVDPM_mf5))^2) / sum(E_cxM + E_cxF)

for(i in 1:(max_age - age_baseline)){
  R_WRMSE_Gom_mf5[i] <- mean(WRMSE_Gom_mf5_i[1:i])
  R_WRMSE_AVDPM_mf5[i] <- mean(WRMSE_AVDPM_mf5_i[1:i])
}

plot(c(age_baseline:(max_age - 1)), R_WRMSE_Gom_mf5, type='l')#, ylim=c(0, 1))
lines(c(age_baseline:(max_age - 1)), R_WRMSE_AVDPM_mf5, col='grey')

# - Plot of number of deaths for each age
plot(c(age_baseline:(max_age - 1)), deaths_mf5, type='l', ylim=c(0, 35))
lines(c(age_baseline:(max_age - 1)), colSums(Int_hf_Gom_mf5), col='grey')
lines(c(age_baseline:(max_age - 1)), colSums(Int_hf_AVDPM_mf5), col='red')



#===================== - WAIC by age - ===========================

WAIC_age_p1_i <- function(unit){
  
  M_theta <- length(seq_bi_thn)
  
  sum_over_m <- 0
  
  for(m in 1:M_theta){
    pf_vector_num <- log(pi_mw_thn[m,]) - exp(alpha_thn[m,1] + beta_thn[m,1] * st_ageM[unit,j] + gamma_star_thn[m,seq(1,M*K-1,M)]) * (exp(beta_thn[m,1] * E_cxM[unit,j]) - 1) / beta_thn[m,1] + d_xM[unit,j] * (alpha_thn[m,1] + beta_thn[m,1] * (st_ageM[unit,j] + E_cxM[unit,j]) + gamma_star_thn[m,seq(1,M*K-1,M)]) - 
                                          exp(alpha_thn[m,2] + beta_thn[m,2] * st_ageF[unit,j] + gamma_star_thn[m,seq(2,M*K  ,M)]) * (exp(beta_thn[m,2] * E_cxF[unit,j]) - 1) / beta_thn[m,2] + d_xF[unit,j] * (alpha_thn[m,2] + beta_thn[m,2] * (st_ageF[unit,j] + E_cxF[unit,j]) + gamma_star_thn[m,seq(2,M*K  ,M)]) + 
      dnorm(Z$agediff[unit], mean=z_star_ad_thn[m,], sd=sqrt(z_ad_v_thn[m]), log = TRUE) + Z$maleelder[unit] * log(z_star_mo_thn[m,]) + (1 - Z$maleelder[unit]) * log(1 - z_star_mo_thn[m,])
    
    max_pf_n <- max(pf_vector_num)
    
    pf_vector_den <- log(pi_mw_thn[m,]) + dnorm(Z$agediff[unit], mean=z_star_ad_thn[m,], sd=sqrt(z_ad_v_thn[m]), log = TRUE) + Z$maleelder[unit] * log(z_star_mo_thn[m,]) + (1 - Z$maleelder[unit]) * log(1 - z_star_mo_thn[m,])
    
    max_pf_d <- max(pf_vector_den)
    
    f_MF <- exp(max_pf_n + log(sum(exp(pf_vector_num - max_pf_n))) - max_pf_d - log(sum(exp(pf_vector_den - max_pf_d))))
    
    sum_over_m <- sum_over_m + f_MF
    
  }
  
  return(-log(M_theta) + log(sum_over_m))
  
}


WAIC_age_p2_i <- function(unit){
  
  M_theta <- length(seq_bi_thn)
  
  sum_over_m <- 0
  
  for(m in 1:M_theta){
    pf_vector_num <- log(pi_mw_thn[m,]) - exp(alpha_thn[m,1] + beta_thn[m,1] * st_ageM[unit,j] + gamma_star_thn[m,seq(1,M*K-1,M)]) * (exp(beta_thn[m,1] * E_cxM[unit,j]) - 1) / beta_thn[m,1] + d_xM[unit,j] * (alpha_thn[m,1] + beta_thn[m,1] * (st_ageM[unit,j] + E_cxM[unit,j]) + gamma_star_thn[m,seq(1,M*K-1,M)]) - 
                                          exp(alpha_thn[m,2] + beta_thn[m,2] * st_ageF[unit,j] + gamma_star_thn[m,seq(2,M*K  ,M)]) * (exp(beta_thn[m,2] * E_cxF[unit,j]) - 1) / beta_thn[m,2] + d_xF[unit,j] * (alpha_thn[m,2] + beta_thn[m,2] * (st_ageF[unit,j] + E_cxF[unit,j]) + gamma_star_thn[m,seq(2,M*K  ,M)]) + 
      dnorm(Z$agediff[unit], mean=z_star_ad_thn[m,], sd=sqrt(z_ad_v_thn[m]), log = TRUE) + Z$maleelder[unit] * log(z_star_mo_thn[m,]) + (1 - Z$maleelder[unit]) * log(1 - z_star_mo_thn[m,])
    
    max_pf_n <- max(pf_vector_num)
    
    pf_vector_den <- log(pi_mw_thn[m,]) + dnorm(Z$agediff[unit], mean=z_star_ad_thn[m,], sd=sqrt(z_ad_v_thn[m]), log = TRUE) + Z$maleelder[unit] * log(z_star_mo_thn[m,]) + (1 - Z$maleelder[unit]) * log(1 - z_star_mo_thn[m,])
    
    max_pf_d <- max(pf_vector_den)
    
    log_f_MF <- max_pf_n + log(sum(exp(pf_vector_num - max_pf_n))) - max_pf_d - log(sum(exp(pf_vector_den - max_pf_d)))
    
    sum_over_m <- sum_over_m + log_f_MF
    
  }
  
  return(sum_over_m/M_theta)
  
}


## - AVM
WAIC_age_p1_vector_os <- matrix(0, sample_size_test, max_age - age_baseline)# sapply(c(1:sample_size_test), WAIC_p1_i)
WAIC_age_p2_vector_os <- matrix(0, sample_size_test, max_age - age_baseline)# #sapply(c(1:sample_size_test), WAIC_p2_i)

for(j in 1:(max_age - age_baseline)){
  WAIC_age_p1_vector_os[,j] <- sapply(c(1:sample_size_test), WAIC_age_p1_i)
  WAIC_age_p2_vector_os[,j] <- sapply(c(1:sample_size_test), WAIC_age_p2_i)
  print(j)
}



p_WAIC_os <- 2 * (sum(WAIC_p1_vector_os) - sum(WAIC_p2_vector_os))
WAIC_os <- 2*(- sum(WAIC_p1_vector_os) + p_WAIC_os)



WAIC_male_p1_i <- function(unit){
  
  M_theta <- length(seq_bi_thn)
  
  sum_over_m <- 0
  
  for(m in 1:M_theta){
    log_d_Z <- dnorm(Z$agediff[unit], mean=z_star_ad_thn[m,], sd=sqrt(z_ad_v_thn[m]), log = TRUE) + Z$maleelder[unit] * log(z_star_mo_thn[m,]) + (1 - Z$maleelder[unit]) * log(1 - z_star_mo_thn[m,])
    log_d_F <-  - exp(alpha_thn[m,2] + beta_thn[m,2] * X$female[unit] + gamma_star_thn[m,seq(2,M*K  ,M)]) * (exp(beta_thn[m,2] * (t2_i[unit] + a_i[unit])) - exp(beta_thn[m,2] * a_i[unit])) / beta_thn[m,2] + d_ci$female[unit] * (alpha_thn[m,2] + beta_thn[m,2] * (t2_i[unit] + X$female[unit] + a_i[unit]) + gamma_star_thn[m,seq(2,M*K  ,M)])
    
    pf_vector_num <- log(pi_mw_thn[m,]) - exp(alpha_thn[m,1] + beta_thn[m,1] * X$male[unit]   + gamma_star_thn[m,seq(1,M*K-1,M)]) * (exp(beta_thn[m,1] * (t1_i[unit] + a_i[unit])) - exp(beta_thn[m,1] * a_i[unit])) / beta_thn[m,1] + d_ci$male[unit]   * (alpha_thn[m,1] + beta_thn[m,1] * (t1_i[unit] + X$male[unit]   + a_i[unit]) + gamma_star_thn[m,seq(1,M*K-1,M)]) + log_d_F + log_d_Z
      
    max_pf_n <- max(pf_vector_num)
    
    pf_vector_den <- log(pi_mw_thn[m,]) + log_d_F + log_d_Z
    
    max_pf_d <- max(pf_vector_den)
    
    f_MF <- exp(max_pf_n + log(sum(exp(pf_vector_num - max_pf_n))) - max_pf_d - log(sum(exp(pf_vector_den - max_pf_d))))
    
    sum_over_m <- sum_over_m + f_MF
    
  }
  
  return(-log(M_theta) + log(sum_over_m))
  
}


WAIC_male_p2_i <- function(unit){
  
  M_theta <- length(seq_bi_thn)
  
  sum_over_m <- 0
  
  for(m in 1:M_theta){
    log_d_Z <- dnorm(Z$agediff[unit], mean=z_star_ad_thn[m,], sd=sqrt(z_ad_v_thn[m]), log = TRUE) + Z$maleelder[unit] * log(z_star_mo_thn[m,]) + (1 - Z$maleelder[unit]) * log(1 - z_star_mo_thn[m,])
    log_d_F <-  - exp(alpha_thn[m,2] + beta_thn[m,2] * X$female[unit] + gamma_star_thn[m,seq(2,M*K  ,M)]) * (exp(beta_thn[m,2] * (t2_i[unit] + a_i[unit])) - exp(beta_thn[m,2] * a_i[unit])) / beta_thn[m,2] + d_ci$female[unit] * (alpha_thn[m,2] + beta_thn[m,2] * (t2_i[unit] + X$female[unit] + a_i[unit]) + gamma_star_thn[m,seq(2,M*K  ,M)])
    
    pf_vector_num <- log(pi_mw_thn[m,]) - exp(alpha_thn[m,1] + beta_thn[m,1] * X$male[unit]   + gamma_star_thn[m,seq(1,M*K-1,M)]) * (exp(beta_thn[m,1] * (t1_i[unit] + a_i[unit])) - exp(beta_thn[m,1] * a_i[unit])) / beta_thn[m,1] + d_ci$male[unit]   * (alpha_thn[m,1] + beta_thn[m,1] * (t1_i[unit] + X$male[unit]   + a_i[unit]) + gamma_star_thn[m,seq(1,M*K-1,M)]) + log_d_F + log_d_Z
    
    max_pf_n <- max(pf_vector_num)
    
    pf_vector_den <- log(pi_mw_thn[m,]) + log_d_F + log_d_Z
    
    max_pf_d <- max(pf_vector_den)
    
    log_f_MF <- max_pf_n + log(sum(exp(pf_vector_num - max_pf_n))) - max_pf_d - log(sum(exp(pf_vector_den - max_pf_d)))
    
    sum_over_m <- sum_over_m + log_f_MF
    
  }
  
  return(sum_over_m/M_theta)
  
}


WAIC_female_p1_i <- function(unit){
  
  M_theta <- length(seq_bi_thn)
  
  sum_over_m <- 0
  
  for(m in 1:M_theta){
    log_d_Z <- dnorm(Z$agediff[unit], mean=z_star_ad_thn[m,], sd=sqrt(z_ad_v_thn[m]), log = TRUE) + Z$maleelder[unit] * log(z_star_mo_thn[m,]) + (1 - Z$maleelder[unit]) * log(1 - z_star_mo_thn[m,])
    log_d_M <-  - exp(alpha_thn[m,1] + beta_thn[m,1] * X$male[unit] + gamma_star_thn[m,seq(1,M*K-1,M)]) * (exp(beta_thn[m,1] * (t1_i[unit] + a_i[unit])) - exp(beta_thn[m,1] * a_i[unit])) / beta_thn[m,1] + d_ci$male[unit] * (alpha_thn[m,1] + beta_thn[m,1] * (t1_i[unit] + X$male[unit] + a_i[unit]) + gamma_star_thn[m,seq(1,M*K-1,M)])
    
    pf_vector_num <- log(pi_mw_thn[m,]) - exp(alpha_thn[m,2] + beta_thn[m,2] * X$female[unit] + gamma_star_thn[m,seq(2,M*K  ,M)]) * (exp(beta_thn[m,2] * (t2_i[unit] + a_i[unit])) - exp(beta_thn[m,2] * a_i[unit])) / beta_thn[m,2] + d_ci$female[unit]   * (alpha_thn[m,2] + beta_thn[m,2] * (t2_i[unit] + X$female[unit] + a_i[unit]) + gamma_star_thn[m,seq(2,M*K  ,M)]) + log_d_M + log_d_Z
    
    max_pf_n <- max(pf_vector_num)
    
    pf_vector_den <- log(pi_mw_thn[m,]) + log_d_M + log_d_Z
    
    max_pf_d <- max(pf_vector_den)
    
    f_MF <- exp(max_pf_n + log(sum(exp(pf_vector_num - max_pf_n))) - max_pf_d - log(sum(exp(pf_vector_den - max_pf_d))))
    
    sum_over_m <- sum_over_m + f_MF
    
  }
  
  return(-log(M_theta) + log(sum_over_m))
  
}

WAIC_female_p2_i <- function(unit){
  
  M_theta <- length(seq_bi_thn)
  
  sum_over_m <- 0
  
  for(m in 1:M_theta){
    log_d_Z <- dnorm(Z$agediff[unit], mean=z_star_ad_thn[m,], sd=sqrt(z_ad_v_thn[m]), log = TRUE) + Z$maleelder[unit] * log(z_star_mo_thn[m,]) + (1 - Z$maleelder[unit]) * log(1 - z_star_mo_thn[m,])
    log_d_M <-  - exp(alpha_thn[m,1] + beta_thn[m,1] * X$male[unit] + gamma_star_thn[m,seq(1,M*K-1,M)]) * (exp(beta_thn[m,1] * (t1_i[unit] + a_i[unit])) - exp(beta_thn[m,1] * a_i[unit])) / beta_thn[m,1] + d_ci$male[unit] * (alpha_thn[m,1] + beta_thn[m,1] * (t1_i[unit] + X$male[unit] + a_i[unit]) + gamma_star_thn[m,seq(1,M*K-1,M)])
    
    pf_vector_num <- log(pi_mw_thn[m,]) - exp(alpha_thn[m,2] + beta_thn[m,2] * X$female[unit] + gamma_star_thn[m,seq(2,M*K  ,M)]) * (exp(beta_thn[m,2] * (t2_i[unit] + a_i[unit])) - exp(beta_thn[m,2] * a_i[unit])) / beta_thn[m,2] + d_ci$female[unit]   * (alpha_thn[m,2] + beta_thn[m,2] * (t2_i[unit] + X$female[unit] + a_i[unit]) + gamma_star_thn[m,seq(2,M*K  ,M)]) + log_d_M + log_d_Z
    
    max_pf_n <- max(pf_vector_num)
    
    pf_vector_den <- log(pi_mw_thn[m,]) + log_d_M + log_d_Z
    
    max_pf_d <- max(pf_vector_den)
    
    log_f_MF <- max_pf_n + log(sum(exp(pf_vector_num - max_pf_n))) - max_pf_d - log(sum(exp(pf_vector_den - max_pf_d)))
    
    sum_over_m <- sum_over_m + log_f_MF
    
  }
  
  return(sum_over_m/M_theta)
  
}




WAIC_m_p1_vector_os <- sapply(c(1:sample_size_test), WAIC_male_p1_i)
WAIC_m_p2_vector_os <- sapply(c(1:sample_size_test), WAIC_male_p2_i)

p_WAIC_m_os <- 2 * (sum(WAIC_m_p1_vector_os) - sum(WAIC_m_p2_vector_os))
WAIC_m_os <- 2*(- sum(WAIC_m_p1_vector_os) + p_WAIC_m_os)


WAIC_f_p1_vector_os <- sapply(c(1:sample_size_test), WAIC_female_p1_i)
WAIC_f_p2_vector_os <- sapply(c(1:sample_size_test), WAIC_female_p2_i)

p_WAIC_f_os <- 2 * (sum(WAIC_f_p1_vector_os) - sum(WAIC_f_p2_vector_os))
WAIC_f_os <- 2*(- sum(WAIC_f_p1_vector_os) + p_WAIC_f_os)






