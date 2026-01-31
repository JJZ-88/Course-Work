library(ggplot2)
library(cmdstanr)
library(bayesplot)
library(posterior)
library(loo)
library(gridExtra)
library(ggpubr)
library(mvtnorm)
library(ggplot2)
library(gtable)
library(grid)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

#### Q1

## Data
y <- c(3.79, 5.8, 12.79, 15.52, 9.98, 18.65, 13.21, 13.91, 8.16, 4.81, 4.59, 2.23)
t <- c(0.083, 0.167, 0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8)
n <- length(y)

t0 <- 0

ggplot(data = data.frame(y = y, t = t), aes(x = t, y = y)) +
  geom_point() + 
  labs(x = "time", y = "drug concentration")

## Fit Basic PCK STAN Model
u0 <- c(1200, 0)

data_pharm <- list(y = y, t = t, n = n, t0 = t0, u0 = u0)

mod_pharm <- cmdstan_model("pharmacokinetic.stan")

init <- function(){
  list(cl = exp(rnorm(1, log(10), 0.5)),
       vc = exp(rnorm(1, log(35), 0.5)),
       ka = exp(rnorm(1, log(2.5), 1)),
       sigma = abs(rnorm(1, 0, 1)))
}

vars_pharm <- c('central clearance', 'compartment volume', 'absorption rate', 'sigma')

fit_pharm <- mod_pharm$sample(data = data_pharm, 
                              seed = 547,
                              chains = 4, 
                              parallel_chains = 4,
                              init = init,
                              save_warmup = TRUE)

## Convergence and Posterior Checks for Basic PCK Model
summary_pharm <- as.matrix(fit_pharm$summary(variables = c('cl', 'vc', 'ka', 'sigma')))
summary_pharm[,-1] <- round(apply(summary_pharm[,-1], 2, as.numeric), 2)

pdf('pck1_summary.pdf', width = 10, height = 1.5)
grid.table(summary_pharm, rows = NULL)
dev.off()


draws_pharm <- as_draws_df(fit_pharm$draws(inc_warmup = TRUE))
names(draws_pharm)[2:5] <- vars_pharm

pdf('pck1_trace.pdf', width = 6, height = 4)
bayesplot::mcmc_trace(draws_pharm, 
                      n_warmup= 1000, 
                      pars = vars_pharm)
dev.off()

pdf('pck1_dens.pdf', width = 6, height = 4)
bayesplot::mcmc_dens_overlay(draws_pharm, pars = vars_pharm)
dev.off()

y_pred_pharm <- as.matrix(as_draws_df(fit_pharm$draws(variables = c("y_pred"))))[,-(13:15)]

ppc_pharma <- bayesplot::ppc_ribbon(y = y, yrep = y_pred_pharm, x = t, y_draw = 'point') +
  labs(x = 'time', y = 'drug concentration')

## Fit Peripheral PCK STAN Model
u0_2 <- c(1200, 0, 0)

data_pharm_periph <- list(y = y, t = t, n = n, t0 = t0, u0 = u0_2)

mod_pharm_periph <- cmdstan_model("pharmacokinetic_peripheral.stan")

init <- function(){
  list(cl = exp(rnorm(1, log(10), 0.5)),
       vc = exp(rnorm(1, log(35), 0.5)),
       ka = exp(rnorm(1, log(2.5), 1)),
       sigma = abs(rnorm(1, 0, 1)),
       q = exp(rnorm(1, log(15), sqrt(0.5))),
       vp = exp(rnorm(1, log(105), sqrt(0.5))))
}

vars_pharm_periph <- c('central clearance', 'central compartment volume', 'absorption rate', 
                'sigma', 'intercompartmental clearance', 'peripheral compartment volume')

fit_pharm_periph <- mod_pharm_periph$sample(data = data_pharm_periph, 
                              seed = 547,
                              chains = 4, 
                              parallel_chains = 4,
                              init = init,
                              save_warmup = TRUE)

## Convergence and Posterior Checks for Peripheral PCK Model
summary_pharm_periph <- as.matrix(fit_pharm_periph$summary(variables = c('cl', 'vc', 'ka', 'sigma', 'q', 'vp')))
summary_pharm_periph[,-1] <- round(apply(summary_pharm_periph[,-1], 2, as.numeric), 2)

pdf('pck2_summary.pdf', width = 10, height = 2)
grid.table(summary_pharm_periph, rows = NULL)
dev.off()

draws_pharm_periph <- as_draws_df(fit_pharm_periph$draws(inc_warmup = TRUE))
names(draws_pharm_periph)[2:7] <- vars_pharm_periph

pdf('pck2_trace.pdf', width = 8, height = 4)
bayesplot::mcmc_trace(draws_pharm_periph, 
                      n_warmup= 1000, 
                      pars = vars_pharm_periph)
dev.off()

pdf('pck2_dens.pdf', width = 8, height = 4)
bayesplot::mcmc_dens_overlay(draws_pharm_periph, pars = vars_pharm_periph)
dev.off()

y_pred_pharm_periph <- as.matrix(as_draws_df(fit_pharm_periph$draws(variables = c("y_pred"))))[,-(13:15)]

ppc_pharma_periph <- bayesplot::ppc_ribbon(y = y, yrep = y_pred_pharm_periph, x = t, y_draw = 'point') +
  labs(x = 'time', y = 'drug concentration')

pdf("pck_ppc.pdf", width = 8, height = 4)
ggarrange(ppc_pharma, ppc_pharma_periph, nrow = 1)
dev.off()

## LOOCV Estimates for Both Models

log_lik_draws <- fit_pharm$draws("loglik")
loo_estimate <- loo(log_lik_draws, r_eff = relative_eff(log_lik_draws))

print(loo_estimate)
print(loo_estimate$diagnostics)


log_lik_draws_periph <- fit_pharm_periph$draws("loglik")
loo_estimate_periph <-
  loo(log_lik_draws_periph, r_eff = relative_eff(log_lik_draws_periph))

print(loo_estimate_periph)
print(loo_estimate_periph$diagnostics)


loo_cent_est <- loo_estimate$estimates
loo_cent_par <- loo_estimate$diagnostics$pareto_k
loo_cent_ess <- loo_estimate$diagnostics$n_eff
pck_central <- c(loo_cent_est[1,1], loo_cent_est[1,2], loo_cent_est[2,1], loo_cent_est[2, 2], max(loo_cent_par), mean(loo_cent_par), mean(loo_cent_ess))

loo_per_est <- loo_estimate_periph$estimates
loo_per_par <- loo_estimate_periph$diagnostics$pareto_k
loo_per_ess <- loo_estimate_periph$diagnostics$n_eff
pck_peripheral <- c(loo_per_est[1,1], loo_per_est[1,2], loo_per_est[2,1], loo_per_est[2, 2], max(loo_per_par), mean(loo_per_par),  mean(loo_per_ess))

loo_mat <- cbind(c('pck_central', 'pck_peripheral'), round(rbind(pck_central, pck_peripheral), 3))
colnames(loo_mat) <- c('model', 'elpd_loo estimate', 'elpd_loo se', 'p_loo estimate', 'p_loo se', 'max pareto k', 'mean pareto k', 'mean ESS')
loo_mat

pdf('pck_loo.pdf', width = 10, height = 1.5)
grid.table(loo_mat, rows = NULL)
dev.off()



###### Q2

##test function and gradient

cov <- matrix(c(2.82, 2.71, 0.18, 2.42, -1.87, 2.71, 7.35, -4.02, -1.17, -1.05, 0.18,
                -4.02, 8.85, 1.64, 3.81, 2.42, -1.17, 1.64, 5.94, -3.95, -1.87, -1.05, 3.81, -3.95, 6.4), 
              nrow = 5, ncol = 5)
dcov <- det(cov)
icov <- solve(cov)

testp <- function(x){
  #cov <- matrix(c(2.82, 2.71, 0.18, 2.42, -1.87, 2.71, 7.35, -4.02, -1.17, -1.05, 0.18,
  #                -4.02, 8.85, 1.64, 3.81, 2.42, -1.17, 1.64, 5.94, -3.95, -1.87, -1.05, 3.81, -3.95, 6.4), 
  #              nrow = 5, ncol = 5)
  return((-5 / 2) * log(2 * pi) - 1 / 2 * log(dcov) - t(x) %*% icov %*% x / 2)
}

testg <- function(x){
  return(-icov %*% x)
}

##HMC with leapfrog proposal

leapfrog_proposal <- function(log_p, grad_log_p, eps, l, pos, m){
  d <- dim(m)[1]
  z <- pos
  rho_0 <- rmvnorm(1, rep(0, d), m)
  m_inv <- solve(m) #using poor method, but fine for diagonal matrices 
  rho <- rho_0 + eps / 2 * as.vector(grad_log_p(z))
  for(i in 1:l){
    z <- z + eps * m_inv %*% t(rho)
    if(i != l){
      rho <- rho + eps * as.vector(grad_log_p(z))
    }
  }
  rho <- rho + eps / 2 * as.vector(grad_log_p(z))
  
  return(list(rho_0 = rho_0, rho = rho, z = z))
}

my_hmc <- function(log_p, grad_log_p, theta0, n_iter, eps, l, m = diag(5)){
  d <- length(theta0)
  chain <- matrix(0, ncol = d, nrow = n_iter + 1)
  chain[1,] <- theta0

  for(i in 1:n_iter){
    lp <- leapfrog_proposal(log_p, grad_log_p, eps, l, chain[i,], m)
    #print(lp)
    h_old <- -log_p(chain[i,]) - log(dmvnorm(lp$rho_0, rep(0, d), m))
    h_new <- -log_p(lp$z) - log(dmvnorm( -lp$rho, rep(0, d), m))

    if(runif(1) < min(1, exp(-h_new + h_old))){
      chain[i + 1,] <- lp$z
    }
    else{
      chain[i + 1,] <- chain[i,]
    }
  }
  return(chain)
}

##Helper Function to Run Multiple Times

run_hmc <- function(c, eps, l, n_iter, testp, testg){
  all_chains <- list()
  #all_chains <- matrix(0, nrow = 5, ncol = n_iter / 2 * c)
  for(j in 1:c){
    theta0 <- rmvnorm(1, rep(0, 5), cov)
    hmc <- my_hmc(testp, testg, theta0, n_iter, eps, l)
    colnames(hmc) <- c('z_1', 'z_2', 'z_3', 'z_4', 'z_5')
    all_chains[[j]] <- hmc
  }
  return(all_chains)
}

##Original parameters - diagnostics

set.seed(547)
hmc_1 <- run_hmc(8, 0.05, 10, 1000, testp, testg)

hmc_combined <- do.call(rbind, lapply(hmc_1, function(h){h[501:1001,]}))
hmc_1_summary <- cbind(c(1, 2, 3, 4, 5), round(as.matrix(summarize_draws(hmc_combined)[,c('mean', 'median', 'sd', 'q5', 'q95', 'rhat', 'ess_bulk')]), 3))

pdf('hmc_b.pdf', width = 10, height = 2)
grid.table(hmc_1_summary)
dev.off()


pdf('hmc_trace.pdf', width = 6, height = 4)
bayesplot::mcmc_trace(hmc_1, 
                      n_warmup= 500)
dev.off()

pdf('hmc_dens.pdf', width = 8, height = 4)
bayesplot::mcmc_dens_overlay(lapply(hmc_1, function(h){h[501:1001,]}))
dev.off()

##Monte Carlo Squared Error
hmc_combined_2 <- do.call(cbind, lapply(hmc_1, function(h){h[501:1001, 1]}))
hmc_means <- rep(0, 500)
hmc_x2 <- rep(0, 500)
for(i in 1:500){
  hmc_means[i] = mean(hmc_combined_2[1:i,])^2
  hmc_x2[i] = mean(hmc_combined_2[1:i,]^2)^2
}

error_x <- ggplot(data.frame(x = seq(1, 500, 1), y = hmc_means), aes(x = x, y = y)) +
  geom_path() +
  labs(x = 'sampling length', y = expression('squared error of ' * x[1]))
  

error_x2 <- ggplot(data.frame(x = seq(1, 500, 1), y = hmc_x2), aes(x = x, y = y)) +
  geom_path() +
  labs(x = 'sampling length', y = expression('squared error of ' * x[1]^2))

pdf('hmc_error.pdf', width = 8, height = 4)
ggarrange(error_x, error_x2)
dev.off()

##Tuning Trajectory
l <- c(10, 80, 120)
hmc_l <- matrix(c(1, 2, 3, 4, 5), nrow = 5)
set.seed(547)
for(i in c(1, 2, 3)){
  hmc <- run_hmc(8, 0.05, l[i], 1000, testp, testg)
  hmc_summary <- summarize_draws(do.call(rbind, lapply(hmc, function(h){h[501:1001,]})))

  hmc_l <- cbind(hmc_l, round(as.matrix(hmc_summary[,c('rhat', 'ess_bulk')]), 3), 
                 round(as.matrix(hmc_summary[,'ess_bulk']) / l[i] / 4000, 3) )
}

colnames(hmc_l) <- c('', 'rhat_10', 'ess_10', 'ess/gradient_10','rhat_80', 'ess_80', 'ess/gradient_80','rhat_120', 'ess_120', 'ess/gradient_120')




pdf('hmc_l.pdf', width = 10, height = 2)
grid.table(hmc_l)
dev.off()

##Jitter Trajectory

my_hmc_jitter <- function(log_p, grad_log_p, theta0, n_iter, eps, m = diag(5)){
  d <- length(theta0)
  chain <- matrix(0, ncol = d, nrow = n_iter + 1)
  chain[1,] <- theta0
  l_iter <- sample(1:120, n_iter, replace=T)
  
  for(i in 1:n_iter){
    lp <- leapfrog_proposal(log_p, grad_log_p, eps, l_iter[i], chain[i,], m)
    #print(lp)
    h_old <- -log_p(chain[i,]) - log(dmvnorm(lp$rho_0, rep(0, d), m))
    h_new <- -log_p(lp$z) - log(dmvnorm( -lp$rho, rep(0, d), m))
    
    if(runif(1) < min(1, exp(-h_new + h_old))){
      chain[i + 1,] <- lp$z
    }
    else{
      chain[i + 1,] <- chain[i,]
    }
  }
  return(list(chain = chain, grad_count = sum(l_iter[501:1000])))
}

run_hmc_jitter <- function(c, eps, n_iter, testp, testg){
  all_chains <- list()
  #all_chains <- matrix(0, nrow = 5, ncol = n_iter / 2 * c)
  total_grad <- 0
  for(j in 1:c){
    theta0 <- rmvnorm(1, rep(0, 5), cov)
    out <- my_hmc_jitter(testp, testg, theta0, n_iter, eps)
    hmc <- out[['chain']]
    colnames(hmc) <- c('z_1', 'z_2', 'z_3', 'z_4', 'z_5')
    all_chains[[j]] <- hmc
    total_grad <- total_grad + out[['grad_count']]
  }
  return(list(all_chains, total_grad))
}

set.seed(547)
ret <- run_hmc_jitter(8, 0.05, 1000, testp, testg)
hmc_jitter <- ret[[1]]

hmc_comb_jitter <- summarize_draws(do.call(rbind, lapply(hmc_jitter, function(h){h[501:1001,]})))
hmc_jitter_summary <- cbind(c(1, 2, 3, 4, 5), round(as.matrix(hmc_comb_jitter[,c('rhat', 'ess_bulk')]), 3),
                            round(as.matrix(hmc_comb_jitter[,'ess_bulk']) / ret[[2]], 3) )

pdf('hmc_jitter.pdf', width = 5, height = 2)
grid.table(hmc_jitter_summary)
dev.off()









