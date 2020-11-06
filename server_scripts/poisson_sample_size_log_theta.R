args <- commandArgs(trailingOnly = TRUE)
n <- as.numeric(args[1])
u_0 <- as.numeric(args[2])
step_size <- ifelse(u_0 == 0.1, 0.2, 0.3)
thin <- ifelse(u_0 == 0.1, 10, 40)
options(warn=-1)
library(loo)
library(cmdstanr)
library(transport)
library(deSolve)

flat_params <- function(ys) {
  # Calculate the parameters for the posterior 
  # distribution for theta in the Poisson model with 
  # a flat prior 
  return(list(alpha = sum(ys) + 1, beta = length(ys)))
}

W2 <- function(draws_1, draws_2) {
  return(transport::wasserstein1d(draws_1, draws_2, 2)^2)
}

prior <- function (t, x, params) {
  e_m_u_theta <- x[1]
  c <- params["c"]
  de_m_u_dt <- sqrt(c / e_m_u_theta - 2 * (1 - log(e_m_u_theta))) * -e_m_u_theta
  list(c(de_m_u_dt))
}

log_p_th <- function(th, y, p) { 
  
  parms <- p['c']
  times <- seq(from=0,to=exp(th),length.out = 2)
  xstart <- c(u_0 = exp(-p['u_0']))
  
  res <- ode(
    func=prior,
    y=xstart,
    times=times,
    parms=parms
  ) 
  log_prior <- log(res[2,2]) + th
  log_lik <- sum(dpois(y,exp(th),log=T))
  if (is.nan(log_prior)) {
    log_prior <- -Inf
  }
  return(log_prior + log_lik)
}

init_fun <- ifelse(u_0 == 0.1, 
                   function() list(theta = runif(1, min = 0.5, max = 2.0)), 
                   function() list(theta = runif(1, min = 0.5, max = 1.0)))

param_guess <- function(dat, p, mod, seed) {
  # Initial guess at mode of posterior 
  # of log(theta) using LBFGS
  N <- length(dat)
  mu <- mean(dat)
  s_dat <- list(y = dat, N = N, c = p[1], u_0 = p[2])
  s_dat$acc <- ifelse(u_0 == 0.1, 1e-8, 1e-11)
  fit <- mod$optimize(data = s_dat, init = init_fun, refresh = 0, algorithm = 'lbfgs', seed = seed)
  
  mu <- log(fit$mle())
  return(
   list(
     mean = mu,
     sd = 0.2
   ) 
  )
}

log_IW_t <- function(theta, dat, guess, p) {
  denom <- dt((theta-guess$mean)/guess$sd,df = 7, log = T) - log(guess$sd) 
  num <- sapply(theta, function(t) log_p_th(t,dat, p))
  return(num - denom)
}

gen_is_samps <- function(n_draws, dat_m, p, mod, seed) {
  par_guess <- param_guess(dat_m, p, mod, seed)
  test_samps <- rt(n_draws, df = 7) * par_guess$sd + par_guess$mean
  m_post_conj_log_IW <- log_IW_t(test_samps, dat = dat_m, par_guess, p)
  finite_draws <- which(is.finite(m_post_conj_log_IW))
  print(length(finite_draws))
  psis_IW <- loo::psis(m_post_conj_log_IW[finite_draws])
  lw <- psis_IW$log_weights[,1]
  weights <- exp(lw)
  m_post_conj <- sample(exp(test_samps[finite_draws]), 
                        size = n_draws, replace = T, prob = weights)
  return(list(samps = m_post_conj, weights = weights))
}

  
refit_both_models_signed_bb <- function(y1ton, 
                                        c_prior,
                                        u_0_prior,
                                        mod,
                                        lower,
                                        upper,
                                        n_draws,
                                        seed) {
  n <- length(y1ton)
  # Fit model that employs an informative prior
  stan_dat <- list(N = n, y = y1ton, c = c_prior, u_0 = u_0_prior)
  stan_dat$acc <- ifelse(u_0 == 0.1, 1e-8, 1e-11)
  n_post_conj_fit <- mod$sample(data = stan_dat, 
                                num_cores = 1, 
                                num_chains = 4,
                                num_samples = floor(n_draws / 4) * thin, 
                                num_warmup = 8000,
                                adapt_delta = 0.99,
                                max_treedepth = 12,
                                inv_metric = as.array(0.0407221),
                                step_size = step_size,
                                adapt_engaged = T,
                                init = init_fun,
                                thin = thin,
                                refresh=0,
                                seed = seed
  )
  s_samps <- drop(n_post_conj_fit$draws(variables = 'theta', inc_warmup=F))
  d_samps <- n_post_conj_fit$sampler_diagnostics()
  n_chains <- dim(s_samps)[2] 
  chains_needed <- ifelse(is.null(n_chains), 3, 4 - n_chains)
  if (chains_needed > 0) {
	  n_post_conj_fit <- mod$sample(data = stan_dat, 
					num_cores = 1, 
					num_chains = 4,
					num_samples = floor(n_draws / 4) * thin, 
					num_warmup = 8000,
					adapt_delta = 0.99,
					max_treedepth = 12,
					inv_metric = as.array(0.0407221),
					step_size = step_size,
					adapt_engaged = T,
					init = init_fun,
					thin = thin,
					refresh=0,
					seed = seed + 1e4
	  )
	  new_samps <- drop(n_post_conj_fit$draws(variables = 'theta', inc_warmup=F))
	  s_samps <- cbind(s_samps,new_samps[,1:chains_needed])
  }
  draws <- s_samps
  n_post_conj <- as.vector(draws)
  
  # Fit model with flat prior to 
  # n observed data points
  n_post_flat <- flat_params(y1ton)
  n_draws_flat <- rgamma(n_draws, shape = n_post_flat$alpha, rate = n_post_flat$beta)
  
  # Generate one draw from the posterior of informative prior model
  theta_now <- sample(n_post_conj,1)
  p_vec <- unlist(stan_dat[c('c','u_0')])
  
  mrange = 1:upper
  dists_nf_mc <- rep(NA_real_,upper+1)
  dists_nc_mf <- rep(NA_real_,upper+1)
  mean_seq <- matrix(NA_real_,upper+1, 2)
  var_weights <- rep(NA_real_,upper)
  conj_samps <- matrix(NA_real_,upper + 1, n_draws)
  
  # Measure W2 between initial posteriors
  dists_nf_mc[1] <- W2(n_post_conj, n_draws_flat)
  
  # Measure W2 between initial posteriors
  dists_nc_mf[1] <- dists_nf_mc[1]
  mean_seq[1,1] <- mean(y1ton)
  mean_seq[1,2] <- mean_seq[1,1]
  conj_samps[1,] <- n_post_conj
  for(m_i in seq_along(mrange)){
    # Generate m_i independent draws of data conditional
    # on a draw from the posterior employing the "informative" prior
    yseq = rpois(m_i, theta_now)
    y1tom = c(y1ton,yseq)
    mean_seq[m_i+1,2] <- sum(y1tom)
    # Refit the model with flat prior to 
    # m data points
    m_post_flat <- flat_params(y1tom)
    m_draws_flat <- rgamma(n_draws, shape = m_post_flat$alpha,
                           rate = m_post_flat$beta)
    
    # Refit the model with conjugate prior to 
    # m data points
    yseq = rpois(m_i, theta_now)
    y1tom = c(y1ton,yseq)
    mean_seq[m_i+1,1] <- sum(y1tom)
    ## Generate samples from the new posterior for the informative model
    is_samps <- gen_is_samps(n_draws, y1tom, p_vec, mod, seed)
    conj_samps[m_i+1,] <- is_samps$samps
    var_weights[m_i] <- sum((is_samps$weights/sum(is_samps$weights))^2)
    # Measure distance from informative prior fitted to expanded dataset
    # to flat fitted to original dataset as the base 
    dists_nf_mc[m_i+1] <- W2(is_samps$samps, n_draws_flat)
    # Measure distance from informative prior fitted to original dataset
    # to new flat prior model fitted to expanded dataset
    dists_nc_mf[m_i+1] <- W2(n_post_conj, m_draws_flat)
  }
 
  dists <- rbind(dists_nf_mc, dists_nc_mf)
  dists_mins <- apply(dists,2,which.min)
  dists_arg_min <- which.min(apply(dists,2,min))
  dists_arg_sign <- dists_mins[dists_arg_min]
  PSS_sign <- ifelse(dists_arg_sign == 2, 1, -1)
  PSS <- (dists_arg_min - 1) * PSS_sign
  return(
    list(
      PSS = PSS,
      seed = seed,
      theta_draw = theta_now,
      y_seq = mean_seq,
      var_weights = var_weights
    )
  )
}

mod <- cmdstan_model('objective-prior-poisson-lik.stan')

seed <- 123 + n
set.seed(seed)

y1ton <- readRDS('pois_draws_n_20_mean_1.35.RDS')

arg_list_refit_both_models_W2_pp_pi_nc_bb_signed <- list(
  fun_to_apply = function(y1ton) {
    return(refit_both_models_signed_bb(
      y1ton
      , 2
      , u_0
      , mod
      , lower = NULL
      , upper = 100
      , n_draws = 1e4
      , seed
    )
    )
  },
  upper = 100
)

run <- arg_list_refit_both_models_W2_pp_pi_nc_bb_signed$fun_to_apply(y1ton)

fl_nm <- paste0('run_',n,'_u0_',u_0,'.RDS')
saveRDS(run, file = fl_nm)
