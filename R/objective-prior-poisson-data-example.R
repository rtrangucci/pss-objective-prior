library(cmdstanr)
library(loo)
mod <- cmdstan_model('stan/objective-prior-poisson-lik.stan')
y <- readRDS('data/pois_draws_n_20_mean_1.35.RDS')
dat <- list(c = 2, u_0 = 0.65, acc = 1e-11)
dat$y <- y
dat$N <- length(y)
fit <- mod$sample(dat, 
                  iter_warmup =  8000,
                  iter_sampling =  2500 * 40,
                  parallel_chains = 4,
                  chains = 4,
                  adapt_delta = 0.99,
                  max_treedepth = 12,
                  inv_metric = as.array(0.0407221),
                  step_size = 0.3,
                  adapt_engaged = T,
                  refresh = 5e3,
                  thin = 40,
                  seed = 123,
                  init = function() list(theta = runif(1,min=0.5,max=1.0))
        )
d <- fit$draws()

dat <- list(c = 2, u_0 = 0.1, acc = 1e-8)
dat$y <- y
dat$N <- length(y) 
fit <- mod$sample(dat, 
                  iter_warmup =  8000,
                  iter_sampling =  2500 * 10,
                  parallel_chains = 4,
                  chains = 4,
                  adapt_delta = 0.99,
                  max_treedepth = 12,
                  inv_metric = as.array(0.0407221),
                  step_size = 0.2,
                  adapt_engaged = T,
                  thin = 10,
                  refresh = 5e3,
                  seed = 123,
                  init = function() list(theta = runif(1,min=0.5,max=2.0))
    )
d_01 <- fit$draws()

pdf('output/post-comparison.pdf',width = 4, height = 6)
par(mfrow=c(2,1))
hist(d_01[,,'theta'],freq = F,breaks=100, main = bquote('p('*theta*'|y), '*u(0)*'=0.10'),
     xlab = bquote(''*theta*''), xlim=range(c(d[,,'theta'],d_01[,,'theta'])))
hist(d[,,'theta'],freq = F,breaks=50, main = bquote('p('*theta*'|y), '*u(0)*'=0.65'),
     xlab = bquote(''*theta*''), xlim=range(c(d[,,'theta'],d_01[,,'theta'])))
dev.off()