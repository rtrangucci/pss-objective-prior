fls <- list.files('data/run_01/',full.names = T)
pss_01 <- rep(NA_real_, length(fls))
min_ess_01 <- rep(NA_real_, length(fls))
for (fl_i in seq_along(fls)) {
  fl <- fls[fl_i] 
  pss_run <- readRDS(fls[fl_i])
  pss_01[fl_i] <- pss_run$PSS
  min_ess_01[fl_i] <- min(1/pss_run$var_weights)
}

fls <- list.files('data/run_65/',full.names = T)
pss_65 <- rep(NA_real_, length(fls))
min_ess_65 <- rep(NA_real_, length(fls))
for (fl_i in seq_along(fls)) {
  pss_run <- readRDS(fls[fl_i])
  pss_65[fl_i] <- pss_run$PSS
  min_ess_65[fl_i] <- min(1/pss_run$var_weights)
}

xlims <- range(c(pss_01,pss_65))

par(mfrow = c(2,1))
pdf('output/opess-comparions.pdf',width = 4, height = 6)
hist(pss_01,breaks=20, main = bquote('OPESS '*u(0)*'=0.10'),
     xlab = 'PSS', xlim=xlims)
hist(pss_65,breaks=50, main = bquote('OPESS '*u(0)*'=0.65'),
     xlab = bquote(''*lambda*''), xlim=xlims)
dev.off()