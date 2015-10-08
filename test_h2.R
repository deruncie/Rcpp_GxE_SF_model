library(mnormt)
Z = BSFG_state$data$Z_1
A = BSFG_state$run_pa$setup$A
n = nrow(Z)

ZAZ = Z %*% A %*% t(Z)
Eye = diag(1,n)

p = 3

resid_prec = rgamma(p,10,1/10)
# E_prec = rgamma(p,10,1/10)
h2 = seq(0,.9,length=p)
E_prec = 1/((h2/resid_prec)/(1-h2))

h2_divisions = 50
h2_priors = rep(1,h2_divisions)

Y = sapply(1:p,function(x) rmnorm(1,rep(0,n),1/E_prec[x]*ZAZ + 1/resid_prec[x]*Eye))

set.seed(1);est = sample_prec_discrete_conditional_c(Y,h2_divisions,h2_priors,BSFG_state$run_var$invert_aI_bZAZ,resid_prec)
set.seed(1);est2 = sample_prec_discrete_conditional(Y,h2_divisions,h2_priors,BSFG_state$run_var$invert_aI_bZAZ,resid_prec)
# plot(h2,est2)
plot(est,E_prec)

ests = sapply(1:100,function(x) sample_prec_discrete_conditional_c(Y,h2_divisions,h2_priors,BSFG_state$run_var$invert_aI_bZAZ,resid_prec))
plot(rowMeans(1/ests),1/E_prec)
# plot(est,E_prec)