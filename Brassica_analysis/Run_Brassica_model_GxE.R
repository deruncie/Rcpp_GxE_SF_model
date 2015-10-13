model_path = '~/Documents/Git_repositories/Rcpp_GxE_SF_model'
source(paste(model_path,'GxE_sampler_init.R',sep='/'))
source(paste(model_path,'GxE_sampler.R',sep='/'))
source(paste(model_path,'GxE_functions.R',sep='/'))
source(paste(model_path,'plotting_diagnostics_GxE.R',sep='/'))

library(Rcpp)
library(RcppArmadillo)

sourceCpp(paste(model_path,'new_BSFG_functions_c.cpp',sep='/'))

seed = 1234
set.seed(seed)

rep = 3
folder = sprintf('R_rep_%d',rep)
try(dir.create(folder))
setwd(folder)

# initialize priors



run_parameters = list(
    b0           = 1,
    b1           = 0.0005,
    epsilon      = 1e-2,
    prop         = 1.00,
    h2_divisions = 50,
    save_freq    = 100,
    burn         = 1000,
    thin         = 10,
    draw_iter    = 200
    )

priors = list(
    k_init             =   20,
    resid_Y_prec_shape =   2,
    resid_Y_prec_rate  =   1/2,
    prec_B_F_shape     =   2,
    prec_B_F_rate      =   1,
    prec_F_resid_shape =   2,
    prec_F_resid_rate  =   1,
    E_a2_prec_shape    =   2,
    E_a2_prec_rate     =   1/2,
    Lambda_df          =   3,
    delta_1_shape      =   2.1,
    delta_1_rate       =   1/5,
    delta_2_shape      =   5,
    delta_2_rate       =   1,
    B_df               =   1,
    B_shape_shape      =   2.1,
    B_shape_rate       =   1/5,
    h2_priors_resids   =   c(0,rep(1,run_parameters$h2_divisions-1)),
    h2_priors_factors  =   dt((1:run_parameters$h2_divisions)/run_parameters$h2_divisions,.1)
)

# state = GxE_state$current_state

print('Initializing')
save(priors,file = 'Priors.RData')
# Initialize Chain, prep runs
GxE_state = GxE_sampler_init(priors,run_parameters)

GxE_state = clear_Posterior(GxE_state)

# Run Gibbs sampler. Run in smallish chunks. Output can be used to re-start chain where it left off.
# c1=fix(clock);
n_samples = 200;
for(i  in 1:10) {
    print(sprintf('Run %d',i))
    GxE_state = GxE_sampler(GxE_state,n_samples)
    print(i)
}


Posterior = GxE_state$Posterior
p = ncol(GxE_state$data_matrices$Y_full)
n = nrow(GxE_state$data_matrices$Y_full)
k = nrow(Posterior$Lambda)/p
I = ncol(Posterior$Lambda)
par(mfrow=c(4,4))
for(i in 1:k){
    j = p*(i-1)+1:p
    mean = rowMeans(Posterior$Lambda[j,])
    cors = c(cor(mean,Posterior$Lambda[j,])^2)
    plot(cors,type='l',ylim=c(0,1),main = i)
}
sizes = apply(Posterior$Lambda,2,function(x) tapply(x,gl(k,p),function(y) mean(y^2)))
par(mfrow=c(1,1))
boxplot(t(sizes))
boxplot(t(Posterior$F_h2))
boxplot(t(Posterior$B_F))
boxplot(t(apply(Posterior$delta,2,function(x) 1/cumprod(x))))
plot(rowMeans(Posterior$F_h2),rowMeans(Posterior$B_F))
image(t(GxE_state$current_state$Lambda),col=palette(viridis(250)))
hist(abs(GxE_state$current_state$Lambda[,1]),breaks=100)
Lambda = matrix(0,p,k)
for(i in 1:k){
    j = p*(i-1)+1:p
    Lambda = Lambda + matrix(Posterior$Lambda,p,k)/I
}
Lambda_thresh = apply(Lambda,2,function(x) {x[abs(x)<quantile(abs(x),.5)]=0;x})
# for(i in k:1){
#     Lambda = Lambda[order(Lambda_thresh[,k]),]
#     Lambda_thresh = Lambda_thresh[order(Lambda_thresh[,k]),]
# }
Lambda = Lambda[order(Lambda_thresh[,1],Lambda_thresh[,2],Lambda_thresh[,3],Lambda_thresh[,4],Lambda_thresh[,5],Lambda_thresh[,6]),]
# Lambda = Lambda[order(lapply(1:k,function(x) Lambda_thresh[,k])),]
image(t((Lambda)),col=viridis(250))


