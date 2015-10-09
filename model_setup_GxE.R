model_path = '~/Documents/Git_repositories/Rcpp_GxE_SF_model'
source(paste(model_path,'GxE_sampler_init.R',sep='/'))
source(paste(model_path,'GxE_sampler.R',sep='/'))
source(paste(model_path,'GxE_functions.R',sep='/'))
source(paste(model_path,'plotting_diagnostics_GxE.R',sep='/'))

library(Rcpp)
library(RcppArmadillo)

sourceCpp(paste(model_path,'new_BSFG_functions_c.cpp',sep='/'))


setwd('Sim_FE4_1')

seed = 1234
set.seed(seed)

rep = 2
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
    burn         = 100,
    thin         = 2,
    draw_iter    = 200
    )

priors = list(
    k_init             =   20,
    resid_Y_prec_shape =   2,
    resid_Y_prec_rate  =   1/2,
    prec_B_F_shape     =   2,
    prec_B_F_rate      =   1/2,
    prec_F_resid_shape =   2,
    prec_F_resid_rate  =   1/2,
    E_a2_prec_shape    =   2,
    E_a2_prec_rate     =   1/2,
    Lambda_df          =   3,
    delta_1_shape      =   2.1,
    delta_1_rate       =   1/20,
    delta_2_shape      =   3,
    delta_2_rate       =   1,
    B_df               =   5,
    B_shape_shape      =   2.1,
    B_shape_rate       =   1/20,
    h2_priors_resids   =   c(0,rep(1,run_parameters$h2_divisions-1)),
    h2_priors_factors  =   dt((1:run_parameters$h2_divisions)/run_parameters$h2_divisions,.1)
)

print('Initializing')
save(priors,file = 'Priors.RData')
# Initialize Chain, prep runs
GxE_state = GxE_sampler_init(priors,run_parameters)

GxE_state = clear_Posterior(GxE_state)


# # optional: To load from end of previous run, run above code, then run these lines:
# load('current_state')
# load('Posterior')
# load('Priors')
# start_i = run_parameters$nrun;

# Run Gibbs sampler. Run in smallish chunks. Output can be used to re-start chain where it left off.
# c1=fix(clock);
n_samples = 200;
for(i  in 1:10) {
    print(sprintf('Run %d',i))
source(paste(model_path,'plotting_diagnostics_GxE.R',sep='/'))
# sourceCpp(paste(model_path,'new_GxE_functions_c.cpp',sep='/'))
source(paste(model_path,'GxE_sampler.R',sep='/'))
    GxE_state = GxE_sampler(GxE_state,n_samples)
    print(i)
}
