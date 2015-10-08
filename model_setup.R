model_path = '~/Documents/Statistics/Sparse_factor_G_matrix/R_version'
source(paste(model_path,'fast_BSFG_sampler_init.R',sep='/'))
source(paste(model_path,'fast_BSFG_sampler.R',sep='/'))
source(paste(model_path,'BSFG_functions.R',sep='/'))
source(paste(model_path,'plotting_diagnostics.R',sep='/'))

library(Rcpp)
library(RcppArmadillo)

sourceCpp(paste(model_path,'BSFG_functions_c.cpp',sep='/'))


setwd('/Users/der7/Documents/Statistics/Sparse_factor_G_matrix/Submission2/Analyses/Website/BSF-G_new/Runcie_Mukherjee_analyses/Simulations/Sim_2')

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
    burn = 100,
    thin = 2,
    draw_iter = 200
    )

priors = list(
    k_init              =   20,
    resid_Y_prec_shape  =   2,
    resid_Y_prec_rate   =   1/2,
    E_a_prec_shape      =   2,
    E_a_prec_rate       =   1/2,
    W_prec_shape        =   2,
    W_prec_rate         =   1/2,
    Lambda_df           =   1.5,
    delta_1_shape       =   2.1,
    delta_1_rate        =   1/20,
    delta_2_shape       =   3,
    delta_2_rate        =   1,
    # h2_priors_factors   =   (1+c(run_parameters$h2_divisions-2,rep(0,run_parameters$h2_divisions-1)/(-2+2*run_parameters$h2_divisions),
    # h2_priors_resids    =   c(0,dt((2:run_parameters$h2_divisions)/run_parameters$h2_divisions,1)),
    h2_priors_resids    =   c(0,rep(1,run_parameters$h2_divisions-1)),
    # h2_priors_resids    =   dt((1:run_parameters$h2_divisions)/run_parameters$h2_divisions,.1),
    h2_priors_factors   =   dt((1:run_parameters$h2_divisions)/run_parameters$h2_divisions,.1)
)

print('Initializing')
save(priors,file = 'Priors.RData')
# Initialize Chain, prep runs
BSFG_state = fast_BSFG_sampler_init(priors,run_parameters)

BSFG_state = clear_Posterior(BSFG_state)


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
    BSFG_state = fast_BSFG_sampler(BSFG_state,n_samples)
    print(i)
}
