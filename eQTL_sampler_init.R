library(R.matlab)
eQTL_sampler_init = function(priors,run_parameters){
    # require(PEIP)
    # require(Matrix)

    # data_matrices,run_parameters,priors,current_state,Posterior,simulation = F)
#   Preps everything for an analysis with the function MA_sampler. 
#       Loads data from setup.mat
#       Decides if the data is from a simulation
#       Uses prior to find initial values for all parameters
#       Initializes Posterior struct
#       Saves run parameters
#
#   Input file setup.mat should contain:
#     Y      gene expression data data: n x p
#     Z_Y    genes-probes incidence matrix: g x p
#     X_E_a2    fixed-effect design matrix for transcripts
#     Z_E_a2    lines-replicates incidence matrix: n x r
#     A      kinship matrix of lines: r x r (should be identify)
#     optional:
#         U_act
#         gen_factor_Lambda
#         error_factor_Lambda
#         h2_act
#         G
#         R
#         B
#         E_a2_act
#         factor_h2s
#         name

# ----------------------- #
# ------read data-------- #
# ----------------------- #
    # load('../setup.RData')
    load('../Sim_data.RData')

    Y            = sim_data$Y
    X            = sim_data$X
    Xf           = sim_data$X
    Z_1          = sim_data$Z1
    Z_2f         = sweep(Z_1,1,Xf,'*')
    cisGenotypes = sim_data$cisGenotypes
    A            = sim_data$A1

    n = nrow(Y)
    p = ncol(Y)
    r = ncol(Z_1)
    r2f = ncol(Z_2f)

    sim_data$gen_factor_Lambda = sim_data$Lambda
    sim_data$error_factor_Lambda = sim_data$Lambda

    #Determine if 'setup.mat' contains output of a simulation, based on if
    #known factor loadings are included. Affects plotting functions
    simulation = F
    if('gen_factor_Lambda' %in% names(sim_data)){
        simulation = T
        print('simulated data')
        run_parameters$setup = sim_data
    }
    run_parameters$simulation = simulation

    #normalize Y to have zero mean and unit variances among observed values,
    #allowing for NaNs.

    Mean_Y = colMeans(Y,na.rm=T)
    VY = apply(Y,2,var,na.rm=T)
    # Don't remove the mean and standardize the variance if it's a
    # simulation because this makes comparing to the simulated values more
    # difficult. 
    # Do we want to do this for real data, or let the user do it?
    if(simulation) {
        VY = rep(1,p)
    } else{
        Y = sweep(Y,2,Mean_Y,'-')
        Y = sweep(Y,2,sqrt(VY),'/')
    }
    
    Y_full = Y
    
    
#     X = zeros(n,0)
#     X = [ones(size(X_f,1),1) X_f]
    
    #determine if a design matrix (X) exists (is loaded from setup.mat). If
    #not, make a dummy X-matrix with no columns.
    if(! 'X' %in% ls()) {
        X=matrix(0,nr = n,nc = 0)
    }
    if(nrow(X) == 1) X = t(X)       # because I used to give X transposed
    stopifnot(nrow(X) == n)
    b = ncol(X)
    
    
    #Determine if a second random effects design matrix exists. If not, make a
    #dummy matrix
    if(! 'Z_2' %in% ls()) {
        Z_2=matrix(0,nr = n,nc = 0)
    }
    stopifnot(nrow(Z_2) == n)
    r2 = ncol(Z_2)

    bf = ncol(Xf)
    r2f = ncol(Z_2f)

    data_matrices = list(
            Y_full       = Y_full,
            Z_1          = Z_1,
            Z_2          = Z_2,
            X            = X,
            Z_2f         = Z_2f,
            Xf           = Xf,
            cisGenotypes = cisGenotypes
            )
    
# ----------------------------- #
# -----Initialize variables---- #
# ----------------------------- # 

   # --- transcript-level model
    # p-vector of probe residual precisions. 
    #  Prior: Gamma distribution for each element
    #       shape = resid_Y_prec_shape
    #       rate = resid_Y_prec_rate
    resid_Y_prec_shape  = priors$resid_Y_prec_shape
    resid_Y_prec_rate   = priors$resid_Y_prec_rate
    resid_Y_prec        = rgamma(p,shape = resid_Y_prec_shape,rate = resid_Y_prec_rate)
   
    # Factors:
    #  initial number of factors
    k = priors$k_init

    # Factor loading precisions (except column penalty tauh).
     #  Prior: Gamma distribution for each element. 
     #       shape = Lambda_df/2
     #       rate = Lambda_df/2
     #    Marginilizes to t-distribution with Lambda_df degrees of freedom
     #    on each factor loading, conditional on tauh
    Lambda_df   = priors$Lambda_df
    Lambda_prec = matrix(rgamma(p*k,shape = Lambda_df/2,rate = Lambda_df/2),nr = p,nc = k)
    
    # Factor penalty. tauh(h) = \prod_{i=1}^h \delta_i
     #  Prior: Gamma distribution for each element of delta
     #     delta_1:
     #       shape = delta_1_shape
     #       rate = delta_1_rate
     #     delta_2 ... delta_m:
     #       shape = delta_2_shape
     #       rate = delta_2_rate
    delta_1_shape  = priors$delta_1_shape
    delta_1_rate   = priors$delta_1_rate
    delta_2_shape  = priors$delta_2_shape
    delta_2_rate   = priors$delta_2_rate
    delta          = c(rgamma(1,shape = delta_1_shape,rate = delta_1_rate),rgamma(k-1,shape = delta_2_shape,rate = delta_2_rate))
    tauh           = cumprod(delta)
    
    # Total Factor loading precisions Lambda_prec * tauh
    Plam = sweep(Lambda_prec,2,tauh,'*')
    
    # Lambda - factor loadings
     #   Prior: Normal distribution for each element.
     #       mu = 0
     #       sd = sqrt(1/Plam)
    Lambda = matrix(rnorm(p*k,0,sqrt(1/Plam)),nr = p,nc = k)

    # g-vector of specific precisions of genetic effects. 
    #  Prior: uniform on 0,1
    resid_h2 = runif(p,0,1)
        
    # Genetic effects not accounted for by factors.
    #   Prior: Normal distribution on each element.
    #       mean = 0
    #       sd = 1./sqrt(E_a_prec)' on each row
    E_a = matrix(rnorm(p*r,0,sqrt(resid_h2/resid_Y_prec)),nr = r,nc = p, byrow = T)
        
    # Latent factor heritabilties. h2 can take h2_divisions values
    #   between 0 and 1.
     #   Prior: 0.5: h2=0, .05: h2 > 0. 
    F_h2 = runif(k)
    
    # Fixed effects on the factor scores.
    #  Prior: Normal distribution for each element
    #       mean = 0
    #       sd = sqrt(F_h2') for each row.
    B_F = matrix(rnorm(bf*k,0,1),bf,k, byrow = T)
    
    # Genetic effects on the factor scores.
    #  Prior: Normal distribution for each element
    #       mean = 0
    #       sd = sqrt(F_h2') for each row.
    F_a = matrix(rnorm(k*r,0,sqrt(F_h2)),nr = r,nc = k, byrow = T)

    prec_F_resid = rgamma(k,shape = priors$prec_F_resid_shape,rate = priors$prec_F_resid_rate)
    prec_B_F = rgamma(k,shape = priors$prec_B_F_shape,rate = priors$prec_B_F_rate)
    prec_E_a2_F = rgamma(k,shape = priors$prec_E_a2_F_shape,rate = priors$prec_E_a2_F_rate)
    
    # Full Factor scores. Combination of genetic and residual variation on
    # each factor.
    #  Prior: Normal distribution for each element
    #       mean = Z_1 * F_a
    #       sd = sqrt(1-F_h2') for each row.
    F = X %*% B_F + Z_1 %*% F_a + matrix(rnorm(k*n,0,sqrt(1-F_h2)),nr = n,nc = k, byrow = T)    
    
    # g-vector of specific precisions of genetic effects. 
    #  Prior: Gamma distribution for each element
    #       shape = E_a_prec_shape
    #       rate = E_a_prec_rate
    E_a2_prec_shape = priors$E_a2_prec_shape
    E_a2_prec_rate  = priors$E_a2_prec_rate
    E_a2_prec       = rgamma(p,shape = E_a2_prec_shape,rate = E_a2_prec_rate)
    if(r2 == 0) {
        E_a2_prec = matrix(0,nr=0,nc=p)
    }
    
    # Genetic effects not accounted for by factors.
    #   Prior: Normal distribution on each element.
    #       mean = 0
    #       sd = 1./sqrt(E_a_prec)' on each row
    E_a2 = matrix(rnorm(p*r2,0,sqrt(1/E_a2_prec)),nr = r2, nc = p, byrow = T)
    E_a2_F = matrix(rnorm(k*r2f,0,sqrt(1/prec_E_a2_F)),nr = r2f, nc = k, byrow = T)

    # Fixed effect coefficients.
    #  Prior: Normal distribution for each element
    #       mean = 0
    #       sd = sqrt(1/fixed_effects_prec)
    B_resid = matrix(rnorm(b*p),nr = b, nc = p)
    cis_effects = matrix(rnorm(p),nr = 1, nc = p)
    mu = rnorm(p)

    B_prec = matrix(rgamma(b*p,3,3),b,p)
    B_shape = rgamma(1,3,3)

    
# ----------------------- #
# -Initialize Posterior-- #
# ----------------------- #
    Posterior = list(
            Lambda       = matrix(0,nr=0,nc=0),
            F            = matrix(0,nr=0,nc=0),
            B_F          = matrix(0,nr=0,nc=0),
            E_a2_F       = matrix(0,nr=0,nc=0),
            F_a          = matrix(0,nr=0,nc=0),
            mu           = matrix(0,nr=p,nc=0),
            B_resid      = matrix(0,nr=b*p,nc=0),
            cis_effects  = matrix(0,nr=p,nc=0),
            delta        = matrix(0,nr=0,nc=0),
            B_shape      = matrix(0,nr=1,nc=0),
            F_h2         = matrix(0,nr=0,nc=0),
            prec_B_F     = matrix(0,nr=0,nc=0),
            prec_E_a2_F  = matrix(0,nr=0,nc=0),
            prec_F_resid = matrix(0,nr=0,nc=0),
            resid_Y_prec = matrix(0,nr=p,nc=0),
            resid_h2     = matrix(0,nr=p,nc=0),
            E_a_prec     = matrix(0,nr=p,nc=0),
            E_a2_prec    = matrix(0,nr=p,nc=0),
            E_a          = matrix(0,nr=r,nc=p),
            E_a2         = matrix(0,nr=r2,nc=p)
        )
# ----------------------- #
# ---Save initial values- #
# ----------------------- #
    current_state = list(
                            resid_Y_prec = resid_Y_prec,
                            B_prec       = B_prec,
                            E_a2_prec    = E_a2_prec,
                            resid_h2     = resid_h2,
                            Lambda_prec  = Lambda_prec,
                            Plam         = Plam,
                            delta        = delta,
                            tauh         = tauh,
                            B_shape      = B_shape,
                            cis_effects  = cis_effects,
                            prec_B_F     = prec_B_F,
                            prec_E_a2_F  = prec_E_a2_F,
                            prec_F_resid = prec_F_resid,
                            F_h2         = F_h2,
                            Lambda       = Lambda,
                            F            = F,
                            mu           = mu,
                            B_resid      = B_resid,
                            E_a2         = E_a2,
                            B_F          = B_F,
                            E_a2_F       = E_a2_F,
                            F_a          = F_a,
                            E_a          = E_a,
                            nrun         = 0
                        ) 


# ------------------------------------ #
# ----Precalculate some matrices------ #
# ------------------------------------ #

    # recover()
    #invert the random effect covariance matrices
    Ainv = solve(A)
    A_2_inv = diag(1,r2) #Z_2 random effects are assumed to have covariance proportional to the identity. Can be modified.

    #pre-calculate transformation parameters to diagonalize aI + bZAZ for fast
    #inversion: inv(aI + bZAZ) = 1/b*U*diag(1./(s+a/b))*U'
    #uses singular value decomposition of ZAZ for stability when ZAZ is low
    
    result = svd(Z_1 %*% A %*% t(Z_1))
    invert_aI_bZAZ = list(
        U = result$u,
        s = result$d
    )

    L = t(chol(A))
    Li = solve(L)
    svd_LZZL = svd(Li %*% t(Z_1) %*% Z_1 %*% t(Li))
    Q = svd_LZZL$u
    sample_a_mats = list(
                            QtLiZt = t(Q) %*% Li %*% t(Z_1),
                            d = svd_LZZL$d,
                            LitQ = t(Li) %*% Q
                        )
 
# ----------------------------- #
# ----Save run parameters------ #
# ----------------------------- #

    run_variables = list(
            p              = p,
            n              = n,
            r              = r,
            r2             = r2,
            b              = b,
            r2f            = r2f,
            bf             = bf,
            Mean_Y         = Mean_Y,
            VY             = VY,
            Ainv           = Ainv,
            A_2_inv        = A_2_inv,
            invert_aI_bZAZ = invert_aI_bZAZ,
            sample_a_mats  = sample_a_mats
    )

    RNG = list(
        Random.seed = .Random.seed,
        RNGkind = RNGkind()
    )

 return(list(
            data_matrices  = data_matrices,
            run_parameters = run_parameters,
            run_variables  = run_variables,
            priors         = priors,
            current_state  = current_state,
            Posterior      = Posterior,
            simulation     = simulation,
            RNG            = RNG
        ))
}
