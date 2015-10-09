GxE_sampler = function(GxE_state,n_samples) {
	# -- Daniel Runcie -- #
	#
	# Gibbs sampler for genetic covariance estimation based on mixed effects
	# model, with missing data
	#
	# Code for:
	# 
	# Runcie and Mukherjee (2013) Dissecting high-dimensional traits
	# with Bayesian sparse factor analysis of genetic covariance matrices.
	# GENETICS.
	#
	# (c) July 30, 2013
	#
	# code based on original provided by Anirban Bhattacharya
	#
	#         
	# This function implements the BSF-G partially collapsed Gibbs sampler.
	# Variable notation follows the Runcie and Mukherjee manuscript as closely
	# as possible.
	#
	# All input and initialization functions are carried out by the function
	# MA_sampler_init. See the documentation of that function for details.
	# 
	# The sampler is designed to do a short-medium length run and then return
	# the state of the chain. After assessing the progress of the chain,
	# (optionally), the Posterior matrices can be reset, and the chain
	# re-started where it left off. Right now, the random seed is not saved. Probably should
	# be if I can figure out how.
	# 
	# This function takes the following inputs:
	#     data_matrices: struct holding the (should be imutable) data, design and incidence matrices:
	#          Y_full: Full probe data, may include NaNs. n x p
	#          Z_Y:    Incidence matrix for transcripts -> probes. g x p
	#          Z_E_a2:    Incidence matrix for lines -> replicates. n x r
	#          X_E_a2:    Fixed effect design matrix. n x p
	#     start_i: iteration number of the end of the last run.
	#     draw_iter: frequency of updating diagnostic plots
	#     burn: number of burnin samples
	#     sp: total number of samples to collect
	#     thin: thinning rate of chain
	#     simulation: boolean. Is this a simulation?
	#     params: struct with chain parameters.
	#     priors: struct with all relevant prior hyperparameters
	#     Posterior: struct with posterior matrices, or running posterior means. 
	#            Note: not sure how well posterior means work after re-starting chain. Probably not well.
	#     current_state: current (initial) conditions of all model parameters
	# 
	# Several diagnostic plots are produced during the run. 
	#     Their interpretation is described within the source codes:
	#         draw_simulation_diagnostics.m: For simulated data with known true values
	#         draw_results_diagnostics.m: Otherwise
	#         
	# This is a new version of the BSFG sampler with the sampling generally going in the opposite direction (across genes instead of individuals)
	# may not ever sample E_a for F_a


	data_matrices  = GxE_state$data_matrices
	params         = GxE_state$params
	priors         = GxE_state$priors
	Posterior      = GxE_state$Posterior
	current_state  = GxE_state$current_state
	run_parameters = GxE_state$run_parameters
	run_variables  = GxE_state$run_variables


	# ----------------------------------------------- #
	# ----------------Load data matrices------------- #
	# ----------------------------------------------- #

	Y_full  = data_matrices$Y_full
	Z_1     = data_matrices$Z_1
	Z_2     = data_matrices$Z_2
	X       = data_matrices$X 			# Note: Do not include global mean in X. Just treatment contrasts

	p   = ncol(Y_full)
	n   = nrow(Y_full)
	r   = ncol(Z_1)
	r2  = ncol(Z_2)
	b   = ncol(X)


	# ----------------------------------------------- #
	# ----------------Load priors-------------------- #
	# ----------------------------------------------- #

	resid_Y_prec_shape = priors$resid_Y_prec_shape
	resid_Y_prec_rate  = priors$resid_Y_prec_rate
	prec_B_F_shape     = priors$prec_B_F_shape
	prec_B_F_rate      = priors$prec_B_F_rate
	prec_F_resid_shape = priors$prec_F_resid_shape
	prec_F_resid_rate  = priors$prec_F_resid_rate
	E_a2_prec_shape    = priors$E_a2_prec_shape
	E_a2_prec_rate     = priors$E_a2_prec_rate
	Lambda_df          = priors$Lambda_df
	delta_1_shape      = priors$delta_1_shape
	delta_1_rate       = priors$delta_1_rate
	delta_2_shape      = priors$delta_2_shape
	delta_2_rate       = priors$delta_2_rate
	B_df               = priors$B_df
	B_shape_shape      = priors$B_shape_shape
	B_shape_rate       = priors$B_shape_rate
	h2_priors_factors  = priors$h2_priors_factors
	h2_priors_resids   = priors$h2_priors_resids


	# ----------------------------------------------- #
	# ----------------Load current state------------- #
	# ----------------------------------------------- #
	resid_Y_prec = current_state$resid_Y_prec
	B_prec       = current_state$B_prec      
	E_a2_prec    = current_state$E_a2_prec   
	resid_h2     = current_state$resid_h2    
	Lambda_prec  = current_state$Lambda_prec 
	Plam         = current_state$Plam        
	delta        = current_state$delta       
	tauh         = current_state$tauh       
	B_shape      = current_state$B_shape        
	prec_B_F     = current_state$prec_B_F    
	prec_F_resid = current_state$prec_F_resid
	F_h2         = current_state$F_h2        
	Lambda       = current_state$Lambda      
	F            = current_state$F           
	mu           = current_state$mu          
	B_resid      = current_state$B_resid           
	E_a2         = current_state$E_a2        
	B_F          = current_state$B_F         
	F_a          = current_state$F_a         
	E_a          = current_state$E_a   
	start_i      = current_state$nrun      

	k = ncol(F)

	# ----------------------------------------------- #
	# -----------Reset Global Random Number Stream--- #
	# ----------------------------------------------- #
	do.call("RNGkind",as.list(GxE_state$RNG$RNGkind))  ## must be first!
	assign(".Random.seed", GxE_state$RNG$Random.seed, .GlobalEnv)

	# ----------------------------------------------- #
	# -----------Load pre-calcualted matrices-------- #
	# ----------------------------------------------- #
	Ainv                             = run_variables$Ainv
	A_2_inv                          = run_variables$A_2_inv
	invert_aI_bZAZ                   = run_variables$invert_aI_bZAZ
	Q = invert_aI_bZAZ$U
	d = invert_aI_bZAZ$s
	Qt = t(Q)
	QtX = Qt %*% X
	QtXZ2 = Qt %*% cbind(1,X,Z_2)
	Y = Y_full
	QtY = Qt %*% Y

	# ----------------------------------------------- #
	# ----------------Set up run--------------------- #
	# ----------------------------------------------- #
	#     b0,b1: parameters controlling rate of adaptation of factor model size
	#     h2_divisions: number of discrete steps for each factor heritability parameter
	#     epsilon: truncation point for factor loadings during adaptation
	b0           = run_parameters$b0
	b1           = run_parameters$b1
	h2_divisions = run_parameters$h2_divisions
	epsilon      = run_parameters$epsilon
	prop         = run_parameters$prop
	save_freq    = run_parameters$save_freq
	burn         = run_parameters$burn
	thin         = run_parameters$thin
	draw_iter    = run_parameters$draw_iter
	simulation   = run_parameters$simulation


	# ----------------------------------------------- #
	# ---Extend posterior matrices for new samples--- #
	# ----------------------------------------------- #

	sp = (start_i + n_samples - burn)/thin - ncol(Posterior$Lambda)
	if(sp > 0){
		Posterior$Lambda       = cbind(Posterior$Lambda,matrix(0,nr = nrow(Posterior$Lambda),nc = sp))
		Posterior$F            = cbind(Posterior$F,matrix(0,nr = nrow(Posterior$F),nc = sp))
		Posterior$B_F          = cbind(Posterior$B_F,matrix(0,nr = nrow(Posterior$B_F),nc = sp))
		Posterior$F_a          = cbind(Posterior$F_a,matrix(0,nr = nrow(Posterior$F_a),nc = sp))
		Posterior$mu           = cbind(Posterior$mu,matrix(0,nr = nrow(Posterior$mu),nc = sp))
		Posterior$B_resid      = cbind(Posterior$B_resid,matrix(0,nr = nrow(Posterior$B_resid),nc = sp))
		Posterior$delta        = cbind(Posterior$delta,matrix(0,nr = nrow(Posterior$delta),nc = sp))
		Posterior$B_shape      = c(Posterior$B_shape,rep(0,nc = sp))		
		Posterior$F_h2         = cbind(Posterior$F_h2,matrix(0,nr = nrow(Posterior$F_h2),nc = sp))
		Posterior$prec_B_F     = cbind(Posterior$prec_B_F,matrix(0,nr = nrow(Posterior$prec_B_F),nc = sp))
		Posterior$prec_F_resid = cbind(Posterior$prec_F_resid,matrix(0,nr = nrow(Posterior$prec_F_resid),nc = sp))
		Posterior$resid_Y_prec = cbind(Posterior$resid_Y_prec,matrix(0,nr = nrow(Posterior$resid_Y_prec),nc = sp))
		Posterior$resid_h2     = cbind(Posterior$resid_h2,matrix(0,nr = nrow(Posterior$resid_h2),nc = sp))
		Posterior$E_a_prec     = cbind(Posterior$E_a_prec,matrix(0,nr = nrow(Posterior$E_a_prec),nc = sp))
		Posterior$E_a2_prec    = cbind(Posterior$E_a2_prec,matrix(0,nr = nrow(Posterior$E_a2_prec),nc = sp))
	}

	# ----------------------------------------------- #
	# --------------start gibbs sampling------------- #
	# ----------------------------------------------- #

	start_time = Sys.time()
	for(i in start_i+(1:n_samples)){
	   
	 # -----fill in missing phenotypes----- #
		#conditioning on everything else
		# phenMissing = is.na(Y_full)
		# if(sum(phenMissing)>0) {
		# 	Y = Y_full
		# 	meanTraits = cbind(1,X) %*% rbind(mu,B) + F %*% t(Lambda) + Z_1 %*% E_a + Z_2 %*% E_a2
		# 	resids = matrix(rnorm(p*n,0,sqrt(1/resid_Y_prec)),nr = n,nc = p,byrow=T)
		# 	Y(phenMissing) = meanTraits(phenMissing) + resids(phenMissing)
		# 	QtY = Qt %*% Y
		# }
		# recover()
		  
	 # -----Sample B, Lambda, E_a2 ------------------ #
		# recover()
		QtW = Qt %*% cbind(1,X,F,Z_2)
		QtY_prec = t(resid_Y_prec/(matrix(resid_h2,nc=1) %*% matrix(d,nr=1) + (1-c(resid_h2))))
		# prior_prec = matrix(rbind(0,B_prec,t(Plam), matrix(E_a2_prec,nr=r2,nc=p)),nr = 1+b+k+r2, ncol = p)
		prior_prec = rbind(0,B_prec,t(Plam), matrix(E_a2_prec,nr=r2,nc=p))
		prior_mean = matrix(0,nr = 1+b+k+r2, ncol = p)
		result_mat = sample_coefs_c(QtY, QtW, QtY_prec, prior_mean,prior_prec)
		mu_row = 1
		b_rows = c()
		lambda_rows = c()
		e_a2_rows = c()
		if(b > 0) b_rows     = 1+c(1:b)
		lambda_rows          = max(c(1,b_rows)) + 1:k
		if(r2 > 0) e_a2_rows = max(lambda_rows) + 1:r2
		mu = result_mat[1,]
		B_resid = matrix(result_mat[b_rows,],nr=b)
		Lambda = t(result_mat[lambda_rows,])
		E_a2 = result_mat[e_a2_rows,]

	
	 # -----Sample F, conditioning on Lambda, B, E_a2 ------------------ #
		QtY_resid_t = t(QtY - QtXZ2 %*% rbind(mu,B_resid,E_a2))
		QtY_prec_t = t(QtY_prec)
		prior_prec = prec_F_resid/(matrix(F_h2,nc = 1) %*% matrix(d,nr = 1) + (1-c(F_h2)))
		# prior_prec = 1/(matrix(F_h2,nc = 1) %*% matrix(d,nr = 1) + (1-c(F_h2)))
		prior_mean = t(QtX %*% B_F)
		FtQ = sample_coefs_c(QtY_resid_t, Lambda, QtY_prec_t,prior_mean, prior_prec)
		QtF = t(FtQ)
		F = Q %*% QtF
		
	 # -----Sample resid_h2 conditioning on B,Lambda,F,E_a2, resid_Y_prec -------------------- #
		QtY_resid = QtY - QtXZ2 %*% rbind(mu,B_resid,E_a2) - t(FtQ) %*% t(Lambda)
		resid_h2 = sample_h2s_discrete_c_new(QtY_resid,resid_Y_prec,h2_divisions,h2_priors_resids,d)
		E_a_prec = (1-resid_h2)/resid_h2 * resid_Y_prec

	 # -----Sample resid_Y_prec------------ #
		QtY_prec = t(1/(matrix(resid_h2,nc=1) %*% matrix(d,nr=1) + (1-c(resid_h2))))
		resid_Y_prec = rgamma(p,shape = resid_Y_prec_shape + n/2, rate = resid_Y_prec_rate+1/2*colSums(QtY_resid^2*QtY_prec)) #model residual precision
	
	 # -----Sample F components------------ #
		# note, have prec_B_F and prec_F_resid both of which are constrained to be >1 and > 0, but are not constrained to sum to 1

	 # -----Sample B_F ------------ #
		# recover()
		QtF_prec = t(prec_F_resid/(matrix(F_h2,nc=1) %*% matrix(d,nr=1) + (1-c(F_h2))))
		# QtF_prec = t(1/(matrix(F_h2,nc=1) %*% matrix(d,nr=1) + (1-c(F_h2))))
		# prior_prec = matrix(0,b,k,byrow=T)
		prior_prec = matrix(prec_B_F,b,k,byrow=T)
		prior_mean = matrix(0,b,k)
		B_F = sample_coefs_c(QtF, QtX, QtF_prec,prior_mean, prior_prec)

	 # -----Sample prec_B_F ------------ #
		# want a truncated gamma
		shape = prec_B_F_shape + (b)/2
		rate  = prec_B_F_rate + colSums(B_F^2)/2
		gam_unif = runif(k,pgamma(1,shape=shape,rate=rate),1)
		prec_B_F = qgamma(gam_unif,shape=shape,rate=rate)
		prec_B_F[gam_unif==1] = 1

	 # -----Sample prec_F_resid ------------ #
		# want a truncated gamma
		QtF_resid = QtF - QtX %*% B_F
		QtF_prec = t(1/(matrix(F_h2,nc=1) %*% matrix(d,nr=1) + (1-c(F_h2))))
		shape = prec_F_resid_shape + (n)/2
		rate  = prec_F_resid_rate+1/2*colSums(QtF_resid^2*QtF_prec)
		gam_unif = runif(k,pgamma(1,shape=shape,rate=rate),1)
		prec_F_resid = qgamma(gam_unif,shape=shape,rate=rate)
		prec_F_resid[gam_unif==1] = 1

		# prec_F_resid = rep(1,k)


	 # -----Sample F_h2 conditioning on F and B_F-------------------- #
		#conditioning on F, marginalizing over F_a
		# F_h2 = sample_h2s_discrete(F - X_f %*% B_F,h2_divisions,h2_priors_factors,invert_aI_bZAZ)
		# F_h2 = sample_h2s_discrete_c(F - X %*% B_F,h2_divisions,h2_priors_factors,invert_aI_bZAZ)
		F_h2 = sample_h2s_discrete_c_new(t(FtQ) - QtX %*% B_F,prec_F_resid,h2_divisions,h2_priors_factors,d)

	 # # -----Sample F_a ------------ #
		# QtF = t(FtQ)
		# QtF_prec = t(prec_F_resid/(matrix(F_h2,nc=1) %*% matrix(d,nr=1) + (1-c(F_h2))))
		# prior_prec = matrix(prec_B_F,b,k,byrow=T)
		# prior_mean = matrix(0,b,k)
		# B_F = sample_coefs_c(QtF, QtX, QtF_prec,prior_mean, prior_prec)
		
 	# -----Sample Lambda_prec------------- #
		Lambda2 = Lambda^2
		Lambda_prec = matrix(rgamma(p*k,shape = (Lambda_df + 1)/2,rate = (Lambda_df + sweep(Lambda2,2,tauh,'*'))/2),nr = p,nc = k)

 	# -----Sample B_prec------------- #
		shape = B_shape_shape + length(B_resid)/2
		rate = B_shape_rate + sum(B_prec * B_resid^2)/2
		B_shape = rgamma(1,shape=shape,rate=rate)
		B_prec = B_shape*matrix(rgamma(b*p,shape = (B_df + b)/2,rate = (B_df + colSums(B_resid^2*B_shape))/2),nr = b,nc = p)

		
	 # -----Sample E_a2_prec------------------ #
		if(ncol(Z_2) > 0) {
			E_a2_prec =  rgamma(p, E_a2_prec_shape + r2/2,rate = E_a2_prec_rate + 1/2*diag(t(E_a2) %*% A_2_inv %*% E_a2))
		}
	 
	 # -----Sample delta, update tauh------ #
		delta = sample_delta( delta,tauh,Lambda_prec,delta_1_shape,delta_1_rate,delta_2_shape,delta_2_rate,Lambda2 )
		tauh  = cumprod(delta)
		
	 # -----Update Plam-------------------- #
		Plam = sweep(Lambda_prec,2,tauh,'*')

	 # -- adapt number of factors to samples ---#
		adapt_result = update_k( 								
								Lambda_prec   = Lambda_prec,
								delta         = delta,
								tauh          = tauh,
								Plam          = Plam,
								Lambda        = Lambda,
								prec_B_F      = prec_B_F,
								prec_F_resid  = prec_F_resid,
								F_h2          = F_h2,
								B_F           = B_F,
								F_a           = F_a,
								F             = F,
								Z_1           = Z_1,
								Lambda_df     = Lambda_df,
								delta_2_shape = delta_2_shape,
								delta_2_rate  = delta_2_rate,
								b0            = b0,
								b1            = b1,
								i             = i,
								epsilon       = epsilon,
								prop          = prop 
								)
		Lambda_prec   = adapt_result$Lambda_prec
		delta         = adapt_result$delta
		tauh          = adapt_result$tauh
		Plam          = adapt_result$Plam
		Lambda        = adapt_result$Lambda
		prec_B_F      = adapt_result$prec_B_F
		prec_F_resid  = adapt_result$prec_F_resid
		F_h2          = adapt_result$F_h2
		B_F           = adapt_result$B_F
		F_a           = adapt_result$F_a
		F             = adapt_result$F
		k = ncol(F)

	 # -- save sampled values (after thinning) -- #
		if( (i-burn) %% thin == 0 && i > burn) {
				
			sp_num = (i-burn)/thin
			
			Posterior = save_posterior_samples( 
												sp_num       = sp_num,
												Posterior    = Posterior,
												Lambda       = Lambda,
												F            = F,
												B_F          = B_F,
												F_a          = F_a,
												mu           = mu,
												B_resid      = B_resid,
												E_a2         = E_a2,
												E_a          = E_a,
												delta        = delta,
												B_shape      = B_shape,
												F_h2         = F_h2,
												prec_B_F     = prec_B_F,
												prec_F_resid = prec_F_resid,
												resid_Y_prec = resid_Y_prec,
												resid_h2     = resid_h2,
												E_a_prec     = E_a_prec,
												E_a2_prec    = E_a2_prec
											)
			
			if(sp_num %% save_freq == 0) save(Posterior,file='Posterior.RData')
		}


	 # -- provide run diagnostics and plots -- #
	   if(i %% draw_iter  == 0) {
			# directory=strread(pwd,'#s','delimiter','/')
			# disp(directory(end))
			# disp(i)
			# elapsed=toc
			# #output some running statistics on the current factors and their genetic variances
			# # each row represents a factor:
			# #    delta[i], [i], F_h2[i], V_P(F_i), V_A(F_i)
			# [delta,[1:size(F,2)]',F_h2,sum(F.^2)'./(size(F,1)-1),sum(F_a.^2)'./(size(F_a,1)-1)]
			# disp(strcat('Time remaining:',num2str((nrun-i) * (elapsed/i) * 1/60)))

			#make some plots of some running statistics
			if(simulation) {
			   gen_factor_Lambda    = run_parameters$gen_factor_Lambda
			   error_factor_Lambda  = run_parameters$Lambda
			   G_act                = run_parameters$G
			   E_act                = run_parameters$R
			   h2                   = run_parameters$h2
			   # draw_simulation_diagnostics(sp_num,run_parameters,run_variables,Posterior,Lambda,F_h2,E_a_prec,resid_Y_prec)
			   draw_simulation_diagnostics_new(
			   									GxE_state = GxE_state,
			   									sp_num = sp_num,
												Posterior    = Posterior,
												Lambda       = Lambda,
												F            = F,
												B_F          = B_F,
												F_a          = F_a,
												mu           = mu,
												B_resid      = B_resid,
												E_a2         = E_a2,
												E_a          = E_a,
												delta        = delta,
												B_shape      = B_shape,
												F_h2         = F_h2,
												prec_B_F     = prec_B_F,
												prec_F_resid = prec_F_resid,
												resid_Y_prec = resid_Y_prec,
												resid_h2     = resid_h2,
												E_a_prec     = E_a_prec,
												E_a2_prec    = E_a2_prec
			   	)
			   # optional to also draw diagnostics for fixed effects. 10: fixed
			   # effects on factors 11: fixed effects on genes.
	#            figure(10)plot(X_f*current_state$B_F*current_state$Lambda',X_f*run_parameters$B_f*run_parameters$Lambda','.')line(xlim,xli
	#            figure(11)plot(X*current_state$B,X*run_parameters$B,'.')line(xlim,xli
			} else {
			   draw_results_diagnostics(sp_num,run_parameters,run_variables,Lambda, F_h2, Posterior)
			}
	   }
	}
	end_time = Sys.time()
	print(end_time - start_time)
	save(Posterior,file = 'Posterior.RData')


	# ----------------------------------------------- #
	# ------------Save state for restart------------- #
	# ----------------------------------------------- #

	
	current_state$resid_Y_prec = resid_Y_prec
	current_state$B_prec       = B_prec
	current_state$E_a2_prec    = E_a2_prec
	current_state$resid_h2     = resid_h2
	current_state$Lambda_prec  = Lambda_prec
	current_state$Plam         = Plam
	current_state$delta        = delta
	current_state$tauh         = tauh
	current_state$B_shape         = B_shape
	current_state$prec_B_F     = prec_B_F
	current_state$prec_F_resid = prec_F_resid
	current_state$F_h2         = F_h2
	current_state$Lambda       = Lambda
	current_state$F            = F
	current_state$mu           = mu
	current_state$B_resid      = B_resid
	current_state$E_a2         = E_a2
	current_state$B_F          = B_F
	current_state$F_a          = F_a
	current_state$E_a          = E_a
	current_state$nrun         = i

	save(current_state,file='current_state.RData')


	GxE_state$current_state = current_state
	GxE_state$Posterior = Posterior
	GxE_state$RNG = list(
		Random.seed = .Random.seed,
		RNGkind = RNGkind()
	)
	return(GxE_state)
}