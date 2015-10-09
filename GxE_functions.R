library(R.matlab)

require(pracma)


cholcov = function(X){
	# calculates a matrix U such that t(U) %*% U == X for X that is not PD
	E = svd(X)
	cols = E$d > 1e-14
	U = E$u[,cols] %*% diag(sqrt(E$d[cols]))
	return(t(U))
}

sample_delta = function( delta,tauh,Lambda_prec,delta_1_shape,delta_1_rate,delta_2_shape,delta_2_rate,Lambda2 ) {
	#sample delta and tauh parameters that control the magnitudes of higher
	#index factor loadings.

	# matlab = readMat('sample_delta_data.mat')
	# for(i in 1:10) names(matlab) = sub('.','_',names(matlab),fixed=T)
	# delta = matlab$delta
	# tauh = matlab$tauh
	# Lambda_prec = matlab$Lambda_prec
	# delta_1_shape = matlab$delta_1_shape
	# delta_1_rate = matlab$delta_1_rate
	# delta_2_shape = matlab$delta_2_shape
	# delta_2_rate = matlab$delta_2_rate
	# Lambda2	 = matlab$Lambda2	

	k = length(tauh)
	mat = Lambda_prec * Lambda2
	n_genes = nrow(mat)

	shape = delta_1_shape + 0.5*n_genes*k
	rate = delta_1_rate + 0.5*(1/delta[1])*sum(tauh*colSums(mat))
	delta[1] = rgamma(1,shape = shape,rate = rate)
	tauh = cumprod(delta)

	for(h in 2:(k-1)) {
		shape = delta_2_shape + 0.5*n_genes*(k-h+1)
		if(h<k){
			rate = delta_2_rate + 0.5*(1/delta[h])*sum(tauh[h:k]*colSums(mat[,h:k]))
		} else{
			rate = delta_2_rate + 0.5*(1/delta[h])*sum(tauh[h:k]*sum(mat[,h:k]))    	
		}
		delta[h] = rgamma(1,shape = shape, rate = rate)
		tauh = cumprod(delta)
	}
	
	return(delta)
}


update_k = function( 
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
				) {
# adapt the number of factors by dropping factors with only small loadings
# if they exist, or adding new factors sampled from the prior if all factors
# appear important. The rate of adaptation decreases through the chain,
# controlled by b0 and b1. Should work correctly over continuations of
# previously stopped chains.

	n = nrow(F)
	k = ncol(F)
	b = nrow(B_F)
	r = nrow(F_a)
	p = nrow(Lambda)
	gene_rows = 1:p

	prob = 1/exp(b0 + b1*i)                # probability of adapting
	uu = runif(1)
	lind = colMeans(abs(Lambda) < epsilon)    # proportion of elements in each column less than eps in magnitude
	vec = lind >= prop
	num = sum(vec)       # number of redundant columns


	if(uu < prob && i>200){
		if(i > 20 && num == 0 && all(lind < 0.995) && k < 2*p) { #add a column
			k           =k+1
			Lambda_prec = cbind(Lambda_prec,rgamma(p,shape = Lambda_df/2, rate = Lambda_df/2))
			delta[k]    = rgamma(1,shape = delta_2_shape,rate = delta_2_rate)
			tauh        = cumprod(delta)
			Plam        = sweep(Lambda_prec,2,tauh,'*')
			Lambda      = cbind(Lambda,rnorm(p,0,sqrt(1/Plam[,k])))
			prec_B_F[k] = runif(1,1,3)
			prec_F_resid[k] = runif(1,1,3)
			F_h2[k]     = runif(1)
			B_F         = cbind(B_F,rnorm(b,0,1/sqrt(prec_B_F[k])))
			F_a         = cbind(F_a,rnorm(r,0,sqrt(F_h2[k]/prec_F_resid[k])))
			F           = cbind(F,rnorm(n,Z_1 %*% F_a[,k],sqrt((1-F_h2[k])/prec_F_resid[k])))
		} else if(num > 0) { # drop redundant columns
			nonred      = which(vec == 0) # non-redundant loadings columns
			Lambda_prec = Lambda_prec[,nonred]
			for(red in which(vec == 1)){
				if(red == k) next
				# combine deltas so that the shrinkage of kept columns doesnt
				# decrease after dropping redundant columns
				delta[red+1] = delta[red+1]*delta[red]
			}
			delta        = delta[nonred]
			tauh         = cumprod(delta)
			Plam         = sweep(Lambda_prec,2,tauh,'*')
			Lambda       = Lambda[,nonred]
			prec_B_F     = prec_B_F[nonred]
			prec_F_resid = prec_F_resid[nonred]
			F_h2         = F_h2[nonred]
			B_F          = B_F[,nonred]
			F_a          = F_a[,nonred]
			F            = F[,nonred]
		}
	}


	return(list(
				Lambda_prec  = Lambda_prec,
				delta        = delta,
				tauh         = tauh,
				Plam         = Plam,
				Lambda       = Lambda,
				prec_B_F     = prec_B_F,
				prec_F_resid = prec_F_resid,
				F_h2         = F_h2,
				B_F          = B_F,
				F_a          = F_a,
				F            = F
		))
}


save_posterior_samples = function( 
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
									prec_B_F         = prec_B_F,
									prec_F_resid         = prec_F_resid,
									resid_Y_prec = resid_Y_prec,
									resid_h2     = resid_h2,
									E_a_prec     = E_a_prec,
									E_a2_prec    = E_a2_prec
								) {
	# save posteriors. Full traces are kept of the more interesting parameters.
	# Only the posterior means are kept of less interesting parameters. These
	# should be correctly calculated over several re-starts of the sampler.

	sp = ncol(Posterior$Lambda)
	size(Posterior$Lambda,2)

	#save factor samples
	if(length(Lambda) > nrow(Posterior$Lambda)){
		# expand factor sample matrix if necessary
		Posterior$Lambda       = rbind(Posterior$Lambda, 	matrix(0,nr = length(Lambda)-nrow(Posterior$Lambda),nc = sp))
		Posterior$F            = rbind(Posterior$F, 	   	matrix(0,nr = length(F)     -nrow(Posterior$F),		nc = sp))
		Posterior$B_F          = rbind(Posterior$B_F,    	matrix(0,nr = length(B_F) 	-nrow(Posterior$B_F),	nc = sp))
		Posterior$F_a          = rbind(Posterior$F_a, 	matrix(0,nr = length(F_a) 	-nrow(Posterior$F_a),	nc = sp))
		Posterior$delta        = rbind(Posterior$delta, 	matrix(0,nr = length(delta) -nrow(Posterior$delta),	nc = sp))
		Posterior$F_h2         = rbind(Posterior$F_h2, 	matrix(0,nr = length(F_h2) 	-nrow(Posterior$F_h2),	nc = sp))
		Posterior$prec_B_F     = rbind(Posterior$prec_B_F, 	matrix(0,nr = length(prec_B_F) 	-nrow(Posterior$prec_B_F),	nc = sp))
		Posterior$prec_F_resid = rbind(Posterior$prec_F_resid, 	matrix(0,nr = length(prec_F_resid) 	-nrow(Posterior$prec_F_resid),	nc = sp))
	}
	Posterior$Lambda[1:length(Lambda),sp_num] = c(Lambda)
	Posterior$F[1:length(F),sp_num]     = c(F)
	Posterior$B_F[1:length(B_F),sp_num] = c(B_F)
	Posterior$F_a[1:length(F_a),sp_num] = c(F_a)
	Posterior$delta[1:length(delta),sp_num] = delta
	Posterior$F_h2[1:length(F_h2),sp_num] = F_h2
	Posterior$prec_B_F[1:length(prec_B_F),sp_num] = prec_B_F
	Posterior$prec_F_resid[1:length(prec_F_resid),sp_num] = prec_F_resid

	Posterior$mu[,sp_num] = mu
	Posterior$B_shape[sp_num] = B_shape
	Posterior$B_resid[,sp_num] = c(B_resid)
	Posterior$resid_Y_prec[,sp_num] = resid_Y_prec
	Posterior$resid_h2[,sp_num]     = resid_h2
	Posterior$E_a_prec[,sp_num]       = E_a_prec
	try({Posterior$E_a2_prec[,sp_num]       = E_a2_prec},silent=T)

	# save B,U,W
	Posterior$E_a = (Posterior$E_a*(sp_num-1) + E_a)/sp_num
	try({Posterior$E_a2   = (Posterior$E_a2*(sp_num-1) + E_a2)/sp_num},silent=T)

	return(Posterior)
}


clear_Posterior = function(GxE_state) {
	# resets Posterior samples if burnin was not sufficient
	Posterior = GxE_state$Posterior
	run_parameters = GxE_state$run_parameters

	if(!is.null(ncol(Posterior$Lambda))) {    
    	run_parameters$burn = run_parameters$burn + run_parameters$thin*ncol(Posterior$Lambda)
    }

    p = nrow(Posterior$resid_Y_prec)
    b = nrow(Posterior$B_resid)/p
    n = nrow(Posterior$W)
    r = nrow(Posterior$E_a)
    r2 = nrow(Posterior$E_a2)
    
    Posterior = list(
            Lambda       = matrix(0,nr=0,nc=0),
            F            = matrix(0,nr=0,nc=0),
            B_F          = matrix(0,nr=0,nc=0),
            F_a          = matrix(0,nr=0,nc=0),
            mu           = matrix(0,nr=p,nc=0),
            B_resid      = matrix(0,nr=b*p,nc=0),
            delta        = matrix(0,nr=0,nc=0),
            B_shape      = matrix(0,nr=1,nc=0),
            F_h2         = matrix(0,nr=0,nc=0),
            prec_B_F     = matrix(0,nr=0,nc=0),
            prec_F_resid = matrix(0,nr=0,nc=0),
            resid_Y_prec = matrix(0,nr=p,nc=0),
            resid_h2     = matrix(0,nr=p,nc=0),
            E_a_prec     = matrix(0,nr=p,nc=0),
            E_a2_prec    = matrix(0,nr=p,nc=0),
            E_a          = matrix(0,nr=r,nc=p),
            E_a2         = matrix(0,nr=r2,nc=p)
    	)

    GxE_state$Posterior = Posterior
    GxE_state$run_parameters = run_parameters
    return(GxE_state)

}