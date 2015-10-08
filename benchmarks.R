benchmark(
sample_Lambda_c( Y_tilde,F,resid_Y_prec, E_a_prec,Plam,invert_aI_bZAZ ),
sample_prec_discrete_conditional_c(Y - X %*% B-F %*% t(Lambda) - Z_2 %*% W,h2_divisions,h2_priors_resids,invert_aI_bZAZ,resid_Y_prec),
sample_means_c( Y_tilde, resid_Y_prec, E_a_prec, invert_aPXA_bDesignDesignT ),
sample_h2s_discrete(F - X_f %*% F_b,h2_divisions,h2_priors_factors,invert_aI_bZAZ),
sample_F_a_c(F,X_f,cbind(X_f,Z_1),F_h2,invert_aPXfA_bDesignDesignT),
sample_factors_scores_c( Y_tilde, X_f,Z_1,Lambda,resid_Y_prec,F_b,F_a,F_h2 ),
matrix(rgamma(p*k,shape = (Lambda_df + 1)/2,rate = (Lambda_df + sweep(Lambda2,2,tauh,'*'))/2),nr = p,nc = k),
rgamma(p,shape = resid_Y_prec_shape + n/2, rate = resid_Y_prec_rate+1/2*colSums(Y_tilde^2)),
sample_delta( delta,tauh,Lambda_prec,delta_1_shape,delta_1_rate,delta_2_shape,delta_2_rate,Lambda2 ),
sweep(Lambda_prec,2,tauh,'*'),
update_k( F,Lambda,F_b,F_a,F_h2,Lambda_prec,Plam,delta,tauh,Z_1,Lambda_df,delta_2_shape,delta_2_rate,b0,b1,i,epsilon,prop ),order='relative')

benchmark(
sample_Lambda( Y_tilde,F,resid_Y_prec, E_a_prec,Plam,invert_aI_bZAZ ),
sample_prec_discrete_conditional(Y - X %*% B-F %*% t(Lambda) - Z_2 %*% W,h2_divisions,h2_priors_resids,invert_aI_bZAZ,resid_Y_prec),
sample_means( Y_tilde, resid_Y_prec, E_a_prec, invert_aPXA_bDesignDesignT ),
sample_h2s_discrete(F - X_f %*% F_b,h2_divisions,h2_priors_factors,invert_aI_bZAZ),
sample_F_a(F,X_f,cbind(X_f,Z_1),F_h2,invert_aPXfA_bDesignDesignT),
sample_factors_scores( Y_tilde, X_f,Z_1,Lambda,resid_Y_prec,F_b,F_a,F_h2 ),
matrix(rgamma(p*k,shape = (Lambda_df + 1)/2,rate = (Lambda_df + sweep(Lambda2,2,tauh,'*'))/2),nr = p,nc = k),
rgamma(p,shape = resid_Y_prec_shape + n/2, rate = resid_Y_prec_rate+1/2*colSums(Y_tilde^2)),
sample_delta( delta,tauh,Lambda_prec,delta_1_shape,delta_1_rate,delta_2_shape,delta_2_rate,Lambda2 ),
sweep(Lambda_prec,2,tauh,'*'),
update_k( F,Lambda,F_b,F_a,F_h2,Lambda_prec,Plam,delta,tauh,Z_1,Lambda_df,delta_2_shape,delta_2_rate,b0,b1,i,epsilon,prop ),order='relative')

benchmark({
sample_Lambda_c( Y_tilde,F,resid_Y_prec, E_a_prec,Plam,invert_aI_bZAZ );
sample_prec_discrete_conditional_c(Y - X %*% B-F %*% t(Lambda) - Z_2 %*% W,h2_divisions,h2_priors_resids,invert_aI_bZAZ,resid_Y_prec);
sample_means_c( Y_tilde, resid_Y_prec, E_a_prec, invert_aPXA_bDesignDesignT );
sample_h2s_discrete(F - X_f %*% F_b,h2_divisions,h2_priors_factors,invert_aI_bZAZ);
sample_F_a_c(F,X_f,cbind(X_f,Z_1),F_h2,invert_aPXfA_bDesignDesignT);
sample_factors_scores_c( Y_tilde, X_f,Z_1,Lambda,resid_Y_prec,F_b,F_a,F_h2 );
matrix(rgamma(p*k,shape = (Lambda_df + 1)/2,rate = (Lambda_df + sweep(Lambda2,2,tauh,'*'))/2),nr = p,nc = k);
rgamma(p,shape = resid_Y_prec_shape + n/2, rate = resid_Y_prec_rate+1/2*colSums(Y_tilde^2));
sample_delta( delta,tauh,Lambda_prec,delta_1_shape,delta_1_rate,delta_2_shape,delta_2_rate,Lambda2 );
sweep(Lambda_prec,2,tauh,'*');
update_k( F,Lambda,F_b,F_a,F_h2,Lambda_prec,Plam,delta,tauh,Z_1,Lambda_df,delta_2_shape,delta_2_rate,b0,b1,i,epsilon,prop );
},{
sample_Lambda( Y_tilde,F,resid_Y_prec, E_a_prec,Plam,invert_aI_bZAZ );
sample_prec_discrete_conditional(Y - X %*% B-F %*% t(Lambda) - Z_2 %*% W,h2_divisions,h2_priors_resids,invert_aI_bZAZ,resid_Y_prec);
sample_means( Y_tilde, resid_Y_prec, E_a_prec, invert_aPXA_bDesignDesignT );
sample_h2s_discrete(F - X_f %*% F_b,h2_divisions,h2_priors_factors,invert_aI_bZAZ);
sample_F_a(F,X_f,cbind(X_f,Z_1),F_h2,invert_aPXfA_bDesignDesignT);
sample_factors_scores( Y_tilde, X_f,Z_1,Lambda,resid_Y_prec,F_b,F_a,F_h2 );
matrix(rgamma(p*k,shape = (Lambda_df + 1)/2,rate = (Lambda_df + sweep(Lambda2,2,tauh,'*'))/2),nr = p,nc = k);
rgamma(p,shape = resid_Y_prec_shape + n/2, rate = resid_Y_prec_rate+1/2*colSums(Y_tilde^2));
sample_delta( delta,tauh,Lambda_prec,delta_1_shape,delta_1_rate,delta_2_shape,delta_2_rate,Lambda2 );
sweep(Lambda_prec,2,tauh,'*');
update_k( F,Lambda,F_b,F_a,F_h2,Lambda_prec,Plam,delta,tauh,Z_1,Lambda_df,delta_2_shape,delta_2_rate,b0,b1,i,epsilon,prop );
},
order='relative')
