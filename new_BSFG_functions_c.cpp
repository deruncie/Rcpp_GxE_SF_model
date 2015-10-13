# include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]


using namespace Rcpp;
using namespace arma;


mat sweep_times(mat x, int MARGIN, vec STATS){
	int m = x.n_rows;
	int n = x.n_cols;
	mat sweep_mat;
	if(MARGIN == 1) sweep_mat = repmat(STATS,1,n);
	if(MARGIN == 2) sweep_mat = repmat(STATS.t(),m,1);

	x = x % sweep_mat;

	return(x);
}

mat dnorm_std_log(mat x){
	return(log(1/sqrt(2*M_PI))-0.5* (x % x));
}


// [[Rcpp::export()]]
mat sample_coefs_c(
					mat QtY,
					mat QtW,
					mat QtY_prec,
					mat prior_mean,
					mat prior_prec) {
	
	// Sample regression coefficients
	// columns of matrices are independent
	// each column conditional posterior is a MVN due to conjugacy

	int b = QtW.n_cols;
	int p = QtY.n_cols;

	mat result = zeros(b,p);
	mat z = randn(b,p);	

	for(int j=0; j < p; j++){
		mat WtQDi = sweep_times(QtW.t(),2,QtY_prec.col(j));
		mat Phi = WtQDi * QtW + diagmat(prior_prec.col(j));
		vec m_rot = WtQDi * QtY.col(j) + diagmat(prior_prec.col(j)) * prior_mean.col(j);

		mat L_Phi = chol(Phi);
		vec mean = solve(L_Phi,solve(L_Phi.t(),m_rot));
		vec e = solve(L_Phi,z.col(j));

		result.col(j) = mean + e;
	}

	return(result);
}

// [[Rcpp::export()]]
mat sample_coefs_cis_effects_c(
					mat QtY,
					mat QtcisGenotypes,  // for genes with no genotype info, code with 0 (not 1)
					mat QtW,
					mat QtY_prec,
					mat prior_mean,
					mat prior_prec) {
	
	// Sample regression coefficients
	// columns of matrices are independent
	// each column conditional posterior is a MVN due to conjugacy

	int b = QtW.n_cols + 1; // +1 for cis effect
	int p = QtY.n_cols;

	mat result = zeros(b,p);
	mat z = randn(b,p);	

	for(int j=0; j < p; j++){
		mat QtcisW = join_rows(QtcisGenotypes.col(j),QtW);
		mat WtQDi = sweep_times(QtcisW.t(),2,QtY_prec.col(j));
		mat Phi = WtQDi * QtcisW + diagmat(prior_prec.col(j));
		vec m_rot = WtQDi * QtY.col(j) + diagmat(prior_prec.col(j)) * prior_mean.col(j);

		mat L_Phi = chol(Phi);
		vec mean = solve(L_Phi,solve(L_Phi.t(),m_rot));
		vec e = solve(L_Phi,z.col(j));

		result.col(j) = mean + e;
	}

	return(result);
}


// [[Rcpp::export()]]
vec sample_h2s_discrete_c_new (mat QtY,
						vec resid_prec,
						int h2_divisions,
						vec h2_priors,
						vec d){
	// sample factor heritibilties from a discrete set on [0,1)
	// samples conditional on F, marginalizes over F_a.
	// uses invert_aI_bZAZ.U and invert_aI_bZAZ.s to not have to invert aI + bZAZ
	// each iteration.

	int p = QtY.n_cols;
	int n = QtY.n_rows;
	vec resid_h2 = zeros(p);

	mat YtQ_std = sweep_times(QtY,2,sqrt(resid_prec)).t();

	mat log_ps = zeros(p,h2_divisions);

	for(double i =0; i < h2_divisions; i+=1){
		double h2 = (i)/(h2_divisions);
		vec D = h2*d + (1-h2);
		mat std_scores = sweep_times(YtQ_std,2,1/sqrt(D));
		vec SS = sum(std_scores % std_scores,1);
		double log_det = sum(log(D));
		log_ps.col(i) = -n/2*log(2*M_PI) -0.5 * log_det - 0.5 * SS + log(h2_priors(i));
	}
	for(int j =0; j < p; j++){
		double norm_factor = max(log_ps.row(j))+log(sum(exp(log_ps.row(j)-max(log_ps.row(j)))));
		mat ps_j = exp(log_ps.row(j) - norm_factor);
		log_ps.row(j) = ps_j;
		vec r = randu(1);
		uvec selected = find(repmat(r,1,h2_divisions)>cumsum(ps_j,1));
		resid_h2(j) = double(selected.n_elem)/(h2_divisions);
	}

	return(resid_h2);
}



// [[Rcpp::export()]]
mat sample_a (	mat Y,
				mat QtLiZt,
				vec resid_Y_prec,
				vec resid_h2,
				vec d,
				mat LitQ

			) {
    // generate sample from the posterior of a given y, s2_a, s2_e: y ~ N(Za,s2e), a~N(0,s2_a*A)
    // uses: QZ1 * diag(d) * QZ1' = 1/s2e*L'Z'ZL + 1/s2a*I
    // with: LL' = A
	int r = QtLiZt.n_rows;
	int p = Y.n_cols;

	mat a = zeros(r,p);
	mat a_tilde = randn(r,p);

	for(int j=0; j < p; j++){
		double s2a = resid_h2(j)/resid_Y_prec(j);
		double s2e = (1-resid_h2(j))/resid_Y_prec(j);
		vec D = 1/s2a + d/s2e;
		vec D_inv = 1/D;
		vec mu = sweep_times(QtLiZt,1,D_inv) * Y.col(j)/s2e;
		a.col(j) = LitQ * (mu + sqrt(D_inv) % a_tilde.col(j));
	}
	return(a);  
}


// // [[Rcpp::export()]]
// vec sample_h2s_discrete_c_new (mat QtY,
// 						vec resid_prec,
// 						int h2_divisions,
// 						vec h2_priors,
// 						vec d){
// 	// sample factor heritibilties from a discrete set on [0,1)
// 	// prior places 50% of the weight at h2=0
// 	// samples conditional on F, marginalizes over F_a.
// 	// uses invert_aI_bZAZ.U and invert_aI_bZAZ.s to not have to invert aI + bZAZ
// 	// each iteration.

// 	int p = QtY.n_cols;
// 	vec resid_h2 = zeros(p);

// 	mat YtQ_std = sweep_times(QtY,2,sqrt(resid_prec)).t();

// 	mat log_ps = zeros(p,h2_divisions);

// 	mat det_mat,std_scores;
// 	mat det;
// 	for(double i =0; i < h2_divisions; i+=1){
// 		double h2 = (i)/(h2_divisions);
// 		if(h2 > 0) {
// 			std_scores = 1/sqrt(h2) * sweep_times(YtQ_std,2,1/sqrt(d+(1-h2)/h2));
// 			det = sum(log((d+(1-h2)/h2)*h2)/2) * ones(1,p);
// 		} else {
// 			std_scores = QtY.t();
// 			det = zeros(1,p);
// 		}
// 		log_ps.col(i) = sum(dnorm_std_log(std_scores),1) - det.t() + log(h2_priors(i));
// 	}
// 	for(int j =0; j < p; j++){
// 		double norm_factor = max(log_ps.row(j))+log(sum(exp(log_ps.row(j)-max(log_ps.row(j)))));
// 		mat ps_j = exp(log_ps.row(j) - norm_factor);
// 		log_ps.row(j) = ps_j;
// 		vec r = randu(1);
// 		uvec selected = find(repmat(r,1,h2_divisions)>cumsum(ps_j,1));
// 		resid_h2(j) = double(selected.n_elem)/(h2_divisions);
// 	}

// 	return(resid_h2);
// }


// [[Rcpp::export()]]
vec sample_h2s_discrete_c (mat F,
						int h2_divisions,
						vec h2_priors,
						List invert_aI_bZAZ){
	// sample factor heritibilties from a discrete set on [0,1)
	// prior places 50% of the weight at h2=0
	// samples conditional on F, marginalizes over F_a.
	// uses invert_aI_bZAZ.U and invert_aI_bZAZ.s to not have to invert aI + bZAZ
	// each iteration.

	mat U = as<mat>(invert_aI_bZAZ["U"]);
	vec s = as<vec>(invert_aI_bZAZ["s"]);

	int k = F.n_cols;
	vec F_h2 = zeros(k);

	mat log_ps = zeros(k,h2_divisions);
	mat std_scores_b = F.t() * U;

	mat det_mat,std_scores;
	mat det;
	for(double i =0; i < h2_divisions; i+=1){
		double h2 = (i)/(h2_divisions);
		if(h2 > 0) {
			std_scores = 1/sqrt(h2) * sweep_times(std_scores_b,2,1/sqrt(s+(1-h2)/h2));
			det = sum(sum(log((s+(1-h2)/h2)*h2)/2)) * ones(1,k);
		} else {
			std_scores = F.t();
			det = zeros(1,k);
		}
		log_ps.col(i) = sum(dnorm_std_log(std_scores),1) - det.t() + log(h2_priors(i));
	}
	for(int j =0; j < k; j++){
		double norm_factor = max(log_ps.row(j))+log(sum(exp(log_ps.row(j)-max(log_ps.row(j)))));
		mat ps_j = exp(log_ps.row(j) - norm_factor);
		log_ps.row(j) = ps_j;
		vec r = randu(1);
		uvec selected = find(repmat(r,1,h2_divisions)>cumsum(ps_j,1));
		F_h2(j) = double(selected.n_elem)/(h2_divisions);
	}

	return(F_h2);
}

