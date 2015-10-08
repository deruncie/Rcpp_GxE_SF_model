# include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]
using namespace Rcpp;

arma::mat sweep_times(arma::mat x, int MARGIN, arma::vec STATS){
	int m = x.n_rows;
	int n = x.n_cols;
	arma::mat sweep_mat;
	if(MARGIN == 1) sweep_mat = repmat(STATS,1,n);
	if(MARGIN == 2) sweep_mat = repmat(STATS.t(),m,1);

	x = x % sweep_mat;

	return(x);
}

arma::mat dnorm_std_log(arma::mat x){
	return(log(1/sqrt(2*M_PI))-0.5* (x % x));
}

// [[Rcpp::export()]]
arma::mat sample_Lambda_c (arma::mat Y_tilde,
						arma::mat F,
						arma::vec resid_Y_prec,
						arma::vec E_a_prec,
						arma::mat Plam,
						List invert_aI_bZAZ ){
	// May want to switch out forwardsolve and backsolve for Armadillo functions


	// Sample factor loadings Lambda while marginalizing over residual
	// genetic effects: Y - Z_2W = F*Lambda' + E, vec(E)~N(0,kron(Psi_E,In) + kron(Psi_U, ZAZ^T))
	// note: conditioning on F, but marginalizing over E_a.
	// sampling is done separately by trait because each column of Lambda is
	// independent in the conditional posterior
	// note: invert_aI_bZAZ has parameters that diagonalize aI + bZAZ for fast
	// inversion: inv(aI + bZAZ) = 1/b*U*diag(1./(s+a/b))*U'
	Environment base("package:base");
	Function forwardsolve = base["forwardsolve"];
	Function backsolve = base["backsolve"];

	int p = resid_Y_prec.n_elem;
	int k = F.n_cols;

	arma::mat U = as<arma::mat>(invert_aI_bZAZ["U"]);
	arma::vec s = as<arma::vec>(invert_aI_bZAZ["s"]);

	arma::mat FtU = F.t() * U;
	arma::mat UtY = U.t() * Y_tilde;

	arma::mat Zlams = arma::randn(k,p);	
	// Environment stats("package:stats");
	// Function rnorm = stats["rnorm"];
	// arma::vec z = as<arma::vec>(rnorm(k*p));
	// arma::mat Zlams = reshape(z,k,p);
	// arma::mat Lambda = arma::zeros(p,k);

	arma::mat FUDi, Qlam, Llam_t,Plam_row;
	arma::vec means,vlam,mlam,ylam,Zcol;
	for(int j =0; j<p; j++){
		FUDi = E_a_prec(j) * sweep_times(FtU, 2,1/(s + E_a_prec(j)/resid_Y_prec(j)));
		means = FUDi * UtY.col(j);
		Qlam = FUDi * FtU.t() + diagmat(Plam.row(j));

		Llam_t = chol(Qlam);
		vlam = as<arma::vec>(forwardsolve(Llam_t.t(),means));
		mlam = as<arma::vec>(backsolve(Llam_t,vlam));
		Zcol = Zlams.col(j);
		ylam = as<arma::vec>(backsolve(Llam_t,Zcol));

		Lambda.row(j) = ylam.t() + mlam.t();

	}

	return(Lambda);
}

// [[Rcpp::export()]]
arma::mat sample_factors_scores_c(arma::mat Y_tilde,
								arma::mat X,
								arma::mat Z_1,
								arma::mat Lambda,
								arma::vec resid_Y_prec,
								arma::mat F_b,
								arma::mat F_a,
								arma::vec F_h2 ) {
//Sample factor scores given factor loadings (F_a), factor heritabilities (F_h2) and
//phenotype residuals

	arma::mat Lmsg = sweep_times(Lambda,1,resid_Y_prec);
	arma::vec tau_e = 1/(1-F_h2);
	arma::mat S = chol(Lambda.t() * Lmsg + diagmat(tau_e)).t();
	arma::mat Meta = trans(solve(S,trans(Y_tilde * Lmsg + sweep_times(X * F_b + Z_1 * F_a,2,tau_e))));

	arma::mat Zlams = arma::randn(Meta.n_rows,Meta.n_cols);	
	// Environment stats("package:stats");
	// Function rnorm = stats["rnorm"];
	// arma::vec z = as<arma::vec>(rnorm(Meta.n_rows*Meta.n_cols));
	// arma::mat Zlams = reshape(z,Meta.n_rows,Meta.n_cols);

	arma::mat F = trans(solve(trans(S),trans(Meta + Zlams)));

	return(F);
}

// [[Rcpp::export()]]
arma::mat sample_means_c(arma::mat Y_tilde,
					   arma::vec resid_Y_prec,
					   arma::vec E_a_prec,
					   List invert_aPXA_bDesignDesignT ) {
	// when used to sample [B;E_a]:
	//  W - F*Lambda' = X*B + Z_1*E_a + E, vec(E)~N(0,kron(Psi_E,In)). 
	//  Note: conditioning on F, Lambda and W.
	// The vector [b_j;E_{a_j}] is sampled simultaneously. Each trait is sampled separately because their
	// conditional posteriors factor into independent MVNs.
	// note:invert_aPXA_bDesignDesignT has parameters to diagonalize mixed model equations for fast inversion: 
	// inv(a*blkdiag(fixed_effects_prec*eye(b),Ainv) + b*[X Z_1]'[X Z_1]) = U*diag(1./(a.*s1+b.*s2))*U'
	// Design_U = [X Z_1]*U, which doesn't change each iteration. 
	
	arma::mat U = as<arma::mat>(invert_aPXA_bDesignDesignT["U"]);
	arma::vec s1 = as<arma::vec>(invert_aPXA_bDesignDesignT["s1"]);
	arma::vec s2 = as<arma::vec>(invert_aPXA_bDesignDesignT["s2"]);
	arma::mat Design_U = as<arma::mat>(invert_aPXA_bDesignDesignT["Design_U"]);

	// int n = Y_tilde.n_rows;
	int p = Y_tilde.n_cols;
	int br = Design_U.n_cols;

	arma::mat means = sweep_times(Design_U.t() * Y_tilde,2,resid_Y_prec);
	arma::mat location_sample = arma::zeros(br,p);


	arma::mat Zlams = arma::randn(br,p);	
	// Environment stats("package:stats");
	// Function rnorm = stats["rnorm"];
	// arma::vec z = as<arma::vec>(rnorm(br*p));
	// arma::mat Zlams = reshape(z,br,p);

	arma::vec d, mlam;
	for(int j =0; j<p; j++) {
		d = s1*E_a_prec(j) + s2*resid_Y_prec(j);
		mlam = means.col(j) /d;
		location_sample.col(j) = U * (mlam + Zlams.col(j)/sqrt(d));
	}

	return(location_sample);
}

// [[Rcpp::export()]]
arma::vec sample_prec_discrete_conditional_c(arma::mat Y,
										   int h2_divisions,
										   arma::vec h2_priors,
										   List invert_aI_bZAZ,
										   arma::vec res_prec) {
	//sample factor heritibilties conditional on a given residual precision
	//(res_precision)
	//prior given as relative weights on each of h2_divisions points. Doesn't
	//have to be normalized
	//samples conditional on F, marginalizes over F_a.
	//uses invert_aI_bZAZ.U and invert_aI_bZAZ.s to not have to invert aI + bZAZ
	//each iteration.
	//ident_prec is a in the above equation.

	// Environment stats("package:stats");
	// Function runif = stats["rnorm"];

	arma::mat U = as<arma::mat>(invert_aI_bZAZ["U"]);
	arma::vec s = as<arma::vec>(invert_aI_bZAZ["s"]);

	int p = Y.n_cols;
	arma::vec Trait_h2 = arma::zeros(p,1);

	arma::mat log_ps = arma::zeros(p,h2_divisions);
	arma::mat std_scores_b = Y.t() * U;
	arma::mat det_mat,std_scores;
	arma::mat det;
	for(double i =0; i < h2_divisions; i+=1){
		double h2 = (i)/(h2_divisions);
		if(h2 > 0) {
			std_scores = sweep_times(sweep_times(std_scores_b,2,1/sqrt(s+(1-h2)/h2)),1,1/sqrt(h2/(res_prec*(1-h2))));
			det_mat = log((s+(1-h2)/h2) * reshape(h2/(res_prec*(1-h2)),1,p))/2;
			det = sum(det_mat,0);
		} else {
			std_scores = Y.t();
			det = arma::zeros(1,p);
		}
		log_ps.col(i) = sum(dnorm_std_log(std_scores),1) - det.t() + log(h2_priors(i));
	}
	for(int j =0; j < p; j++){
		double norm_factor = max(log_ps.row(j))+log(sum(exp(log_ps.row(j)-max(log_ps.row(j)))));
		arma::mat ps_j = exp(log_ps.row(j) - norm_factor);
		log_ps.row(j) = ps_j;
		arma::vec r = arma::randu(1);
		arma::uvec selected = arma::find(repmat(r,1,h2_divisions)>cumsum(ps_j,1));
		Trait_h2(j) = double(selected.n_elem)/(h2_divisions);
	}
	arma::vec Prec = (res_prec % (1-Trait_h2))/Trait_h2;

	return(Prec);
}