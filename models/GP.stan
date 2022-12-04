functions {
	//GP covariance function
	vector gp(vector[] x, real sdgp, real lscale, vector zgp) { 
		matrix[size(x), size(x)] cov;
		cov = cov_exp_quad(x, sdgp, lscale);
		for (n in 1:size(x)) {
			cov[n, n] = cov[n, n] + 1e-12;
		}
		return cholesky_decompose(cov) * zgp;
	}
}
data {
 int<lower=1> K;  //features dimension
 int<lower=1> D;
 int<lower=1> N1; //samples for training
 int<lower=1> N2; //samples for test
 int<lower=1> NN; //samples for prediction
 int idx1[N1]; //indices for training observations
 int idx2[N2]; //indices for test observations
 int<lower=1> M; //samples for marginal effects est
 int<lower=0,upper=1> y[N1+N2, D];
 int<lower=0,upper=1> Y[N1+N2];
 int trials[N1+N2];
 vector[K] x[N1+N2];
 vector[K] X[NN];
 matrix[NN,K] x_predict;
 matrix[M,K] x_m_1; //marginal effect AREA
 matrix[M,K] x_m_2; //marginal effect CROP
 matrix[M,K] x_m_3; //marginal effect SEP
 }
parameters {
 real<lower=0> alpha[K]; // lengthscale of f
 real<lower=0> rho[K];       // scale of f
 real<lower=0> sigma;         // noise sigma
 matrix[N1+N2,K] eta;
 }
transformed parameters{
	vector[N1+N2] f;
	
	f = gp(x, alpha[1], rho[1], eta[,1]);

}
model {
  // priors
 alpha ~ normal(0, 1);
 rho ~ normal(0, 1);
 sigma ~ normal(0,5);
 to_vector(eta) ~ normal(0,1);
 Y[idx1] ~ binomial_logit(trials[idx1], f[idx1]);
 }    
generated quantities{
	vector[N1] y_predict;
	vector[N1] f_invlogit;
	vector[N1] log_lik;
	// vector[NN] log_y_predict;
	// vector[N] y_new;
	vector[N2] log_y_new;
	
	for(i in 1:(N1)){
		y_predict[i] = binomial_rng(1, inv_logit(f[i]));
		f_invlogit[i] = inv_logit(f[i]);  //bernoulli_logit_rng(f[i])
		log_lik[i] = bernoulli_logit_lpmf(Y[i] | inv_logit(f[i]));
	}
	
	for(i in 1:N2){
		log_y_new[i] = binomial_logit_lpmf(Y[idx2[i]] | trials[idx2[i]], f[idx2[i]]);
	}
	
}
