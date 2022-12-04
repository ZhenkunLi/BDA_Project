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
matrix[N1+N2,K]  x; //xxxxxx
vector[K] X[NN];
matrix[NN,K] x_predict;
matrix[M,K] x_m_1; //marginal effect AREA
matrix[M,K] x_m_2; //marginal effect CROP
matrix[M,K] x_m_3; //marginal effect SEP
}
parameters {
real beta0;
vector[K] betak; //beta k
real mu_beta0;
real <lower=0> sigma_beta0;
real mu_betak;
real <lower=0> sigma_betak;
}
transformed parameters{
vector[N1+N2] f;
f = betak[1] * x[,1] + betak[2] * x[,2] + betak[3] * x[,3] + beta0;
}
model {
// priors
mu_beta0 ~ normal(0, 1);               // mu_beta0 ~ normal(0, 10); //in Section 10
sigma_beta0 ~ gamma(2, 2);             // sigma_beta0 ~ gamma(1, 0.5); //in Section 10
beta0 ~ normal(mu_beta0, sigma_beta0); 
mu_betak ~ normal(0, 3);               // mu_betak ~ normal(0, 30); //in Section 10
sigma_betak ~ gamma(3, 3);             // sigma_betak ~ gamma(3, 0.5); //in Section 10

for (i in 1:K)
betak[i] ~ normal(mu_betak, sigma_betak);
//likelihood
Y[idx1] ~ binomial_logit(trials[idx1], f[idx1]);
}    
generated quantities{
	vector[N1+N2] y_predict;
	vector[N1+N2] f_invlogit;
	vector[N1+N2] log_lik;
	vector[N2] log_y_new;

for(i in 1:(N1+N2)){
y_predict[i] = binomial_rng(1, inv_logit(f[i]));
f_invlogit[i] = inv_logit(f[i]);
log_lik[i] = bernoulli_logit_lpmf(Y[i] | inv_logit(f[i]));
}

for(i in 1:N2){
log_y_new[i] = binomial_logit_lpmf(Y[idx2[i]] | trials[idx2[i]], f[idx2[i]]);
}
}


