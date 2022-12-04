
data {
 int<lower=1> K;
 int<lower=1> D;
 int<lower=1> N;
 int<lower=1> NN;
 int<lower=1> M; //samples for marginal effects est
 int<lower=0,upper=1> y[N, D];
 vector[K] x[N];
 matrix[N,K] x_data;
 matrix[NN,K] x_predict;
 matrix[M,K] x_m_1; //marginal effect AREA
 matrix[M,K] x_m_2; //marginal effect CROP
 matrix[M,K] x_m_3; //marginal effect SEP
 }
parameters {
 vector[D] beta0;
 matrix[D,K] beta; 
 cholesky_factor_corr[D] L_Omega;
 
 corr_matrix[K] Sigma;        // prior correlation
 vector<lower=0>[K] tau;
 
 // matrix[K, K] Sigma_beta;
 //model 3 vector[K] 
 real<lower=0, upper=1> u[N,D]; // nuisance that absorbs inequality constraints
 
real mu_beta0;
real <lower=0> sigma_beta0;

vector[K] mu_betak;
 }
model {
 L_Omega ~ lkj_corr_cholesky(2);
 mu_beta0 ~ normal(0,1); // weakly informative prior
sigma_beta0 ~ gamma(1,1); // weakly informative prior

for (i in 1:D) {
beta0[i] ~ normal(mu_beta0, sigma_beta0);
}
 // beta0 ~ normal(0,5);
 // to_vector(beta) ~ student_t(1, 0, 1.25);
 
 // has to be [D, K+1], but Sigma is K+1 x K+1\
 // so it generates a 1 x K sample, need to make it 
 // D x K, where each d=1,...,D is independent
 tau ~ cauchy(0, 1);
 Sigma ~ lkj_corr(2);
mu_betak ~ normal(0,5);
 // for (i in 1:D)
 //   beta[i] ~ multi_normal(rep_vector(0, K), quad_form_diag(Sigma, tau));
 for (i in 1:D)
   beta[i] ~ multi_normal(mu_betak, quad_form_diag(Sigma, tau));
 
 { // likelihood
  for (n in 1:N){
   vector[D] mu;
   vector[D] z;
   real prev;
   // vector[D] beta0;
   // TODO assign beta0 to D-dim array of last column 
   // beta0 = beta[:,K+1]; // should be D x 1
   mu = beta0 + (beta * x[n]); // make beta a D x K matrix 
   prev = 0;
   for (d in 1:D) { // Phi and inv_Phi may overflow and / or be numerically inaccurate
    real bound; // threshold at which utility = 0
    bound = Phi( -(mu[d] + prev) / L_Omega[d, d] );
    if (y[n, d] == 1) {
     real t;
     t = bound + (1 - bound) * u[n, d];
     z[d] = inv_Phi(t); // implies utility is positive
     target += log1m(bound); // Jacobian adjustment
     }
    else {
     real t;
     t = bound * u[n, d];
     z[d] = inv_Phi(t); // implies utility is negative
     target += log(bound); // Jacobian adjustment
     }
    if (d < D) prev = L_Omega[d+1, 1:d] * head(z, d);
    // Jacobian adjustments imply z is truncated standard normal
    // thus utility --- mu + L_Omega * z --- is truncated multivariate normal
    }
   }
  }
 }    
generated quantities {
 matrix[D,N] z_new;
 matrix[D,NN] z_predict;
 matrix[D,M] z_m_1;
 matrix[D,M] z_m_2;
 matrix[D,M] z_m_3;
 corr_matrix[D] Omega;

 Omega = multiply_lower_tri_self_transpose(L_Omega);

 {
  matrix[K,D] t_beta;
  matrix[D,N] beta_x;
  matrix[D,NN] beta_x_predict;
  matrix[D,M] beta_x_m1;
  matrix[D,M] beta_x_m2;
  matrix[D,M] beta_x_m3;
   t_beta = beta';
   for (d in 1:D){
    beta_x[d,] = beta0[d] + (x_data * t_beta[,d])';
    beta_x_predict[d,] = beta0[d] + (x_predict * t_beta[,d])';
    beta_x_m1[d,] = beta0[d] + (x_m_1 * t_beta[,d])';
    beta_x_m2[d,] = beta0[d] + (x_m_2 * t_beta[,d])';
    beta_x_m3[d,] = beta0[d] + (x_m_3 * t_beta[,d])';
    }
   for (n in 1:N)
    z_new[,n] = multi_normal_cholesky_rng(beta_x[,n], L_Omega);
   for (nn in 1:NN){
    z_predict[,nn] = multi_normal_cholesky_rng(beta_x_predict[,nn], L_Omega);
    }
   for (m in 1:M) {
    z_m_1[,m] = multi_normal_cholesky_rng(beta_x_m1[,m], L_Omega);
    z_m_2[,m] = multi_normal_cholesky_rng(beta_x_m2[,m], L_Omega);
    z_m_3[,m] = multi_normal_cholesky_rng(beta_x_m3[,m], L_Omega);
    }
  }
 }
