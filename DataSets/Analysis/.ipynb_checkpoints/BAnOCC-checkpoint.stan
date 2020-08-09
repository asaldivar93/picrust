data {
  int<lower=0> P; // Number of elements in one composition
  int<lower=0> N; // Number of samples
  matrix<lower=0, upper=1>[N,P] C; // Composition matrix
  vector[P] n; // the mean prior parameter for m (one mean prior per P element in sampleś composition)
  cov_matrix[P] L; // the scale prior parameter for m
  real<lower=0> a; // the shape prior parameter for lambda
  real<lower=0> b; // the rate prior parameter for lambda
}

transformed data{
  matrix[P, P] L_chol; //cholesky decomposition of covariance matrix
  matrix[P, N] log_C_t; //log transformed composition matrix
  matrix[N, P] one_mat; //matrix of ones
  
  L_chol = cholesky_decompose(L);
  log_C_t = log(C)';
  one_mat = rep_matrix(1, N, P);
}

parameters {
  vector[P] m_raw; // initial value for m
  cov_matrix[P] O; // log-basis precision matrix
  real<lower=0> lamb; // shrinkage prior lambda
}

transformed parameters{
  vector[P] m; // lognormal centrality parameter
  
  m = n + L_chol*m_raw;
}

model {
  matrix[P,N] alpha_star;
  real s_sq_star;
  vector[N] m_star;
  vector[P] O_diag; // diagonal values of O
  vector[4] lik; //vector for the four parts of the likelihooh
  
  m_raw ~ normal(0, 1); // implies: m ~ multi_normal(n, L)
  lamb ~ gamma(a, b);
  
  O_diag = diagonal(O); // Slice of the diagonal elements of O
  O_diag ~ exponential(lamb/2); // ¿? Updates the diagonal elements with a draw form the exp distribution
  
  for (k in 1:(P-1)){
    vector[P- k] O_tri_row; // row k of the uuper-triangular elements of O
    O_tri_row = sub_col(O, k +1, k, P - k);
    O_tri_row ~ double_exponential(0, lamb); // ¿? Updates the row of the upper-triangular elements of O with a draw from doub-exp distribution
  }
  
  lik[1] = 0.5*N*log_determinant(O);
  
  s_sq_star = inv(sum(O));
  lik[2] = 0.5*N*log(s_sq_star);
  
  alpha_star = rep_matrix(m, N) - log_C_t;
  lik[3] = -0.5*trace_quad_form(O, alpha_star);
  
  m_star = rows_dot_product(alpha_star'*O, one_mat);
  lik[4] = 0.5*dot_self(m_star)*s_sq_star;
  
  target += sum(lik);
}
