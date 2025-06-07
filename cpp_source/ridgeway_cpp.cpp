#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;



arma::vec update_logw(double lb, double ub, int index, arma::vec logw){
  
  double logpli = R::pnorm(lb, 0.0, 1.0, 1, 1);
  double logpui = R::pnorm(ub, 0.0, 1.0, 1, 1);

  if(index == -1){
    if( lb > -arma::datum::inf && ub < arma::datum::inf) {
      logw.fill( logpui + log(1.0 - exp(logpli - logpui)) );
    } else if(ub < arma::datum::inf) {
      logw.fill( logpui ); 
    } else {
      logw.fill(log(1.0 - exp(logpli))); 
      Rcout << "...";
    }
    
  }else{  

  if( lb > -arma::datum::inf && ub < arma::datum::inf) {
    logw[index] += logpui + log(1.0 - exp(logpli - logpui));
  } else if(ub < arma::datum::inf) {
    logw[index] += logpui; 
  } else {
    Rcout << logpli << "\n" << lb << "\n";
    logw[index] += log(1.0 - exp(logpli)); 
    Rcout << "...";
  }
  }
return(logw);
}

double logSumExp_cpp(const arma::vec logv) {
  double maxv = logv.max();
  return maxv + log(arma::accu(exp(logv - maxv)));
}

// [[Rcpp::export]]
double rtruncnorm_cpp(double a, 
                      double b,
                      double mu = 0.0, 
                      double sigma = 1.0) {
  
  if(a>b){
    Rcout << "Error: lower limit is "<<a<< " and upper limit is "<<b<<"\n"; 
  }
  
  double sample = 0.0;
  // Standardize the truncation points
  //double beta  = (b - mu) / sigma;
  // Calculate standard normal CDF at truncation points
  double phi_beta  = R::pnorm(b,  mu, sigma, true, false);
  
  if( a == -arma::datum::inf){
    // Generate uniform random numbers in [0, phi_beta]
    double u = R::runif(0,phi_beta);
    
    if( u == 0){
      
      Rcout << "Warning: Underflow occurred";
      sample = b - 1e-3;
      
    }else{
    
    // Apply inverse CDF transform
    sample = R::qnorm(u, mu, sigma, true, false);
    }
    // Transform back to desired location and scale
    
  }else{
    
    // Standardize the truncation points
    //double alpha = (a - mu) / sigma;
    // Calculate standard normal CDF at truncation points
    double phi_alpha = R::pnorm(a, mu, sigma, true, false);
    
    // Generate uniform random numbers in [phi_alpha, phi_beta]
    double u = R::runif(phi_alpha,phi_beta);
    //u = phi_alpha + u * (phi_beta - phi_alpha);
    
    // Apply inverse CDF transform
    sample = R::qnorm(u, mu, sigma, true, false);
    
    // Transform back to desired location and scale
  }
  
  return sample;
}


arma::rowvec rtruncnorm_cpp_row(int N,
                                double a, 
                                double b,
                                double mu = 0.0, 
                                double sigma = 1.0) {
  
  arma::rowvec sample(N);
  for(int n = 0; n < N; n ++){
    sample[n] = rtruncnorm_cpp(a, b, mu, sigma);
  }
  
  return sample;
}

arma::rowvec rtruncnorm_cpp_row2( arma::rowvec a, 
                                  arma::rowvec b,
                                  double mu = 0.0, 
                                  double sigma = 1.0) {
  
  int N = a.n_cols;
  arma::rowvec sample(N);
  
  for(int qqq = 0; qqq < N; qqq ++){
    sample[qqq] = rtruncnorm_cpp( a[qqq], b[qqq], mu, sigma);
  }
  
  return sample;
}

arma::vec pnorm_cpp(arma::vec x){
  arma::vec p = x; //p.zeros();
  for(int i=0; i< x.n_elem; i++){
    p[i] = R::pnorm5(x[i], 0.0, 1.0, true, false);
  }
  return(p);
}

arma::mat ESS_low(int i, int M, arma::colvec logw, double logZ, 
                  arma::mat eta_particle_mat, 
                  arma::mat C, arma::colvec b){ 
  
  arma::uvec index_0_M = arma::linspace<arma::uvec>(0, M-1, M);
  
  // resampling step
  // a1: update normalizing constant
  // a2: resample particles and set all weights equal to 1 (log = 0)      
  arma::vec wprob = exp(logw - logSumExp_cpp(logw));
  arma::uvec resampledIndices = 
    Rcpp::RcppArmadillo::sample_main(index_0_M, M, true,  wprob);
  
  arma::mat sm0 = eta_particle_mat.rows(0,i);
  eta_particle_mat.rows(0,i) = sm0.cols(resampledIndices);
  
  // logw.zeros();
  
  // a3: Gibbs move from full conditionals
  for (int j = 0; j < (i+1); j++) {
    
    int n_rows = i+1 - j; // rows from j to i-1 in C++
    arma::mat C_sub = C.submat(j, 0, i, i); // C[j:(i), 0:(i)]
    C_sub.shed_col(j);
    
    arma::mat eta_sub = eta_particle_mat.rows(0, i);
    eta_sub.shed_row(j); // remove j-th row
    
    arma::mat mat_b = arma::repmat(b.subvec(j, i), 1, M);
    
    arma::colvec C_jij = C.submat(j, j, i, j);  // C[j:(i-1), j]
    arma::mat denom = arma::repmat(arma::abs(C_jij), 1, M); // replicate abs col
    
    arma::mat limits = (mat_b - C_sub * eta_sub) / denom;
    
    arma::uvec pos_idx = arma::find(C_jij > 0);
    arma::uvec neg_idx = arma::find(C_jij < 0);
    
    arma::rowvec ub3(M);
    arma::rowvec lb3(M);
    
    
    if (!pos_idx.is_empty()) {
      arma::mat limits_pos = limits.rows(pos_idx);
      ub3 = arma::min(limits_pos, 0);  // column-wise min
    } else {
      ub3.fill(arma::datum::inf);
    }
    
    if (!neg_idx.is_empty()) {
      arma::mat limits_neg = -limits.rows(neg_idx);  // negate values
      lb3 = arma::max(limits_neg, 0); // column-wise max
    } else {
      lb3.fill(-arma::datum::inf);
    }
    
    // Sample from truncated normal
    eta_particle_mat.row(j) = rtruncnorm_cpp_row2(lb3, ub3, 0.0, 1.0 ); 
  }
  return(eta_particle_mat);
  //return List::create(Named("eta_particle_mat") = eta_particle_mat, 
  //                     Named("logw") = logw);
}


// [[Rcpp::export]]
List reorder_Sigma_cpp(arma::vec x, arma::mat Sigma) {
  
  int m = x.n_rows;
  
  arma::vec a(m); 
  a.fill(-arma::datum::inf);
  arma::vec b = x;
  
  arma::uvec sorted_indices(m);
  arma::vec eta(m); eta.zeros();
  arma::mat C(m,m); C.zeros();
  
  // Initial step
  // Compute standardized bounds
  arma::vec denom = arma::sqrt( Sigma.diag() );
  
  arma::vec lb = a / denom;
  arma::vec ub = b / denom;
  
  // Compute the range of truncated normal probabilities
  arma::vec probs = pnorm_cpp(ub) - pnorm_cpp(lb);
  
  // Find index with smallest probability range
  sorted_indices[0] = probs.index_min();
  
  // Swap rows and columns in Sigma
  arma::rowvec Sigma_1_old = Sigma.row(0);
  arma::rowvec Sigma_1_new = Sigma.row(sorted_indices[0]);
  
  double var_1_old = Sigma_1_old[0];
  double var_1_new = Sigma_1_new[sorted_indices[0]];
  double cov_old_new = Sigma_1_old[sorted_indices[0]];
  
  Sigma_1_old(0) = cov_old_new;
  Sigma_1_old(sorted_indices[0]) = var_1_old;
  
  Sigma_1_new(0) = var_1_new;
  Sigma_1_new(sorted_indices[0]) = cov_old_new;
  
  // Apply swaps
  Sigma.row(0) = Sigma_1_new;
  Sigma.col(0) = Sigma_1_new.t();
  
  Sigma.row(sorted_indices[0]) = Sigma_1_old;
  Sigma.col(sorted_indices[0]) = Sigma_1_old.t();
  
  // Swap a and b
  double a_old = a(0); a(0) = a(sorted_indices[0]); a(sorted_indices[0]) = a_old;
  
  double b_old = b(0); b(0) = b(sorted_indices[0]); b(sorted_indices[0]) = b_old;
  
  C(0,0)  = sqrt(Sigma(0,0));
  C(arma::span(1,m-1),0) = Sigma(arma::span(1,m-1),0) / C(0,0);
  
  // compute mean eta_1 (mu_1 in Trinh and Genz)
  double hat_a = a[0]/C(0,0);
  double hat_b = b[0]/C(0,0);
  eta(0) = (R::dnorm4(hat_a, 0.0, 1.0, false) - R::dnorm4(hat_b, 0.0, 1.0, false)) /
    (R::pnorm(hat_b, 0.0, 1.0, true, false) -  R::pnorm5(hat_a, 0.0, 1.0, true, false));
  
  // ---------------------------------------------------------------------------
  
  for(int i=1; i<m; i++) {
    
    arma::vec sigma = Sigma.diag();
    
    // find optimal index and eta[i]
    arma::vec num0 = C.cols( 0, i-1 ) * eta.rows( 0, i-1);
    arma::vec num = num0.rows(i,m-1);
    
    arma::vec den0 = arma::sum(C.cols( 0, i-1 )%C.cols( 0, i-1 ), 1);
    arma::vec den  = sqrt(sigma.rows(i,m-1) - den0.rows(i,m-1));  
    
    lb = ((a.rows(i,m-1) - num)/den);
    ub = ((b.rows(i,m-1) - num)/den);
    
    
    // Compute the range of truncated normal probabilities
    arma::vec probs2 =  pnorm_cpp(ub) - pnorm_cpp(lb);
    
    // Find index with smallest probability range
    
    // switch indeces of variables in Sigma
    int opt_i   = (i) + probs2.index_min(); //# sum (i-1) to get an index on the original 1,...,m scale
    
    arma::rowvec Sigma_i_old = Sigma.row(i);
    arma::rowvec Sigma_i_new = Sigma.row(opt_i);
    
    
    double var_i_old      = Sigma_i_old[i];
    double var_i_new      = Sigma_i_new[opt_i];
    double cov_old_new    = Sigma_i_old[opt_i];
    Sigma_i_old[i] = cov_old_new;
    Sigma_i_old[opt_i] = var_i_old;
    Sigma_i_new[i] = var_i_new;
    Sigma_i_new[opt_i] = cov_old_new;
    
    Sigma.row(i) = Sigma_i_new;
    Sigma.col(i) = Sigma_i_new.t();
    Sigma.row(opt_i) = Sigma_i_old;
    Sigma.col(opt_i) = Sigma_i_old.t();
    
    // switch indeces in a and b
    a_old     = a[i]; a[i]  = a[opt_i]; a[opt_i]  = a_old;
    b_old     = b[i]; b[i]  = b[opt_i]; b[opt_i]  = b_old;
    
    // switch rows in C
    arma::rowvec c_old    = C.row(i); C.row(i)    = C.row(opt_i);  C.row(opt_i) = c_old;
    
    // complete the i-th column of C
    C(i,i) = sqrt( Sigma(i,i) - arma::accu( pow( C(i,arma::span(0,i-1) ),2) ) );
    if(i < m-1){
      C(arma::span(i+1,m-1),i) = (Sigma(arma::span(i+1,m-1),i) - 
        C(arma::span(i+1,m-1), arma::span(0,i-1)) * C(i, arma::span(0,i-1)).t() )/C(i,i);
    }
    
    arma::mat Q = (C(i,arma::span(0,i-1)) * eta.rows(0,i-1) ) /C(i,i);
    
    // update eta_i
    hat_a = (a[i]/C(i,i) - Q[0]);
    hat_b = (b[i]/C(i,i) - Q[0]);
    eta(i) = (R::dnorm4(hat_a, 0.0, 1.0, false) - R::dnorm4(hat_b, 0.0, 1.0, false)) /
      (R::pnorm5(hat_b, 0.0, 1.0, true, false) -  R::pnorm5(hat_a, 0.0, 1.0, true, false));
    
  }
  
  return List::create(Named("a") = a, 
                      Named("b") = b,
                      Named("C") = C);
}



// [[Rcpp::export]]
List ridgway_smc_cpp(arma::vec x,
                     arma::mat Sigma, 
                     double ESS_fraction = 0.5,
                     int M = 1e4, bool verbose = false) {
  
  int m = x.n_rows;
  
  List L = reorder_Sigma_cpp(x, Sigma);
  
  arma::vec a = L["a"];
  arma::vec b = L["b"];
  arma::mat C = L["C"];
  
  arma::mat eta_particle_mat(m, M); 
  eta_particle_mat.zeros();
  arma::vec logw(M);
  logw.zeros();
  double logZ = 0.0;
  
  // Step 1
  double lb = a[0] / C(0,0);
  double ub = b[0] / C(0,0);
  eta_particle_mat.row(0) = rtruncnorm_cpp_row(M, lb, ub, 0.0, 1.0); 

  logw =  update_logw(lb, ub, -1, logw);
  
  for (int i = 1; i < m; i++) {
    R_CheckUserInterrupt();
    
    arma::rowvec C_i = C(i, arma::span(0, i-1));
    
    for (int part = 0; part < M; part++) {
      
      arma::vec eta_prev = eta_particle_mat( arma::span(0, i-1), part );
      double Ci_eta = arma::dot(C_i, eta_prev);
      double lbi = (a[i] - Ci_eta) / C(i,i);
      double ubi = (b[i] - Ci_eta) / C(i,i);
      eta_particle_mat(i, part) = rtruncnorm_cpp(lbi, ubi, 0.0, 1.0);
      
      logw =  update_logw(lbi, ubi, part, logw);
    
    }
    double ESS = exp(2 * logSumExp_cpp(logw) - logSumExp_cpp(2*logw));
    if(verbose){Rcout << "ESS: "<<ESS << "-- Obs:" << i << """\n";}
    if (ESS < M * ESS_fraction) {
      logZ += (-log(M) + logSumExp_cpp(logw));
      
      eta_particle_mat = ESS_low(i,
                                 M, 
                                 logw, 
                                 logZ, 
                                 eta_particle_mat, 
                                 C, 
                                 b);
      logw.zeros();
    }
  }
  double logCDF = logZ - log(M) + logSumExp_cpp(logw);
  return List::create(Named("eta") = eta_particle_mat, 
                      Named("logCDF") = logCDF,
                      Named("logw") = logw);
}
