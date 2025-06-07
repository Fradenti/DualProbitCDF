#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
#include "A_aux.h"

// [[Rcpp::export]]
Rcpp::List getParamsEP_priorMean_algo1_cpp(
                                arma::mat X,
                                arma::vec y,
                                arma::vec priorMean,
                                double nu2,
                                double tolerance,
                                int maxIter=1e6){
 
  int n = X.n_rows;
  int p = X.n_cols;
  
  double log_tolerance = log(tolerance);
  
  arma::mat Xt = X.t(); 
  double diff = 1;
  
  int count =0;
  arma::vec r = priorMean/nu2;
  double logDetinvQ = -p*log(nu2);
  
  // Containers
  arma::colvec logZ(n);
  arma::colvec k(n);
  arma::mat repl_k(n,p);
  arma::colvec m(n);
  
  arma::vec meanBeta(p);
  
  
  arma::mat invQ = arma::eye(p,p) * nu2 ; 
  
  int nnIter = 0;
  
  for(int nIter = 0; nIter<maxIter; nIter++){
    Rcpp::checkUserInterrupt();
    if(diff < log_tolerance){nnIter = nIter; 
      break; }
    
    diff = - arma::datum::inf;
    count = 0;
    
    for(int i=0; i<n; i++){
      
      arma::vec xi = (Xt.col(i));
      arma::vec r_i = r - m[i]*xi;
      
      arma::mat Oxi = invQ * xi;
      double xitOxi = arma::dot(xi,Oxi);
      /////////////////////////////////////////
      arma::mat Oi = invQ + (Oxi*Oxi.t())*k[i] / (1.0 - k[i]*xitOxi);
      /////////////////////////////////////////
      arma::mat Oixi = Oi*xi;
      double xiOixi = arma::dot(xi, Oixi);
        
      
      if(xiOixi>0){
        
        double r_iOixi = arma::dot(r_i,Oixi);
        
        double s = (2.0 * y[i] - 1.0) / sqrt(1.0 + xiOixi);
        double tau = s * r_iOixi;
        
        double z1 = zeta1(tau);
        double z2 = zeta2(tau,z1);
        
        double kNew = -z2/(1.0 + xiOixi + z2 * xiOixi);
        double delta_k = kNew - k[i];
        double mNew = (z1*s + kNew * r_iOixi + kNew * z1 * s * xiOixi);
        double delta_m = mNew - m[i];
        
        double rel_logdelta_k = log(std::fabs(delta_k)) - log(std::fabs(k[i]));
        double rel_logdelta_m = log(std::fabs(delta_m)) - log(std::fabs(m[i]));
        
        
        k[i] = kNew;
        m[i] = mNew;
        
        double prev_logZ = logZ[i];
        logZ[i] = (( 2 * m[i] * r_iOixi + 
          pow(m[i],2) * xiOixi - k[i]*pow(r_iOixi,2))/
               ( 1 + k[i] * xiOixi ) - 
              log(1.0 + k[i] * xiOixi ) )*.5 - 
              log(arma::normcdf(tau));
        double delta_Z = exp(logZ[i]) - exp(prev_logZ);
        double rel_logdelta_Z = log(std::fabs(delta_Z)) - logZ[i];
        
        
        r = r_i + m[i]*xi;
        
        double denominator = (1.+(delta_k) * quadform(xi,invQ) );
        logDetinvQ += log(denominator);
        invQ = Oi + z2 * pow(s,2)* (Oixi*Oixi.t());
        
        diff = check_convergence_diffs(rel_logdelta_k,
                                       rel_logdelta_m,
                                       rel_logdelta_Z,
                                       diff);
        
      }else{
        count = count+1;
        Rcpp::Rcout << count << " units skipped\n";;
      }
     
     
      
    }
    

  }

  meanBeta  = invQ * r;


  double logML = (arma::dot(r,meanBeta) - logDetinvQ - p*log(nu2) - 
                  arma::dot(priorMean,priorMean)/nu2) *.5 - arma::accu(logZ);
  
  arma::colvec id = invQ.diag();
  List results = List::create(_["meanBeta"] = meanBeta,
                              _["diagOmega"] = id,
                              _["logML"] = logML,
                              _["nIter"] = nnIter,
                              _["kEP"] = k,
                              _["mEP"] = m);
  return(results);
  
}