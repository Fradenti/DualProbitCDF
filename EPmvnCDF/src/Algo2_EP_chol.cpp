#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
#include "A_aux.h"

// [[Rcpp::export]]
Rcpp::List getParamsEP_priorMean_cpp(
                          arma::mat X,
                          arma::vec y,
                          arma::vec priorMean,
                          double nu2,
                          double tolerance,
                          int maxIter=1e6){
  
  int n = X.n_rows;
  int p = X.n_cols;
  
  arma::mat Xt = X.t(); 
  double diff = 1;
  
  double log_tolerance = log(tolerance);
  arma::vec Log_rel_differences(3);
  
  int count =0;
  arma::vec r = priorMean/nu2;
  double logDetinvQ = -p*log(nu2);
  
  // Containers
  arma::colvec logZ(n);
  arma::colvec k(n);
  arma::mat repl_k(n,p);
  arma::colvec m(n);
  

  arma::mat invQ = arma::eye(p,p) * nu2 ; 
  arma::mat V = nu2*Xt;
  arma::vec diagOmega(p);
  arma::vec meanBeta(p);
  int nnIter = 0;
  
  for(int nIter = 0; nIter<maxIter; nIter++){
    Rcpp::checkUserInterrupt();
    if(diff < log_tolerance){nnIter = nIter; 
      break; }

      diff = - arma::datum::inf;
      count = 0;
      
      for(int i=0; i<n; i++){
        
        arma::colvec v = V.col(i);
        arma::colvec xi = (X.row(i)).t();
        double xTv = arma::dot(xi,v);
        
        double d = 1.0 - k[i]*xTv;
        arma::colvec w = v/d;
        double xTw = xTv/d;

        if(xTw>0){
          
          double r_iTw = arma::dot(r,w) - m[i]*xTw;
          
          double s = (2 * y[i] - 1.0) / sqrt(1.0 + xTw);
          double tau = s*r_iTw;

          double z1 = zeta1(tau);
          double z2 = zeta2(tau,z1);
          
          double kNew = -z2/(1.0 + xTw + z2*xTw);
          double delta_k = kNew - k[i];
          double mNew = (z1*s + kNew*r_iTw + kNew*z1*s*xTw);
          double delta_m = mNew - m[i];
          
          Log_rel_differences(0) = log(std::fabs(delta_k)) - log(std::fabs(k[i]));
          Log_rel_differences(1) = log(std::fabs(delta_m)) - log(std::fabs(m[i]));
          
          
          k[i] = kNew;
          m[i] = mNew;
          
          double prev_logZ = logZ[i];
          logZ[i] = (( 2 * m[i] * r_iTw + 
                        pow(m[i],2) * xTw - k[i]*pow(r_iTw,2))/
                     ( 1 + k[i] * xTw ) - log(1.0 + k[i] * xTw ) )*.5 - log(arma::normcdf(tau));
          double delta_Z = exp(logZ[i]) - exp(prev_logZ);
          Log_rel_differences(2) = log(std::fabs(delta_Z)) - logZ[i];
          
          r = r + delta_m * xi;
          
          double denominator = (1.+(delta_k)*xTv);
          double ratio = -(delta_k)/denominator ;
          logDetinvQ += log(denominator);
          V += (v*ratio)*crossprod_vectomat(xi,V);
          
          diff = check_convergence_diffs(Log_rel_differences(0),
                                         Log_rel_differences(1),
                                         Log_rel_differences(2),
                                         diff);

        }else{
          count = count+1;
          Rcpp::Rcout << count << " units skipped\n";;
        }
        
      }
  
  arma::mat repl_k = arma::repelem(k,1,p);
  meanBeta  = nu2 * (r - V*( repl_k % X * r ));
  diagOmega = nu2 * ( 1.0 - (arma::sum(V %  (repl_k % X).t() ,1)));
   
  }
  
  
  double logML = (arma::dot(r,meanBeta) - logDetinvQ - p*log(nu2) - 
                  arma::dot(priorMean,priorMean)/nu2) *.5 - arma::accu(logZ);
  
  List results = List::create(_["meanBeta"] = meanBeta,
                                _["diagOmega"] = diagOmega,
                                _["logML"] = logML,
                                _["nIter"] = nnIter,
                                _["kEP"] = k,
                                _["mEP"] = m);
  return(results);
    
}


//////--------------------------------------------------------------------------
//////--------------------------------------------------------------------------
//////--------------------------------------------------------------------------
//////--------------------------------------------------------------------------
//////--------------------------------------------------------------------------

