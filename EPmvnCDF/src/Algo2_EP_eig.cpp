#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
#include "A_aux.h"

// [[Rcpp::export]]
Rcpp::List getParamsEP_priorMean_priorVar_cpp(arma::mat X,
                                              arma::vec y,
                                              arma::vec priorMean,
                                              arma::vec priorVariances, 
                                              double tolerance,
                                              int maxIter=1e6){
  
  int n = X.n_rows;
  int p = X.n_cols;
  arma::mat Xt = X.t(); 
  double diff = 1;

  double log_tolerance = log(tolerance);

  int count = 0;
  
  arma::vec r_prior = priorMean / priorVariances;
  arma::vec r = r_prior;
  
  double logDetInvPrior = -arma::accu(log(priorVariances));
  double logDetinvQ = logDetInvPrior;
  
  // Containers
  arma::colvec logZ(n);
  arma::colvec k(n);
  arma::mat repl_k(n,p);
  arma::colvec m(n);
  
  arma::mat OmXt = arma::diagmat(priorVariances)*Xt; 
  
  arma::mat V = OmXt;
  arma::vec diagOmega(p);
  arma::vec meanBeta(p);
  int nnIter = 0;
  
  //Rcpp::Rcout << 1;
  // Algorithm
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
        
        double rel_logdelta_k = log(std::fabs(delta_k)) - log(std::fabs(k[i]));
        double rel_logdelta_m = log(std::fabs(delta_m)) - log(std::fabs(m[i]));
        
        k[i] = kNew;
        m[i] = mNew;
        
        double prev_logZ = logZ[i];
        logZ[i] = (( 2 * m[i] * r_iTw + 
          pow(m[i],2) * xTw - k[i]*pow(r_iTw,2))/
            ( 1 + k[i] * xTw ) - log(1.0 + k[i] * xTw ) )*.5 - log(arma::normcdf(tau));
        double delta_Z = exp(logZ[i]) - exp(prev_logZ);
        double rel_logdelta_Z = log(std::fabs(delta_Z)) - logZ[i];
        
        r = r + delta_m * xi;
        
        double denominator = (1.+(delta_k)*xTv);
        double ratio = -(delta_k)/denominator ;
        logDetinvQ += log(denominator);
        V += (v*ratio)*crossprod_vectomat(xi,V);
        
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
  
  /// Posterior Approximate Moments
  arma::mat XOm = OmXt.t();
  //Rcpp::Rcout << 2;
  diagOmega = priorVariances - arma::sum( V % (OmXt * arma::diagmat(k) ), 1);
  //Rcpp::Rcout << arma::sum( V * (OmXt * arma::diagmat(k)), 1);
  //Rcpp::Rcout << 3;
  arma::vec Om_r = priorVariances % r;
  //Rcpp::Rcout << 4;
  meanBeta = Om_r - V * (arma::diagmat(k) * X * Om_r);
  
  double logML =  (arma::dot(r,meanBeta) - logDetinvQ + logDetInvPrior - 
                   arma::dot(r_prior,priorMean)) *.5 - arma::accu(logZ);
  
  List results = List::create(_["meanBeta"] = meanBeta,
                              _["diagOmega"] = diagOmega,
                              _["logML"] = logML,
                              _["nIter"] = nnIter,
                              _["kEP"] = k,
                              _["mEP"] = m);
  
  return(results);
  
}
