#include "A_aux.h"

int mod(int a, int n){
  return a - floor(a/n)*n;
}

arma::vec zeta1_v(arma::vec x){
  arma::vec y = arma::log_normpdf(x) - log(arma::normcdf(x));
  return( exp(y) );
}

arma::vec zeta2_v(arma::vec x , arma::vec z1){
  arma::vec y = - (z1 + x)%z1;;
  return( y );
}

double zeta1(double x){
  double y = arma::log_normpdf(x) - log(arma::normcdf(x));
  return( exp(y) );
}

double zeta2(double x , double z1){
  double y = - (z1 + x)*z1;;
  return( y );
}

// return rowvec
arma::rowvec crossprod_vectomat(arma::vec x, arma::mat y){
  arma::mat C = x.t()*y;
  return(C);
}

arma::mat tcrossprod(arma::mat x){
  return(x * x.t());
}

double quadform(arma::vec x, arma::mat A){
  arma::mat B = x.t() * A * x;
  return( B(0,0) );
}


double check_convergence_diffs(double rel_logdelta_k,
                         double rel_logdelta_m,
                         double rel_logdelta_Z,
                         double diffs){
  
  double maxDiff = 0;
  
  arma::vec D(3);
  D(0) = rel_logdelta_k;
  D(1) = rel_logdelta_m;
  D(2) = rel_logdelta_Z;
  maxDiff = D.max();
  
  if( maxDiff > diffs ){diffs = maxDiff;}
  
  return(diffs);
  
  
}

