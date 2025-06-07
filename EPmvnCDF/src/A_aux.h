#ifndef A_aux
#define A_aux
#include <RcppArmadillo.h>

int mod(int a, int n);

arma::vec zeta1_v(arma::vec x);

arma::vec zeta2_v(arma::vec x , arma::vec z1);

double zeta1(double x);

double zeta2(double x , double z1);

arma::rowvec crossprod_vectomat(arma::vec x, arma::mat y);

arma::mat tcrossprod(arma::mat x);

double quadform(arma::vec x, arma::mat A);

double check_convergence_diffs(double rel_logdelta_k,
                               double rel_logdelta_m,
                               double rel_logdelta_Z,
                               double diff);

#endif
