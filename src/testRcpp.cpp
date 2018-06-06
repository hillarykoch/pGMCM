#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace RcppArmadillo;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

// estimation and model selection of penalized GMM
// [[Rcpp::export]]
arma::mat C_fpGMM(List sigma) {
    arma::mat tmp_sigma = sigma[1];
    return tmp_sigma;
}

// [[Rcpp::export]]
arma::vec C_oneL(arma::vec prob0) {
    arma::vec h_est = prob0/(sum(prob0));
    return h_est;
}

// [[Rcpp::export]]
List C_slice(int n, int d, int k, arma::mat x, arma::mat mu_new, arma::mat h_est, List sigma_old){
    // update sigma
    arma::cube dist(n, d, k, arma::fill::none);
    List sigma_new = sigma_old;
    arma::mat ones(n, 1, arma::fill::ones);
    for(int i = 0; i < k; ++i) {
         dist.slice(i) = x - ones * mu_new.row(i);
         sigma_new[i] = trans(dist.slice(i)) * diagmat(h_est.col(i)) * dist.slice(i) / (sum(h_est.col(i))  * 1.0L);
     }
    return sigma_new;
}

// [[Rcpp::export]]
List choose_list(List sigma, arma::uvec& idx) {
    int n_keep = idx.size();
    List out(n_keep);
    for(int i = 0; i < n_keep; ++i) {
        out[i] = sigma[idx(i)];
    }
    return out;
}

// [[Rcpp::export]]
List getuvec(arma::vec input, List sigma){
    arma::uvec idx = find(input > 0);
    sigma = choose_list(sigma, idx);
    return sigma;
}
