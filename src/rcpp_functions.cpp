// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace RcppArmadillo;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

// choose slices from a cube given index
// [[Rcpp::export]]
arma::cube choose_slice(arma::cube& Q, arma::uvec idx, int d, int k) {
    arma::cube y(d, d, k, arma::fill::none);
    for(int i = 0; i < k; ++i) {
        y.slice(i) = Q.slice(idx(i));
    }
    return y;
}

// compute Mahalanobis distance for multivariate normal
// [[Rcpp::export]]
arma::vec Mahalanobis(arma::mat x, arma::rowvec mu, arma::mat sigma){
    const int n = x.n_rows;
    arma::mat x_cen;
    x_cen.copy_size(x);
    for (int i=0; i < n; i++) {
        x_cen.row(i) = x.row(i) - mu;
    }
    return sum((x_cen*sigma.i()) % x_cen, 1);
}

// Compute density of multivariate normal
// [[Rcpp::export]]
arma::vec cdmvnorm(arma::mat x, arma::rowvec mu, arma::mat sigma){
    arma::vec mdist = Mahalanobis(x, mu, sigma);
    double logdet = log(arma::det(sigma));
    const double log2pi = std::log(2.0 * M_PI);
    arma::vec logretval = -(x.n_cols * log2pi + logdet + mdist)/2;

    return exp(logretval);
}

// estimation and model selection of penalized GMM
// [[Rcpp::export]]
Rcpp::List cfpGMM(arma::mat& x,
                arma::rowvec prop,
                arma::mat& mu,
                arma::cube& sigma,
                int k,
                double df,
                int citermax,
                double lambda) {

    const int n = x.n_rows;
    const int d = x.n_cols;
    double delta = 1;
    double tol = 1E-06;
    arma::rowvec prop_old = prop;
    arma::mat mu_old = mu;
    arma::cube sigma_old = sigma;
    arma::uvec tag(n, arma::fill::none);
    arma::mat pdf_est(n, k, arma::fill::none);
    arma::mat prob0(n, k, arma::fill::none);
    arma::mat h_est(n, k, arma::fill::none);

    for(int step = 0; step < citermax; ++step) {
        // E step
        for(int i = 0; i < k; ++i) {
            arma::rowvec tmp_mu = mu_old.row(i);
            arma::mat tmp_sigma = sigma_old.slice(i);
            pdf_est.col(i) = cdmvnorm(x, tmp_mu, tmp_sigma);
            prob0.col(i) = pdf_est.col(i) * prop_old(i);
        }

        arma::mat h_est(n, k, arma::fill::none);
        for(int i = 0; i < n; ++i) {
            h_est.row(i) = prob0.row(i)/(sum(prob0.row(i)) * 1.0L);
        }

        // M step
        // update mean
        arma::mat mu_new(k, d, arma::fill::none);
        for(int i = 0; i < k; ++i) {
            mu_new.row(i) = trans(h_est.col(i))*x/(sum(h_est.col(i)) * 1.0L);
        }
        // update sigma
        arma::cube dist(n, d, k, arma::fill::none);
        arma::cube sigma_new = sigma_old;
        arma::mat ones(n, 1, arma::fill::ones);
        for(int i = 0; i < k; ++i) {
            dist.slice(i) = x - ones * mu_new.row(i);
            sigma_new.slice(i) = trans(dist.slice(i)) * diagmat(h_est.col(i)) * dist.slice(i) / (sum(h_est.col(i)) * 1.0L);
        }
        // update proportion
        arma::rowvec prop_new(k, arma::fill::none);
        for(int i = 0; i < k; ++i){
            prop_new(i) = (sum(h_est.col(i)) - lambda * df) / (n-k*lambda*df) * 1.0L;
            if(prop_new(i) < 1E-04) // tolerance greater than 0 for numerical stability (Huang2013)
                prop_new(i) = 0;
        }
        prop_new = prop_new/(sum(prop_new) * 1.0L);

        // calculate difference between two iterations
        delta = sum(abs(prop_new - prop_old));
        
        // eliminate small clusters
        if(sum(prop_new == 0) > 0) {
            arma::uvec idx = find(prop_new > 0);
            k = idx.size();
            prop_old.set_size(k);
            prop_old = trans(prop_new.elem(idx));
            mu_old.set_size(k,d);
            mu_old = mu_new.rows(idx);
            sigma_old.set_size(d,d,k);
            sigma_old = choose_slice(sigma_new, idx, d, k);
            pdf_est = pdf_est.cols(idx);
            prob0 = prob0.cols(idx);
            h_est = h_est.cols(idx);
            delta = 1;
        }
        else{
            prop_old = prop_new;
            mu_old = mu_new;
            sigma_old = sigma_new;
        }
        //calculate cluster with maximum posterior probability
        tag = index_max(h_est, 1);
        
        if(delta < tol)
            break;

        if(k <= 1)
            break;

    }

    return Rcpp::List::create(Rcpp::Named("k") = k,
                              Rcpp::Named("prop") = prop_old,
                              Rcpp::Named("mu") = mu_old,
                              Rcpp::Named("sigma") = sigma_old,
                              Rcpp::Named("pdf_est") = pdf_est,
                              Rcpp::Named("ll") = sum(prob0),
                              Rcpp::Named("cluster") = tag+1);
}

// for testing stuff
//arma::rowvec testfun(arma::mat& x,
//                  arma::rowvec prop,
//                  arma::mat& mu,
//                  arma::cube& sigma,
//                  int k, double df, int citermax,
//                  double lambda) {
//
//    const int n = x.n_rows;
//    const int d = x.n_cols;
//    double delta = 1;
//    double tol = 1E-06;
//    arma::rowvec prop_old = prop;
//    arma::mat mu_old = mu;
//    arma::cube sigma_old = sigma;
//    arma::uvec tag(n, arma::fill::none);
//    arma::mat pdf_est(n, k, arma::fill::none);
//    arma::mat prob0(n, k, arma::fill::none);
//    arma::mat h_est(n, k, arma::fill::none);
//
//    // E step
//    pdf_est.set_size(n, k);
//    h_est.set_size(n, k);
//    h_est = Estep_pGMM(x, n, k, mu_old, sigma_old, prop_old, pdf_est);
//
//    // M step
//    // update mean
//    arma::mat mu_new(k, d, arma::fill::none);
//    for(int i = 0; i < k; ++i) {
//        mu_new.row(i) = trans(h_est.col(i))*x/(sum(h_est.col(i)));
//    }
//    // update sigma
//    arma::cube dist(n, d, k, arma::fill::none);
//    arma::cube sigma_new = sigma_old;
//    arma::mat ones(n, 1, arma::fill::ones);
//    for(int i = 0; i < k; ++i) {
//        dist.slice(i) = x - ones * mu_new.row(i);
//        sigma_new.slice(i) = trans(dist.slice(i)) * diagmat(h_est.col(i)) * dist.slice(i) / (sum(h_est.col(i)));
//    }
//    // update proportion
//    // This is *not* the SCAD penalty
//    arma::rowvec prop_new(k, arma::fill::none);
//    for(int i = 0; i < k; ++i){
//        prop_new(i) = (mean(h_est.col(i)) - lambda * df) / (1-k*lambda*df);
//        if(prop_new(i) < 0)
//            prop_new(i) = 0;
//    }
//    prop_new = prop_new/(sum(prop_new));
//
//    return prop_new;
//
//
//}
