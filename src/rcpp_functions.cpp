// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace RcppArmadillo;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]
// [[Rcpp::plugins(cpp11)]]

//
// Stuff for general penalized Gaussian mixture model
//

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
                double lambda,
                int citermax,
                double tol
                ) {

    const int n = x.n_rows;
    const int d = x.n_cols;
    double delta = 1;
    arma::rowvec prop_old = prop;
    arma::mat mu_old = mu;
    arma::cube sigma_old = sigma;
    arma::uvec tag(n, arma::fill::none);
    arma::mat pdf_est(n, k, arma::fill::none);
    arma::mat prob0(n, k, arma::fill::none);

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

    // update likelihood for output
    for(int i = 0; i < k; ++i) {
        arma::rowvec tmp_mu = mu_old.row(i);
        arma::mat tmp_sigma = sigma_old.slice(i);
        pdf_est.col(i) = cdmvnorm(x, tmp_mu, tmp_sigma);
        prob0.col(i) = pdf_est.col(i) * prop_old(i);
    }

    return Rcpp::List::create(Rcpp::Named("k") = k,
                              Rcpp::Named("prop") = prop_old,
                              Rcpp::Named("mu") = mu_old,
                              Rcpp::Named("sigma") = sigma_old,
                              Rcpp::Named("pdf_est") = pdf_est,
                              Rcpp::Named("ll") = sum(log(sum(prob0,1))),
                              Rcpp::Named("cluster") = tag+1);
}

//
// Stuff for constrained penalized Gaussian mixture model
//

// Made this so that absolute value of double returns double, not integer
// [[Rcpp::export]]
double abs3(double val){
    return std::abs(val);
}

// Calculate variance covariance matrix for constrained pGMM
// [[Rcpp::export]]
arma::mat cget_constr_sigma(arma::rowvec sigma, double rho, arma::rowvec combos, int d){
    arma::mat Sigma = diagmat(sigma);
     for(int i = 0; i < d-1; ++i){
         for(int j = i+1; j < d; ++j){
            if(combos(i) == combos(j) & combos(i) != 0){
                Sigma(i,j) = rho;
                Sigma(j,i) = rho;
            } else if(combos(i) == -combos(j) & combos(i) != 0){
                Sigma(i,j) = -rho;
                Sigma(j,i) = -rho;
            }
        }
    }
    return Sigma;
}

// Transform rho to constrain to be in [-1,1]
// [[Rcpp::export]]
double trans_rho(double rho) {
    double out = log(rho + 1) - log(1 - rho);
    return out;

}

// Transform rho to constrain to be in [-1,1]
// [[Rcpp::export]]
double trans_rho_inv(double rho) {
    double out = (exp(rho) - 1) / (exp(rho) + 1);
    return out;
}

// objective function to be optimized
// [[Rcpp::export]]
double func_to_optim(const arma::colvec& init_val,
                     arma::mat& x,
                     arma::mat& h_est,
                     arma::mat& combos) {

    double mu = init_val(0);
    double sigma = exp(init_val(1));
    double rho = trans_rho_inv(init_val(2));
    int n = x.n_rows;
    int d = x.n_cols;
    int k = h_est.n_cols;
    double nll;

    arma::mat tmp_sigma(d, d, arma::fill::none);
    arma::rowvec tmp_mu(d, arma::fill::none);
    arma::mat pdf_est(n, k, arma::fill::none);

    for(int i = 0; i < k; ++i) {
        // get sigma_in to pass to cget_constr_sigma
        // This amount to finding the diagonal of Sigma
        arma::rowvec sigma_in(abs(sigma*combos.row(i)));
        arma::uvec zeroidx = find(combos.row(i) == 0);
        sigma_in.elem(zeroidx).ones();

        tmp_sigma = cget_constr_sigma(sigma_in, rho, combos.row(i), d);
        tmp_mu = mu*combos.row(i);
        pdf_est.col(i) = cdmvnorm(x, tmp_mu, tmp_sigma);
    }

    nll = -accu(h_est % log(pdf_est));
    return nll;
}

// optimize objective function using 'optim' is R-package 'stats'
// [[Rcpp::export]]
arma::vec optim_rcpp(const arma::vec& init_val,
                     arma::mat& x,
                     arma::mat& h_est,
                     arma::mat& combos){

    Rcpp::Environment stats("package:stats");
    Rcpp::Function optim = stats["optim"];

    Rcpp::List opt = optim(Rcpp::_["par"] = init_val,
                           Rcpp::_["fn"] = Rcpp::InternalFunction(&func_to_optim),
                           Rcpp::_["method"] = "Nelder-Mead",
                           Rcpp::_["x"] = x,
                           Rcpp::_["h_est"] = h_est,
                           Rcpp::_["combos"] = combos);

    arma::vec mles = Rcpp::as<arma::vec>(opt["par"]);

    return mles;
}

// estimation and model selection of constrained penalized GMM
// [[Rcpp::export]]
Rcpp::List cfconstr_pGMM(arma::mat& x,
                         arma::rowvec prop,
                         arma::mat mu,
                         arma::mat sigma,
                         double rho,
                         arma::mat combos,
                         int k,
                         arma::rowvec df,
                         int lambda,
                         int citermax,
                         double tol) {

    const int n = x.n_rows;
    const int d = x.n_cols;
    double delta = 1;
    arma::rowvec prop_old = prop;
    arma::mat mu_old = mu;
    arma::mat sigma_old = sigma;
    double rho_old = rho;
    arma::uvec tag(n, arma::fill::none);
    arma::mat pdf_est(n, k, arma::fill::none);
    arma::mat prob0(n, k, arma::fill::none);
    arma::mat tmp_sigma(d,d,arma::fill::none);

    for(int step = 0; step < citermax; ++step) {
        // E step
        for(int i = 0; i < k; ++i) {
            arma::rowvec tmp_mu = mu_old.row(i);
            tmp_sigma = cget_constr_sigma(sigma_old.row(i), rho_old, combos.row(i), d);

            pdf_est.col(i) = cdmvnorm(x, tmp_mu, tmp_sigma);
            prob0.col(i) = pdf_est.col(i) * prop_old(i);
        }

        arma::mat h_est(n, k, arma::fill::none);
        for(int i = 0; i < n; ++i) {
            h_est.row(i) = prob0.row(i)/(sum(prob0.row(i)) * 1.0L);
        }

        // M step
        // update mean and variance covariance with numerical optimization

        // Select the mean and variance associated with reproducibility
        arma::uvec repidx = find(combos, 0);
        int idx = repidx(0);
        double mu_in = abs3(mu_old(idx));
        double sigma_in = sigma_old(idx);
        arma::colvec init_val = arma::colvec({mu_in, log(sigma_in), trans_rho(rho_old)});

        // Optimize using optim (for now)
        arma::colvec param_new(3, arma::fill::none);
        param_new = optim_rcpp(init_val, x, h_est, combos);

        // transform sigma, rho back
        param_new(1) = exp(param_new(1));
        param_new(2) = trans_rho_inv(param_new(2));

        // update proportion
        arma::rowvec prop_new(k, arma::fill::none);
        for(int i = 0; i < k; ++i){
            prop_new(i) = (sum(h_est.col(i)) - lambda * df(i)) / (n-lambda*sum(df)) * 1.0L;
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
            prop_old = trans(prop_new.elem(idx));
            combos = combos.rows(idx);
            df = trans(df.elem(idx));

            mu_old.set_size(k,d);
            mu_old = combos*param_new(0);

            sigma_old.set_size(k,d);
            sigma_old = abs(combos*param_new(1));
            arma::uvec zeroidx2 = find(sigma_old == 0);
            sigma_old.elem(zeroidx2).ones();

            rho_old = param_new(2);

            pdf_est = pdf_est.cols(idx);
            prob0 = prob0.cols(idx);
            h_est = h_est.cols(idx);
            delta = 1;
        }
        else{
            prop_old = prop_new;
            mu_old = combos*param_new(0);

            sigma_old = abs(combos*param_new(1));
            arma::uvec zeroidx2 = find(sigma_old == 0);
            sigma_old.elem(zeroidx2).ones();

            rho_old = param_new(2);
        }

        //calculate cluster with maximum posterior probability
        tag = index_max(h_est, 1);

        if(delta < tol)
            break;

        if(k <= 1)
            break;

    }

    // update the likelihood for output
    for(int i = 0; i < k; ++i) {
        arma::rowvec tmp_mu = mu_old.row(i);
        tmp_sigma = cget_constr_sigma(sigma_old.row(i), rho_old, combos.row(i), d);
        pdf_est.col(i) = cdmvnorm(x, tmp_mu, tmp_sigma);
        prob0.col(i) = pdf_est.col(i) * prop_old(i);
    }

    return Rcpp::List::create(Rcpp::Named("k") = k,
                              Rcpp::Named("prop") = prop_old,
                              Rcpp::Named("mu") = mu_old,
                              Rcpp::Named("sigma") = sigma_old,
                              Rcpp::Named("rho") = rho_old,
                              Rcpp::Named("df") = df,
                              Rcpp::Named("pdf_est") = pdf_est,
                              Rcpp::Named("ll") = sum(log(sum(prob0,1))),
                              Rcpp::Named("cluster") = tag+1);
}

// test stuff
// [[Rcpp::export]]
arma::mat teststuff(arma::mat combos, double sigma0) {
    arma::uvec zeroidx = find(combos == 0);
    arma::mat sigma = abs(sigma0*combos);
    sigma.elem(zeroidx).ones();
    return sigma;
}

