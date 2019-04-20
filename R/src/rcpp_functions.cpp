// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include <RcppArmadillo.h>
#include <algorithm>

using namespace Rcpp;
using namespace RcppArmadillo;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

//
// Stuff for general penalized Gaussian mixture model
//

// Made this so that absolute value of double returns double, not integer
// [[Rcpp::export]]
double abs3(double val){
    return std::abs(val);
}

// First derivative of SCAD penalty
// NOTE THAT THIS PENALTY LOOKS A LIL DIFFERENT THAN SCAD, BECAUSE HUANG PUTS LAMBDA OUTSIDE THE PENALIZATION TERM
// [[Rcpp::export]]
arma::rowvec SCAD_1d(arma::rowvec prop, double lambda, int k, double a = 3.7) {
    arma::colvec term2(k, arma::fill::zeros);
    arma::rowvec out(k, arma::fill::none);

    for(int i = 0; i < k; ++i) {
        if(a*lambda - prop(i) > 0) {
            term2(i) = (a*lambda - prop(i))/(a*lambda-lambda);
        }
        out(i) = ((prop(i) <= lambda) + term2(i)*(prop(i) > lambda));
    }
    return out;
}

// First derivative of SCAD penalty, when only passing an integer and not a vector
// NOTE THAT THIS PENALTY LOOKS A LIL DIFFERENT THAN SCAD, BECAUSE HUANG PUTS LAMBDA OUTSIDE THE PENALIZATION TERM
// [[Rcpp::export]]
double double_SCAD_1d(double prop, double lambda, double a = 3.7) {
    double term2 = 0.0;
    double out;

    if(a*lambda - prop > 0) {
        term2 = (a*lambda - prop)/(a*lambda-lambda);
    }
    out = ((prop <= lambda) + term2*(prop > lambda));
    return out;
}

// SCAD penalty
// NOTE THAT THIS PENALTY LOOKS A LIL DIFFERENT THAN SCAD, BECAUSE HUANG PUTS LAMBDA OUTSIDE THE PENALIZATION TERM
// [[Rcpp::export]]
arma::rowvec SCAD(arma::rowvec prop, double lambda, int k, double a = 3.7) {
    arma::rowvec out(k, arma::fill::none);
    double val;

    for(int i = 0; i < k; ++i) {
        val = abs3(prop(i));
        if(val <= lambda) {
            out(i) = lambda;
        } else if(lambda < val & val <= a * lambda) {
            out(i) = -((pow(val,2)-2*a*lambda*val+pow(lambda,2)) / (2 * (a-1) * lambda));
        } else {
            out(i) = ((a+1) * lambda) / 2;
        }
    }
    return out;
}

// SCAD penalty, when only passing an integer and not a vector
// NOTE THAT THIS PENALTY LOOKS A LIL DIFFERENT THAN SCAD, BECAUSE HUANG PUTS LAMBDA OUTSIDE THE PENALIZATION TERM
// [[Rcpp::export]]
double double_SCAD(double prop, double lambda, double a = 3.7) {
    double out;
    double val;

    val = abs3(prop);
    if(val <= lambda) {
        out = lambda;
    } else if(lambda < val & val <= a * lambda) {
        out = -((pow(val,2)-2*a*lambda*val+pow(lambda,2)) / (2 * (a-1) * lambda));
    } else {
        out = ((a+1) * lambda) / 2;
    }
    return out;
}

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

    return sum((x_cen * sigma.i()) % x_cen, 1); // can probably change sigma.i() to inversion for positive semi-definite matrices (for speed!)
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
    arma::mat h_est(n, k, arma::fill::none);
    //double thresh = 1/(log(n) * sqrt(n));
    double thresh = 1E-03;

    for(int step = 0; step < citermax; ++step) {
        // E step
        for(int i = 0; i < k; ++i) {
            arma::rowvec tmp_mu = mu_old.row(i);
            arma::mat tmp_sigma = sigma_old.slice(i);
            pdf_est.col(i) = cdmvnorm(x, tmp_mu, tmp_sigma);
            prob0.col(i) = pdf_est.col(i) * prop_old(i);
        }

        h_est.set_size(n, k);
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
            if(prop_new(i) < thresh) // tolerance greater than 0 for numerical stability (Huang2013)
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
                              Rcpp::Named("cluster") = tag+1,
                              Rcpp::Named("post_prob") = h_est);
}

//
// Stuff for constrained penalized Gaussian mixture model
//

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

// objective function to be optimized
// [[Rcpp::export]]
double func_to_optim(const arma::colvec& init_val,
                     arma::mat& x,
                     arma::mat& h_est,
                     arma::mat& combos) {

    double mu = init_val(0);
    double sigma = exp(init_val(1));
    double rho = init_val(2);
    int n = x.n_rows;
    int d = x.n_cols;
    int k = h_est.n_cols;
    double nll;

    arma::mat tmp_sigma(d, d, arma::fill::none);
    arma::rowvec tmp_mu(d, arma::fill::none);
    arma::mat pdf_est(n, k, arma::fill::none);

    if( (abs3(rho) >= sigma)) {// || (abs3(rho) >= 1) ) {
        return std::numeric_limits<double>::infinity();
    }

    for(int i = 0; i < k; ++i) {
        // get sigma_in to pass to cget_constr_sigma
        // This amount to finding the diagonal of Sigma
        arma::rowvec sigma_in(abs(sigma*combos.row(i)));
        arma::uvec zeroidx = find(combos.row(i) == 0);
        sigma_in.elem(zeroidx).ones();

        // This min part accounts for the possibility that sigma is actually bigger than 1
        // Need to enforce strict inequality between rho and sigma
        tmp_sigma = cget_constr_sigma(sigma_in, rho, combos.row(i), d); // rho * sigma should constrain rho to be less than sigma in the optimization

        tmp_mu = mu*combos.row(i);
        pdf_est.col(i) = cdmvnorm(x, tmp_mu, tmp_sigma);
    }

    // Adjust for numerically zero values
    arma::uvec zeroidx = find(pdf_est == 0);
    if(zeroidx.size() > 0) {
        arma::uvec posidx = find(pdf_est != 0);
        double clamper = pdf_est.elem(posidx).min();

        arma::colvec populate_vec = arma::ones(zeroidx.size());
        populate_vec = populate_vec * clamper;
        pdf_est.elem(zeroidx) = populate_vec;
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

    try{
        Rcpp::List opt = optim(Rcpp::_["par"] = init_val,
                               Rcpp::_["fn"] = Rcpp::InternalFunction(&func_to_optim),
                               Rcpp::_["method"] = "Nelder-Mead",
                               Rcpp::_["x"] = x,
                               Rcpp::_["h_est"] = h_est,
                               Rcpp::_["combos"] = combos);
        arma::vec mles = Rcpp::as<arma::vec>(opt["par"]);

        return mles;
    }
    catch(...){
        arma::colvec err = { NA_REAL, NA_REAL, NA_REAL };
        return err;
    }
}

// // estimation and model selection of constrained penalized GMM
// // [[Rcpp::export]]
// Rcpp::List cfconstr_pGMM(arma::mat& x,
//                          arma::rowvec prop,
//                          arma::mat mu,
//                          arma::mat sigma,
//                          double rho,
//                          arma::mat combos,
//                          int k,
//                          arma::rowvec df,
//                          int lambda,
//                          int citermax,
//                          double tol,
//                          unsigned int LASSO) {
//
//     const int n = x.n_rows;
//     const int d = x.n_cols;
//     double delta = 1;
//     arma::rowvec prop_old = prop;
//     arma::rowvec prop_new;
//     arma::mat mu_old = mu;
//     arma::mat sigma_old = sigma;
//     double rho_old = rho;
//     arma::uvec tag(n, arma::fill::none);
//     arma::mat pdf_est(n, k, arma::fill::none);
//     arma::mat prob0(n, k, arma::fill::none);
//     arma::mat tmp_sigma(d,d,arma::fill::none);
//     arma::mat h_est(n, k, arma::fill::none);
//     double term; // for SCAD
//     arma::colvec err_test =  { NA_REAL };
//
//     double thresh = 1E-03;
//
//     for(int step = 0; step < citermax; ++step) {
//         // E step
//         for(int i = 0; i < k; ++i) {
//             arma::rowvec tmp_mu = mu_old.row(i);
//             tmp_sigma = cget_constr_sigma(sigma_old.row(i), rho_old, combos.row(i), d);
//
//             try {
//                 pdf_est.col(i) = cdmvnorm(x, tmp_mu, tmp_sigma);
//                 prob0.col(i) = pdf_est.col(i) * prop_old(i);
//             }
//             catch(...){
//                 arma::colvec err = { NA_REAL, NA_REAL, NA_REAL };
//                 return Rcpp::List::create(Rcpp::Named("optim_err") = NA_REAL);
//             }
//         }
//
//         h_est.set_size(n, k);
//         for(int i = 0; i < n; ++i) {
//             h_est.row(i) = prob0.row(i)/(sum(prob0.row(i)) * 1.0L);
//         }
//
//         // M step
//         // update mean and variance covariance with numerical optimization
//
//         // Select the mean and variance associated with reproducibility
//         arma::uvec repidx = find(combos, 0);
//         int idx = repidx(0);
//         double mu_in = abs3(mu_old(idx));
//         double sigma_in = sigma_old(idx);
//         arma::colvec init_val = arma::colvec( { mu_in, log(sigma_in), rho_old } );
//
//         // Optimize using optim (for now)
//         arma::colvec param_new(3, arma::fill::none);
//         param_new = optim_rcpp(init_val, x, h_est, combos);
//
//         if(param_new(0) == err_test(0)) {
//             return Rcpp::List::create(Rcpp::Named("optim_err") = NA_REAL);
//         }
//
//         // transform sigma back
//         param_new(1) = exp(param_new(1));
//
//         prop_new.set_size(k);
//         if(LASSO == 1) {
//             // update proportion via LASSO penalty
//             for(int i = 0; i < k; ++i){
//                 prop_new(i) = (sum(h_est.col(i)) - lambda * df(i)) / (n-lambda*sum(df)) * 1.0L;
//                 if(prop_new(i) < thresh) // tolerance greater than 0 for numerical stability (Huang2013)
//                     prop_new(i) = 0;
//             }
//         } else {
//             // proportion update via SCAD penalty
//             term = accu((SCAD_1d(prop, lambda, k) % prop_old) % (1 / (1E-06 + SCAD(prop, lambda, k))));
//             for(int i = 0; i < k; ++i) {
//                 prop_new(i) = sum(h_est.col(i)) /
//                     (n - (double_SCAD_1d(prop_old(i), lambda) / (1E-06 + double_SCAD(prop_old(i), lambda)) +
//                     term) * lambda * df(i)) * 1.0L;
//                 if(prop_new(i) < thresh) // tolerance greater than 0 for numerical stability (Huang2013)
//                     prop_new(i) = 0;
//             }
//         }
//         prop_new = prop_new/(sum(prop_new) * 1.0L); // renormalize weights
//
//         // calculate difference between two iterations
//         delta = sum(abs(prop_new - prop_old));
//
//         // eliminate small clusters
//         if(sum(prop_new == 0) > 0) {
//             arma::uvec idx = find(prop_new > 0);
//             k = idx.size();
//             prop_old = trans(prop_new.elem(idx));
//             combos = combos.rows(idx);
//             df = trans(df.elem(idx));
//
//             mu_old.set_size(k,d);
//             mu_old = combos*param_new(0);
//
//             sigma_old.set_size(k,d);
//             sigma_old = abs(combos*param_new(1));
//             arma::uvec zeroidx2 = find(sigma_old == 0);
//             sigma_old.elem(zeroidx2).ones();
//
//             rho_old = param_new(2);
//
//             pdf_est = pdf_est.cols(idx);
//             prob0 = prob0.cols(idx);
//             h_est = h_est.cols(idx);
//             delta = 1;
//         }
//         else{
//             prop_old = prop_new;
//             mu_old = combos*param_new(0);
//
//             sigma_old = abs(combos*param_new(1));
//             arma::uvec zeroidx2 = find(sigma_old == 0);
//             sigma_old.elem(zeroidx2).ones();
//
//             rho_old = param_new(2);
//         }
//
//         //calculate cluster with maximum posterior probability
//         tag = index_max(h_est, 1);
//
//         if(delta < tol)
//             break;
//
//         if(k <= 1)
//             break;
//
//     }
//
//     // update the likelihood for output
//     for(int i = 0; i < k; ++i) {
//         arma::rowvec tmp_mu = mu_old.row(i);
//         tmp_sigma = cget_constr_sigma(sigma_old.row(i), rho_old, combos.row(i), d);
//         pdf_est.col(i) = cdmvnorm(x, tmp_mu, tmp_sigma);
//         prob0.col(i) = pdf_est.col(i) * prop_old(i);
//     }
//
//     return Rcpp::List::create(Rcpp::Named("k") = k,
//                               Rcpp::Named("prop") = prop_old,
//                               Rcpp::Named("mu") = mu_old,
//                               Rcpp::Named("sigma") = sigma_old,
//                               Rcpp::Named("rho") = rho_old,
//                               Rcpp::Named("df") = df,
//                               Rcpp::Named("pdf_est") = pdf_est,
//                               Rcpp::Named("ll") = sum(log(sum(prob0,1))),
//                               Rcpp::Named("cluster") = tag+1,
//                               Rcpp::Named("post_prob") = h_est,
//                               Rcpp::Named("combos") = combos);
// }



// objective function to be optimized
// [[Rcpp::export]]
double func_to_optim_bound(const arma::colvec& init_val,
                     arma::mat& x,
                     arma::mat& h_est,
                     arma::mat& combos,
                     double& bound) {

    double mu = init_val(0);
    double sigma = exp(init_val(1));
    double rho = init_val(2);
    int n = x.n_rows;
    int d = x.n_cols;
    int k = h_est.n_cols;
    double nll;

    if( (abs3(rho) >= sigma)) {
        return std::numeric_limits<double>::infinity();
    }

    if( (abs3(mu) < bound)) {
        return std::numeric_limits<double>::infinity();
    }

    arma::mat tmp_sigma(d, d, arma::fill::none);
    arma::rowvec tmp_mu(d, arma::fill::none);
    arma::mat pdf_est(n, k, arma::fill::none);

    for(int i = 0; i < k; ++i) {
        // get sigma_in to pass to cget_constr_sigma
        // This amount to finding the diagonal of Sigma
        arma::rowvec sigma_in(abs(sigma*combos.row(i)));
        arma::uvec zeroidx = find(combos.row(i) == 0);
        sigma_in.elem(zeroidx).ones();

        // This min part accounts for the possibility that sigma is actually bigger than 1
        // Need to enforce strict inequality between rho and sigma
        tmp_sigma = cget_constr_sigma(sigma_in, rho, combos.row(i), d); // rho * sigma should constrain rho to be less than sigma in the optimization
        tmp_mu = mu*combos.row(i);
        pdf_est.col(i) = cdmvnorm(x, tmp_mu, tmp_sigma);
    }

    // Adjust for numerically zero values
    arma::uvec zeroidx = find(pdf_est == 0);
    if(zeroidx.size() > 0) {
        arma::uvec posidx = find(pdf_est != 0);
        double clamper = pdf_est.elem(posidx).min();

        arma::colvec populate_vec = arma::ones(zeroidx.size());
        populate_vec = populate_vec * clamper;
        pdf_est.elem(zeroidx) = populate_vec;
    }

    nll = -accu(h_est % log(pdf_est));
    return nll;
}

// optimize objective function using 'optim' is R-package 'stats'
// [[Rcpp::export]]
arma::vec optim_rcpp_bound(const arma::vec& init_val,
                     arma::mat& x,
                     arma::mat& h_est,
                     arma::mat& combos,
                     double& bound){

    Rcpp::Environment stats("package:stats");
    Rcpp::Function optim = stats["optim"];

    try{
        Rcpp::List opt = optim(Rcpp::_["par"] = init_val,
                               Rcpp::_["fn"] = Rcpp::InternalFunction(&func_to_optim_bound),
                               Rcpp::_["method"] = "Nelder-Mead",
                               Rcpp::_["x"] = x,
                               Rcpp::_["h_est"] = h_est,
                               Rcpp::_["combos"] = combos,
                               Rcpp::_["bound"] = bound);
        arma::vec mles = Rcpp::as<arma::vec>(opt["par"]);

        return mles;
    }
    catch(...){
        arma::colvec err = { NA_REAL, NA_REAL, NA_REAL };
        return err;
    }
}



// estimation and model selection of constrained penalized GMM
// [[Rcpp::export]]
Rcpp::List cfconstr_pgmm(arma::mat& x,
                         arma::rowvec prop,
                         arma::mat mu,
                         arma::mat sigma,
                         double rho,
                         arma::mat combos,
                         int k,
                         arma::rowvec df,
                         int lambda,
                         int citermax,
                         double tol,
                         unsigned int LASSO,
                         double bound = 0.0) {

    const int n = x.n_rows;
    const int d = x.n_cols;
    double delta = 1;
    arma::rowvec prop_old = prop;
    arma::rowvec prop_new;
    arma::mat mu_old = mu;
    arma::mat sigma_old = sigma;
    double rho_old = rho;
    arma::uvec tag(n, arma::fill::none);
    arma::mat pdf_est(n, k, arma::fill::none);
    arma::mat prob0(n, k, arma::fill::none);
    arma::mat tmp_sigma(d,d,arma::fill::none);
    arma::mat h_est(n, k, arma::fill::none);
    double term; // for SCAD
    arma::colvec err_test =  { NA_REAL };

    double thresh = 1E-03;

    for(int step = 0; step < citermax; ++step) {
        // E step
        for(int i = 0; i < k; ++i) {
            arma::rowvec tmp_mu = mu_old.row(i);
            tmp_sigma = cget_constr_sigma(sigma_old.row(i), rho_old, combos.row(i), d);

            try {
                pdf_est.col(i) = cdmvnorm(x, tmp_mu, tmp_sigma);
                prob0.col(i) = pdf_est.col(i) * prop_old(i);
            }
            catch(...){
                arma::colvec err = { NA_REAL, NA_REAL, NA_REAL };
                return Rcpp::List::create(Rcpp::Named("optim_err") = NA_REAL);
            }
        }

        h_est.set_size(n, k);
        for(int i = 0; i < n; ++i) {
            h_est.row(i) = prob0.row(i)/(sum(prob0.row(i)) * 1.0L);
        }

        // M step
        // update mean and variance covariance with numerical optimization

        // Select the mean and variance associated with reproducibility
        arma::uvec repidx = find(combos, 0);
        int idx = repidx(0);
        double mu_in = std::max(abs3(mu_old(idx)), bound);

        double sigma_in = sigma_old(idx);
        arma::colvec init_val = arma::colvec( { mu_in, log(sigma_in), rho_old } );

        // Optimize using optim (for now)
        arma::colvec param_new(3, arma::fill::none);

        if(bound == 0.0) {
            param_new = optim_rcpp(init_val, x, h_est, combos);
        } else {
            param_new = optim_rcpp_bound(init_val, x, h_est, combos, bound);
        }


        if(param_new(0) == err_test(0)) {
            return Rcpp::List::create(Rcpp::Named("optim_err") = NA_REAL);
        }

        // transform sigma back
        param_new(1) = exp(param_new(1));

        prop_new.set_size(k);
        if(LASSO == 1) {
            // update proportion via LASSO penalty
            for(int i = 0; i < k; ++i){
                prop_new(i) = (sum(h_est.col(i)) - lambda * df(i)) / (n-lambda*sum(df)) * 1.0L;
                if(prop_new(i) < thresh) // tolerance greater than 0 for numerical stability (Huang2013)
                    prop_new(i) = 0;
            }
        } else {
            // proportion update via SCAD penalty
            term = accu((SCAD_1d(prop, lambda, k) % prop_old) % (1 / (1E-06 + SCAD(prop, lambda, k))));
            for(int i = 0; i < k; ++i) {
                prop_new(i) = sum(h_est.col(i)) /
                    (n - (double_SCAD_1d(prop_old(i), lambda) / (1E-06 + double_SCAD(prop_old(i), lambda)) +
                        term) * lambda * df(i)) * 1.0L;
                if(prop_new(i) < thresh) // tolerance greater than 0 for numerical stability (Huang2013)
                    prop_new(i) = 0;
            }
        }
        prop_new = prop_new/(sum(prop_new) * 1.0L); // renormalize weights

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
                              Rcpp::Named("cluster") = tag+1,
                              Rcpp::Named("post_prob") = h_est,
                              Rcpp::Named("combos") = combos);
}


// //
// // Stuff for the less constrained (constr0) penalized Gaussian mixture model
// //
//
// // Bind diagonals from covariances into one matrix (for use by optim)
// // [[Rcpp::export]]
// arma::mat bind_diags(arma::cube Sigma_in) {
//     int n = Sigma_in.n_slices;
//     arma::mat diagbind(n, 2, arma::fill::none);
//     for(int i = 0; i < n; ++i) {
//         diagbind.row(i) = Sigma_in.slice(i).diag().t();
//     }
//     return diagbind;
// }
//
// // Bind off-diagonals from covariances into one matrix (for use by optim)
// // [[Rcpp::export]]
// arma::colvec bind_offdiags(arma::cube Sigma_in) {
//     arma::colvec rhobind = Sigma_in.tube(1,0);
//     return rhobind;
// }
//
// // get vector of association-related variances for optim input
// // [[Rcpp::export]]
// arma::colvec get_sigma_optim(arma::mat diagbind, arma::mat combos) {
//     arma::uvec neg = find(combos == -1);
//     arma::uvec pos = find(combos == 1);
//     arma::colvec sig_out;
//
//     if(neg.n_rows > 0 & pos.n_rows > 0){
//         sig_out = { diagbind(neg(0)), diagbind(pos(0)) };
//     } else if(neg.n_rows > 0 & pos.n_rows == 0){
//         sig_out = { diagbind(neg(0)) };
//     } else if(neg.n_rows == 0 & pos.n_rows > 0){
//         sig_out = { diagbind(pos(0)) };
//     }
//     return sig_out;
// }
//
// // get vector of association-related means for optim input
// // [[Rcpp::export]]
// arma::colvec get_mu_optim(arma::mat mu_in, arma::mat combos) {
//     arma::uvec neg = find(combos == -1);
//     arma::uvec pos = find(combos == 1);
//     arma::colvec mu_out;
//
//     if(neg.n_rows > 0 & pos.n_rows > 0){
//         mu_out = { mu_in(neg(0)), mu_in(pos(0)) };
//     } else if(neg.n_rows > 0 & pos.n_rows == 0){
//         mu_out = { mu_in(neg(0)) };
//     } else if(neg.n_rows == 0 & pos.n_rows > 0){
//         mu_out = { mu_in(pos(0)) };
//     }
//     return mu_out;
// }
//
// // get rho for optim
// // [[Rcpp::export]]
// arma::colvec get_rho_optim(arma::colvec rhobind, arma::mat combos) {
//     int n = combos.n_rows;
//     int len = 0;
//     arma::colvec rho_out(3);
//
//     arma::colvec twoneg = { -1, -1 };
//     arma::colvec twopos = { 1, 1 };
//     arma::colvec cross1 = { -1, 1 };
//     arma::colvec cross2 = { 1, -1 };
//
//     for(int i = 0; i < n; ++i){
//         if(combos(i,0) == twoneg(0) & combos(i,1) == twoneg(1)) {
//             rho_out(len) = rhobind(i);
//             ++len;
//             break;
//         }
//     }
//
//     for(int i = 0; i < n; ++i){
//         if((combos(i,0) == cross1(0) & combos(i,1) == cross1(1)) ||
//                 (combos(i,0) == cross2(0) & combos(i,1) == cross2(1))){
//             rho_out(len) = rhobind(i);
//             ++len;
//             break;
//         }
//     }
//
//     for(int i = 0; i < n; ++i){
//         if(combos(i,0) == twopos(0) & combos(i,1) == twopos(1)){
//             rho_out(len) = rhobind(i);
//             ++len;
//             break;
//         }
//     }
//     return rho_out.head(len);
// }
//
// // objective function to be optimized
// // [[Rcpp::export]]
// double func_to_optim0(const arma::colvec& init_val,
//                      const arma::mat& x,
//                      const arma::mat& h_est,
//                      const arma::mat& combos,
//                      const int& a, const int& b, const int& c,
//                      const arma::uvec& negidx) {
//     arma::colvec mu(a);
//     arma::colvec sigma(b);
//     arma::colvec rho(c);
//
//     mu.insert_rows(0, exp(init_val.subvec(0,a-1)));
//     mu.elem(negidx) = -mu.elem(negidx);
//     sigma.insert_rows(0, exp(init_val.subvec(a,a+b-1)));
//     rho.insert_rows(0, init_val.subvec(a+b, a+b+c-1));
//     for(int i = 0; i < c; ++i){
//         rho(i) = trans_rho_inv(rho(i));
//     }
//
//     int n = x.n_rows;
//     int d = x.n_cols;
//     int k = h_est.n_cols;
//     double nll;
//
//     // get appropriate mean vector and covariance matrix to pass to cdmvnorm
//     arma::mat tmp_sigma(d, d, arma::fill::zeros);
//     arma::rowvec tmp_mu(d, arma::fill::zeros);
//     arma::mat pdf_est(n, k, arma::fill::none);
//
//     for(int i = 0; i < k; ++i) {
//         for(int j = 0; j < d; ++j){
//             if(combos(i,j) == -1){
//                 tmp_mu(j) = mu(0);
//                 tmp_sigma(j,j) = sigma(0);
//             } else if(combos(i,j) == 0){
//                 tmp_sigma(j,j) = 1;
//             } else{
//                 tmp_mu(j) = mu(a-1);
//                 tmp_sigma(j,j) = sigma(b-1);
//             }
//         }
//
//         if(combos(i,0) == -1 & combos(i,1) == -1) {
//             tmp_sigma(0,1) = rho(0);
//             tmp_sigma(1,0) = rho(0);
//         } else if(combos(i,0) == 1 & combos(i,1) == 1){
//              tmp_sigma(0,1) = rho(c-1);
//              tmp_sigma(1,0) = rho(c-1);
//         } else if((combos(i,0) == -1 & combos(i,1) == 1) ||
//              (combos(i,0) == 1 & combos(i,1) == -1)){
//             if(rho(0) < 0){
//                 tmp_sigma(0,1) = rho(0);
//                 tmp_sigma(1,0) = rho(0);
//             } else if(rho(1) < 0){
//                 tmp_sigma(0,1) = rho(1);
//                 tmp_sigma(1,0) = rho(1);
//             } else{
//                 tmp_sigma(0,1) = rho(2);
//                 tmp_sigma(1,0) = rho(2);
//             }
//         }
//         pdf_est.col(i) = cdmvnorm(x, tmp_mu, tmp_sigma);
//     }
//     nll = -accu(h_est % log(pdf_est));
//     return nll;
// }
//
// // optimize objective function using 'optim' is R-package 'stats'
// // [[Rcpp::export]]
// arma::vec optim0_rcpp(const arma::vec& init_val,
//                      arma::mat& x,
//                      arma::mat& h_est,
//                      arma::mat& combos,
//                      int& a, int& b, int& c,
//                      arma::uvec& negidx){
//
//     Rcpp::Environment stats("package:stats");
//     Rcpp::Function optim = stats["optim"];
//
//     Rcpp::List opt = optim(Rcpp::_["par"] = init_val,
//                            Rcpp::_["fn"] = Rcpp::InternalFunction(&func_to_optim0),
//                            Rcpp::_["method"] = "Nelder-Mead",
//                            Rcpp::_["x"] = x,
//                            Rcpp::_["h_est"] = h_est,
//                            Rcpp::_["combos"] = combos,
//                            Rcpp::_["a"] = a,
//                            Rcpp::_["b"] = b,
//                            Rcpp::_["c"] = c,
//                            Rcpp::_["negidx"] = negidx);
//     arma::vec mles = Rcpp::as<arma::vec>(opt["par"]);
//
//     return mles;
// }
//
// // estimation and model selection of constrained penalized GMM
// // [[Rcpp::export]]
// Rcpp::List cfconstr0_pGMM(arma::mat& x,
//                          arma::rowvec prop,
//                          arma::mat mu,
//                          arma::cube Sigma,
//                          arma::mat combos,
//                          int k,
//                          arma::rowvec df,
//                          int lambda,
//                          int citermax,
//                          double tol) {
//
//     const int n = x.n_rows;
//     const int d = x.n_cols;
//     double delta = 1;
//     arma::rowvec prop_old = prop;
//     arma::mat mu_old = mu;
//     arma::cube Sigma_old = Sigma;
//     arma::uvec tag(n, arma::fill::none);
//     arma::mat pdf_est(n, k, arma::fill::none);
//     arma::mat prob0(n, k, arma::fill::none);
//     arma::mat tmp_Sigma(d,d,arma::fill::none);
//     arma::mat h_est(n, k, arma::fill::none);
//     arma::colvec init_val;
//     arma::colvec param_new;
//     arma::rowvec tmp_mu(d, arma::fill::none);
//     arma::rowvec prop_new;
//     double thresh = 1E-03;
//
//     for(int step = 0; step < citermax; ++step) {
//         // E step
//         for(int i = 0; i < k; ++i) {
//             tmp_mu = mu_old.row(i);
//             tmp_Sigma = Sigma_old.slice(i);
//
//             pdf_est.col(i) = cdmvnorm(x, tmp_mu, tmp_Sigma);
//             prob0.col(i) = pdf_est.col(i) * prop_old(i);
//         }
//
//         h_est.set_size(n, k);
//         for(int i = 0; i < n; ++i) {
//             h_est.row(i) = prob0.row(i)/(sum(prob0.row(i)) * 1.0L);
//         }
//
//         // M step
//         // update mean and variance covariance with numerical optimization
//         arma::colvec sgn_mu_in = get_mu_optim(mu_old, combos);
//         arma::uvec negidx = find(sgn_mu_in < 0);
//         arma::colvec sigma_in = get_sigma_optim(bind_diags(Sigma_old), combos);
//         arma::colvec rho_in = get_rho_optim(bind_offdiags(Sigma_old), combos);
//         int a = sgn_mu_in.n_rows;
//         int b = sigma_in.n_rows;
//         int c = rho_in.n_rows;
//
//         for(int i = 0; i < c; ++i){
//             rho_in(i) = trans_rho(rho_in(i));
//         }
//
//         init_val.set_size(a+b+c);
//         init_val.subvec(0,a-1) = log(abs(sgn_mu_in));
//         init_val.subvec(a,a+b-1) = log(sigma_in);
//         init_val.subvec(a+b, a+b+c-1) = rho_in;
//
//         // Optimize using optim
//         param_new.copy_size(init_val);
//         param_new = optim0_rcpp(init_val, x, h_est, combos, a, b, c, negidx);
//
//         // Transform parameters back
//         param_new.subvec(0,a-1) = exp(param_new.subvec(0,a-1));
//         param_new.elem(negidx) = -param_new.elem(negidx);
//         param_new.subvec(a,a+b-1) = exp(param_new.subvec(a,a+b-1));
//         for(int i = 0; i < c; ++i){
//             param_new(a+b+i) = trans_rho_inv(param_new(a+b+i));
//         }
//
//         // Make Sigma cube, mu matrix again
//         for(int i = 0; i < k; ++i) {
//             for(int j = 0; j < d; ++j){
//                 if(combos(i,j) == -1){
//                     mu_old(i,j) = param_new(0);
//                     Sigma_old(j,j,i) = param_new(a);
//                 } else if(combos(i,j) == 0){
//                     mu_old(i,j) = 0;
//                     Sigma_old(j,j,i) = 1;
//                 } else{
//                     mu_old(i,j) = param_new(a-1);
//                     Sigma_old(j,j,i) = param_new(a+b-1);
//                 }
//             }
//             if(combos(i,0) == -1 & combos(i,1) == -1) {
//                 Sigma_old(0,1,i) = param_new(a+b);
//                 Sigma_old(1,0,i) = param_new(a+b);
//             } else if(combos(i,0) == 1 & combos(i,1) == 1){
//                 Sigma_old(0,1,i) = param_new(a+b+c-1);
//                 Sigma_old(1,0,i) = param_new(a+b+c-1);
//             } else if((combos(i,0) == -1 & combos(i,1) == 1) ||
//                 (combos(i,0) == 1 & combos(i,1) == -1)){
//                 if(param_new(a+b) < 0){
//                     Sigma_old(0,1,i) = param_new(a+b);
//                     Sigma_old(1,0,i) = param_new(a+b);
//                 } else if(param_new(a+b+1) < 0){
//                     Sigma_old(0,1,i) = param_new(a+b+1);
//                     Sigma_old(1,0,i) = param_new(a+b+1);
//                 } else{
//                     Sigma_old(0,1,i) = param_new(a+b+c-1);
//                     Sigma_old(1,0,i) = param_new(a+b+c-1);
//                 }
//             }
//         }
//
//         // update proportion
//         prop_new.set_size(k);
//         for(int i = 0; i < k; ++i){
//             prop_new(i) = (sum(h_est.col(i)) - lambda * df(i)) / (n-lambda*sum(df)) * 1.0L;
//             if(prop_new(i) < thresh) // tolerance greater than 0 for numerical stability (Huang2013)
//                 prop_new(i) = 0;
//         }
//         prop_new = prop_new/(sum(prop_new) * 1.0L);
//
//         // calculate difference between two iterations
//         delta = sum(abs(prop_new - prop_old));
//
//         // eliminate small clusters
//         if(sum(prop_new == 0) > 0) {
//             arma::uvec idx = find(prop_new > 0);
//             k = idx.size();
//             prop_old = trans(prop_new.elem(idx));
//             combos = combos.rows(idx);
//             df = trans(df.elem(idx));
//
//             mu_old = mu_old.rows(idx);
//             Sigma_old = choose_slice(Sigma_old, idx, d, k);
//
//             pdf_est = pdf_est.cols(idx);
//             prob0 = prob0.cols(idx);
//             h_est = h_est.cols(idx);
//             delta = 1;
//         }
//         else{
//             prop_old = prop_new;
//         }
//
//         //calculate cluster with maximum posterior probability
//         tag = index_max(h_est, 1);
//
//     if(delta < tol)
//         break;
//     if(k <= 1)
//         break;
// }
//
//     // update the likelihood for output
//     for(int i = 0; i < k; ++i) {
//         arma::rowvec tmp_mu = mu_old.row(i);
//         tmp_Sigma = Sigma_old.slice(i);
//         pdf_est.col(i) = cdmvnorm(x, tmp_mu, tmp_Sigma);
//         prob0.col(i) = pdf_est.col(i) * prop_old(i);
//     }
//
//     return Rcpp::List::create(Rcpp::Named("k") = k,
//                               Rcpp::Named("prop") = prop_old,
//                               Rcpp::Named("mu") = mu_old,
//                               Rcpp::Named("Sigma") = Sigma_old,
//                               Rcpp::Named("df") = df,
//                               Rcpp::Named("pdf_est") = pdf_est,
//                               Rcpp::Named("ll") = sum(log(sum(prob0,1))),
//                               Rcpp::Named("cluster") = tag+1,
//                               Rcpp::Named("post_prob") = h_est);
// }

// Compute density of univariate normal
// [[Rcpp::export]]
arma::colvec cduvnorm(arma::colvec x, double mu, double sigma){
    arma::colvec mdist = ((x-mu) % (x-mu))/(sigma);
    const double log2pi = std::log(2.0 * M_PI);
    double logcoef = -(log2pi + log(sigma));
    arma::colvec logretval = (logcoef - mdist)/2;

    return exp(logretval);
}

// computing marginal likelihood of constrained gmm (for copula likelihood)
// [[Rcpp::export]]
double cmarg_ll_gmm(arma::mat& z,
                    arma::mat mu,
                    arma::mat sigma,
                    arma::rowvec prop,
                    arma::mat combos,
                    int k) {
    const int n = z.n_rows;
    const int d = z.n_cols;
    arma::mat mll(n,d,arma::fill::none);
    arma::mat pdf_est(n,k,arma::fill::none);
    double tmp_sigma;
    double tmp_mu;

    for(int i = 0; i < d; ++i){
        for(int j = 0; j < k; ++j) {
            tmp_mu = mu(j,i);
            tmp_sigma = sigma(j,i);
            pdf_est.col(j) = prop(j) * cduvnorm(z.col(i), tmp_mu, tmp_sigma);
        }
        mll.col(i) = log(sum(pdf_est,1));
    }
    return accu(mll);
}

// computing joint likelihood of constrained gmm (for copula likelihood)
// [[Rcpp::export]]
double cll_gmm(arma::mat& z,
               arma::mat mu,
               arma::mat sigma,
               double rho,
               arma::rowvec prop,
               arma::mat combos,
               int k) {
    const int n = z.n_rows;
    const int d = z.n_cols;
    arma::rowvec tmp_mu(d, arma::fill::none);
    arma::mat tmp_sigma(d, d, arma::fill::none);
    arma::mat pdf_est(n,k,arma::fill::none);

    for(int i = 0; i < k; ++i) {
        tmp_mu = mu.row(i);
        tmp_sigma = cget_constr_sigma(sigma.row(i), rho, combos.row(i), d);
        pdf_est.col(i) = prop(i) * cdmvnorm(z, tmp_mu, tmp_sigma);
    }

    return accu(log(sum(pdf_est,1)));
}

// // computing marginal likelihood of LESS constrained (0) gmm (for copula likelihood)
// // [[Rcpp::export]]
// double cmarg0_ll_gmm(arma::mat& z,
//                     arma::mat mu,
//                     arma::cube Sigma,
//                     arma::rowvec prop,
//                     int k) {
//     const int n = z.n_rows;
//     const int d = z.n_cols;
//     arma::mat mll(n,d,arma::fill::none);
//     arma::mat pdf_est(n,k,arma::fill::none);
//     double tmp_sigma;
//     double tmp_mu;
//     arma::mat tmp_Sigma(d,d,arma::fill::none);
//
//     for(int i = 0; i < d; ++i){
//         for(int j = 0; j < k; ++j) {
//             tmp_mu = mu(j,i);
//             tmp_Sigma = Sigma.slice(j);
//             tmp_sigma = tmp_Sigma(i,i);
//             pdf_est.col(j) = prop(j) * cduvnorm(z.col(i), tmp_mu, tmp_sigma);
//         }
//         mll.col(i) = log(sum(pdf_est,1));
//     }
//     return accu(mll);
// }
//
// // computing joint likelihood of LESS xconstrained (0) gmm (for copula likelihood)
// // [[Rcpp::export]]
// double cll0_gmm(arma::mat& z,
//                 arma::mat mu,
//                 arma::cube Sigma,
//                 arma::rowvec prop,
//                 int k) {
//     const int n = z.n_rows;
//     const int d = z.n_cols;
//     arma::rowvec tmp_mu(d, arma::fill::none);
//     arma::mat tmp_sigma(d, d, arma::fill::none);
//     arma::mat pdf_est(n,k,arma::fill::none);
//
//     for(int i = 0; i < k; ++i) {
//         tmp_mu = mu.row(i);
//         tmp_sigma = Sigma.slice(i);
//         pdf_est.col(i) = prop(i) * cdmvnorm(z, tmp_mu, tmp_sigma);
//     }
//     return accu(log(sum(pdf_est,1)));
// }
//
//
