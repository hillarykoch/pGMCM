// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// cgetPaths
arma::mat cgetPaths(std::string filepath, int len_filt_h, Rcpp::List nonconsec, Rcpp::List mus, arma::mat labels, int n, int dist_tol);
RcppExport SEXP _pGMCM_cgetPaths(SEXP filepathSEXP, SEXP len_filt_hSEXP, SEXP nonconsecSEXP, SEXP musSEXP, SEXP labelsSEXP, SEXP nSEXP, SEXP dist_tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type filepath(filepathSEXP);
    Rcpp::traits::input_parameter< int >::type len_filt_h(len_filt_hSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type nonconsec(nonconsecSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type mus(musSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type labels(labelsSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type dist_tol(dist_tolSEXP);
    rcpp_result_gen = Rcpp::wrap(cgetPaths(filepath, len_filt_h, nonconsec, mus, labels, n, dist_tol));
    return rcpp_result_gen;
END_RCPP
}
// crowMatch
arma::uvec crowMatch(arma::mat assoc, arma::mat nonconsec);
RcppExport SEXP _pGMCM_crowMatch(SEXP assocSEXP, SEXP nonconsecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type assoc(assocSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type nonconsec(nonconsecSEXP);
    rcpp_result_gen = Rcpp::wrap(crowMatch(assoc, nonconsec));
    return rcpp_result_gen;
END_RCPP
}
// cget_prior_count
arma::colvec cget_prior_count(arma::mat red_class, Rcpp::List mus, arma::mat labels, int d, int n, int dist_tol);
RcppExport SEXP _pGMCM_cget_prior_count(SEXP red_classSEXP, SEXP musSEXP, SEXP labelsSEXP, SEXP dSEXP, SEXP nSEXP, SEXP dist_tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type red_class(red_classSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type mus(musSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type labels(labelsSEXP);
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type dist_tol(dist_tolSEXP);
    rcpp_result_gen = Rcpp::wrap(cget_prior_count(red_class, mus, labels, d, n, dist_tol));
    return rcpp_result_gen;
END_RCPP
}
// cget_true_assoc_idx
arma::colvec cget_true_assoc_idx(arma::mat red_class, arma::mat true_assoc);
RcppExport SEXP _pGMCM_cget_true_assoc_idx(SEXP red_classSEXP, SEXP true_assocSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type red_class(red_classSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type true_assoc(true_assocSEXP);
    rcpp_result_gen = Rcpp::wrap(cget_true_assoc_idx(red_class, true_assoc));
    return rcpp_result_gen;
END_RCPP
}
// trans_func
int trans_func(double& x);
RcppExport SEXP _pGMCM_trans_func(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(trans_func(x));
    return rcpp_result_gen;
END_RCPP
}
// get_list_names
std::vector<std::string> get_list_names(Rcpp::List L);
RcppExport SEXP _pGMCM_get_list_names(SEXP LSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type L(LSEXP);
    rcpp_result_gen = Rcpp::wrap(get_list_names(L));
    return rcpp_result_gen;
END_RCPP
}
// cpaste0
std::string cpaste0(std::vector<std::string> str1);
RcppExport SEXP _pGMCM_cpaste0(SEXP str1SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<std::string> >::type str1(str1SEXP);
    rcpp_result_gen = Rcpp::wrap(cpaste0(str1));
    return rcpp_result_gen;
END_RCPP
}
// cstr_split
arma::mat cstr_split(std::vector<std::string> strings, std::string split);
RcppExport SEXP _pGMCM_cstr_split(SEXP stringsSEXP, SEXP splitSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<std::string> >::type strings(stringsSEXP);
    Rcpp::traits::input_parameter< std::string >::type split(splitSEXP);
    rcpp_result_gen = Rcpp::wrap(cstr_split(strings, split));
    return rcpp_result_gen;
END_RCPP
}
// caccept
arma::vec caccept(arma::mat x, arma::colvec y);
RcppExport SEXP _pGMCM_caccept(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(caccept(x, y));
    return rcpp_result_gen;
END_RCPP
}
// caccept2
arma::vec caccept2(arma::mat x, double y1, double y2);
RcppExport SEXP _pGMCM_caccept2(SEXP xSEXP, SEXP y1SEXP, SEXP y2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type y1(y1SEXP);
    Rcpp::traits::input_parameter< double >::type y2(y2SEXP);
    rcpp_result_gen = Rcpp::wrap(caccept2(x, y1, y2));
    return rcpp_result_gen;
END_RCPP
}
// cprune_path
bool cprune_path(Rcpp::List nonconsec, std::vector<int> assoc);
RcppExport SEXP _pGMCM_cprune_path(SEXP nonconsecSEXP, SEXP assocSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type nonconsec(nonconsecSEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type assoc(assocSEXP);
    rcpp_result_gen = Rcpp::wrap(cprune_path(nonconsec, assoc));
    return rcpp_result_gen;
END_RCPP
}
// cassociate
std::vector<int> cassociate(std::vector<int> path, std::string filepath, int len_filt_h);
RcppExport SEXP _pGMCM_cassociate(SEXP pathSEXP, SEXP filepathSEXP, SEXP len_filt_hSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<int> >::type path(pathSEXP);
    Rcpp::traits::input_parameter< std::string >::type filepath(filepathSEXP);
    Rcpp::traits::input_parameter< int >::type len_filt_h(len_filt_hSEXP);
    rcpp_result_gen = Rcpp::wrap(cassociate(path, filepath, len_filt_h));
    return rcpp_result_gen;
END_RCPP
}
// cprune_path2
double cprune_path2(std::vector<int> assoc, Rcpp::List mus, arma::mat labels, int d, int n, int dist_tol);
RcppExport SEXP _pGMCM_cprune_path2(SEXP assocSEXP, SEXP musSEXP, SEXP labelsSEXP, SEXP dSEXP, SEXP nSEXP, SEXP dist_tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<int> >::type assoc(assocSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type mus(musSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type labels(labelsSEXP);
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type dist_tol(dist_tolSEXP);
    rcpp_result_gen = Rcpp::wrap(cprune_path2(assoc, mus, labels, d, n, dist_tol));
    return rcpp_result_gen;
END_RCPP
}
// abs3
double abs3(double val);
RcppExport SEXP _pGMCM_abs3(SEXP valSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type val(valSEXP);
    rcpp_result_gen = Rcpp::wrap(abs3(val));
    return rcpp_result_gen;
END_RCPP
}
// SCAD_1d
arma::rowvec SCAD_1d(arma::rowvec prop, double lambda, int k, double a);
RcppExport SEXP _pGMCM_SCAD_1d(SEXP propSEXP, SEXP lambdaSEXP, SEXP kSEXP, SEXP aSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::rowvec >::type prop(propSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    rcpp_result_gen = Rcpp::wrap(SCAD_1d(prop, lambda, k, a));
    return rcpp_result_gen;
END_RCPP
}
// double_SCAD_1d
double double_SCAD_1d(double prop, double lambda, double a);
RcppExport SEXP _pGMCM_double_SCAD_1d(SEXP propSEXP, SEXP lambdaSEXP, SEXP aSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type prop(propSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    rcpp_result_gen = Rcpp::wrap(double_SCAD_1d(prop, lambda, a));
    return rcpp_result_gen;
END_RCPP
}
// SCAD
arma::rowvec SCAD(arma::rowvec prop, double lambda, int k, double a);
RcppExport SEXP _pGMCM_SCAD(SEXP propSEXP, SEXP lambdaSEXP, SEXP kSEXP, SEXP aSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::rowvec >::type prop(propSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    rcpp_result_gen = Rcpp::wrap(SCAD(prop, lambda, k, a));
    return rcpp_result_gen;
END_RCPP
}
// double_SCAD
double double_SCAD(double prop, double lambda, double a);
RcppExport SEXP _pGMCM_double_SCAD(SEXP propSEXP, SEXP lambdaSEXP, SEXP aSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type prop(propSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    rcpp_result_gen = Rcpp::wrap(double_SCAD(prop, lambda, a));
    return rcpp_result_gen;
END_RCPP
}
// choose_slice
arma::cube choose_slice(arma::cube& Q, arma::uvec idx, int d, int k);
RcppExport SEXP _pGMCM_choose_slice(SEXP QSEXP, SEXP idxSEXP, SEXP dSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cube& >::type Q(QSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type idx(idxSEXP);
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(choose_slice(Q, idx, d, k));
    return rcpp_result_gen;
END_RCPP
}
// Mahalanobis
arma::vec Mahalanobis(arma::mat x, arma::rowvec mu, arma::mat sigma);
RcppExport SEXP _pGMCM_Mahalanobis(SEXP xSEXP, SEXP muSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::rowvec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma(sigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(Mahalanobis(x, mu, sigma));
    return rcpp_result_gen;
END_RCPP
}
// cdmvnorm
arma::vec cdmvnorm(arma::mat x, arma::rowvec mu, arma::mat sigma);
RcppExport SEXP _pGMCM_cdmvnorm(SEXP xSEXP, SEXP muSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::rowvec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma(sigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(cdmvnorm(x, mu, sigma));
    return rcpp_result_gen;
END_RCPP
}
// cfpGMM
Rcpp::List cfpGMM(arma::mat& x, arma::rowvec prop, arma::mat& mu, arma::cube& sigma, int k, double df, double lambda, int citermax, double tol);
RcppExport SEXP _pGMCM_cfpGMM(SEXP xSEXP, SEXP propSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP kSEXP, SEXP dfSEXP, SEXP lambdaSEXP, SEXP citermaxSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::rowvec >::type prop(propSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::cube& >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< double >::type df(dfSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< int >::type citermax(citermaxSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(cfpGMM(x, prop, mu, sigma, k, df, lambda, citermax, tol));
    return rcpp_result_gen;
END_RCPP
}
// cget_constr_sigma
arma::mat cget_constr_sigma(arma::rowvec sigma, double rho, arma::rowvec combos, int d);
RcppExport SEXP _pGMCM_cget_constr_sigma(SEXP sigmaSEXP, SEXP rhoSEXP, SEXP combosSEXP, SEXP dSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::rowvec >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< arma::rowvec >::type combos(combosSEXP);
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    rcpp_result_gen = Rcpp::wrap(cget_constr_sigma(sigma, rho, combos, d));
    return rcpp_result_gen;
END_RCPP
}
// trans_rho
double trans_rho(double rho);
RcppExport SEXP _pGMCM_trans_rho(SEXP rhoSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    rcpp_result_gen = Rcpp::wrap(trans_rho(rho));
    return rcpp_result_gen;
END_RCPP
}
// trans_rho_inv
double trans_rho_inv(double rho);
RcppExport SEXP _pGMCM_trans_rho_inv(SEXP rhoSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    rcpp_result_gen = Rcpp::wrap(trans_rho_inv(rho));
    return rcpp_result_gen;
END_RCPP
}
// func_to_optim
double func_to_optim(const arma::colvec& init_val, arma::mat& x, arma::mat& h_est, arma::mat& combos);
RcppExport SEXP _pGMCM_func_to_optim(SEXP init_valSEXP, SEXP xSEXP, SEXP h_estSEXP, SEXP combosSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type init_val(init_valSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type h_est(h_estSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type combos(combosSEXP);
    rcpp_result_gen = Rcpp::wrap(func_to_optim(init_val, x, h_est, combos));
    return rcpp_result_gen;
END_RCPP
}
// optim_rcpp
arma::vec optim_rcpp(const arma::vec& init_val, arma::mat& x, arma::mat& h_est, arma::mat& combos);
RcppExport SEXP _pGMCM_optim_rcpp(SEXP init_valSEXP, SEXP xSEXP, SEXP h_estSEXP, SEXP combosSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type init_val(init_valSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type h_est(h_estSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type combos(combosSEXP);
    rcpp_result_gen = Rcpp::wrap(optim_rcpp(init_val, x, h_est, combos));
    return rcpp_result_gen;
END_RCPP
}
// cfconstr_pGMM
Rcpp::List cfconstr_pGMM(arma::mat& x, arma::rowvec prop, arma::mat mu, arma::mat sigma, double rho, arma::mat combos, int k, arma::rowvec df, int lambda, int citermax, double tol, unsigned int LASSO);
RcppExport SEXP _pGMCM_cfconstr_pGMM(SEXP xSEXP, SEXP propSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP rhoSEXP, SEXP combosSEXP, SEXP kSEXP, SEXP dfSEXP, SEXP lambdaSEXP, SEXP citermaxSEXP, SEXP tolSEXP, SEXP LASSOSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::rowvec >::type prop(propSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type combos(combosSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< arma::rowvec >::type df(dfSEXP);
    Rcpp::traits::input_parameter< int >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< int >::type citermax(citermaxSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type LASSO(LASSOSEXP);
    rcpp_result_gen = Rcpp::wrap(cfconstr_pGMM(x, prop, mu, sigma, rho, combos, k, df, lambda, citermax, tol, LASSO));
    return rcpp_result_gen;
END_RCPP
}
// bind_diags
arma::mat bind_diags(arma::cube Sigma_in);
RcppExport SEXP _pGMCM_bind_diags(SEXP Sigma_inSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cube >::type Sigma_in(Sigma_inSEXP);
    rcpp_result_gen = Rcpp::wrap(bind_diags(Sigma_in));
    return rcpp_result_gen;
END_RCPP
}
// bind_offdiags
arma::colvec bind_offdiags(arma::cube Sigma_in);
RcppExport SEXP _pGMCM_bind_offdiags(SEXP Sigma_inSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cube >::type Sigma_in(Sigma_inSEXP);
    rcpp_result_gen = Rcpp::wrap(bind_offdiags(Sigma_in));
    return rcpp_result_gen;
END_RCPP
}
// get_sigma_optim
arma::colvec get_sigma_optim(arma::mat diagbind, arma::mat combos);
RcppExport SEXP _pGMCM_get_sigma_optim(SEXP diagbindSEXP, SEXP combosSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type diagbind(diagbindSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type combos(combosSEXP);
    rcpp_result_gen = Rcpp::wrap(get_sigma_optim(diagbind, combos));
    return rcpp_result_gen;
END_RCPP
}
// get_mu_optim
arma::colvec get_mu_optim(arma::mat mu_in, arma::mat combos);
RcppExport SEXP _pGMCM_get_mu_optim(SEXP mu_inSEXP, SEXP combosSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type mu_in(mu_inSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type combos(combosSEXP);
    rcpp_result_gen = Rcpp::wrap(get_mu_optim(mu_in, combos));
    return rcpp_result_gen;
END_RCPP
}
// get_rho_optim
arma::colvec get_rho_optim(arma::colvec rhobind, arma::mat combos);
RcppExport SEXP _pGMCM_get_rho_optim(SEXP rhobindSEXP, SEXP combosSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::colvec >::type rhobind(rhobindSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type combos(combosSEXP);
    rcpp_result_gen = Rcpp::wrap(get_rho_optim(rhobind, combos));
    return rcpp_result_gen;
END_RCPP
}
// func_to_optim0
double func_to_optim0(const arma::colvec& init_val, const arma::mat& x, const arma::mat& h_est, const arma::mat& combos, const int& a, const int& b, const int& c, const arma::uvec& negidx);
RcppExport SEXP _pGMCM_func_to_optim0(SEXP init_valSEXP, SEXP xSEXP, SEXP h_estSEXP, SEXP combosSEXP, SEXP aSEXP, SEXP bSEXP, SEXP cSEXP, SEXP negidxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type init_val(init_valSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type h_est(h_estSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type combos(combosSEXP);
    Rcpp::traits::input_parameter< const int& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const int& >::type b(bSEXP);
    Rcpp::traits::input_parameter< const int& >::type c(cSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type negidx(negidxSEXP);
    rcpp_result_gen = Rcpp::wrap(func_to_optim0(init_val, x, h_est, combos, a, b, c, negidx));
    return rcpp_result_gen;
END_RCPP
}
// optim0_rcpp
arma::vec optim0_rcpp(const arma::vec& init_val, arma::mat& x, arma::mat& h_est, arma::mat& combos, int& a, int& b, int& c, arma::uvec& negidx);
RcppExport SEXP _pGMCM_optim0_rcpp(SEXP init_valSEXP, SEXP xSEXP, SEXP h_estSEXP, SEXP combosSEXP, SEXP aSEXP, SEXP bSEXP, SEXP cSEXP, SEXP negidxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type init_val(init_valSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type h_est(h_estSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type combos(combosSEXP);
    Rcpp::traits::input_parameter< int& >::type a(aSEXP);
    Rcpp::traits::input_parameter< int& >::type b(bSEXP);
    Rcpp::traits::input_parameter< int& >::type c(cSEXP);
    Rcpp::traits::input_parameter< arma::uvec& >::type negidx(negidxSEXP);
    rcpp_result_gen = Rcpp::wrap(optim0_rcpp(init_val, x, h_est, combos, a, b, c, negidx));
    return rcpp_result_gen;
END_RCPP
}
// cfconstr0_pGMM
Rcpp::List cfconstr0_pGMM(arma::mat& x, arma::rowvec prop, arma::mat mu, arma::cube Sigma, arma::mat combos, int k, arma::rowvec df, int lambda, int citermax, double tol);
RcppExport SEXP _pGMCM_cfconstr0_pGMM(SEXP xSEXP, SEXP propSEXP, SEXP muSEXP, SEXP SigmaSEXP, SEXP combosSEXP, SEXP kSEXP, SEXP dfSEXP, SEXP lambdaSEXP, SEXP citermaxSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::rowvec >::type prop(propSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type Sigma(SigmaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type combos(combosSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< arma::rowvec >::type df(dfSEXP);
    Rcpp::traits::input_parameter< int >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< int >::type citermax(citermaxSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(cfconstr0_pGMM(x, prop, mu, Sigma, combos, k, df, lambda, citermax, tol));
    return rcpp_result_gen;
END_RCPP
}
// cduvnorm
arma::colvec cduvnorm(arma::colvec x, double mu, double sigma);
RcppExport SEXP _pGMCM_cduvnorm(SEXP xSEXP, SEXP muSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::colvec >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(cduvnorm(x, mu, sigma));
    return rcpp_result_gen;
END_RCPP
}
// cmarg_ll_gmm
double cmarg_ll_gmm(arma::mat& z, arma::mat mu, arma::mat sigma, arma::rowvec prop, arma::mat combos, int k);
RcppExport SEXP _pGMCM_cmarg_ll_gmm(SEXP zSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP propSEXP, SEXP combosSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< arma::rowvec >::type prop(propSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type combos(combosSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(cmarg_ll_gmm(z, mu, sigma, prop, combos, k));
    return rcpp_result_gen;
END_RCPP
}
// cll_gmm
double cll_gmm(arma::mat& z, arma::mat mu, arma::mat sigma, double rho, arma::rowvec prop, arma::mat combos, int k);
RcppExport SEXP _pGMCM_cll_gmm(SEXP zSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP rhoSEXP, SEXP propSEXP, SEXP combosSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< arma::rowvec >::type prop(propSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type combos(combosSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(cll_gmm(z, mu, sigma, rho, prop, combos, k));
    return rcpp_result_gen;
END_RCPP
}
// cmarg0_ll_gmm
double cmarg0_ll_gmm(arma::mat& z, arma::mat mu, arma::cube Sigma, arma::rowvec prop, int k);
RcppExport SEXP _pGMCM_cmarg0_ll_gmm(SEXP zSEXP, SEXP muSEXP, SEXP SigmaSEXP, SEXP propSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type Sigma(SigmaSEXP);
    Rcpp::traits::input_parameter< arma::rowvec >::type prop(propSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(cmarg0_ll_gmm(z, mu, Sigma, prop, k));
    return rcpp_result_gen;
END_RCPP
}
// cll0_gmm
double cll0_gmm(arma::mat& z, arma::mat mu, arma::cube Sigma, arma::rowvec prop, int k);
RcppExport SEXP _pGMCM_cll0_gmm(SEXP zSEXP, SEXP muSEXP, SEXP SigmaSEXP, SEXP propSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type Sigma(SigmaSEXP);
    Rcpp::traits::input_parameter< arma::rowvec >::type prop(propSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(cll0_gmm(z, mu, Sigma, prop, k));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_pGMCM_cgetPaths", (DL_FUNC) &_pGMCM_cgetPaths, 7},
    {"_pGMCM_crowMatch", (DL_FUNC) &_pGMCM_crowMatch, 2},
    {"_pGMCM_cget_prior_count", (DL_FUNC) &_pGMCM_cget_prior_count, 6},
    {"_pGMCM_cget_true_assoc_idx", (DL_FUNC) &_pGMCM_cget_true_assoc_idx, 2},
    {"_pGMCM_trans_func", (DL_FUNC) &_pGMCM_trans_func, 1},
    {"_pGMCM_get_list_names", (DL_FUNC) &_pGMCM_get_list_names, 1},
    {"_pGMCM_cpaste0", (DL_FUNC) &_pGMCM_cpaste0, 1},
    {"_pGMCM_cstr_split", (DL_FUNC) &_pGMCM_cstr_split, 2},
    {"_pGMCM_caccept", (DL_FUNC) &_pGMCM_caccept, 2},
    {"_pGMCM_caccept2", (DL_FUNC) &_pGMCM_caccept2, 3},
    {"_pGMCM_cprune_path", (DL_FUNC) &_pGMCM_cprune_path, 2},
    {"_pGMCM_cassociate", (DL_FUNC) &_pGMCM_cassociate, 3},
    {"_pGMCM_cprune_path2", (DL_FUNC) &_pGMCM_cprune_path2, 6},
    {"_pGMCM_abs3", (DL_FUNC) &_pGMCM_abs3, 1},
    {"_pGMCM_SCAD_1d", (DL_FUNC) &_pGMCM_SCAD_1d, 4},
    {"_pGMCM_double_SCAD_1d", (DL_FUNC) &_pGMCM_double_SCAD_1d, 3},
    {"_pGMCM_SCAD", (DL_FUNC) &_pGMCM_SCAD, 4},
    {"_pGMCM_double_SCAD", (DL_FUNC) &_pGMCM_double_SCAD, 3},
    {"_pGMCM_choose_slice", (DL_FUNC) &_pGMCM_choose_slice, 4},
    {"_pGMCM_Mahalanobis", (DL_FUNC) &_pGMCM_Mahalanobis, 3},
    {"_pGMCM_cdmvnorm", (DL_FUNC) &_pGMCM_cdmvnorm, 3},
    {"_pGMCM_cfpGMM", (DL_FUNC) &_pGMCM_cfpGMM, 9},
    {"_pGMCM_cget_constr_sigma", (DL_FUNC) &_pGMCM_cget_constr_sigma, 4},
    {"_pGMCM_trans_rho", (DL_FUNC) &_pGMCM_trans_rho, 1},
    {"_pGMCM_trans_rho_inv", (DL_FUNC) &_pGMCM_trans_rho_inv, 1},
    {"_pGMCM_func_to_optim", (DL_FUNC) &_pGMCM_func_to_optim, 4},
    {"_pGMCM_optim_rcpp", (DL_FUNC) &_pGMCM_optim_rcpp, 4},
    {"_pGMCM_cfconstr_pGMM", (DL_FUNC) &_pGMCM_cfconstr_pGMM, 12},
    {"_pGMCM_bind_diags", (DL_FUNC) &_pGMCM_bind_diags, 1},
    {"_pGMCM_bind_offdiags", (DL_FUNC) &_pGMCM_bind_offdiags, 1},
    {"_pGMCM_get_sigma_optim", (DL_FUNC) &_pGMCM_get_sigma_optim, 2},
    {"_pGMCM_get_mu_optim", (DL_FUNC) &_pGMCM_get_mu_optim, 2},
    {"_pGMCM_get_rho_optim", (DL_FUNC) &_pGMCM_get_rho_optim, 2},
    {"_pGMCM_func_to_optim0", (DL_FUNC) &_pGMCM_func_to_optim0, 8},
    {"_pGMCM_optim0_rcpp", (DL_FUNC) &_pGMCM_optim0_rcpp, 8},
    {"_pGMCM_cfconstr0_pGMM", (DL_FUNC) &_pGMCM_cfconstr0_pGMM, 10},
    {"_pGMCM_cduvnorm", (DL_FUNC) &_pGMCM_cduvnorm, 3},
    {"_pGMCM_cmarg_ll_gmm", (DL_FUNC) &_pGMCM_cmarg_ll_gmm, 6},
    {"_pGMCM_cll_gmm", (DL_FUNC) &_pGMCM_cll_gmm, 7},
    {"_pGMCM_cmarg0_ll_gmm", (DL_FUNC) &_pGMCM_cmarg0_ll_gmm, 5},
    {"_pGMCM_cll0_gmm", (DL_FUNC) &_pGMCM_cll0_gmm, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_pGMCM(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
