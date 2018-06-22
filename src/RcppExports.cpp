// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

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
// teststuff
arma::mat teststuff(arma::mat combos);
RcppExport SEXP _pGMCM_teststuff(SEXP combosSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type combos(combosSEXP);
    rcpp_result_gen = Rcpp::wrap(teststuff(combos));
    return rcpp_result_gen;
END_RCPP
}
// cfconstr_pGMM
Rcpp::List cfconstr_pGMM(arma::mat& x, arma::rowvec prop, arma::mat mu, arma::mat sigma, double rho, arma::mat combos, int k, arma::rowvec df, int lambda, int citermax, double tol);
RcppExport SEXP _pGMCM_cfconstr_pGMM(SEXP xSEXP, SEXP propSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP rhoSEXP, SEXP combosSEXP, SEXP kSEXP, SEXP dfSEXP, SEXP lambdaSEXP, SEXP citermaxSEXP, SEXP tolSEXP) {
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
    rcpp_result_gen = Rcpp::wrap(cfconstr_pGMM(x, prop, mu, sigma, rho, combos, k, df, lambda, citermax, tol));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_pGMCM_choose_slice", (DL_FUNC) &_pGMCM_choose_slice, 4},
    {"_pGMCM_Mahalanobis", (DL_FUNC) &_pGMCM_Mahalanobis, 3},
    {"_pGMCM_cdmvnorm", (DL_FUNC) &_pGMCM_cdmvnorm, 3},
    {"_pGMCM_cfpGMM", (DL_FUNC) &_pGMCM_cfpGMM, 9},
    {"_pGMCM_cget_constr_sigma", (DL_FUNC) &_pGMCM_cget_constr_sigma, 4},
    {"_pGMCM_trans_rho", (DL_FUNC) &_pGMCM_trans_rho, 1},
    {"_pGMCM_trans_rho_inv", (DL_FUNC) &_pGMCM_trans_rho_inv, 1},
    {"_pGMCM_func_to_optim", (DL_FUNC) &_pGMCM_func_to_optim, 4},
    {"_pGMCM_optim_rcpp", (DL_FUNC) &_pGMCM_optim_rcpp, 4},
    {"_pGMCM_teststuff", (DL_FUNC) &_pGMCM_teststuff, 1},
    {"_pGMCM_cfconstr_pGMM", (DL_FUNC) &_pGMCM_cfconstr_pGMM, 11},
    {NULL, NULL, 0}
};

RcppExport void R_init_pGMCM(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
