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
// teststuff
arma::colvec teststuff(arma::mat mu_old, arma::cube Sigma_old, arma::mat combos);
RcppExport SEXP _pGMCM_teststuff(SEXP mu_oldSEXP, SEXP Sigma_oldSEXP, SEXP combosSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type mu_old(mu_oldSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type Sigma_old(Sigma_oldSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type combos(combosSEXP);
    rcpp_result_gen = Rcpp::wrap(teststuff(mu_old, Sigma_old, combos));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_pGMCM_choose_slice", (DL_FUNC) &_pGMCM_choose_slice, 4},
    {"_pGMCM_Mahalanobis", (DL_FUNC) &_pGMCM_Mahalanobis, 3},
    {"_pGMCM_cdmvnorm", (DL_FUNC) &_pGMCM_cdmvnorm, 3},
    {"_pGMCM_cfpGMM", (DL_FUNC) &_pGMCM_cfpGMM, 9},
    {"_pGMCM_abs3", (DL_FUNC) &_pGMCM_abs3, 1},
    {"_pGMCM_cget_constr_sigma", (DL_FUNC) &_pGMCM_cget_constr_sigma, 4},
    {"_pGMCM_trans_rho", (DL_FUNC) &_pGMCM_trans_rho, 1},
    {"_pGMCM_trans_rho_inv", (DL_FUNC) &_pGMCM_trans_rho_inv, 1},
    {"_pGMCM_func_to_optim", (DL_FUNC) &_pGMCM_func_to_optim, 4},
    {"_pGMCM_optim_rcpp", (DL_FUNC) &_pGMCM_optim_rcpp, 4},
    {"_pGMCM_cfconstr_pGMM", (DL_FUNC) &_pGMCM_cfconstr_pGMM, 11},
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
    {"_pGMCM_teststuff", (DL_FUNC) &_pGMCM_teststuff, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_pGMCM(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
