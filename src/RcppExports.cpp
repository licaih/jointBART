// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// JointBart
void JointBart(const IntegerVector& n, const size_t& p, const IntegerVector& np, const List& x, const List& y, const List& xp, const size_t& m, IntegerVector& nc, const size_t& nd, const size_t& burn, const double& mybeta, const double& alpha, const NumericVector& tau, const NumericVector& nu, const NumericVector& lambda, const NumericVector& sigest, const List& w, const bool& dart, const double& theta, const double& omega, IntegerVector& igrp, const double& a, const double& b, const double& rho, const bool& aug, const size_t& keeptrain, const size_t& keeptest, const size_t& keeptestme, const size_t& keeptreedraws, const size_t& printevery, const List& iXinfo);
RcppExport SEXP _jointBART_JointBart(SEXP nSEXP, SEXP pSEXP, SEXP npSEXP, SEXP xSEXP, SEXP ySEXP, SEXP xpSEXP, SEXP mSEXP, SEXP ncSEXP, SEXP ndSEXP, SEXP burnSEXP, SEXP mybetaSEXP, SEXP alphaSEXP, SEXP tauSEXP, SEXP nuSEXP, SEXP lambdaSEXP, SEXP sigestSEXP, SEXP wSEXP, SEXP dartSEXP, SEXP thetaSEXP, SEXP omegaSEXP, SEXP igrpSEXP, SEXP aSEXP, SEXP bSEXP, SEXP rhoSEXP, SEXP augSEXP, SEXP keeptrainSEXP, SEXP keeptestSEXP, SEXP keeptestmeSEXP, SEXP keeptreedrawsSEXP, SEXP printeverySEXP, SEXP iXinfoSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerVector& >::type n(nSEXP);
    Rcpp::traits::input_parameter< const size_t& >::type p(pSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type np(npSEXP);
    Rcpp::traits::input_parameter< const List& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const List& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const List& >::type xp(xpSEXP);
    Rcpp::traits::input_parameter< const size_t& >::type m(mSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type nc(ncSEXP);
    Rcpp::traits::input_parameter< const size_t& >::type nd(ndSEXP);
    Rcpp::traits::input_parameter< const size_t& >::type burn(burnSEXP);
    Rcpp::traits::input_parameter< const double& >::type mybeta(mybetaSEXP);
    Rcpp::traits::input_parameter< const double& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type sigest(sigestSEXP);
    Rcpp::traits::input_parameter< const List& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const bool& >::type dart(dartSEXP);
    Rcpp::traits::input_parameter< const double& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const double& >::type omega(omegaSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type igrp(igrpSEXP);
    Rcpp::traits::input_parameter< const double& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const double& >::type b(bSEXP);
    Rcpp::traits::input_parameter< const double& >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< const bool& >::type aug(augSEXP);
    Rcpp::traits::input_parameter< const size_t& >::type keeptrain(keeptrainSEXP);
    Rcpp::traits::input_parameter< const size_t& >::type keeptest(keeptestSEXP);
    Rcpp::traits::input_parameter< const size_t& >::type keeptestme(keeptestmeSEXP);
    Rcpp::traits::input_parameter< const size_t& >::type keeptreedraws(keeptreedrawsSEXP);
    Rcpp::traits::input_parameter< const size_t& >::type printevery(printeverySEXP);
    Rcpp::traits::input_parameter< const List& >::type iXinfo(iXinfoSEXP);
    JointBart(n, p, np, x, y, xp, m, nc, nd, burn, mybeta, alpha, tau, nu, lambda, sigest, w, dart, theta, omega, igrp, a, b, rho, aug, keeptrain, keeptest, keeptestme, keeptreedraws, printevery, iXinfo);
    return R_NilValue;
END_RCPP
}
// rcpp_hello
List rcpp_hello();
RcppExport SEXP _jointBART_rcpp_hello() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(rcpp_hello());
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP cwbart(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_jointBART_JointBart", (DL_FUNC) &_jointBART_JointBart, 31},
    {"_jointBART_rcpp_hello", (DL_FUNC) &_jointBART_rcpp_hello, 0},
    {"cwbart", (DL_FUNC) &cwbart, 31},
    {NULL, NULL, 0}
};

RcppExport void R_init_jointBART(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}