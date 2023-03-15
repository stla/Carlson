// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// Carlson_RF_
Rcomplex Carlson_RF_(Rcomplex xr, Rcomplex yr, Rcomplex zr, double err);
RcppExport SEXP _Carlson_Carlson_RF_(SEXP xrSEXP, SEXP yrSEXP, SEXP zrSEXP, SEXP errSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcomplex >::type xr(xrSEXP);
    Rcpp::traits::input_parameter< Rcomplex >::type yr(yrSEXP);
    Rcpp::traits::input_parameter< Rcomplex >::type zr(zrSEXP);
    Rcpp::traits::input_parameter< double >::type err(errSEXP);
    rcpp_result_gen = Rcpp::wrap(Carlson_RF_(xr, yr, zr, err));
    return rcpp_result_gen;
END_RCPP
}
// Carlson_RD_
Rcomplex Carlson_RD_(Rcomplex xr, Rcomplex yr, Rcomplex zr, double err);
RcppExport SEXP _Carlson_Carlson_RD_(SEXP xrSEXP, SEXP yrSEXP, SEXP zrSEXP, SEXP errSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcomplex >::type xr(xrSEXP);
    Rcpp::traits::input_parameter< Rcomplex >::type yr(yrSEXP);
    Rcpp::traits::input_parameter< Rcomplex >::type zr(zrSEXP);
    Rcpp::traits::input_parameter< double >::type err(errSEXP);
    rcpp_result_gen = Rcpp::wrap(Carlson_RD_(xr, yr, zr, err));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_Carlson_Carlson_RF_", (DL_FUNC) &_Carlson_Carlson_RF_, 4},
    {"_Carlson_Carlson_RD_", (DL_FUNC) &_Carlson_Carlson_RD_, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_Carlson(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
