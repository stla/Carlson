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
// Carlson_RJ_
Rcomplex Carlson_RJ_(Rcomplex xr, Rcomplex yr, Rcomplex zr, Rcomplex pr, double err);
RcppExport SEXP _Carlson_Carlson_RJ_(SEXP xrSEXP, SEXP yrSEXP, SEXP zrSEXP, SEXP prSEXP, SEXP errSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcomplex >::type xr(xrSEXP);
    Rcpp::traits::input_parameter< Rcomplex >::type yr(yrSEXP);
    Rcpp::traits::input_parameter< Rcomplex >::type zr(zrSEXP);
    Rcpp::traits::input_parameter< Rcomplex >::type pr(prSEXP);
    Rcpp::traits::input_parameter< double >::type err(errSEXP);
    rcpp_result_gen = Rcpp::wrap(Carlson_RJ_(xr, yr, zr, pr, err));
    return rcpp_result_gen;
END_RCPP
}
// ellEcpp
Rcpp::ComplexVector ellEcpp(Rcpp::ComplexVector phi_, Rcpp::ComplexVector m_, double err);
RcppExport SEXP _Carlson_ellEcpp(SEXP phi_SEXP, SEXP m_SEXP, SEXP errSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::ComplexVector >::type phi_(phi_SEXP);
    Rcpp::traits::input_parameter< Rcpp::ComplexVector >::type m_(m_SEXP);
    Rcpp::traits::input_parameter< double >::type err(errSEXP);
    rcpp_result_gen = Rcpp::wrap(ellEcpp(phi_, m_, err));
    return rcpp_result_gen;
END_RCPP
}
// ellFcpp
Rcpp::ComplexVector ellFcpp(Rcpp::ComplexVector phi_, Rcpp::ComplexVector m_, double err);
RcppExport SEXP _Carlson_ellFcpp(SEXP phi_SEXP, SEXP m_SEXP, SEXP errSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::ComplexVector >::type phi_(phi_SEXP);
    Rcpp::traits::input_parameter< Rcpp::ComplexVector >::type m_(m_SEXP);
    Rcpp::traits::input_parameter< double >::type err(errSEXP);
    rcpp_result_gen = Rcpp::wrap(ellFcpp(phi_, m_, err));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_Carlson_Carlson_RF_", (DL_FUNC) &_Carlson_Carlson_RF_, 4},
    {"_Carlson_Carlson_RD_", (DL_FUNC) &_Carlson_Carlson_RD_, 4},
    {"_Carlson_Carlson_RJ_", (DL_FUNC) &_Carlson_Carlson_RJ_, 5},
    {"_Carlson_ellEcpp", (DL_FUNC) &_Carlson_ellEcpp, 3},
    {"_Carlson_ellFcpp", (DL_FUNC) &_Carlson_ellFcpp, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_Carlson(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
