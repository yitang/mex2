// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// rcpp_hello_world
List rcpp_hello_world();
RcppExport SEXP HTModel_rcpp_hello_world() {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        List __result = rcpp_hello_world();
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// obj_cpp
double obj_cpp(NumericVector para, NumericVector yi, NumericVector yj);
RcppExport SEXP HTModel_obj_cpp(SEXP paraSEXP, SEXP yiSEXP, SEXP yjSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericVector >::type para(paraSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type yi(yiSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type yj(yjSEXP );
        double __result = obj_cpp(para, yi, yj);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
