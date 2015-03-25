// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// irtpp
List irtpp(IntegerMatrix data, CharacterVector nameOfModel, IntegerVector dim, CharacterVector nameOfInitVal, NumericVector vEpsilonConv, IntegerVector maxIt, LogicalVector vVerbose);
RcppExport SEXP IRTpp_irtpp(SEXP dataSEXP, SEXP nameOfModelSEXP, SEXP dimSEXP, SEXP nameOfInitValSEXP, SEXP vEpsilonConvSEXP, SEXP maxItSEXP, SEXP vVerboseSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< IntegerMatrix >::type data(dataSEXP );
        Rcpp::traits::input_parameter< CharacterVector >::type nameOfModel(nameOfModelSEXP );
        Rcpp::traits::input_parameter< IntegerVector >::type dim(dimSEXP );
        Rcpp::traits::input_parameter< CharacterVector >::type nameOfInitVal(nameOfInitValSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type vEpsilonConv(vEpsilonConvSEXP );
        Rcpp::traits::input_parameter< IntegerVector >::type maxIt(maxItSEXP );
        Rcpp::traits::input_parameter< LogicalVector >::type vVerbose(vVerboseSEXP );
        List __result = irtpp(data, nameOfModel, dim, nameOfInitVal, vEpsilonConv, maxIt, vVerbose);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// itemfitpp
List itemfitpp(NumericVector parametersO, IntegerMatrix data, CharacterVector nameOfModel, LogicalVector vVerbose);
RcppExport SEXP IRTpp_itemfitpp(SEXP parametersOSEXP, SEXP dataSEXP, SEXP nameOfModelSEXP, SEXP vVerboseSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericVector >::type parametersO(parametersOSEXP );
        Rcpp::traits::input_parameter< IntegerMatrix >::type data(dataSEXP );
        Rcpp::traits::input_parameter< CharacterVector >::type nameOfModel(nameOfModelSEXP );
        Rcpp::traits::input_parameter< LogicalVector >::type vVerbose(vVerboseSEXP );
        List __result = itemfitpp(parametersO, data, nameOfModel, vVerbose);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
