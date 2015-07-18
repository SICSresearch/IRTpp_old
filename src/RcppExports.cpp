// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// irtppinterfacevalues
Rcpp::List irtppinterfacevalues(Rcpp::NumericMatrix dat, int e_model, Rcpp::NumericMatrix quads, Rcpp::NumericMatrix init_val);
RcppExport SEXP IRTpp_irtppinterfacevalues(SEXP datSEXP, SEXP e_modelSEXP, SEXP quadsSEXP, SEXP init_valSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type dat(datSEXP);
    Rcpp::traits::input_parameter< int >::type e_model(e_modelSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type quads(quadsSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type init_val(init_valSEXP);
    __result = Rcpp::wrap(irtppinterfacevalues(dat, e_model, quads, init_val));
    return __result;
END_RCPP
}
// irtppinterface
Rcpp::List irtppinterface(Rcpp::NumericMatrix dat, int e_model, Rcpp::NumericMatrix quads);
RcppExport SEXP IRTpp_irtppinterface(SEXP datSEXP, SEXP e_modelSEXP, SEXP quadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type dat(datSEXP);
    Rcpp::traits::input_parameter< int >::type e_model(e_modelSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type quads(quadsSEXP);
    __result = Rcpp::wrap(irtppinterface(dat, e_model, quads));
    return __result;
END_RCPP
}
// eapinterface
Rcpp::List eapinterface(Rcpp::NumericMatrix zita_par, Rcpp::NumericMatrix dat, int e_model, Rcpp::NumericMatrix quads);
RcppExport SEXP IRTpp_eapinterface(SEXP zita_parSEXP, SEXP datSEXP, SEXP e_modelSEXP, SEXP quadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type zita_par(zita_parSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type dat(datSEXP);
    Rcpp::traits::input_parameter< int >::type e_model(e_modelSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type quads(quadsSEXP);
    __result = Rcpp::wrap(eapinterface(zita_par, dat, e_model, quads));
    return __result;
END_RCPP
}
// mapinterface
Rcpp::List mapinterface(Rcpp::NumericMatrix zita_par, Rcpp::NumericMatrix dat, int e_model, Rcpp::NumericMatrix quads);
RcppExport SEXP IRTpp_mapinterface(SEXP zita_parSEXP, SEXP datSEXP, SEXP e_modelSEXP, SEXP quadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type zita_par(zita_parSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type dat(datSEXP);
    Rcpp::traits::input_parameter< int >::type e_model(e_modelSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type quads(quadsSEXP);
    __result = Rcpp::wrap(mapinterface(zita_par, dat, e_model, quads));
    return __result;
END_RCPP
}
