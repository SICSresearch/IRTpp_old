#include <Rcpp.h>
#include <string.h>

#include <util_package.hpp>

Rcpp::List irtpp_file(std::string dat, int e_model, Rcpp::NumericMatrix quads,
                           Rcpp::NumericMatrix init_val, bool init_val_flag)
{
  PatternMatrix *datSet = getPatternMatrix(dat);
  
  Rcpp::List result = irtpp_aux(datSet, e_model, quads, init_val, init_val_flag);
  
  delete datSet;

  return result;
}

Rcpp::List irtpp_r(Rcpp::NumericMatrix dat, int e_model, Rcpp::NumericMatrix quads,
                        Rcpp::NumericMatrix init_val, bool init_val_flag)
{
  PatternMatrix *datSet = getPatternMatrix(dat);

  Rcpp::List result = irtpp_aux(datSet, e_model, quads, init_val, init_val_flag);

  delete datSet;

  return result;
}

//[[Rcpp::export]]
Rcpp::List irtppinterfacevalues(Rcpp::NumericMatrix dat, int e_model, Rcpp::NumericMatrix quads,
                                Rcpp::NumericMatrix init_val)
{ return irtpp_r(dat, e_model, quads, init_val, true); }

//[[Rcpp::export]]
Rcpp::List irtppinterface(Rcpp::NumericMatrix dat, int e_model, Rcpp::NumericMatrix quads)
{
  Rcpp::NumericMatrix init_val(1,1);
  return irtpp_r(dat, e_model, quads, init_val, false);
}

//[[Rcpp::export]]
Rcpp::List irtppinterfacefilevalues(std::string dat, int e_model, Rcpp::NumericMatrix quads,
                                    Rcpp::NumericMatrix init_val)
{ return irtpp_file(dat, e_model, quads, init_val, true); }

//[[Rcpp::export]]
Rcpp::List irtppinterfacefile(std::string dat, int e_model, Rcpp::NumericMatrix quads)
{
  Rcpp::NumericMatrix init_val(1,1);
  return irtpp_file(dat, e_model, quads, init_val, false);
}

//[[Rcpp::export]]
Rcpp::List mapinterfacefile(Rcpp::NumericMatrix zita_par, std::string dat,
                                int e_model, Rcpp::NumericMatrix quads)
{
  PatternMatrix *datSet = getPatternMatrix(dat);

  Rcpp::List result = abilityinterface(zita_par, datSet, e_model, quads, 1);

  //delete datSet;

  return result;
}

//[[Rcpp::export]]
Rcpp::List mapinterface(Rcpp::NumericMatrix zita_par, Rcpp::NumericMatrix dat,
                        int e_model, Rcpp::NumericMatrix quads)
{
  PatternMatrix *datSet = getPatternMatrix(dat);

  Rcpp::List result = abilityinterface(zita_par, datSet, e_model, quads, 1);

  //delete datSet;

  return result;
}

//[[Rcpp::export]]
Rcpp::List eapinterfacefile(Rcpp::NumericMatrix zita_par, std::string dat,
                            int e_model, Rcpp::NumericMatrix quads)
{
  PatternMatrix *datSet = getPatternMatrix(dat);

  Rcpp::List result = abilityinterface(zita_par, datSet, e_model, quads, 0);

  //delete datSet;

  return result;
}

//[[Rcpp::export]]
Rcpp::List eapinterface(Rcpp::NumericMatrix zita_par, Rcpp::NumericMatrix dat,
                        int e_model, Rcpp::NumericMatrix quads)
{
  PatternMatrix *datSet = getPatternMatrix(dat);

  Rcpp::List result = abilityinterface(zita_par, datSet, e_model, quads, 0);

  //delete datSet;

  return result;
}