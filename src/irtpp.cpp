#include <Rcpp.h>
#include <string.h>

#include <util_package.hpp>

// Pasar todos los argumentos como una lista
// controlar todo con la menor cantidad de funciones posibles

Rcpp::List irtpp_file(std::string dat, int e_model, Rcpp::NumericMatrix quads,
                      Rcpp::NumericMatrix init_val, bool init_val_flag,
                      bool to_flag_file, std::string output_path)
{
  PatternMatrix *datSet = getPatternMatrix(dat);
  
  Rcpp::List result = irtpp_aux(datSet, e_model, quads, init_val, init_val_flag, to_flag_file, output_path);
  
  delete datSet;

  return result;
}

Rcpp::List irtpp_r(Rcpp::NumericMatrix dat, int e_model, Rcpp::NumericMatrix quads,
                   Rcpp::NumericMatrix init_val, bool init_val_flag,
                   bool to_flag_file, std::string output_path)
{
  PatternMatrix *datSet = getPatternMatrix(dat);

  Rcpp::List result = irtpp_aux(datSet, e_model, quads, init_val, init_val_flag, to_flag_file, output_path);

  delete datSet;

  return result;
}

//[[Rcpp::export]]
Rcpp::List irtppinterfacevalues(Rcpp::NumericMatrix dat, int e_model, Rcpp::NumericMatrix quads,
                                Rcpp::NumericMatrix init_val, bool to_flag_file, std::string output_path)
{ return irtpp_r(dat, e_model, quads, init_val, true, to_flag_file, output_path); }

//[[Rcpp::export]]
Rcpp::List irtppinterface(Rcpp::NumericMatrix dat, int e_model, Rcpp::NumericMatrix quads,
                          bool to_flag_file, std::string output_path)
{
  Rcpp::NumericMatrix init_val(1,1);
  return irtpp_r(dat, e_model, quads, init_val, false, to_flag_file, output_path);
}

//[[Rcpp::export]]
Rcpp::List irtppinterfacefilevalues(std::string dat, int e_model, Rcpp::NumericMatrix quads,
                                    Rcpp::NumericMatrix init_val, bool to_flag_file, std::string output_path)
{ return irtpp_file(dat, e_model, quads, init_val, true, to_flag_file, output_path); }

//[[Rcpp::export]]
Rcpp::List irtppinterfacefile(std::string dat, int e_model, Rcpp::NumericMatrix quads,
                              bool to_flag_file, std::string output_path)
{
  Rcpp::NumericMatrix init_val(1,1);
  return irtpp_file(dat, e_model, quads, init_val, false, to_flag_file, output_path);
}

//[[Rcpp::export]]
Rcpp::List mapinterfacefile(Rcpp::NumericMatrix zita_par, std::string dat,
                            int e_model, Rcpp::NumericMatrix quads,
                            bool to_flag_file, std::string output_path)
{
  PatternMatrix *datSet = getPatternMatrix(dat);

  Rcpp::List result = abilityinterface(zita_par, datSet, e_model, quads, 1, to_flag_file, output_path);

  //delete datSet;

  return result;
}

//[[Rcpp::export]]
Rcpp::List mapinterface(Rcpp::NumericMatrix zita_par, Rcpp::NumericMatrix dat,
                        int e_model, Rcpp::NumericMatrix quads,
                        bool to_flag_file, std::string output_path)
{
  PatternMatrix *datSet = getPatternMatrix(dat);

  Rcpp::List result = abilityinterface(zita_par, datSet, e_model, quads, 1, to_flag_file, output_path);

  //delete datSet;

  return result;
}

//[[Rcpp::export]]
Rcpp::List eapinterfacefile(Rcpp::NumericMatrix zita_par, std::string dat,
                            int e_model, Rcpp::NumericMatrix quads,
                            bool to_flag_file, std::string output_path)
{
  PatternMatrix *datSet = getPatternMatrix(dat);

  Rcpp::List result = abilityinterface(zita_par, datSet, e_model, quads, 0, to_flag_file, output_path);

  //delete datSet;

  return result;
}

//[[Rcpp::export]]
Rcpp::List eapinterface(Rcpp::NumericMatrix zita_par, Rcpp::NumericMatrix dat,
                        int e_model, Rcpp::NumericMatrix quads,
                        bool to_flag_file, std::string output_path)
{
  PatternMatrix *datSet = getPatternMatrix(dat);

  Rcpp::List result = abilityinterface(zita_par, datSet, e_model, quads, 0, to_flag_file, output_path);

  //delete datSet;

  return result;
}