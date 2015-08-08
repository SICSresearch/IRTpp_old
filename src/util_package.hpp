#ifndef UTILP_HPP
#define UTILP_HPP

#include <estimation/bayesian/LatentTraitEstimation.h>
#include <estimation/classical/EMEstimation.h>
#include <model/SICSGeneralModel.h>
#include <model/ModelFactory.h>
#include <type/PatternMatrix.h>
#include <type/LatentTraits.h>
#include <type/Matrix.h>
#include <input/Input.h>
#include <model/Model.h>
#include <string.h>
#include <fstream>
#include <Rcpp.h>

using namespace std;

PatternMatrix * getPatternMatrix(string r_path);

Rcpp::List transformParameterOutput(void *);

Rcpp::List transformAbilityOutput(void *);

PatternMatrix * getPatternMatrix(Rcpp::NumericMatrix r_dataSet);

Rcpp::List irtpp_aux(PatternMatrix *datSet, int e_model, Rcpp::NumericMatrix quads,
                     Rcpp::NumericMatrix init_val, bool init_val_flag,
                     bool to_flag_file, string output_path);

Rcpp::List abilityinterface(Rcpp::NumericMatrix zita_par, PatternMatrix * datSet,
                            int e_model, Rcpp::NumericMatrix quads, int method,
                            bool to_flag_file, string output_path);

#endif