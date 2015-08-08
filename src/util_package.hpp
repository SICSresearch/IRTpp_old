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
#include <Rcpp.h>
#include <string.h>

using namespace std;

PatternMatrix * getPatternMatrix(string r_path);

PatternMatrix * getPatternMatrix(Rcpp::NumericMatrix r_dataSet);

Rcpp::List irtpp_aux(PatternMatrix *datSet, int e_model, Rcpp::NumericMatrix quads,
                     Rcpp::NumericMatrix init_val, bool init_val_flag);

Rcpp::List abilityinterface(Rcpp::NumericMatrix zita_par, PatternMatrix * datSet,
                            int e_model, Rcpp::NumericMatrix quads, int method);

#endif