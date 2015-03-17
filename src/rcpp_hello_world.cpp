#include <Rcpp.h>
#include <string.h>

#include <model/parameter/OnePLACModel.h>
#include <model/parameter/OnePLACModel.cpp>
#include <model/parameter/OnePLModel.h>
#include <model/parameter/OnePLModel.cpp>
#include <model/parameter/ParameterModel.h>
#include <model/parameter/ParameterModel.cpp>
#include <model/parameter/RaschModel.h>
#include <model/parameter/RaschModel.cpp>
#include <model/parameter/TwoPLModel.h>
#include <model/parameter/TwoPLModel.cpp>
#include <model/parameter/ThreePLModel.h>
#include <model/parameter/ThreePLModel.cpp>

//util
#include <util/asa111.hpp>
#include <util/asa111.cpp>
#include <util/util.h>

//type
#include <type/Constant.h>
#include <type/Constant.cpp>
#include <type/DataSet.h>
#include <type/DataSet.cpp>
#include <type/Matrix.h>
#include <type/PatternMatrix.h>
#include <type/PatternMatrix.cpp>
#include <type/QuadratureNodes.h>
#include <type/QuadratureNodes.cpp>

//trace
#include <trace/Timer.h>
#include <trace/Trace.h>

//output
#include <output/Output.h>
#include <output/Output.cpp> 

//optimizer
#include <optimizer/BFGSOptimizer.h>
#include <optimizer/Optimizer.h>
#include <optimizer/Optimizer.cpp>

//input
#include <input/Input.h>
#include <input/Input.cpp>

//estimation
#include <estimation/Estimation.h>
#include <estimation/Estimation.cpp>
#include <estimation/classical/ClassicalEstimation.h>
#include <estimation/classical/ClassicalEstimation.cpp>
#include <estimation/classical/EMEstimation.h>

#include <estimation/classical/EMEstimation.cpp>
#include <estimation/classical/EMEstimators/EM1PL.h>
#include <estimation/classical/EMEstimators/EM1PLAC.h>
#include <estimation/classical/EMEstimators/EM2PL.h>
#include <estimation/classical/EMEstimators/EM3PL.h>
#include <estimation/classical/EMEstimators/EMEstimator.h>


//model
#include <model/Model.h>
#include <model/Model.cpp>
#include <model/ModelFactory.h>
#include <model/ModelFactory.cpp>
#include <model/SICSGeneralModel.h>
#include <model/SICSGeneralModel.cpp>
#include <model/dimension/DimensionModel.h>
#include <model/dimension/DimensionModel.cpp>
#include <model/dimension/MultidimensionalModel.h>
#include <model/dimension/MultidimensionalModel.cpp>
#include <model/dimension/MultiUniDimModel.h>
#include <model/dimension/MultiUniDimModel.cpp>
#include <model/dimension/UnidimensionalModel.h>
#include <model/dimension/UnidimensionalModel.cpp>
#include <model/item/DichotomousModel.h>
#include <model/item/DichotomousModel.cpp>
#include <model/item/ItemModel.h>
#include <model/item/ItemModel.cpp>
#include <model/item/PolytomousModel.h>
#include <model/item/PolytomousModel.cpp>

#include "interface.h"

using namespace Rcpp;

// [[Rcpp::export]]
List irtpp( IntegerMatrix data,  CharacterVector nameOfModel, IntegerVector dim, CharacterVector nameOfInitVal, NumericVector vEpsilonConv, IntegerVector maxIt, LogicalVector vVerbose) {
    Rcout<<"model "<< ", Column = "<<data.ncol()<<", Row = "<<data.nrow()<<endl;
    Rcout<<"name of Model: "<<nameOfModel[0]<<endl;
    Rcout<<"number of dimensions: ";
    Rcout<<dim[0]<<endl;
    Rcout<<"name of initials values: ";
    Rcout<<nameOfInitVal[0]<<endl;
    Rcout<<"Epsilon for convergence: ";
    Rcout<<vEpsilonConv[0]<<endl;
    Rcout<<"maximum number of iterations: ";
    Rcout<<maxIt[0]<<endl;
    Rcout<<"value of verbose: ";
    Rcout<<vVerbose[0]<<endl;
    int nColumn = data.ncol(), nRow = data.nrow(); 
    int **DataI, nuM, nPar;
    DataI = new int *[nRow];
  	for ( int i = 0; i < nRow; i++ ) DataI[i] = new int[nColumn];
  	char *model, *initValues;    
  	for ( int j = 0; j<  nColumn; j++ ) for ( int i = 0; i < nRow; i++ ) DataI[i][j] = data[i+j*nRow];
  	model = nameOfModel[0];
  	initValues = "ANDRADE";
  	double epsilon = vEpsilonConv[0];
  	int maxNIteration = maxIt[0];
  	bool verbose = true;
  	double *parameters;
  	
    int numberOfCycles = -1;
    double logLik = -1;
	  double convEp = -1;
    Rcout<<model<<endl;
    if ( strcmp(model , "RASCH_A1" ) == 0)
    {
      nuM = Constant::RASCH_A1;
      nPar = nColumn;
    }
    else if ( strcmp(model , "RASCH_A_CONSTANT" ) == 0)
    {
      nuM = Constant::RASCH_A_CONSTANT;
      nPar = nColumn+1;
    }
    else if( strcmp(model , "TWO_PL" ) == 0)
    {
      nuM = Constant::TWO_PL;
      nPar = 2*nColumn;
    }
    else
    {
      nuM = Constant::THREE_PL;
      nPar = 3*nColumn;
    }
    parameters = new double[nPar];
  	estimatingParameterse(DataI, data.nrow() , data.ncol(), nuM , 1, initValues, epsilon, maxNIteration, verbose, parameters, numberOfCycles, logLik, convEp);   
    NumericVector parametersA(nPar);
    for ( int i = 0;i  < nPar; i++ )
    {
  		parametersA[i] = parameters[i];
  	}
    NumericVector numberOfCyclesA(1);
    numberOfCyclesA[0] = numberOfCycles;
    NumericVector logLikA(1);
    logLikA[0] = logLik;
    NumericVector convEpA(1);
    convEpA[0] = convEp;
    Rcout<<"output :\n1) parameters\n2) number of cycles\n3) loglik\n4) convEp"<<endl;
    List z = List::create( parametersA, numberOfCyclesA, logLikA, convEpA ) ;
    return z ;
}
