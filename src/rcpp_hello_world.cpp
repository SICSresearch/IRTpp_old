
#include <Rcpp.h>
#include "interface.h"
#include "interface.cpp"

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
#include <estimation/bayesian/BayesianEstimation.h>
//#include <estimation/bayesian/BayesianEstimation.cpp>
#include <estimation/bayesian/LatentTraitEstimation.h>
//#include <estimation/bayesian/LatentTraitEstimation.cpp>
#include <estimation/classical/ClassicalEstimation.h>
//#include <estimation/classical/ClassicalEstimation.cpp>
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
    double b  = randomd();
    const int items = 5, peoples = 5;
    int **DataI;
  	DataI = new int *[peoples];
  	for ( int i = 0; i < peoples; i++ ) DataI[i] = new int[items];
  	char *model, *initValues;
    int nColumn = 4, nRow = 4;    
  	for ( int i = 0; i < nRow; i++ ) for ( int j = 0; j<  nColumn; j++ ) DataI[i][j] = data[j+i*nColumn];
  	model = "RASCH_A_CONSTANT";
  	initValues = "ANDRADE";
  	double epsilon = 0.001;
  	int maxNIteration = 200;
  	bool verbose = true;
  	double *parameters;
  	parameters = new double[items+1];
    int numberOfCycles = -1;
    double logLik = -1;
	  double convEp = -1;
  	estimatingParameters(DataI, data.ncol(),data.nrow(), nameOfModel[0], 1, initValues, epsilon, maxNIteration, verbose, parameters);
  	for ( int i = 0;i  <= items; i++ )
  	{
  		std::cout<<parameters[i] <<" ";
  	}
    std::cout<<endl;
    NumericVector parametersA(items+1); // items + 1 is for RASCH_A_CONSTANT
    for ( int i = 0;i  <= items; i++ )
    {
  		parametersA[i] = parameters[i];
  	}
    NumericVector numberOfCyclesA(1);
    numberOfCyclesA[0] = numberOfCycles;
    NumericVector logLikA(1);
    logLikA[0] = logLik;
    NumericVector convEpA(1);
    convEpA[0] = convEp;
    List z = List::create( parametersA, numberOfCyclesA, logLikA, convEpA ) ;
    return z ;
}
