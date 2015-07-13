#include <Rcpp.h>
#include <string.h>

#include <model/parameter/OnePLACModel.h>
#include <model/parameter/OnePLACModel.cpp>
#include <model/parameter/OnePLModel.h>
#include <model/parameter/OnePLModel.cpp>
#include <model/parameter/ParameterModel.h>
#include <model/parameter/TwoPLModel.h>
#include <model/parameter/TwoPLModel.cpp>
#include <model/parameter/ThreePLModel.h>
#include <model/parameter/ThreePLModel.cpp>

//util
#include <util/asa111.h>
#include <util/asa111.cpp>
#include <util/util.h>

//type
#include <type/Constant.h>
#include <type/Constant.cpp>
#include <type/DataSet.h>
#include <type/Matrix.h>
#include <type/PatternMatrix.h>
#include <type/PatternMatrix.cpp>
#include <type/QuadratureNodes.h>
#include <type/QuadratureNodes.cpp>

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
#include <estimation/classical/ClassicalEstimation.h>
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


//[[Rcpp::export]]
Rcpp::List irtppinterface(Rcpp::NumericMatrix dat, int e_model, Rcpp::NumericMatrix quads){
  //Load the model

  Model *model = new Model();
  ModelFactory *modelFactory;
  modelFactory = new SICSGeneralModel();
  model->setModel(modelFactory, e_model);
  delete modelFactory;

  //Cast the datset
  PatternMatrix *datSet = new PatternMatrix(0);
  for (int i = 0;  i<dat.nrow(); i++) {
    vector <char> dset(dat.ncol());
    for (int j = 0; j < dat.ncol(); j++) {
      dset[j] = dat[j*dat.nrow()+i];
    }
    datSet->size = dat.ncol();
    datSet->push(dset);
  }
  datSet->size = dat.ncol();

  //Matrices for thetas and weights
  Matrix<double> *theta;
  Matrix<double> *weight;
  theta = new Matrix<double>(1,41);
  weight = new Matrix<double>(1,41);

  //Cast the quadrature matrices
  for (int k = 0; k < quads.nrow(); k++) {
    (*theta)(0,k)=quads[k];
  }
  for (int k = 0; k < quads.nrow(); k++) {
    (*weight)(0,k)=quads[k+quads.nrow()];
  }
  QuadratureNodes nodes(theta,weight);
  //Set datset to model
  model->getItemModel()->setDataset(datSet);

  //Build parameter set
  model->getParameterModel()->buildParameterSet(model->getItemModel(),model->getDimensionModel());

  //create estimation
  EMEstimation em;

  em.setQuadratureNodes(&nodes);

  em.setModel(model);
  em.setInitialValues(Constant::ANDRADE);

  //We estimate here
  em.estimate();
  double* returnpars;
  //TODO size of model.
  returnpars = new double[3*dat.ncol()];
  model->parameterModel->getParameters(returnpars);

  //Return in list
  Rcpp::NumericVector pars(3*dat.ncol());
  for (int i = 0;i < 3*dat.ncol();i++) {
    pars[i] = returnpars[i];
  }
  Rcpp::List z = Rcpp::List::create(pars);
  delete model;
  delete datSet;
  delete theta;
  delete weight;

  return z;
}
