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
Rcpp::List irtppinterface(Rcpp::IntegerMatrix data, int e_model, Rcpp::NumericMatrix quads){
  //Load the model

  Model *model = new Model();
  ModelFactory *modelFactory;
  modelFactory = new SICSGeneralModel();
  model->setModel(modelFactory, e_model);
  delete modelFactory;

  //Cast the dataset
  cout<<"data rows"<<data.nrow()<<endl;
  cout<<"data cols"<<data.ncol()<<endl;
  PatternMatrix *dataSet = new PatternMatrix(data.ncol());
  for (int i = 0;  i<data.nrow(); i++) {
    vector <char> dset(data.ncol());
    for (int j = 0; j < data.ncol(); j++) {
      dset[j] = data[i*data.ncol()+j];
    }
    dataSet->push(dset);
  }

  //Matrices for thetas and weights
  Matrix<double> *theta;
  Matrix<double> *weight;
  theta = new Matrix<double>(1,41);
  weight = new Matrix<double>(1,41);

  //Cast the quadrature matrices
  for (int k = 0; k < quads.nrow(); k++) {
    (*theta)(0,k)=quads[k*2];
    (*weight)(0,k)=quads[k*2+1];
  }

  QuadratureNodes nodes(theta,weight);
  //Set dataset to model
  model->getItemModel()->setDataset(dataSet);

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
  returnpars = new double[3*data.ncol()];
  cout<<returnpars<<endl<<"address"<<endl;
  model->parameterModel->getParameters(returnpars);

  //Return in list
  Rcpp::NumericVector pars(3*data.ncol());
  for (int i = 0;i < 3*data.ncol();i++) {
    pars[i] = returnpars[i];
  }
  Rcpp::List z = Rcpp::List::create(pars);

  return z;
}


/*
Attaching package: ‘IRTpp’

The following objects are masked _by_ ‘.GlobalEnv’:

    checkModel, cronbach, print.sentence, probability.3pl, simulateItemParameters, simulateTest, ua

The following object is masked from ‘package:graphics’:

    curve

  */
