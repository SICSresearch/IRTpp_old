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
#include <type/LatentTraits.h>

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
#include <estimation/bayesian/LatentTraitEstimation.h>

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

Rcpp::List irtpp_aux(PatternMatrix *datSet, int e_model, Rcpp::NumericMatrix quads,
                     Rcpp::NumericMatrix init_val, bool init_val_flag)
{
  EMEstimation em;
  Model *model;
  ModelFactory *modelFactory;
  Matrix<double> *theta;
  Matrix<double> *weight;
  double *** zita_set;
  int items;

  model = new Model();
  modelFactory = new SICSGeneralModel();

  //Matrices for thetas and weights
  theta = new Matrix<double>(1,41);
  weight = new Matrix<double>(1,41);

  //Cast the quadrature matrices
  for (int k = 0; k < quads.nrow(); k++)
    (*theta)(0,k)=quads[k];
  
  for (int k = 0; k < quads.nrow(); k++)
    (*weight)(0,k)=quads[k+quads.nrow()];

  model->setModel(modelFactory, e_model);

  delete modelFactory;
  
  QuadratureNodes nodes(theta,weight);
  //Set datset to model
  model->getItemModel()->setDataset(datSet);

  //Build parameter set
  model->getParameterModel()->buildParameterSet(model->getItemModel(),model->getDimensionModel());

  //create estimation
  em.setQuadratureNodes(&nodes);

  em.setModel(model);
  
  if(!init_val_flag)
    em.setInitialValues(Constant::ANDRADE);
  else
  {
    // Cast the initial values;
    zita_set = new double**[3];

    for(int i = 0; i < 3; i++)
      zita_set[i] = new double *[1];

    items = model->getItemModel()->countItems();

    for(int i = 0; i < 3; i++)
      zita_set[i][0] = new double[items + 1];

    for (int i = 0; i < init_val.ncol(); i++)
      for (int j = 0; j < init_val.nrow(); j++)
        zita_set[i][0][j] = init_val[i * init_val.nrow() + j];

    em.setInitialValues(zita_set);
  }


  //We estimate here
  em.estimate();
  double* returnpars;
  //TODO size of model.
  returnpars = new double[3*datSet->size];
  model->parameterModel->getParameters(returnpars);
  //For 2pl & 1pl
  if(e_model < 3)
    for (int i = 2*datSet->size;i < 3*datSet->size;i++)
      returnpars[i]=0;

  //Return in list
  Rcpp::NumericVector pars(3*datSet->size);

  for (int i = 0;i < 3*datSet->size;i++)
    pars[i] = returnpars[i];

  Rcpp::List z = Rcpp::List::create(pars);

  delete model;
  delete theta;
  delete weight;

  return z;
}

Rcpp::List irtpp_from_file(std::string dat, int e_model, Rcpp::NumericMatrix quads,
                           Rcpp::NumericMatrix init_val, bool init_val_flag)
{
  char * path = new char[dat.size() + 1];
  std::copy(dat.begin(), dat.end(), path);
  path[dat.size()] = '\0';
  Input input;
  PatternMatrix *datSet;
  datSet = new PatternMatrix(0);
  input.importCSV(path, *datSet, 1, 0);

  Rcpp::List result = irtpp_aux(datSet, e_model, quads, init_val, init_val_flag);

  delete datSet;

  return result;
}

Rcpp::List irtpp_from_r(Rcpp::NumericMatrix dat, int e_model, Rcpp::NumericMatrix quads,
                         Rcpp::NumericMatrix init_val, bool init_val_flag)
{
  PatternMatrix *datSet;
  datSet = new PatternMatrix(0);

  //Cast the datset
  for (int i = 0;  i<dat.nrow(); i++)
  {
    vector <char> dset(dat.ncol());

    for (int j = 0; j < dat.ncol(); j++)
      dset[j] = dat[j*dat.nrow()+i];
    
    datSet->size = dat.ncol();
    datSet->push(dset);
  }

  datSet->size = dat.ncol();

  Rcpp::List result = irtpp_aux(datSet, e_model, quads, init_val, init_val_flag);

  delete datSet;

  return result;
}

//[[Rcpp::export]]
Rcpp::List irtppinterfacevalues(Rcpp::NumericMatrix dat, int e_model, Rcpp::NumericMatrix quads, Rcpp::NumericMatrix init_val)
{
  // .-.
  return irtpp_from_r(dat, e_model, quads, init_val, true);
}

//[[Rcpp::export]]
Rcpp::List irtppinterface(Rcpp::NumericMatrix dat, int e_model, Rcpp::NumericMatrix quads)
{
  Rcpp::NumericMatrix init_val(1,1);
  return irtpp_from_r(dat, e_model, quads, init_val, false);
}

//[[Rcpp::export]]
Rcpp::List irtppinterfacefile(std::string dat, int e_model, Rcpp::NumericMatrix quads)
{
  Rcpp::NumericMatrix init_val(1,1);
  return irtpp_from_file(dat, e_model, quads, init_val, false);
}

//[[Rcpp::export]]
Rcpp::List irtppinterfacefilevalues(std::string dat, int e_model, Rcpp::NumericMatrix quads, Rcpp::NumericMatrix init_val)
{
  // .-.
  return irtpp_from_file(dat, e_model, quads, init_val, true);
}

Rcpp::List abilityinterface(Rcpp::NumericMatrix zita_par, Rcpp::NumericMatrix dat, int e_model, Rcpp::NumericMatrix quads, int method)
{
  //Load the model
  Model *model;
  ModelFactory *modelFactory;
  PatternMatrix *datSet;
  Matrix<double> *theta;
  Matrix<double> *weight;
  double *** zita_set;
  double ** result;
  int items;

  model = new Model();
  modelFactory = new SICSGeneralModel();
  model->setModel(modelFactory, e_model);

  delete modelFactory;

  //Cast the datset
  datSet = new PatternMatrix(0);

  for (int i = 0; i < dat.nrow(); i++)
  {
    vector <char> dset(dat.ncol());
    for (int j = 0; j < dat.ncol(); j++)
      dset[j] = dat[j*dat.nrow()+i];
    datSet->size = dat.ncol();
    datSet->push(dset);
  }

  datSet->size = dat.ncol();

  //Matrices for thetas and weights
  theta = new Matrix<double>(1,41);
  weight = new Matrix<double>(1,41);

  //Cast the quadrature matrices
  for (int k = 0; k < quads.nrow(); k++)
    (*theta)(0,k)=quads[k];

  for (int k = 0; k < quads.nrow(); k++)
    (*weight)(0,k)=quads[k+quads.nrow()];

  QuadratureNodes nodes(theta,weight);
  //Set datset to model
  model->getItemModel()->setDataset(datSet);

  //Build parameter set
  model->getParameterModel()->buildParameterSet(model->getItemModel(),model->getDimensionModel());

  /*
   * Now we will run the estimation of individual parameter
   */

  //Cast the zita set
  zita_set = new double**[3];

  for(int i = 0; i < 3; i++)
    zita_set[i] = new double *[1];

  items = model->getItemModel()->countItems();

  for(int i = 0; i < 3; i++)
    zita_set[i][0] = new double[items + 1];

  for (int i = 0; i < zita_par.ncol(); i++)
    for (int j = 0; j < zita_par.nrow(); j++)
      zita_set[i][0][j] = zita_par[i * zita_par.nrow() + j];

  //Now create the estimation
  LatentTraitEstimation lte(datSet);
  //Pass the model
  lte.setModel(model);
  //Pass the quadrature nodes
  lte.setQuadratureNodes(&nodes);

  //Ready to estimate
  if(method == 0)
    lte.estimateLatentTraitsEAP(zita_set);
  else
    lte.estimateLatentTraitsMAP(zita_set);

  result = lte.lt->getListPatternTheta();

  //Return in list
  Rcpp::NumericVector pars1(lte.lt->pm->countItems() * lte.lt->pm->matrix.size());
  Rcpp::NumericVector pars2(lte.lt->pm->matrix.size());

  for(unsigned int i = 0; i < lte.lt->pm->matrix.size(); i++)
  {
      for(int j = 0; j < items; j++)
        pars1[i * items + j] = result[i][j];

      pars2[i] = result[i][items];
  }

  Rcpp::List z = Rcpp::List::create(pars1, pars2);

  delete model;
  delete datSet;
  delete theta;
  delete weight;

  return z;
}

//[[Rcpp::export]]
Rcpp::List eapinterface(Rcpp::NumericMatrix zita_par, Rcpp::NumericMatrix dat, int e_model, Rcpp::NumericMatrix quads)
{
  return abilityinterface(zita_par, dat, e_model, quads, 0);
}

//[[Rcpp::export]]
Rcpp::List mapinterface(Rcpp::NumericMatrix zita_par, Rcpp::NumericMatrix dat, int e_model, Rcpp::NumericMatrix quads)
{
  return abilityinterface(zita_par, dat, e_model, quads, 1);
}
