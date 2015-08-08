#include <util_package.hpp>

PatternMatrix * getPatternMatrix(string r_path)
{
    PatternMatrix *dataSet;
    Input input;

    char * path = new char[r_path.size() + 1];
    std::copy(r_path.begin(), r_path.end(), path);
    path[r_path.size()] = '\0';
    
    dataSet = new PatternMatrix(0);
    
    input.importCSV(path, *dataSet, 1, 0);

    return dataSet;
}

PatternMatrix * getPatternMatrix(Rcpp::NumericMatrix r_dataSet)
{
    PatternMatrix *dataSet;

    dataSet = new PatternMatrix(0);

    //Cast the datset
    for (int i = 0; i < r_dataSet.nrow(); i++)
    {
        vector <char> dset(r_dataSet.ncol());

        for (int j = 0; j < r_dataSet.ncol(); j++)
            dset[j] = r_dataSet[j*r_dataSet.nrow()+i];

        dataSet->size = r_dataSet.ncol();
        dataSet->push(dset);
    }

    dataSet->size = r_dataSet.ncol();

    return dataSet;
}

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
    void ** status_list;

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
    status_list = em.estimate();
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

    Rcpp::List z = Rcpp::List::create(Rcpp::_["zita"] = pars,
                                      Rcpp::_["iterations"] = (*(int*)status_list[0]),
                                      Rcpp::_["convergence"] = (*(bool*)status_list[1]));

    delete model;
    delete theta;
    delete weight;

    delete (int*)status_list[0];
    delete (bool*)status_list[1];

    return z;
}

Rcpp::List abilityinterface(Rcpp::NumericMatrix zita_par, PatternMatrix * datSet,
                            int e_model, Rcpp::NumericMatrix quads, int method)
{
  //Load the model
  Model *model;
  ModelFactory *modelFactory;
  Matrix<double> *theta;
  Matrix<double> *weight;
  double *** zita_set;
  double ** result;
  int items;

  model = new Model();
  modelFactory = new SICSGeneralModel();
  model->setModel(modelFactory, e_model);

  delete modelFactory;

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