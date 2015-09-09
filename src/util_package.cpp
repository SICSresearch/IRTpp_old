#include <util_package.h>

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
                     Rcpp::NumericMatrix init_val, bool init_val_flag,
                     bool to_file_flag, string output_path)
{
    EMEstimation em;
    Model *model;
    ModelFactory *modelFactory;
    Matrix<double> *theta;
    Matrix<double> *weight;
    double *** zita_set;
    unsigned int items, d = 1, c = 41;
    void ** status_list;

    model = new Model();
    modelFactory = new SICSGeneralModel();

    // Matrices for thetas and weights
    theta = new Matrix<double>(d,c);
    weight = new Matrix<double>(d,c);

    // Cast the quadrature matrices
    for (unsigned int k = 0; k < quads.nrow(); k++)
        (*theta)(0,k)=quads[k];

    for (unsigned int k = 0; k < quads.nrow(); k++)
        (*weight)(0,k)=quads[k+quads.nrow()];

    model->setModel(modelFactory, e_model);


    delete modelFactory;

    QuadratureNodes nodes(theta,weight);
    // Set datset to model
    model->getItemModel()->setDataset(datSet);

    // Build parameter set
    std::cout << "parameterSet is being built, look at the parameters" << std::endl;
    std::cout << model->getDimensionModel() << std::endl;
    model->getParameterModel()->buildParameterSet(model->getItemModel(),model->getDimensionModel());

    // Create estimation
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

    // We estimate here
    status_list = em.estimate();
    double* returnpars;

    returnpars = new double[3*datSet->size];
    model->parameterModel->getParameters(returnpars);
    // For 2pl & 1pl
    if(e_model < 3)
        for (int i = 2*datSet->size;i < 3*datSet->size;i++)
            returnpars[i]=0;

    Rcpp::NumericVector pars(3*datSet->size);
    Rcpp::NumericVector pars_aux(1);

    if(to_file_flag)
    {
        // Return in file
        std::ofstream output;

        output.open(output_path.c_str(), std::ios::app);

        for (int i = 0; i < 3 * datSet->size - 1; i++)
            output << returnpars[i] << " ";

        output << returnpars[3 * datSet->size - 1] << endl;

        output.close();
    }
    else
    {
        // Return in list
        for (int i = 0; i < 3*datSet->size; i++)
            pars[i] = returnpars[i];
    }

    Rcpp::XPtr< Matrix<double> > p_matrix((Matrix<double>*)status_list[2], false);

    Rcpp::List z = Rcpp::List::create(Rcpp::_["zita"] = to_file_flag ? pars_aux : pars,
                                      Rcpp::_["path"] = to_file_flag ? output_path : "No path",
                                      Rcpp::_["iterations"] = (*(int*)status_list[0]),
                                      Rcpp::_["convergence"] = (*(bool*)status_list[1]),
                                      Rcpp::_["probability_matrix"] = p_matrix);

    delete model;
    delete theta;
    delete weight;

    delete (int*)status_list[0];
    delete (bool*)status_list[1];
    delete [] status_list;
    delete [] returnpars;

    return z;
}

Rcpp::List abilityinterface(Rcpp::NumericMatrix zita_par, PatternMatrix * datSet,
                            int e_model, Rcpp::NumericMatrix quads, int method,
                            bool to_file_flag, string output_path,
                            bool matrix_flag, SEXP prob_matrix)
{
    // Load the model
    Model *model;
    ModelFactory *modelFactory;
    Matrix<double> *theta;
    Matrix<double> *weight;
    double *** zita_set;
    double ** result;
    unsigned int items, d = 1, c = 41;

    model = new Model();
    modelFactory = new SICSGeneralModel();
    model->setModel(modelFactory, e_model);

    delete modelFactory;

    // Matrices for thetas and weights
    theta = new Matrix<double>(d,c);
    weight = new Matrix<double>(d,c);

    // Cast the quadrature matrices
    for (int k = 0; k < quads.nrow(); k++)
        (*theta)(0,k)=quads[k];

    for (int k = 0; k < quads.nrow(); k++)
        (*weight)(0,k)=quads[k+quads.nrow()];

    QuadratureNodes nodes(theta,weight);
    // Set datset to model
    model->getItemModel()->setDataset(datSet);

    // Build parameter set
    model->getParameterModel()->buildParameterSet(model->getItemModel(),model->getDimensionModel());

    // Cast the zita set
    zita_set = new double**[3];

    for(int i = 0; i < 3; i++)
        zita_set[i] = new double *[1];

    items = model->getItemModel()->countItems();

    for(int i = 0; i < 3; i++)
        zita_set[i][0] = new double[items + 1];

    for (int i = 0; i < zita_par.ncol(); i++)
        for (int j = 0; j < zita_par.nrow(); j++)
            zita_set[i][0][j] = zita_par[i * zita_par.nrow() + j];

    // Now create the estimation
    LatentTraitEstimation lte(datSet);
    // Pass the model
    lte.setModel(model);
    // Pass the quadrature nodes
    lte.setQuadratureNodes(&nodes);

    // Ready to estimate
    if(matrix_flag)
    {
        Rcpp::XPtr< Matrix<double> > matrix_ref(prob_matrix);
        model->getParameterModel()->probabilityMatrix = matrix_ref;

        if(method == 0)
            lte.estimateLatentTraitsEAP();
        else
            lte.estimateLatentTraitsMAP(zita_set);

        delete model->getParameterModel()->probabilityMatrix;
    }
    else
    {
        if(method == 0)
            lte.estimateLatentTraitsEAP(zita_set);
        else
            lte.estimateLatentTraitsMAP(zita_set);
    }

    result = lte.lt->getListPatternTheta();

    // Return in list
    Rcpp::NumericVector pars1(lte.lt->pm->countItems() * lte.lt->pm->matrix.size());
    Rcpp::NumericVector pars_aux(1);
    Rcpp::NumericVector pars2(lte.lt->pm->matrix.size());

    if(to_file_flag)
    {
        // Return in file
        std::ofstream output;

        output.open(output_path.c_str(), std::ios::app);

        for(int j = 0; j < items; j++)
            output << "V" << j << " ";
        output << "trait" << endl;

        for(unsigned int i = 0; i < lte.lt->pm->matrix.size(); i++)
        {
            for(int j = 0; j < items; j++)
                output << result[i][j] << " ";

            output << result[i][items] << endl;
        }

        output.close();
    }
    else
    {
        for(unsigned int i = 0; i < lte.lt->pm->matrix.size(); i++)
        {
            for(int j = 0; j < items; j++)
                pars1[i * items + j] = result[i][j];

            pars2[i] = result[i][items];
        }
    }

    lte.lt->deleteListPatternTheta(result);

    Rcpp::List z = Rcpp::List::create(Rcpp::_["patterns"] = to_file_flag ? pars_aux : pars1,
                                      Rcpp::_["trait"] = to_file_flag ? pars_aux : pars2,
                                      Rcpp::_["path"] = to_file_flag ? output_path : "No path");

    for(int i = 0; i < 3; i++) delete [] zita_set[i][0];
    for(int i = 0; i < 3; i++) delete [] zita_set[i];
    delete [] zita_set;
    delete model;
    delete theta;
    delete weight;

    return z;
}
