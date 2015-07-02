#ifndef _INTERFACE_H
#include <iostream>
#include <type/Matrix.h>
#include <boost/dynamic_bitset.hpp>
#include <type/PatternMatrix.h>
#include <model/Model.h>
#include <model/ModelFactory.h>
#include <model/SICSGeneralModel.h>
#include <estimation/classical/EMEstimation.h>
#include <input/Input.h>
#include <time.h>
#include <trace/Trace.h>


#define _INTERFACE_H
double loglikitem3pl(double*,double*,int,int);
double loglik3pl(double*,double*,int,int);
double printFunc(double, double);
double loglik3pl(double* args, double* pars, int a, int b){
  double loglik = 0;
  double (*fptr)(double*, double*, int, int);
  fptr = &ThreePLModel::logLikelihood;
  loglik = (*fptr)(args, pars, a,b);
  return loglik;
}

double loglikitem3pl(double* args, double* pars, int a, int b){
  double loglik = 0;
  double (*fptr)(double*, double*, int, int);
  fptr = &ThreePLModel::itemLogLik;
  loglik = (*fptr)(args, pars, a,b);
  return loglik;
}



void estimatingParameters(int**,int,int,int,int,char*,double,int,bool,double*,int&,double&,double&);
//printFunc();
double printFunc(double xa , double xb){
  //call the banana function
  //banana mark : double* args, double* pars, int nargs, int npars
  //call with two arguments and two values
  double* p = new double [2];
  double* a = new double [2];
  int nargs = 2;
  int npars = 2;
  a[0] = xa;
  a[1] = xb;
  double (*fptr)(double*, double*, int, int);
  fptr = &ThreePLModel::banana;
  Rcpp::Rcout<<"Loading the banana"<<endl;
  double result;
  result = (*fptr)(a, p, nargs,npars);
  Rcpp::Rcout << result << " : "<< "Resultado "<<endl;
  return result;
}

void EAP(double* zeta, int model, int** dataset){
  Input input;
  Matrix<double> cuad(41, 2);
  Model *model = new Model();
  ModelFactory *modelFactory = new SICSGeneralModel();
  model->setModel(modelFactory, ESTIMATION_MODEL);
  model->getParameterModel()->buildParameterSet(model->getItemModel(),
  model->getDimensionModel());
  LatentTraits * latentTraits;
  latentTraits = new LatentTraits(dataSet);
  LatentTraitEstimation * lte = new LatentTraitEstimation();
	lte->setModel(model);
	lte->setLatentTraits(latentTraits);
	lte->setQuadratureNodes(&nodes);
	lte->estimateLatentTraitsEAP();
}

// remember dimI, initValI,
void estimatingParameters(int ** dataI, int nRowsDataI, int nColumnsDataI, int modelI, int dimI, char * initValI, double epsilonConvI, int maxIterI, bool verboseI, double *parametersO, int & numberOfCyclesO, double & logLikO, double & convEp) {
  int ESTIMATION_MODEL = modelI;  //*** to Inteface
	Input input;
	Constant::CONVERGENCE_DELTA = epsilonConvI; //*** to Interface
	Constant::MAX_EM_ITERS = maxIterI; //*** to Interface
	Matrix<double> cuad(41, 2);
	//Create the profiler to profile the program
	Trace* profiler = new Trace("Profile.log");
	profiler->resetTimer("total");
	profiler->startTimer("total");
	profiler->resetTimer("input");
	profiler->startTimer("input");
	profiler->resetTimer("for1");
	profiler->resetTimer("for2");
	profiler->resetTimer("fyr");
	profiler->resetTimer("optim");
	input.importCSV((char *) "../IRTpp/inst/SICSRepository/SICS/Cuads.csv", cuad, 1, 0);
	// **** **** Run model complete and ordered process **** ****
	// Create general pars
	Model *model = new Model();
	// Create general model
	ModelFactory *modelFactory = new SICSGeneralModel();
	PatternMatrix *dataSet = new PatternMatrix(0);
	// Load matrix
	//input.importCSV(args, *dataSet, 1, 0);

	//<*** to Interface
	//copy matrix
	for( int _i_ = 0; _i_ < nRowsDataI; _i_++ )
	{
		vector<char> dset(nColumnsDataI);
		for ( int _j_ = 0; _j_ < nColumnsDataI; _j_++ )
		{
			if ( dataI[_i_][_j_] == 1) dset[_j_] = true;  // taken of Input.cpp
		}
		dataSet->size = nColumnsDataI;
		dataSet->push(dset);
	}
	//***> to Interface


	// set dataset
	//RASCH_A1, RASCH_A_CONSTANT, TWO_PL, THREE_PL
	model->setModel(modelFactory, ESTIMATION_MODEL);
	//This is where it is decided what model is the test to make
	model->getItemModel()->setDataset(dataSet);		//Sets the dataset.
	// set Theta and weight for the EM Estimation
  cout<<"Still not declared thetas"<<endl;
  double theta2 [] = {-5.1225, -4.86637, -4.61025, -4.35412, -4.098, -3.84187, -3.58575,-3.32962, -3.0735, -2.81737, -2.56125, -2.30512, -2.049, -1.79287, -1.53675, -1.28062, -1.0245, -0.768375, -0.51225, -0.256125, 0,0.256125, 0.51225, 0.768375, 1.0245, 1.28062, 1.53675, 1.79287, 2.049, 2.30512, 2.56125, 2.81737, 3.0735, 3.32962, 3.58575, 3.84187,4.098, 4.35412, 4.61025, 4.86637, 5.1225};
  double weight2 [] = {0.000000204842,0.000000736153,0.00000247758,0.00000780904,0.0000230504,0.000063719,0.000164957,0.000399927,0.000908035,0.00193079,0.00384482,0.00717015,0.0125225,0.0204816,0.0313724,0.0450029,0.0604567,0.0760603,0.0896154,0.098882,0.102179,0.098882,0.0896154,0.0760603,0.0604567,0.0450029,0.0313724,0.0204816,0.0125225,0.00717015,0.00384482,0.00193079,0.000908035,0.000399927,0.000164957,0.000063719,0.0000230504,0.00000780904,0.00000247758,0.000000736153,0.000000204842};
	cout<<"Declared and set but unused"<<endl;
  Matrix<double> *theta = new Matrix<double>(1, 41);
	Matrix<double> *weight = new Matrix<double>(1, 41);
  cout<<"Ready to copy"<<endl;
	for (int k = 0; k < cuad.nR(); k++) {
		(*theta)(0, k) = cuad(k, 0);
    (*theta)(0,k) = theta2[k];
		(*weight)(0, k) = cuad(k, 1);
    (*weight)(0,k) = weight2[k];
	}
  cout<<"Copied the thetas , lets out the matrices"<<endl;
  cout<<(*theta);
  cout<<(*weight);
  cout<<"Outed the matrices"<<endl;
	// build parameter set
	model->getParameterModel()->buildParameterSet(model->getItemModel(),
			model->getDimensionModel());

	// Create estimation
	profiler->stopTimer("input");
	profiler->resetTimer("initial");
	profiler->startTimer("initial");
	EMEstimation *em = new EMEstimation();
	//Here is where quadratures must be set.
	//create the quad nodes
	em->setProfiler(profiler);
	QuadratureNodes nodes(theta, weight);
	em->setQuadratureNodes(&nodes);
	em->setModel(model);
	em->setInitialValues(Constant::ANDRADE); // initvalI here!
	profiler->stopTimer("initial");
	//Pass the profiler to the estimation object so it can be used to profile each step
	em->setProfiler(profiler);
	//Run the estimation
	em->estimate();
	model->parameterModel->getParameters(parametersO); //*** to Interface
	numberOfCyclesO = Constant::ITER; //*** to Interface
	convEp = Constant::EPSILONC; //*** to Inferface
	logLikO = Constant::LOGLIKO;
	delete modelFactory;
	delete dataSet;
	delete em;
	delete model;
	profiler->stopTimer("total");
	//Out the profiler here
	//profilerOut(profiler, 4);
	delete profiler;
}
#endif
