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
void estimatingParameters(int **, int, int, int, int , char *, double, int, bool, double *, int &, double &, double &);

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
	input.importCSV((char *) "Cuads.csv", cuad, 1, 0);
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
	Matrix<double> *theta = new Matrix<double>(1, 41);
	Matrix<double> *weight = new Matrix<double>(1, 41);
	for (int k = 0; k < cuad.nR(); k++) {
		(*theta)(0, k) = cuad(k, 0);
		(*weight)(0, k) = cuad(k, 1);
	}

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
