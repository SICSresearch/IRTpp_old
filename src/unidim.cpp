#include <iostream>

#include <type/Matrix.h>
#include <type/parameter.h>
#include <type/dataset.h>

#include <utils/andrade.h>
#include <utils/asa111.h>
#include <utils/Input.h>
#include <utils/ramsay.h>

#include <model/model.h>
#include <model/onepl.h>
#include <model/twopl.h>
#include <model/threepl.h>

#include <estimation/emestimation.h>
#include <estimation/estep.h>
#include <estimation/mstep.h>

#include <Rcpp.h>

using std::cout;
using std::endl;

//Please reference dat memory before this or set pointer to NULL.
irtpp::dataset* mat2dat(Rcpp::NumericMatrix mat)
{
  irtpp::dataset* dat;
  
  dat = new irtpp::dataset(0);
  
  for (int r = 0; r < mat.nrow(); r++)
  {
    std::vector<char> v(mat.ncol());
    
    for (int c = 0; c < mat.ncol(); c++)
    {
      v[c] = mat[c*mat.nrow()+r];
    }
    
    dat->size = mat.ncol();
    dat->push(v);
  }
  
  dat->size = mat.ncol();
  
  return dat;
}

Rcpp::NumericMatrix mat2rcpp(Matrix<double>* mat)
{
  Rcpp::NumericMatrix x(mat->nR(),mat->nC());
  
  for (size_t i = 0; i < x.nrow() ; i++)
  {
    for (size_t j = 0; j < x.ncol(); j++)
    {
      x[i*x.nrow()+j*x.ncol()] = (*mat)(i,j);
    }
  }
  
  return x;
}

//' uirtestimate
//' @export
// [[Rcpp::export]]
Rcpp::List uirtestimate(Rcpp::NumericMatrix data , int model)
{
  // Create Dataset.
  void**   status_list;
  irtpp::dataset* d = NULL;
  d = mat2dat(data);
  irtpp::model* m = NULL;
  // d->print();
  
  irtpp::emestimation em(new irtpp::threepl(), d);
  status_list = em.estimate();
  cout<<"seis1"<<endl;
  Matrix<double>* args = em.coef();
  cout<<"seis2"<<endl;
  Rcpp::List pars = Rcpp::List::create(
    Rcpp::_["args"] = mat2rcpp(args),
    Rcpp::_["iterations"] = *((int*)status_list[0])
  );
  
  cout << "Iterations: " << *((int*)status_list[0]) << endl;
  delete (int*)status_list[0];
  delete (bool*)status_list[1];
  delete [] status_list;
  return pars;
  
}