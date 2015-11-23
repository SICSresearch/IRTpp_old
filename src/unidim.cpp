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
  
  for (int i = 0; i < x.nrow() ; i++)
  {
    for (int j = 0; j < x.ncol(); j++)
    {
      x[i+j*x.nrow()] = (*mat)(i,j);
    }
  }
  
  return x;
}

//' uirtestimate
//' @export
// [[Rcpp::export]]
Rcpp::List uirtestimate(Rcpp::NumericMatrix data , int model_)
{
  void**   status_list;
  Matrix<double>* args;
  irtpp::dataset* d = mat2dat(data);

  irtpp::model* m;

  if(model_ == 1)
  {
    irtpp::emestimation em(new irtpp::onepl(), d);
    status_list = em.estimate();
    args = em.coef();
  }
  else if(model_ == 2)
  {
    irtpp::emestimation em(new irtpp::twopl(), d);
    status_list = em.estimate();
    args = em.coef();
  }
  else
  {
    irtpp::emestimation em(new irtpp::threepl(), d);
    status_list = em.estimate();
    args = em.coef();
  }

  double* returnpars;
  Rcpp::NumericVector pars(3*d->size);
  returnpars = new double[3*d->size];

  for(int i = 0; i < args->nR(); i++)
  {
    for(int j = 0; j < args->nC(); j++)
    {
      returnpars[i*args->nC() + j] = (*args)(i, j);
    }
  }

  if(model_ < 3)
    for (int i = 2*d->size;i < 3*d->size;i++)
      returnpars[i]=0;

  for (int i = 0; i < 3*d->size; i++)
    pars[i] = returnpars[i];

  Rcpp::List result = Rcpp::List::create(
    Rcpp::_["z"] = pars,
    Rcpp::_["iterations"] = *((int*)status_list[0])
  );

  delete (int*)status_list[0];
  delete (bool*)status_list[1];
  delete [] status_list;
  delete [] returnpars;

  return result;
}