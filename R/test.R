#######################################################################
# This file contains the testing procedures to generate test data for IRTpp package.
#
#Generate test algorithm
#Generate items
#Generate individuals
#Simulate responses

#' Simulate a test according to a model.
#' 
#' Example \eqn{a + b^2}
#' @author Juan Liberato
#' @return A data frame with the simulated dataset
#' @param model A string with the model to simulate, please refer to the model documentation in irtpp documentation.
#' @param items the number of items to simulate
#' @param individuals the number of individuals to simulate
#' @param reps The number of tests to generate with this settings
#' @param independent Set this to false if all the individuals used to simulate the test must be the same
#' @param dims Optional. The number of dimensions to simulate in the test if the model is multidimensional
#' @param boundaries Optional. The kind of boundaries that are specified for the parameters. 
#' @param itempars Optional. Item parameters to be used in the simulation. When the parameters are not generated, the item parameters must be specified.
#' @param seed Optional. Seed to use to generate all the data
simulateTest<-function(model,items,individuals,independent=TRUE,reps=1,dims=1,boundaries=NULL,generated=TRUE,itempars=NULL,seed=NULL)
{
  #set the seed if not set
  seed = ua(seed,floor(runif(1)*10000000))
  set.seed(seed);
  smodel = 0;
  if(model == "1PLAD") smodel=1
  if(model == "1PL") smodel=2
  if(model == "2PL") smodel=3
  if(model == "3PL") smodel=4
  if(smodel==0) stop("No valid model selected")
  #Generate the persons parameters (or read)
  itempars = ua(itempars,simulateItemParameters(items,model,dims,boundaries));
  #Generate the individual parameters (assume normal for now, change later)
  th = rnorm(individuals);
  #Generate the tests
}

#' Undefined assignment, Helper function
#' @param var , The variable to test
#' @param val , The value to return if the tested variables is NULL
#' @return Returns the value of the tested variable if it is not NULL, otherwise, returns the default
ua<-function(var,val){
  if(is.null(var)){val}
  else{var}
}

#' Simulates item parameters depending on a model
#' @param items , Number of items to generate
#' @param model A string with the model to simulate, please refer to the model documentation in irtpp documentation.
#' @param dims Optional. The number of dimensions to simulate in the test if the model is multidimensional
#' @param boundaries Optional. The kind of boundaries that are specified for the parameters.
simulateItemParameters<- function(items, model, dims=1, boundaries=NULL){
  bd = boundaries;
  bd$b_lower = ua(bd$b_lower,-4); 
  bd$b_upper = ua(bd$b_upper,4); 
  bd$a_upper = ua(bd$a_upper,5); 
  bd$a_lower = ua(bd$a_lower,0.0001); 
  bd$c_upper = ua(bd$c_upper,0.35); 
  bd$c_lower = ua(bd$c_lower,0);
  print("juan david")
  b = rnorm(items);
  if(model == "3PL"){
    a = rlnorm(items,meanlog=0,sdlog=1/4)
    c = runif(items,min=bd$c_lower,max=bd$c_upper)
  }
  if(model == "2PL"){
    a = rlnorm(items,meanlog=0,sdlog=1/4)
    c = rep(0,items)
  }
  if(model == "1PLAD"){
    temp = rlnorm(1,meanlog=0,sdlog=1/4)
    a = rep(temp,items)
    c = rep(0,items) 
  }
  if(model == "1PL"){
    a = rep(1,items)
    c = rep(0,items)
  }
  ret = list(a=a,b=b,c=c);
  ret
}
#######################################################################
# Esquema de pruebas
#
# Tiempo
# Convergencia vs Mirt
# COnvergencia Z
# Convergencia theta
# Loglikelihood
