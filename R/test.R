#' Undefined assignment, Helper function
#' @param var , The variable to test
#' @param val , The value to return if the tested variables is NULL
#' @return Returns the value of the tested variable if it is not NULL, otherwise, returns the default
#' @examples
#' ua(itempars,simulateItemParameters(items,model,dims,boundaries));
ua<-function(var,val){
  if(is.null(var)){val}
  else{var}
}
#' Print Sentence, Helper function
#' Prints different strings concatenating them with a " "
#' @param ... The strings to print.
#' @param sep The separator to pass
print.sentence<-function(...,sep=" "){
  print(paste(sep=sep,...))
}

#' Check Model
#' Checks a test according to the model library to validate if it can be estimated or abort the current procedure
#' @param model The model to check
#' @param msg Additional error message
#' @param stop Optional, If false, this function wont throw a stop. 
#' @return No return, in case that the model is not valid throws a stop, if error is false, Only prints a message
checkModel<-function(model,msg="",error=T){
  #check if a model is valid according to the model list.
  if(!(model=="1PL"||model=="2PL"||model=="3PL"||model=="1PLAD"))
    if(error){
     stop(model," Is not a valid model",call.=F,msg)
    }
  else {
    print.sentence(model,"Is not a valid model",msg,sep=" ");
  }
}

#' SimulateTest.
#' Simulates a test according to a model
#' Example \eqn{a + b^2}
#' @author Juan Liberato
#' @return A List with the model, the seed , itempars the item parameters
#' @param model A string with the model to simulate, please refer to the model documentation in irtpp documentation.
#' @param items the number of items to simulate
#' @param individuals the number of individuals to simulate
#' @param reps The number of tests to generate with this settings
#' @param independent Set this to false if all the individuals used to simulate the test must be the same
#' @param dims Optional. The number of dimensions to simulate in the test if the model is multidimensional TODO (Untested in multidimensional, please do not use this parameter for now)
#' @param boundaries Optional. The kind of boundaries that are specified for the parameters. 
#' @param itempars Optional. Item parameters to be used in the simulation. When the parameters are not generated, the item parameters must be specified.
#' @param seed Optional. Seed to use to generate all the data
#' @param cores Optional. If set to a number set those cores to perform simulations, Only available on UNIX Systems
simulateTest<-function(model,items,individuals,independent=TRUE,reps=1,dims=1,boundaries=NULL,generated=TRUE,itempars=NULL,seed=NULL,cores=NULL)
{
  ret = NULL;
  ret$model = model;
  #set the seed if not set
  seed = ua(seed,floor(runif(1)*10000000))
  set.seed(seed);
  ret$seed = seed;
  checkModel(model);
  #Generate the persons parameters (or read)
  z = ua(itempars,simulateItemParameters(items,model,dims,boundaries));
  ret$itempars = z;
  #Generate the individual parameters (assume normal for now, change later)
  th = matrix(rnorm(individuals*dims),ncol=dims);
  ret$latentTraits = th
  #Generate the tests
  if(!is.null(cores)){
    library(parallel);
    ret$test=replicate(reps, do.call(rbind,mclapply(th,function(x,z) ifelse(runif(1)>probability.3pl(theta=x,z=z),1,0),z=z,mc.cores=cores)),simplify=F)
  }
  else{
    ret$test=replicate(reps, do.call(rbind,lapply(th,function(x,z) ifelse(runif(1)>probability.3pl(theta=x,z=z),1,0),z=z)),simplify=F)
  }
  ret
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
