#' SimulateTest.
#' Simulates a test according to a model
#' 
#' @description This function simulates tests according to a IRT model.
#' @author Juan Liberato
#' @return A List with the model, the seed , itempars the item parameters
#' @param model A string with the model to simulate, please refer to the model documentation in irtpp documentation.
#' @param items the number of items to simulate
#' @param individuals the number of individuals to simulate
#' @param reps The number of tests to generate with this settings
#' @param latentTraits
#' @param independent Set this to false if all the individuals used to simulate the test must be the same
#' @param dims Optional. The number of dimensions to simulate in the test if the model is multidimensional TODO (Untested in multidimensional, please do not use this parameter for now)
#' @param boundaries Optional. The kind of boundaries that are specified for the parameters. 
#' @param itempars Optional. Item parameters to be used in the simulation. When the parameters are not generated, the item parameters must be specified.
#' @param seed Optional. Seed to use to generate all the data
#' @param verbose Optional. If true, output is made to know the status of the algorithm
#' @examples
#' k=simulateTest(items=20,individuals=2000,threshold=0.01,dims=1,reps=3,model="3PL")
simulateTest<-function(model="2PL",items=10,individuals=1000,reps=1,dims=1,boundaries=NULL,generated=TRUE,itempars=NULL,latentTraits=NULL,seed=NULL,verbose=F,threshold=0)
{ 
  model=irtpp.model(model)
  #TODO Implement multidimensional test simulation
  dims=1
  ret = NULL;
  ret$model = model;
  #set the seed if not set
  seed = ua(seed,floor(runif(1)*10000000))
  set.seed(seed);
  ret$seed = seed;
  check.model(model);
  #Generate the persons parameters (or read)
  z = ua(itempars,simulateItemParameters(items,model,dims,boundaries));
  ret$itempars = z;
  #Generate the individual parameters (assume normal for now, change later)
  th=rnorm(individuals*dims,0,1)
  th=(th-mean(th))/sd(th)
  th = matrix(th,ncol=dims);
  ret$latentTraits = th
  #Generate the tests
  if(verbose){print("Starting simulation")}
  ret$prob=replicate(reps,do.call(rbind,lapply(th,function(x,z) probability.3pl(theta=x,z=z),z=z)),simplify=F)
  coins=replicate(reps,matrix(runif(items*individuals),ncol=items),simplify=F);
  ret$test=mapply(function(x,y){ifelse(x>y,1,0)},ret$prob,coins,SIMPLIFY=F);
  if(verbose){print("Simulation finished ... ")}
  #Impute the test to exclude individuals in the threshold
  repeat{
    #scores of individuals and items
    if(verbose){print("")}
    individual.scores = lapply(ret$test,function(x) rowSums(x)/items);
    #outlier scores
    outliers.flags = lapply(individual.scores,function(x) ifelse(x<threshold | x>(1-threshold),T,F))
    outliers.indices = lapply(outliers.flags,function(x) as.list(which(x)))
    outliers.missing = lapply(outliers.indices,length)
    outliers.total = Reduce(sum,outliers.missing)
    if(outliers.total<2){
      print.sentence("No outliers left",verbose=verbose)
      break
    }
    else{
      print.sentence("Outliers left",outliers.total,verbose=verbose)
    }
    #resimulate the coins
    if(verbose){print("Resimulating coins ...")}
    newcoins = sapply(outliers.missing,function(x){matrix(runif(x*items),ncol=items)},simplify=F)
    probs = mapply(function(x,y){x[as.numeric(y),]},ret$prob,outliers.indices,SIMPLIFY=F)
    if(verbose){print("Calculating new scores ...")}
    newscores=mapply(function(x,y){ifelse(x>y,1,0)},probs,newcoins,SIMPLIFY=F);
    #assign the new scores in the the old test
    if(verbose){print("Re-assigning new scores ...")}
    mapply(function(x,y,z){
      if(outliers.missing[[z]]>0){
        mapply(function(a,b){
          ret$test[[z]][a,]<<-y[b,]
        },x,seq(length(x)),SIMPLIFY=F);
      }
    },outliers.indices,newscores,seq(length(outliers.indices)),SIMPLIFY=F);
  }
  ret
}

#' Simulates item parameters depending on a model
#' @param items , Number of items to generate
#' @param model A string with the model to simulate, please refer to the model documentation in irtpp documentation.
#' @param dims Optional. The number of dimensions to simulate in the test if the model is multidimensional
#' @param boundaries Optional. The kind of boundaries that are specified for the parameters.
simulateItemParameters<- function(items, model, dims=1, boundaries=NULL){
  model = irtpp.model(model);
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
  if(model == "Rasch"){
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