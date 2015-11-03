## Function for initial parameterization :
#V* Input Dataset
#v* Output Some sort of cluster or items for each dimension and num dimensions.
# Find out the clusters of items for the proyection of dem all
# To each cluster assign the principal vector of the transformm yea!

#' @name ParameterInitialization
#' @title parameters.initialize
#' @author SICS Research Group
#' @param dataset The dataset to calculate the initial parameters. A matrix of 0's  and 1's
#' @param dims To set the dimensions of use in the test.
#' @param model The model with dimensions.
#' @param method. Optional, "PCA" for multidimensional and "ANDRADE" for unidimensional are teh current implementations
#' @param red.items. Optional, default true. Reduces the dataset to the dataset without trash items.s
#' @return Initial values for a estimation on the dataset and the model.
#' @export
parameters.initialize<-function(dataset, model , dims , method = "DEFAULT" , red.items=T){
  check.model(model)
  if(model == "3PL"){
    
    items = ncol(dataset);
    individuals = nrow(dataset)
    ## Generate the identity for the first a's
    ## For the rest 0.851
    a = rbind(diag(x = 1, 3, 3),matrix(seq(0.851),items-dims,dims))
    ### For all the b's 0's
    b = seq(0,items);
    ### For all the c's 0.2
    c = seq(0.2,items)
  }
  
  ## Here we return the parameter list for a , b and c.
  
  #for the a's is dims'
  
  ## Here there must be a list of items or clusters.
  
  ## A GLM is made to extract parameters b and c.
  
  # Returns the parameters for a, b and C .
  ret =list("a"=a,"b"=b,"c"=c,"dims"=dims,"model"="3PL")
  ret
  }



