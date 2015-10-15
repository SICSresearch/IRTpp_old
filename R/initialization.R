## Function for initial parameterization :
#V* Input Dataset
#v* Output Some sort of cluster or items for each dimension and num dimensions.
# Find out the clusters of items for the proyection of dem all
# To each cluster assign the principal vector of the transformm yea!

#' @name ParameterInitialization
#' @title parameters.initialize
#' @author SICS Research Group
#' @param dataset The dataset to calculate the initial parameters. A matrix of 0's  and 1's
#' @param model The model with dimensions.
#' @param method. Optional, "PCA" for multidimensional and "ANDRADE" for unidimensional are teh current implementations
#' @param red.items. Optional, default true. Reduces the dataset to the dataset without trash items.s
#' @return Initial values for a estimation on the dataset and the model.
#' @export
parameters.initialize<-function(dataset, model , method = "DEFAULT", red.items=T){
  check.model(model)
  if(model == "3PL"){
    
  }
  
  
  ## Here there must be a list of items or clusters.
  
  ## A GLM is made to extract parameters b and c.
  
  # Returns the parameters for a, b and C .
}



