#' Probability function for all models.
#' @param model The model to calculate the probability to
irtpp.p<- function(model){
  model = irtpp.model(model)
  if(model=="3PL") probability.3pl
  if(model=="2PL") probability.2pl
  if(model=="1PL") probability.1pl
  if(model=="Rasch") probability.2pl
}

#'3PL probability function
#'The probability function in the 3PL model.
#'@param z Optional. A list with the parameters a b and c specified by keys.
#'@param a The discrimination parameter
#'@param b The difficulty parameter. (Optional if d is specified)
#'@param c The guessing parameter
#'@param theta The subject's latent trait.
#'@param d Optional. Overrides the b parameter, it is equal to -a*b. Used in some functions.
#'@param cp Optional. Overrides the c parameter, it is logit(c)
probability.3pl = function(z,a=z$a,b=z$b,c=z$c, theta, d=-a*b,cp=NULL){
  if(is.null(cp)){
    c+((1-c)/(1+exp(-a*(theta-b))))
  }
  else{
    exp(cp)/(1+exp(cp))+ (1-(exp(cp)/(1+exp(cp))))*(1 + exp(-(a*theta+d)))^(-1)
  }
}

#'2PL probability function
#'The probability function in the 3PL model.
#'@param z Optional. A list with the parameters a and b specified by keys.
#'@param a The discrimination parameter
#'@param b The difficulty parameter. (Optional if d is specified)
#'@param theta The subject's latent trait.
#'@param d Optional. Overrides the b parameter, it is equal to -a*b. Used in some functions.
probability.2pl = function(z,a=z$a,b=z$b,theta, d=-a*b)((1)/(1+exp(-a*(theta-b))))

#'2PL probability function
#'The probability function in the 2PL model.
#'@param z Optional. A list with the parameter b specified by keys.
#'@param b The difficulty parameter. (Optional if d is specified)
#'@param theta The subject's latent trait.
probability.1pl = function(z,b=z$b,theta)((1)/(1+exp(-(theta-b))))

#'LogLikelihood of a IRT model
#'@param test A matrix of 0's and 1's
#'@param traits A vector with each individual parameter.
#'@param z A list with the parameters a b and c specified by keys.
#' Each key must contain a vector of the item parameters for each parameter
loglik<- function(test,traits,z){
  #prob matrix
  pm = lapply(traits,function(x)probability.3pl(z=z,theta=x))
  #flatten it and flatten the test, then do the if and Sumall with the reduce
  Reduce("+",mapply(function(x,y){
    ifelse(y,log(x),log(1-x))
  },unlist(pm),c(test),SIMPLIFY=F))
}