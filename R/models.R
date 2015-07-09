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

#'LogLikelihood of a IRT model
#'@param model The model to get the probability from
#'@param est The estimation containing the item parameters, test and individual parameters
loglik<- function(model,est){
  
}