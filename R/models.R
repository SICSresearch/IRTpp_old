#' Probability function for all models.
#' @param model The model to calculate the probability to
irtpp.p<- function(model){
  model = irtpp.model(model)
  if(model=="3PL") probability.3pl
  if(model=="2PL") probability.2pl
  if(model=="1PL") probability.1pl
  if(model=="Rasch") probability.2pl
}


####Fast probability functions

#'3PL fast probability function
#'The probability function in the 3PL models
#'@param z. List containing the item parameters a, d and c.
#'@param theta. Vector, contains the latent traits of the individual
#'@export
prob.3pl<- function(z, theta){
  prob <- z$c + (1 - z$c)/(1 + exp(-(sum(theta*z$a)+z$d)))
  return(prob)
}


#'2PL fast probability function
#'The probability function in the 2PL model
#'@param z. List containing the item parameters a and d
#'@param theta. Vector, contains the latent traits of the individual
#'@export
prob.2pl<- function(z,theta){
  prob <- (1)/(1 + exp(-(sum(theta*z$a)+z$d)))
  return(prob)
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
#'@export
probability.3pl<-function (z = NULL, a = z$a, b = z$b, c = NULL, d = z$d, cp = z$cp, 
                           theta) 
{
  if("c" %in% names(z) && is.null(c))
  {
    c = z$c;
  }
  if( length(a)>1 && length(a) == length(theta)){
    ##Multidimensional case.
    if( is.null(c) && !is.null(cp)){
      c = plogis(cp)
    }
    if( is.null(c) && is.null(cp)){
      c = 0;
    }
    if( is.null(d)){
      stop("d parameter must be specified in multidim case.")
    }
    prob <- c + (1 - c)/(1 + exp(-(sum(theta*a)+d)))
    return(prob)
  }
  else{
  if (is.null(d)) {
    d = -a * b
  }
  if (is.null(b)) {
    b = -d/a
  }
  if (is.null(cp)) {
    c + ((1 - c)/(1 + exp(-a * (theta - b))))
  }
  else {
    exp(cp)/(1 + exp(cp)) + (1 - (exp(cp)/(1 + exp(cp)))) * 
      (1 + exp(-(a * theta + d)))^(-1)
    }
  }
}

#'2PL probability function
#'The probability function in the 3PL model.
#'@param z Optional. A list with the parameters a and b specified by keys.
#'@param a The discrimination parameter
#'@param b The difficulty parameter. (Optional if d is specified)
#'@param theta The subject's latent trait.
#'@param d Optional. Overrides the b parameter, it is equal to -a*b. Used in some functions.
#'@export
probability.2pl = function(z,a=z$a,b=z$b,theta, d=-a*b)((1)/(1+exp(-a*(theta-b))))

#'2PL probability function
#'The probability function in the 2PL model.
#'@param z Optional. A list with the parameter b specified by keys.
#'@param b The difficulty parameter. (Optional if d is specified)
#'@param theta The subject's latent trait.
#'@export
probability.1pl = function(z,b=z$b,theta)((1)/(1+exp(-(theta-b))))

#'LogLikelihood of a IRT model for UIRT
#'@param test A matrix of 0's and 1's
#'@param traits A vector with each individual parameter, or list of 
#'@param z A list with the parameters a b and c specified by keys.
#' Each key must contain a vector of the item parameters for each parameter
#' @export
loglik<-function (test, traits, z) 
{
  pm = lapply(traits, function(x) probability.3pl(z = z, theta = x));
  Reduce("+", mapply(function(x, y) {
    ifelse(y, log(x), log(1 - x))
  }, unlist(pm), c(t(test)), SIMPLIFY = F))
}

#'LogLikelihood of a IRT model for UIRT (Fast version)
#'@param test A matrix of 0's and 1's
#'@param traits A vector with each individual parameter, or list of 
#'@param z A list with the parameters a b and c specified by keys.
#' Each key must contain a vector of the item parameters for each parameter
#' @export
loglik.f <- function(test, traits, z){
  pm = lapply(traits, function(x) prob.3pl(z = z, theta = x));
  Reduce("+", mapply(function(x, y) {
    ifelse(y, log(x), log(1 - x))
  }, unlist(pm), c(t(test)), SIMPLIFY = F))
}


