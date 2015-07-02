#######################################################################
#' @name Cronbach's alpha
#' @description Cronbach's alpha measures how correlated are the items in a test
#' License: \tab GPL (>= 2) Extacted from multilevel_2.5 package \cr
#' @param items Dataframe that holds the test response data
#' @return Cronbach's alpha for the test.
cronbach<-function(items)
{
  items<-na.exclude(items)	
  N <- ncol(items)
  TOTVAR <- var(apply(items, 1, sum))
  SUMVAR <- sum(apply(items, 2, var))
  ALPHA <- (N/(N - 1)) * (1 - (SUMVAR/TOTVAR))
  OUT<-list(Alpha=ALPHA,N=nrow(items))        
  return(OUT)
}
