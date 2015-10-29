#######################################################################
## cronbach - funcion del paquete multilevel
## Package: multilevel
## Version: 2.5
## License: GPL (>= 2)
## <environment: namespace:multilevel>

#######################################################################
#' @name Cronbachs alpha
#' @title Alpha de Cronbach
#' @description Cronbach's alpha measures how correlated are the items in a test
#' License :GPL (>= 2) Extacted from multilevel_2.5 package
#' @param test a matrix or a Dataframe that holds the test response data
#' @return Cronbach's alpha for the test.

cronbach <- function (test) 
{
  data <- na.exclude(test)
  N <- ncol(test)
  TOTVAR <- var(apply(test, 1, sum)) #varianza de la suma sobre las filas
  SUMVAR <- sum(apply(test, 2, var)) #suma de las varianzas sobre columnas
  ALPHA <- (N/(N - 1)) * (1 - (SUMVAR/TOTVAR))
  OUT <- list(Alpha = ALPHA, N = nrow(test))
  return(OUT)
}

#######################################################################
#### Coeficiente de correlaci?n biserial

## biserial.cor funcion del paquete ltm
## Package: ltm 
## Version:   1.0-0
## License:	GPL
## <environment: namespace:ltm>  

#######################################################################
#' @name Biserial coefficient 
#' @title Point-Biserial correlation coefficient
#' @description Point-Biserial correlation coefficient is a correlation coefficient
#' used when one variable is continuous and the other variable is dichotomous.
#' License :GPL Adapted from ltm_1.0 package
#' @param x a numeric vector representing the continuous variable. 
#' @param y a numeric vector representing the dichotomous variable.
#' @param use Is a option for the use of missing values. 
#' @param level which level of y to use.
#' @return the value of the point-biserial correlation.


biserial.cor <- function (x, y, use = c("all.obs", "complete.obs"), level = 1) 
{
  if (!is.numeric(x)) 
    stop("'x' must be a numeric variable.\n")
  y <- as.factor(y)
  if (length(levs <- levels(y)) > 2) 
    stop("'y' must be a dichotomous variable.\n")
  if (length(x) != length(y)) 
    stop("'x' and 'y' do not have the same length")
  use <- match.arg(use)
  if (use == "complete.obs") {
    cc.ind <- complete.cases(x, y)
    x <- x[cc.ind]
    y <- y[cc.ind]
  }
  ind <- y == levs[level]
  diff.mu <- mean(x[ind]) - mean(x[!ind])
  prob <- mean(ind)
  diff.mu * sqrt(prob * (1 - prob))/sd(x)
}


#######################################################################
## Curva de Cronbach-Mesbah
## Package:	 CMC
## Version:	 1.0
## License:	 GPL (>= 2)
## <environment: namespace:CMC>

#######################################################################
#' @name Cronbach-Mesbah Curve 
#' @title Cronbach-Mesbah Curve
#' @description Point-Biserial correlation coefficient is a correlation coefficient
#' used when one variable is continuous and the other variable is dichotomous.
#' License :GPL (>= 2) Adapted from CMC_1.0 package
#' @param test a Dataframe that holds the test response data
#' @return The number of items used to calculate the coefficient. 
#' @return The maximum value of the alpha coefficient calculated at each step. 
#' @return The item removed at each step. 
#' @return Tue Cronbach-Mesbah curve plot.


alpha.curve <- function (test) 
{
  data = test
  n = nrow(data)
  k = ncol(data)
  max.vec = c()
  which.max.vec = c()
  complete.alpha = cronbach(data)
  j = 1
  elements.to.remove.1 = seq(1:k)
  alpha.1 = c()
  for (i in 1:length(elements.to.remove.1)) {
    data.reduced = data[, -elements.to.remove.1[i]]
    alpha.1[i] = cronbach(data.reduced)
  }
  max.vec[j] = max(alpha.1)
  which.max.vec[j] = which.max(alpha.1)
  for (j in 2:(k - 2)) {
    elements.to.remove = matrix(0, nrow = (k - j + 1), ncol = (j - 
                                                                 1))
    for (r in 1:length(which.max.vec)) {
      elements.to.remove[, r] = matrix(rep(which.max.vec[r], 
                                           k - j + 1), ncol = 1)
    }
    elements.to.remove = cbind(elements.to.remove, seq(1:k)[-which.max.vec[1:(j - 1)]])
    alpha = c()
    for (i in 1:nrow(elements.to.remove)) {
      data.reduced = data[, -elements.to.remove[i, ]]
      alpha[i] = cronbach(data.reduced)
    }
    max.vec[j] = max(alpha)
    which.max.vec[j] = elements.to.remove[which.max(alpha), j]
  }
  output = data.frame(N.Item = seq(2, k), `Alpha Max` = c(rev(max.vec),
                                                          complete.alpha), `Removed Item` = c(colnames(data)[rev(which.max.vec)], "--"))
  plot(output[, 1], output[, 2], t = "b", xlab = "N.Item", 
       ylab = "Alpha Max")
  text(seq(2, k), output[, 2], c(colnames(data)[rev(which.max.vec)], 
                                 ""), pos = 3, cex = 0.6)
  return(output)
}


################################################################################
## Informaci?n del test
## information funcion del paquete ltm
## Package: ltm 
## Version:	 1.0-0
## License:	GPL
## <environment: namespace:ltm>  


#######################################################################
#' @name Test or Item information 
#' @title Test or Item information
#' @description Computes the amount of test or item information. 
#' License :GPL (>= 2) Adapted from CMC_1.0 package
#' @param object a matrix that contain the coefficient estimates 
#' @param range a interval for which the test information should be computed.
#' @param items the items for which the information shoulb be computed.
#' @return TotalInfo the total amount of information.
#' @return RangeInfo the amount of information in the specified interval.
#' @return RangeProp the proportion of information in the specified interval.
#' @return Range the value of range argument
#' @return items the value of items argument

information <- function (object, range, items = NULL, ...) {
  p <- nrow(object) #numero de items
  
  itms <- if (!is.null(items)) {
    if (!is.numeric(items) && length(items) > p) 
      stop("'items' should be a numeric vector of maximum length ", p,
           ".\n")
    if (any(!items %in% 1:p)) 
      stop("'items' should contain numbers in: ", paste(1:p, collapse = ", "), " indicating the items.\n")
    items
  }
  else 1:p
  
  f <- function(z){
    thetas <- object 
    Z <- cbind(z, 1)
    betas <- thetas[, 1:2] 
    cs <- plogis(thetas[, 3]) 
    pi. <- plogis(Z %*% t(betas)) 
    cs <- matrix(cs, nrow(Z), p, TRUE) 
    pi <- cs + (1 - cs) * pi.
    pqr <- pi * (1 - pi) * (pi./pi)^2
    mat <- t(t(pqr) * (betas[, 1]^2))
    rowSums(mat[, itms, drop = FALSE])
  }
  
  
  I0 <- integrate(f, -10, 10, ...)$value
  I1 <- integrate(f, range[1], range[2], ...)$value
  
  out <- list(InfoRange = I1, InfoTotal = I0, PropRange = I1/I0, 
              range = range, items = items)
  class(out) <- "information"
  out
}