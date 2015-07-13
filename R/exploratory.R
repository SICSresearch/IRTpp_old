## cronbach - funci?n del paquete multilevel
## Package: multilevel
## Version: 2.5
## License: GPL (>= 2)
## <environment: namespace:multilevel>

#######################################################################
#' @name Cronbachs alpha
#' @title Alpha de Cronbach
#' @description Cronbach's alpha measures how correlated are the items in a test
#' License :GPL (>= 2) Extacted from multilevel_2.5 package
#' @param items Dataframe that holds the test response data
#' @return Cronbach's alpha for the test.

cronbach <- function (items) 
{
  items <- na.exclude(items)
  N <- ncol(items)
  TOTVAR <- var(apply(items, 1, sum)) #varianza de la suma sobre las filas
  SUMVAR <- sum(apply(items, 2, var)) #suma de las varianzas sobre columnas
  ALPHA <- (N/(N - 1)) * (1 - (SUMVAR/TOTVAR))
  OUT <- list(Alpha = ALPHA, N = nrow(items), S.e = SEM)
  return(OUT)
}

data(bhr2000)
head(bhr2000)
cronbach(bhr2000[,2:11])

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
#' License :GPL Extracted from ltm_1.0 package
#' @param x a numeric vector representing the continuous variable. 
#' @param y a numeric vector representing the dichotomous variable.
#' @param use Is a option for the use of missing values. 
#' @param level which level of y to use.
#' @details 
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

#### Curva de Cronbach-Mesbah
## Package:	 CMC
## Version:	 1.0
## License:	 GPL (>= 2)
## <environment: namespace:CMC>

#######################################################################
#' @name Cronbach-Mesbah Curve 
#' @title Cronbach-Mesbah Curve
#' @description Point-Biserial correlation coefficient is a correlation coefficient
#' used when one variable is continuous and the other variable is dichotomous.
#' License :GPL (>= 2) Extracted from CMC_1.0 package
#' @param x a Dataframe that holds the test response data
#' @details 
#' @return The number of items used to calculate the coefficient. 
#' @return The maximum value of the alpha coefficient calculated at each step. 
#' @return The item removed at each step. 
#' @return Tue Cronbach-Mesbah curve plot.


alpha.curve <- function (x) 
{
  data = x
  n = nrow(data)
  k = ncol(data)
  max.vec = c()
  which.max.vec = c()
  complete.alpha = alpha.cronbach(data)
  j = 1
  elements.to.remove.1 = seq(1:k)
  alpha.1 = c()
  for (i in 1:length(elements.to.remove.1)) {
    data.reduced = data[, -elements.to.remove.1[i]]
    alpha.1[i] = alpha.cronbach(data.reduced)
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
    elements.to.remove = cbind(elements.to.remove, seq(1:k)[-which.max.vec[1:(j - 
                                                                                1)]])
    alpha = c()
    for (i in 1:nrow(elements.to.remove)) {
      data.reduced = data[, -elements.to.remove[i, ]]
      alpha[i] = alpha.cronbach(data.reduced)
    }
    max.vec[j] = max(alpha)
    which.max.vec[j] = elements.to.remove[which.max(alpha), 
                                          j]
  }
  output = data.frame(N.Item = seq(2, k), `Alpha Max` = c(rev(max.vec), 
                                                          complete.alpha), `Removed Item` = c(colnames(data)[rev(which.max.vec)], 
                                                                                              "--"))
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
#' @description 
#' 
#' License :GPL (>= 2) Extracted from CMC_1.0 package
#' @param object an object from either class ltm, rasch or tpm.
#' @param range a interval for which the test information should be computed.
#' @param items the items for which the information shoulb be computed.
#' @details 
#' @return TotalInfo the total amount of information.
#' @return RangeInfo the amount of information in the specified interval.
#' @return RangeProp the proportion of information in the specified interval.
#' @return Range the value of range argument
#' @return items the value of items argument

information <- function (object, range, items = NULL, ...) 
{
  if (!class(object) %in% c("grm", "gpcm", "ltm", "rasch", 
                            "tpm")) 
    stop("'object' must inherit from either class 'grm', class 'gpcm', class 'ltm', class 'rasch' or class 'tpm'.\n")
  p <- ncol(object$X)
  itms <- if (!is.null(items)) {
    if (!is.numeric(items) && length(items) > p) 
      stop("'items' should be a numeric vector of maximum length ", 
           p, ".\n")
    if (any(!items %in% 1:p)) 
      stop("'items' should contain numbers in: ", paste(1:p, 
                                                        collapse = ", "), " indicating the items.\n")
    items
  }
  else 1:p
  if (class(object) == "ltm" && (object$ltst$factors > 1 | 
                                   any(unlist(object$ltst[c("inter", "quad.z1", "quad.z2")])))) 
    stop("Information is currently computed only for the two-parameter logistic model.\n")
  f <- function(z) {
    switch(class(object), grm = rowSums(infoprobs(object$coefficients,z)[, itms, drop = FALSE]), 
           gpcm = rowSums(infoGPCM(object$coefficients,z, object$IRT.param)[, itms, drop = FALSE]), 
           ltm = { betas <- object$coefficients
                   Z <- cbind(1, z)
                   mat <- t(t(plogis(Z %*% t(betas)) * (1 - plogis(Z %*%t(betas)))) * betas[, 2]^2)
                   rowSums(mat[, itms, drop = FALSE]) },
           rasch = { betas <- object$coefficients
                     Z <- cbind(1, z)
                     mat <- betas[1, 2]^2 * plogis(Z %*% t(betas)) * (1 - plogis(Z %*% t(betas)))
                     rowSums(mat[, itms, drop = FALSE]) },
           tpm = { thetas <- object$coefficients
                   Z <- cbind(1, z)
                   betas <- thetas[, 2:3]
                   cs <- plogis(thetas[, 1]) * object$max.guessing
                   pi. <- plogis(Z %*% t(betas))
                   cs <- matrix(cs, nrow(Z), p, TRUE)
                   pi <- cs + (1 - cs) * pi.
                   pqr <- pi * (1 - pi) * (pi./pi)^2
                   mat <- t(t(pqr) * betas[, 2]^2)
                   rowSums(mat[, itms, drop = FALSE]) })
  }
  I0 <- integrate(f, -10, 10, ...)$value
  I1 <- integrate(f, range[1], range[2], ...)$value
  out <- list(InfoRange = I1, InfoTotal = I0, PropRange = I1/I0, 
              range = range, items = items, call = object$call)
  class(out) <- "information"
  out
}


################################################################################
## Prueba de paralelismo
## parallel funcion del paquete nFActors
## Package: nFactors
## Version:	 2.3.2
## License:	GPL
## <environment: namespace:nFactors>


#######################################################################
#' @name Parallel Analysis 
#' @title Parallel Analysis of a Correlation Matrix
#' @description 
#' License :GPL Extracted from nFactors_2.3.2 package
#' @param subject a number of subjects (default is 100)
#' @param var a number of variables (default is 10)
#' @param rep a number of replications of the correlation matrix (default is 100)
#' @param quantile a quantile of the distribution on which the decision is made
#' @param model a character "components" or "factors"
#' @param sd a vector of standard deviations of the simulated variables
#' @param ... other parameters for the "mvtnorm" or "corr" functions
#' @details 
#' @return eigen a Data frame containing the information on the distribution of the eigenvalues 
#' @return eigen$mevpea Mean of the eigenvalues distribution
#' @return eigen$sevpea Standard deviation of the eigenvalues distribution
#' @return eigen$qevpea Quantile of the eigenvalues distribution
#' @return eigen$sqevpea Standard erro of the quantile of the eigenvalues distribution
#' @return subject a number of subjects
#' @return variables a number of variables
#' @return centile Selected quantile

parallel <- function (subject = 100, var = 10, rep = 100, cent = 0.05, quantile = cent, 
                      model = "components", sd = diag(1, var), ...) 
{
  r <- subject
  c <- var
  y <- matrix(c(1:r * c), nrow = r, ncol = c)
  ycor <- matrix(c(1:c * c), nrow = c, ncol = c)
  evpea <- NULL
  leg.txt <- "Pearson"
  for (k in c(1:rep)) {
    y <- mvrnorm(n = r, mu = rep(0, var), Sigma = sd, empirical = FALSE)
    corY <- cov(y, ...)
    if (model == "components") 
      diag(corY) <- diag(sd)
    if (model == "factors") 
      corY <- corY - ginv(diag(diag(ginv(corY))))
    evpea <- rbind(evpea, eigen(corY)[[1]])
  }
  SEcentile <- function(sd, n = 100, p = 0.95) {
    return(sd/sqrt(n) * sqrt(p * (1 - p))/dnorm(qnorm(p)))
  }
  sprob <- c(cent)
  mevpea <- sapply(as.data.frame(evpea), mean)
  sevpea <- sapply(as.data.frame(evpea), sd)
  qevpea <- moreStats(evpea, quantile = quantile)[3, ]
  sqevpea <- sevpea
  sqevpea <- sapply(as.data.frame(sqevpea), SEcentile, n = rep, 
                    p = cent)
  result <- list(eigen = data.frame(mevpea, sevpea, qevpea, 
                                    sqevpea), subject = r, variables = c, centile = cent)
  class(result) <- "parallel"
  return(result)
}

################################################################################
## Prueba de paralelismo - plot
## plotParallel funcion del paquete nFActors
## Package: nFactors
## Version:	 2.3.2
## License:	GPL
## <environment: namespace:nFactors>


#######################################################################
#' @name Plot Parallel  
#' @title Plot a Parallel Analysis 
#' @description 
#' License: GPL Extracted from nFactors_2.3.2 package
#' @param parallel a vector of the results of parallel analysis 
#' @param eig eigenvalues to analyse
#' @param x a vector of eigenvalues, a matrix of correlations or a data.frame
#' @param model a character "components" or "factors"
#' @param ... other graphics parameters 
#' @details 

plotParallel <- function (parallel, eig = NA, x = eig, model = "components", 
                          legend = TRUE, ylab = "Eigenvalues", xlab = "Components", 
                          main = "Parallel Analysis", ...) 
{
  if (any(!is.na(x))) 
    eig <- eigenComputes(x, ...)
  if (!inherits(parallel, "parallel")) 
    stop("Method is only for parallel objects")
  if (model == "factors") 
    xlab <- "Factors"
  var <- length(parallel$eigen$qevpea)
  if (length(eig) == 1) {
    Component <- var:1
    Location <- seq(from = 0, to = max(parallel$eigen$qevpea) * 
                      3, length.out = var)
    plot.default(as.numeric(Component), as.numeric(Location), 
                 type = "n", main = main, xlab = xlab, ylab = ylab)
  }
  if (length(eig) > 1) {
    plotuScree(eig, main = main, xlab = xlab, ylab = ylab)
  }
  lines(1:var, parallel$eigen$qevpea, col = "green", type = "p", 
        pch = 2)
  lines(1:var, parallel$eigen$mevpea, col = "red")
  if (legend == TRUE) {
    if (length(eig) == 1) {
      leg <- c("Mean Eigenvalues", "Centiles of the Eigenvalues")
      tco <- c("red", "green")
      co <- c("red", "green")
      pc <- c(NA, 2)
    }
    if (length(eig) > 1) {
      leg <- c("Eigenvalues", "Mean Eigenvalues", "Centiles of the Eigenvalues")
      tco <- c("black", "red", "green")
      co <- c("black", "red", "green")
      pc <- c(1, NA, 2)
    }
    legend("topright", legend = leg, text.col = tco, col = co, 
           pch = pc)
  }
}

################################################################################
## Test unidimensonal 
## unidimTest funci?n del paquete ltm
## Package: ltm
## Version:	 1.0-0
## License:	GPL
## <environment: namespace:ltm>  

#######################################################################
#' @name Unidimensional test  
#' @title Test of Unidimensionality using Modified Parallel Analysis 
#' @description 
#' License: GPL Extracted from ltm_1.0. package
#' @param object an object from either class ltm, rasch or tpm.
#' @param data a matrix or a data.frame; used if object is missing.
#' @param thetas a matrix with IRT model parameter values to be used in rmvlogis; used if object is missing.
#' @param IRT
#' @param z.vals
#' @param B a number of samples for the Monte Carlo procedure
#' @details 


unidimTest <- function (object, data, thetas, IRT = TRUE, z.vals = NULL, B = 100, 
                        ...) 
{
  if (missing(object) && (missing(data) & missing(thetas))) 
    stop("either 'object' or both 'data' and 'thetas' must be supplied.\n")
  if (missing(object)) {
    if (!inherits(data, "matrix") && !inherits(data, "data.frame")) 
      stop("'data' must be either a data.frame or a matrix")
    data <- data.matrix(data)
    if (any(its <- apply(data, 2, function(x) {
      x <- x[!is.na(x)]
      length(unique(x))
    }) > 2)) 
      stop("'data' contain more that 2 distinct values for item(s): ", 
           paste(which(its), collapse = ", "))
    data <- apply(data, 2, function(x) if (all(unique(x) %in% 
                                                 c(1, 0, NA))) 
      x
      else x - 1)
    data <- data[complete.cases(data), ]
    n <- nrow(data)
    if (n == 0) 
      stop("zero rows after case-wise missing values deletion.\n")
    p <- ncol(data)
    if (nrow(thetas) != p) 
      stop("the dimensions of 'data' and 'thetas' do not much.\n")
    parms <- thetas
  }
  else {
    if (!class(object) %in% c("ltm", "rasch", "tpm")) 
      stop("Use only with 'ltm', 'rasch' or 'tpm' objects.\n")
    if (inherits(object, "ltm") && any(c(object$ltst$factors > 
                                           1, object$ltst$quad.z1))) 
      stop("\nfor 'ltm' objects it is assumed that the two-parameter logistic model has been fitted\n\t(i.e., one latent variable and no nonlinear terms).")
    data <- object$X
    data <- data.matrix(data)
    data <- data[complete.cases(data), ]
    n <- nrow(data)
    if (n == 0) 
      stop("\nzero rows after case-wise missing values deletion.\n")
    p <- ncol(data)
    parms <- if (inherits(object, "tpm")) 
      cbind(object$coef[, 2:3], plogis(object$coef[, 1]))
    else object$coef
    fsc <- factor.scores(object, resp.patterns = data)$score.dat
    ablts <- fsc$z1
    se.ablts <- fsc$se.z1
    IRT <- FALSE
  }
  ind <- t(combn(p, 2))
  n.ind <- nrow(ind)
  eigenRho <- function(data, ...) {
    rho <- diag(p)
    for (i in 1:n.ind) {
      r <- ind[i, ]
      rho[rbind(r, rev(r))] <- polychor(data[, r[1]], data[, 
                                                           r[2]], ...)
    }
    rho. <- rho
    diag(rho.) <- rep(max(rho[ind]), p)
    list(Rho = rho, ev = eigen(rho., symmetric = TRUE, only.values = TRUE)$values)
  }
  eigR <- eigenRho(data)
  rho <- eigR$Rho
  Tobs <- eigR$ev
  T.boot <- matrix(0, B, length(Tobs))
  for (b in 1:B) {
    if (!missing(object)) 
      z.vals <- rnorm(n, ablts, se.ablts)
    data.new <- rmvlogis(n, parms, IRT = IRT, z.vals = z.vals)
    T.boot[b, ] <- eigenRho(data.new)$ev
  }
  pval <- (sum(T.boot[, 2] >= Tobs[2], na.rm = TRUE) + 1)/(B + 
                                                             1)
  if (!is.null(cnams <- colnames(data))) 
    dimnames(rho) <- list(cnams, cnams)
  out <- list(Tobs = Tobs, T.boot = T.boot, p.value = pval, 
              Rho = rho, call = if (missing(object)) NULL else object$call)
  class(out) <- "unidimTest"
  out
}

