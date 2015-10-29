###Multidimensional simulation of tests

#reference
#Robert, C. P. Simulation of truncated normal variables. Statistics and Computing (1995) 5, 121?125



rtnorm <-function (n, mean = 0, sd = 1, lower = -Inf, upper = Inf)
{
  if (length(n) > 1)
    n <- length(n)
  mean <- rep(mean, length = n)
  sd <- rep(sd, length = n)
  lower <- rep(lower, length = n)
  upper <- rep(upper, length = n)
  lower <- (lower - mean)/sd
  upper <- (upper - mean)/sd
  ind <- seq(length = n)
  ret <- numeric(n)
  alg <- ifelse(lower > upper, -1, ifelse(((lower < 0 & upper ==
                                              Inf) | (lower == -Inf & upper > 0) | (is.finite(lower) &
                                                                                      is.finite(upper) & (lower < 0) & (upper > 0) & (upper -
                                                                                                                                        lower > sqrt(2 * pi)))), 0, ifelse((lower >= 0 & (upper >
                                                                                                                                                                                            lower + 2 * sqrt(exp(1))/(lower + sqrt(lower^2 + 4)) *
                                                                                                                                                                                            exp((lower * 2 - lower * sqrt(lower^2 + 4))/4))),
                                                                                                                                                                           1, ifelse(upper <= 0 & (-lower > -upper + 2 * sqrt(exp(1))/(-upper +
                                                                                                                                                                                                                                         sqrt(upper^2 + 4)) * exp((upper * 2 - -upper * sqrt(upper^2 +
                                                                                                                                                                                                                                                                                               4))/4)), 2, 3))))
  ind.nan <- ind[alg == -1]
  ind.no <- ind[alg == 0]
  ind.expl <- ind[alg == 1]
  ind.expu <- ind[alg == 2]
  ind.u <- ind[alg == 3]
  ret[ind.nan] <- NaN
  while (length(ind.no) > 0) {
    y <- rnorm(length(ind.no))
    done <- which(y >= lower[ind.no] & y <= upper[ind.no])
    ret[ind.no[done]] <- y[done]
    ind.no <- setdiff(ind.no, ind.no[done])
  }
  stopifnot(length(ind.no) == 0)
  while (length(ind.expl) > 0) {
    a <- (lower[ind.expl] + sqrt(lower[ind.expl]^2 + 4))/2
    z <- rexp(length(ind.expl), a) + lower[ind.expl]
    u <- runif(length(ind.expl))
    done <- which((u <= exp(-(z - a)^2/2)) & (z <= upper[ind.expl]))
    ret[ind.expl[done]] <- z[done]
    ind.expl <- setdiff(ind.expl, ind.expl[done])
  }
  stopifnot(length(ind.expl) == 0)
  while (length(ind.expu) > 0) {
    a <- (-upper[ind.expu] + sqrt(upper[ind.expu]^2 + 4))/2
    z <- rexp(length(ind.expu), a) - upper[ind.expu]
    u <- runif(length(ind.expu))
    done <- which((u <= exp(-(z - a)^2/2)) & (z <= -lower[ind.expu]))
    ret[ind.expu[done]] <- -z[done]
    ind.expu <- setdiff(ind.expu, ind.expu[done])
  }
  stopifnot(length(ind.expu) == 0)
  while (length(ind.u) > 0) {
    z <- runif(length(ind.u), lower[ind.u], upper[ind.u])
    rho <- ifelse(lower[ind.u] > 0, exp((lower[ind.u]^2 -
                                           z^2)/2), ifelse(upper[ind.u] < 0, exp((upper[ind.u]^2 -
                                                                                    z^2)/2), exp(-z^2/2)))
    u <- runif(length(ind.u))
    done <- which(u <= rho)
    ret[ind.u[done]] <- z[done]
    ind.u <- setdiff(ind.u, ind.u[done])
  }
  stopifnot(length(ind.u) == 0)
  ret * sd + mean
}

##################################################################
#  Utilitary Functions
# 1. Normalize a set of vectors contained in a matrix
# the vectors are the rows of the matrix
# a simple vector is possible
# using the metric space (R^n,Sigma)
##################################################################
normalize<-function(x){#normaliza un vector(divide por la norma)
  #control section
  if(!is.numeric(x))
    stop("'x' must be numeric")
  #work section
  if(!is.matrix(x))
    return(x/sqrt(sum(x^2)))
  else
    return(x/matrix(sqrt(apply(x*x,1,sum)),dim(x)[1],dim(x)[2]))
} # end normalize



simulateTestMD <- function(items = 10, individuals = 1000, dims = 3, clusters = 4 , seed = 10)
  {

####Start here
### Decide number of items per cluster.
rem = items%%clusters
nitems = items - rem
itemlist = rep(nitems/clusters,clusters)
itemlist[[clusters]] = itemlist[[clusters]] + rem
##split
print(itemlist)

##determinar direcciones principales
idnoisy = diag(dims)+matrix(rnorm(dims*dims,0.15,0.05),nrow=dims,ncol=dims);
idnoisy = idnoisy * ifelse(idnoisy < 1 , 1, 0) + diag(dims)
idnoisy = normalize(idnoisy)
beta_w =rbind(idnoisy,matrix(rnorm(dims*(clusters-dims),0.3,0.05),nrow = (clusters-dims), dims))
beta_w

noise <- seq(0.1,by=.25 / clusters, length.out=clusters)


# Matrix of the beta directions
dir_beta = matrix(NA,sum(itemlist),dims)

# seed for reproducible experiments
set.seed(seed)

#Perturb the dims
##List item limits
ends = cumsum(itemlist)
inits = (c(0,cumsum(itemlist))+1)[1:length(itemlist)]
dir_beta[items,]
i = 1
for (i in 1:clusters) {
  dir_beta[inits[i]:ends[i],] = matrix(beta_w[i,], itemlist[i], dims, byrow=TRUE) + matrix(runif(itemlist[i]*dims,-noise[i],noise[i]), itemlist[i],dims)
}

dir_beta = normalize(dir_beta)

## Quitar negativos
dir_beta = ifelse(dir_beta<0,0,dir_beta)

## Definir vectores de identificacion.
dir_beta[inits[1:dims],] = diag(dims)


##True reference directions
true_w = matrix(NA,clusters,dims)
for (i in 1:clusters) {
  true_w[i,] <- abs(eigen(t(dir_beta[inits[i]:ends[i],])%*%dir_beta[inits[i]:ends[i],])$vectors)[,1] 
}

### Now simulate itempars

l_a=0.25
u_a = Inf
Alpha <- rtnorm(items, mean = 0, sd = 1.0, lower = l_a,  upper = u_a)#genera los alphas
length(Alpha)

# clasical a-parameters
print(dim(dir_beta))
a = dir_beta * matrix(Alpha, items,dims, byrow=FALSE)

# B parametters

sd.gamma <-1
Gamma <- rnorm(items,0,sd.gamma)
Gamma[inits[1:dims]] = 0;

## C parameters

guessing=runif(items,0,0.25)

theta <- matrix(rnorm(individuals*dims,0,1),individuals,dims)
# this line is to garantee the the covariance matrix is the identity
theta <- theta %*% solve((chol(cov(theta))))
cov(theta)
theta.true <- theta

## Setting prob matrix
eta  <- theta%*% t(a) -  matrix(Gamma,individuals,items,byrow=TRUE)
P = guessing + (1-guessing)/(1+exp(-eta))
## Coins
U <- matrix(runif(items*individuals),individuals,items);

responses <- ifelse(P>U,1,0)

t.score = rowSums(P);
c.score = rowSums(responses);

cor(t.score,c.score)

##Return everything

return (list("test"= responses,"z" = list("a"=a,"b"=Gamma,"c"=guessing),"clusters"=itemlist,"direction" = beta_w))
}

