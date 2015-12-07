
##### Simulación de Test Politomicos usando listas

##Probabilidad para el modelo gpcm

prob_gpcm <-  function (thetas, z, eps = .Machine$double.eps^(1/2)) {
  pred_gpcm <- function (thetas, z) {
    lapply(thetas, function (x) {
      n_x <- length(x)
      t(x[n_x] * outer(z, x[-n_x], "-"))
    })
  }
  
  lapply(pred_gpcm(thetas, z), function (x) {
    num <- exp(apply(x, 2, cumsum)) #sumas acumuladas de las exponenciales por columnas
    if (!is.matrix(num))
      num <- t(num)
    den <- 1 + colSums(num) #Calcula el denominador para la probabilidad
    out <- rbind(1/den, num/rep(den, each = nrow(x))) #pone el denominador en la primera fila,
    # y en las siguienes los numeros sobre denominador 
    ## Arregla la probabilidad redondeando
    if (any(ind <- out == 1))
      out[ind] <- 1 - eps
    if (any(ind <- out == 0))
      out[ind] <- eps
    out
  })
}

#############################################
#######################################################################
#' @name Simulation of test polytomous
#' @title Simulate Test polytomous
#' @description Politomous test simulates according to models GRM and GPCM
#' @param n a value indicating the number of patterns to simulate
#' @param thetas a list containing the parameters of the items where the former 
#' are the values for the difficulty and the last value is that of discrimination
#' @param model a value indicating the model to use
#' @return data the polytomous test according to the parameters established


SimulateTestPoly <- function (n, thetas, model = c("gpcm", "grm")) 
{
  ##Vector de trazos latentes segun la distribución
  z <- rnorm(n)
  
  # el valor espeficicado por el usuario
  
  p <- length(thetas)  ##longitud de los parametros en forma de lista
  nk <- sapply(thetas, length) ##numero de categorias es la longitud de los thetas
  
  ###
  
  ## vector de probabilidades 
  prob <- if (model == "grm") { ##si el modelo es grm funcion de probabilidad
    thetas <- lapply(thetas, function(x) {
      n_x <- length(x)
      cbind(plogis(matrix(x[-n_x], n, n_x - 1, TRUE) - x[n_x] * z), 1)
    })
    #a la función gammas la pone como matriz
    lapply(thetas, function(x) {
      nc <- ncol(x)
      cbind(x[, 1], x[, 2:nc] - x[, 1:(nc - 1)])
    })
  }
  else { ##Si el modelo es gpcm entonces 
          lapply(prob_gpcm(thetas, z), t)
  }
  
  data <- matrix(0, n, p)
  for (j in 1:p) {
    ##Extrae las muestras con las probabilidades halladas anteriormente segun el modelo
    for (i in 1:n) data[i, j] <- sample(nk[j], 1, prob = prob[[j]][i, ])
  }
    data
}




######################## Ejemplos ########################
set.seed(13)
thetas <- lapply(1:5, function(u) c(seq(-1, 1, len = 2), 1.2))

thetas <- lapply(1:5, function(u) c(seq(-1, 2, len = 4), 1.2))

SimulateTestPoly(n=100,thetas,model="grm") 

SimulateTestPoly(n=100,thetas,model="gpcm") 

#### Ejemplo para diferente numero de categorias por item

thetas <- lapply(1:5, function(u) c(seq(-1, 1, len = 2), 1.2))
pru <- list(c(seq(-1, 1, len = 3), 1.5))

thetas1 <- c(thetas,pru)

SimulateTestPoly(n=100,thetas,model="grm") 
