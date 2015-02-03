######################################################
# Función principal para correr el algoritmo de SICS #
######################################################

IRTpp = function(data, ini = c("glm","andrade"),item = c("3PL","2PL","Rasch","1PL"),epsilon = 10^(-4),maxit = 500,verbose = TRUE,export = FALSE,
                 export.path = NULL,dims = 1){
  cl <- match.call()
  ######################
  # Seccion de control #
  ######################
  
  difEnt = unique(as.vector(data))
  if(!all(difEnt %in% c(0,1,NA))){
    stop("Las entradas del dataset solo pueden tomar valores 0, 1 y NA")
  }
  data = ifelse(is.na(data),0,data)
  
  ini = match.arg(ini)
  item = match.arg(item)
  
  if(!is.numeric(epsilon)){
    stop("epsilon debe ser un número")
  }
    
  if(!(is.numeric(maxit))){
    stop("maxit debe ser un entero mayor a 1")
  }
  maxit = ceiling(maxit)
  if(maxit < 1){
    stop("el valor de maxit debe ser un entero mayor a 1")
  }
  
  if((!verbose %in% c(TRUE,FALSE,T,F))){
    stop("verbose debe ser un booleano")
  }
  
  if((!export %in% c(TRUE,FALSE,T,F))){
    stop("export debe ser un booleano")
  }
  
  if(!(is.null(export.path) || class(export.path) == "character")){
    stop("export.path debe ser la ruta para exportar reportes")
  }
  
  ####################################
  # Acá se invoca al código de SICS #
  ###################################
  
  fit = list(coefs = 0,conv = 0,LL = 0,cycles = 0,pats = 0,call=cl,numCuads =0,AIC=0,BIC=0)
  class(fit) = "irtClass"
  fit
}


data=matrix(sapply(X = 1:10000,FUN = function(X)ifelse(runif(1) > .5,1,0)),ncol=10)
data[3,7] = NA
#data[5,8] = 2
IRTpp(data,ini = "glm",item = "3PL",epsilon = 10^(-4),export.path = NULL)
