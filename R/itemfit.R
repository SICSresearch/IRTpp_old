
#' itemfit x2
#' evalua que tan bien se ajusta el modelo a los datos.
#' @param model el modelo implementado.
#' @param z Parametros de los items, el c esta en (0,1) y recibe en el siguiente orden a,d,c OJO:d
#' @param patterns Patrones, frecuencias y trazos
#' @param pval.sim si es verdadero se simularara un valor p con bootstrap
#' @param G es el numero de grupos
#' @param FUN es la funcion con la que se calcula la probabilidad esperada en cada grupo
#' @param B es el numero de iteraciones bootstrap, por defecto es 100
#' @return est Estadistica X?? y P-valor
#' @export
#' 
x2_itemf<-function(model,z,patterns,pval.sim,G,FUN,B=NULL){
  
  check.model(model)
  if(is.null(B)){B=100}
  
  x2obs=x2(model=model,z=z,patterns=patterns,G=G,FUN=FUN)
  
    if (pval.sim==F) {  
      df=G-irtpp.model(model,asnumber=T);
      pvals <- pchisq(x2obs, df = df, lower.tail = FALSE)
    }
  
  else {
    T.boot <- matrix(0, ncol(patterns[,-c(ncol(patterns)-1,ncol(patterns))]), B)
      for (b in 1:B) {
      X.new=simulateTest(model,ncol(patterns)-2,individuals = sum(patterns[,ncol(patterns)-1]),itempars = z)$test
      #X.new <- rmvlogis(n, parms, IRT = FALSE) #genera un test
      #object.new <- if (class(object) %in% c("rasch", "tpm")) {
        #update(, data = X.new)
      #}
      object.new=irtpp(X.new,model = model)
      z=parameter.matrix(object.new$zita)
      z=model.transform(z,model=model,"b","d")
      z =parameter.list(z,dpar = T)
      #else {
       # update(object, formula = X.new ~ z1)
      #}
      #parms.new <- object.new$coefficients
      #X.new <- object.new$patterns$X #patrones de respuesta
      pat.new=individual.traits(model=model,parameter.matrix(object.new$zita),method = "EAP",dataset = X.new,probability_matrix = object.new$probability_matrix)
      #obs.new <- object.new$patterns$obs #frecuencias de los patrones
      #z1.new <- factor.scores(object.new, resp.patterns = X.new)$score.dat$z1c#trazos por patron
      freqs=cbind(pattern.freqs(X.new,pat.new),pat.new$trait)
      #T.boot[, b] <- itmFit(X.new, z1.new, parms.new, obs.new) #matriz con el num de filas igual al num de items y el numero de columnas es el num de iteraciones
      T.boot[, b]=x2(model=model,z=z,patterns=freqs,G=G,FUN=FUN)
      }
    pvals <- (rowSums(T.boot >= rep(x2obs, B), na.rm = TRUE) + 
                1)/(B + 1) #compara las filas(para cada item)
  }
  
  return(cbind(x2obs,pvals))
  
  
}


#' Z3 PERSONFIT
#' Calcula la estad??stica Z3
#' @param data: datos
#' @param zita: Parametros de los items (con el c en todo R y el d en lugar del b)
#' @param patterns: Matriz con los patrones, las frecuencias y los trazos.
#' @return Z3_personfit
#' @export 



z3_personf = function(data,zita,patterns){
  #zita  = est$zita #est. de los parametros de items
  #zita[,3] = qlogis(zita[,3]) # c en todo R
  scores = patterns[,-(ncol(patterns)-1)]
  nitems = ncol(data) #numero de items
  nscores = nrow(patterns) #numero de patrones (scores distintos)
  ninds = nrow(data) #numero de individuos
  
  #Expansi??n de patrones sobre los datos originales
  index = indexPat(data,patterns[,-ncol(patterns)])  #que individuos corresponden a que patron
  scoresTot = numeric(nrow(data)) #define vector donde van a ir los trazos por individuo
  for(mm in 1:nrow(scores)){
    scoresTot[index[[mm]]] = scores[mm,ncol(data) +1]  ##trazo por individuo (en el vector anterior)  
  }
  
  #Matriz de probabilidad
  P = lapply(scoresTot,FUN=function(x){probability.3pl(theta=x,z=zita)}) #p_i (probabilidad de contestar correctamente al item i)
  P=matrix(unlist(P),ncol=nitems,byrow=T)
  
  #Calculo de logverosimilitud
  LL = matrix(0,ncol = ncol(P),nrow = nrow(P)) #matriz de tama??o ninds*nitems
  LL[data == 1] = P[data == 1] #p_{i}(theta estimado) para todos los individuos
  LL[data == 0] = 1 - P[data == 0] #q_{i}(theta estimado) para todos los individuos
  LL = rowSums(log(LL)) #log-verosimilitud (7) del articulo
  
  #Calculo de estimado Z3
  mu = sigmaCuad = rep(0,ninds)  
  for( i in 1:nitems){
    Pi = cbind(P[,i],1 - P[,i]) 
    logPi = log(Pi)
    mu = mu+rowSums(Pi * logPi)    
    #sigmaCuad = sigmaCuad + Pi[,1] * Pi[,2] * (log(Pi[,1]/Pi[,2])^2)
    sigmaCuad = sigmaCuad + (Pi[,1] * Pi[,2] * (log(Pi[,1]/Pi[,2])^2))
  }
  
  # print(dim(LL))
  # print(dim(mu))
  # print(dim(sigmaCuad))
  
  Z3 = (LL - mu) / sqrt(sigmaCuad)
  Z3
  
}


#' Z3 ITEMFIT
#' Calcula la estad??stica Z3
#' @param data: datos
#' @param zita: Parametros de los items (con el c en todo R y el d en lugar del b)
#' @param patterns: Matriz con los patrones, las frecuencias y los trazos.
#' @return Z3_itemfit
#' @export 


#funci??n  que calcula item fit basado en Z3 
z3_itemf = function(data,zita,patterns){
  #zita  = est$zita #est. de los parametros de items
  #zita[,3] = qlogis(zita[,3]) # c en todo R
  scores = patterns[,-(ncol(patterns)-1)]
  nitems = ncol(data) #numero de items
  nscores = nrow(patterns) #numero de patrones (scores distintos)
  ninds = nrow(data) #numero de individuos
  
  #Expansi??n de patrones sobre los datos originales
  index = indexPat(data,patterns[,-ncol(patterns)])  #que individuos corresponden a que patron
  scoresTot = numeric(nrow(data)) #define vector donde van a ir los trazos por individuo
  for(mm in 1:nrow(scores)){
    scoresTot[index[[mm]]] = scores[mm,ncol(data) +1]  ##trazo por individuo (en el vector anterior)  
  }
  
  #Matriz de probabilidad
    P = lapply(scoresTot,FUN=function(x){probability.3pl(theta=x,z=zita)}) #p_i (probabilidad de contestar correctamente al item i)
    P=matrix(unlist(P),ncol=nitems,byrow=T)
    
  #Calculo de logverosimilitud
  LL = matrix(0,ncol = ncol(P),nrow = nrow(P)) #matriz de tama??o ninds*nitems
  LL[data == 1] = P[data == 1] #p_{i}(theta estimado) para todos los individuos
  LL[data == 0] = 1 - P[data == 0] #q_{i}(theta estimado) para todos los individuos
  LL = colSums(log(LL)) #log-verosimilitud (7) del articulo
  
  #Calculo de estimado Z3
  mu = sigmaCuad = rep(0,nitems)  #vector donde va a ir E_3 y SIGMA_3
  for( i in 1:nitems){
    Pi = cbind(P[,i],1 - P[,i]) #define una matriz: en la primer columna p_i y en la otra q_i
    logPi = log(Pi) #log de la matriz anterior
    mu[i] = sum(Pi * logPi)    #(13) del articulo
    sigmaCuad[i] = sum(Pi[,1] * Pi[,2] * (log(Pi[,1]/Pi[,2])^2))  #14 del articulo
    
  }
  Z3 = (LL - mu) / sqrt(sigmaCuad) #z3
  Z3
}


#' Orlando's statistic
#' 
#' calculates the S-X2 values from Maria Orlando and David Thisen.
#' @param patterns: matrix of patterns response, the frequency of each pattern and the latent traits
#' @param G number of quadrature Points
#' @param zita: matrix of estimations of the parameters of the items (discrimination,difficulty, guessing)
#' @param model type of model ( "1PL", 2PL", "3PL" )
#' @return Orlando's statistic, degrees of freedom and pvalue for each item
#' @author SICS Research, National University of Colombia \email{ammontenegrod@@unal.edu.co}
#' @export 
#'
#' @seealso
#' \code{\link{z3_itemf}}, \code{\link{x2_itemf}}
#'
#' @references
#'
#' Orlando, M. & Thissen, D. (2000). Likelihood-based item fit indices for dichotomous item
#' response theory models. \emph{Applied Psychological Measurement, 24}, 50-64.
#'
#' @examples
#'
#'irtpp()



orlando_itemf=function(patterns,G,zita,model){
  
  if(model=="3PL"){mod=3}
  if(model=="2PL"){mod=2}
  if(model=="1PL"){mod=1}
  
  pats = as.matrix(patterns[,-ncol(patterns)])   #patrones sin frecuencia
  frec = as.vector(patterns[,ncol(patterns)-1])     #fr. de los patrones de respuesta
  patsSinFrec = as.matrix(pats[,-ncol(pats)])   #patrones sin frecuencia
  nitems=ncol(patsSinFrec)
  
  
  seq=seq(-6,6,length=61)
  pesos=dnorm(seq)/sum(dnorm(seq))
  Cuad=matrix(c(seq,pesos),byrow=F,ncol=2)# cuadraturas
  
  theta=Cuad
  w.cuad = theta[,2] #pesos
  thetaG = theta[,1] #nodos
  
  zita=parameter.list(zita,"3PL")  
  pr = lapply(thetaG,FUN=function(x){probability.3pl(theta=x,z=zita)}) #probabilidad para cada punto de cuadratura
  pr=matrix(unlist(pr),ncol=nitems,byrow=T)
  
  ### TOTALES POR SCORE
  
  score = rowSums(patsSinFrec )  #suma de los patrones(scores sin agrupar)
  Nk=NULL
  for(i in 1:nitems - 1){  #recorriendo los scores(fijando un score)
    inds = print(which(score == i)) #para agrupar llos patrones q determinen un mismo score (i)
    patsCoin = pats[inds,]   #patrones q determinan el score "i": agrupados
    if(class(patsCoin) == "matrix"){  
      if(dim(patsCoin)[1] != 0){     #(si hay 1 patron o mas)
        Nk[i] = sum(patsCoin[,ncol(patsCoin)])
      }else{
        Nk[i] = 0
      }        
    }else{
      if(class(patsCoin) == "numeric"){
        Nk[i] = patsCoin[length(patsCoin)]
      }
    }
  }
  
  
  ### FRECUENCIAS OBSERVADAS
  
  
  O=list()
  print("nitems")
  print(nitems)
  Oik = matrix(0,ncol = nitems -1 ,nrow = nitems)  # 4 scores(columnas) y nitems(filas)
  for(i in 1:nitems - 1){  #recorriendo los scores(fijando un score)
    inds = print(which(score == i)); #para agrupar los patrones q determinen un mismo score (i)
    for(j in 1:(nitems)){  #recorriendo los items, para llenar la matriz por items fiijado un score (i)
      patsCoin = pats[inds,]   #patrones q determinan el score "i": agrupados
      if(class(patsCoin) == "matrix"){  
        if(dim(patsCoin)[1] != 0){     #(si hay 1 patron o mas)
          Oik[j,i] = sum(apply(X = patsCoin,MARGIN = 1,FUN = function(x){ifelse(x[j] == 1,yes = x[nitems + 1],0)})) 
        }else{
          Oik[j,i] = 0
        }        
      }else{
        if(class(patsCoin) == "numeric"){
          Oik[j,i] = ifelse(patsCoin[j] == 1,yes = patsCoin[nitems + 1],0)
        }
      }
      O[[j]]=cbind(Nk-Oik[j,],Oik[j,])
    }
  }
  
  ### S SIN 0
  
  sact=s_ss(pr,nitems=nitems,G=G)
  Denom = colSums(matrix(rep(w.cuad,nitems -1 ),ncol = nitems - 1) * sact[,-c(1,ncol(sact))])  
  
  
  ### S 0
  
  nitems = nitems - 1
  smoo=list();
  
  for(p in 1:(nitems+1)){
    smoo[[p]]=s_ss(pr[,-p],nitems=nitems,G=G) 
  }
  
  ### Eik   
  
  nitems=ncol(patsSinFrec)
  E=list()
  for(i in 1:length(smoo)){
    E[[i]]=cbind(1-colSums(smoo[[i]][,-ncol(smoo[[i]])]*(pr[,i]*w.cuad))/Denom,colSums(smoo[[i]][,-ncol(smoo[[i]])]*(pr[,i]*w.cuad))/Denom)
  }
  
  
  
  #Estad??stica
  
  S_X2=NULL
  df.S_X2=NULL
  for (i in 1:nitems) E[[i]] <- E[[i]] * Nk
  coll <- collapseCells(O, E, mincell = 1)
  O <- coll$O
  E <- coll$E
  for (i in 1:length(O)) {
    S_X2[i] <- sum((O[[i]] - E[[i]])^2/E[[i]], na.rm = TRUE)
    df.S_X2[i] <- sum(!is.na(E[[i]])) - nrow(E[[i]]) - mod
  }
  
  pval=pchisq(S_X2,df.S_X2,lower.tail = F)
  
    lista=cbind("S_X2"=S_X2,"df.SX2"=df.S_X2,"p.val.S_X2"=pval)
  
  return(lista)
  
}


#' @name test.hessian
#' @title test.hessian
#' @param test The original test
#' @param z. The item parameter list
#' @param th The expanded latent traits
#' @param byitems Calculate the hessian by items
#' @param byinds. Calculate by individuals
#' @return A list with $itemhessian and $indhessian
#' @export
test.hessian<-function(test, z , th , byitems = T , byinds = T){
  itemhessian = NULL;
  i=1
  if(byitems){
    for (i in 1:ncol(test)) {
      zi = lapply(z, function(x){x[[i]]})
      itemhessian[[i]] = item.hessian(test, zi, th);
    }
  }
  indhessian = NULL;
  i=1;
  if(byinds){
    for (i in 1:nrow(test)) {
      indhessian[[i]] = ind.hessian(test, z, th[i]);
    }
    summary(unlist(indhessian))
  }
  list(itemhessian,indhessian)
}

#' @name iitem.hessian
#' @title item.hessian
#' @param test The original test
#' @param zitem. The item
#' @param th The expanded latent traits
#' @return A matrix, the hessian
#' @export
item.hessian<-function(test , zitem , th){
  hs = hessian(function(abc){
    z = list("a"=abc[1],"b"=abc[2],"c"=abc[3]) ;
    loglik(test = test, z = z , traits = th)
  },
  x = unlist(zitem))
  hs
}

#' @name ind.hessian
#' @title ind.hessian
#' @param pattern The original individual pattern
#' @param z. The item parameter list
#' @param trait the individual trait
#' @return A numeric hessian
#' @export
ind.hessian<-function(pattern , z , trait){
  hs = hessian(function(th){
    loglik(test = pattern, z = z , traits = trait)
  },
  x = full.th[i])
  hs
}
