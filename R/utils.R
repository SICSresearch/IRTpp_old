

x2<-function(model,z,patterns,G,FUN){
  ##Expandir los patrones
  #fp = full.pattern.expand(patterns,ncol(patterns)-1);
  #theta = fp[,ncol(patterns)]
  nitems=ncol(patterns)-2
  theta = patterns[,ncol(patterns)]
  frec = patterns[,ncol(patterns)-1]
  groups = quantile(theta,seq(0, 1, length.out = G + 1))
  groups[1] = groups[1] - 0.1
  #summary(theta)
  groups[G + 1] = groups[G + 1] + 0.1
  groups.Ind = findInterval(theta,groups)  #que theta pertenece a que intervalo (grupos)
  groups.Ind = factor(groups.Ind, levels = sort(unique(groups.Ind)))  #volverla un factor
  thetaG = tapply(rep(theta, frec), rep(groups.Ind, frec), FUN = FUN) #por grupos calcule la mediana
  #thetaG es los trazos latentes representantes de los grupos
  prs = matrix(unlist(lapply(thetaG,function(x){probability.3pl(z,theta=x)})),ncol=nitems,byrow = T)
  z = model.transform(z,"3PL","cp","c")
  
  Njs = as.vector(tapply(frec, groups.Ind, sum))
  Obss2 = rowsum(frec * patterns[,-c(ncol(patterns)-1,ncol(patterns))], groups.Ind, reorder = T)/Njs
  
  chi.square = Njs * (Obss2 - prs)^2/(prs * (1 - prs))   #matriz pre estadistica  
  x2 = colSums(chi.square, na.rm = TRUE)        #estadistica
  return(x2)
}

#' Autoapply a lapply to a function if the parameter is a list instead of a parameter.
#' @param clist Argument that is possibly a list
#' @param fun The function
#' @param ... additional arguments
autoapply<-function(clist,fun,...){
  if(!is.list(clist)) r=fun(clist,...) else {
    if(length(clist)==1) r=fun(clist[[1]],...) else r=lapply(clist,fun,...)
  }
  r
}

#' Undefined assignment, Helper function 
#' @param var , The variable to test
#' @param val , The value to return if the tested variables is NULL
#' @return Returns the value of the tested variable if it is not NULL, otherwise, returns the default
ua<-function(var,val){
  if(is.null(var)){val}
  else{var}
}
#' Print Sentence, Helper function
#' Prints different strings concatenating them with a " "
#' @param ... The strings to print.
#' @param sep The separator to pass
#' @param verbose Pass the verbose parameter of a function to make this print invisible when verbose is false.
print.sentence<-function(...,sep=" ",verbose=T){
  if(verbose){print(paste(sep=sep,...))}
}

#' Returns a parameter matrix from a parameter list.
#' @param pars The parameter list
#' @param model The model whose the list refers to.
#' @export
parameter.matrix<-function(pars, model="3PL", dims = 1 , byrow = F){
  cols=0;
  if(model=="1PL"){cols=3}
  if(model=="2PL"){cols=3}
  if(model=="3PL"){cols=3}
  if(model=="Rasch"){cols=3}
  if(model=="mirt"){cols=4}
  if(cols==0){stop("No valid model provided")}
  if(dims>1){
    cols = cols + dims - 1 ;
  }
  if(model=="mirt"){
    pars = unlist(pars)
    matrix(pars[1:(length(pars)-2)],ncol=4,byrow = T)[,1:(dims+2)]
  }else{ matrix(unlist(pars,use.names=F,recursive=T),ncol=cols,byrow = byrow)}
  
}

#' Returns a parameter list from a parameter matrix.
#' @param pars The parameter matrix
#' @param model The model whose the matrix refers to.
#' @export
parameter.list<-function (pars, model = "3PL", dpar = NULL, cp = NULL) 
{
  if (is.null(dpar)) {
    b = "b"
  }
  else {
    b = "d"
  }
  if (is.null(cp)) {
    c = "c"
  }
  else {
    c = "cp"
  }
  names = NULL
  if (model == "1PL") {
    names = c("a", b, c)
  }
  if (model == "2PL") {
    names = c("a", b, c)
  }
  if (model == "3PL") {
    names = c("a", b, c)
  }
  if (model == "Rasch") {
    names = c("a", b, c)
  }
  if (is.null(names)) {
    stop("No Valid model provided")
  }
  pars = list(pars[, 1], pars[, 2], pars[, 3])
  names(pars) <- names
  pars
}




#' Return a valid model string for a representation of a model.
#' @param model Model representation in a integer, string, character, or list with a model element
#' @param asnumber Boolean. Set to true if you need the model as a integer (i.e. To interface with cpp)
#' @return model The valid string model for this model.
#' @export
irtpp.model<-function(model,asnumber=F){
  if(typeof(model)=="list"){
    model = model$model;
  }
  if(typeof(model)=="double"){
    if(model==1){model="1PL"}
    if(model==2){model="2PL"}
    if(model==3){model="3PL"}
  }
  if(typeof(model)=="character"){
    if(model=="1"){model="1PL"}
    if(model=="2"){model="2PL"}
    if(model=="3"){model="3PL"}
  }
  model = toupper(model)
  checkModel(model);
  
  #If model return needs to be an integer.
  
  if(asnumber){
    if(model=="1PL"){model=1}
    if(model=="2PL"){model=2}
    if(model=="3PL"){model=3}
    if(model=="RASCH"){model=4}
  }
  
  if(model=="RASCH"){model="Rasch"}
  model
}

#' Check Model
#' Checks a test according to the model library to validate if it can be estimated or abort the current procedure
#' @param model The model to check
#' @param msg Additional error message
#' @param error Optional, If false, this function wont throw a stop.
#' @return No return, in case that the model is not valid throws a stop, if error is false, Only prints a message
#' @export
check.model<-function(model,msg="",error=T){
  checkModel(model,msg,error)
}
#' @export
checkModel<-function(model,msg="",error=T){
  #check if a model is valid according to the model list.
  if(!(model=="1PL"||model=="2PL"||model=="3PL"||model=="Rasch"))
    if(error){
      stop(model," Is not a valid model",call.=F,msg)
    }
  else {
    print.sentence(model,"Is not a valid model",msg,sep=" ");
  }
}
#' Lists all available models
#' @export
irtpp.models<-function(){
  c("1PL","2PL","3PL","Rasch")
}


#' Model transformations
#' Implements parameter transforms from one parameter to others dependent on the model
#' Currently only 3PL supported with passings between b <-> d and c <-> cp'
#' @param z The parameter matrix or named list.
#' @param model The model to transform
#' @param src The original parameter to transform
#' @param target The target parameter of the transform
#' @return z The parameters with the transform
#' @export
model.transform<-function(z,model,src,target, byrow = F){
  z = parameter.matrix(z, byrow = byrow);
  if(irtpp.model(model)=="3PL"){
    #b = -d/a
    #d = -ab
    if(src == "b" && target == "d"){
      z[,2] <- -z[,1]*z[,2];
    }
    else if(src == "d" && target == "b"){
      z[,2] <- -z[,2]/z[,1];
    }
    else if(src == "c" && target == "cp"){
      z[,3] <- qlogis(z[,3])
    }
    else if(src == "cp" && target == "c"){
      z[,3] <- plogis(z[,3])
    }
    else{
      stop("Source and target not recognized")
    }
  }
  else{
    stop("Model not recognized");
  }
  z
}

#' @export
#' @title pattern.expand
#' @description Expands a pattern frequency matrix by the frequencies in the last column. Inverse of pattern.freqs
#' @param pts The pattern matrix with frequencies
#' 
#' @return test. The original test of the patterns. 
pattern.expand = function(pts){
  pts[rep(1:nrow(pts),pts[,ncol(pts)]),-ncol(pts)]
}




#' @export
#' @title full.pattern.expand
#' @description Parameter expansion by a given column.
#' @param pts A matrix with frequencies in some column
#' @param expandcol A given frequency column, defaults to the last column
full.pattern.expand = function(pts,expandcol = ncol(pts)){
  pts[rep(1:nrow(pts),pts[,expandcol]),]
}


#' @name pattern.freqs
#' @title pattern.freqs
#' @description Calculates the pattern frequencies given a dataset.
#' @param test The data matrix
#' @return Patterns and frequencies in a matrix , the last column are the frequencies
#' @export
pattern.freqs = function(data, traitobj=NULL, onlyTheta = T){
  ret = NULL;
  if(is.null(traitobj)){
  d <- data.frame(data)
  ones <- rep(1,nrow(d))
  pts = aggregate(ones,by=as.list(d),FUN=sum)
  plist = lapply(names(pts)[-length(pts)],function(x){c(pts[x])}[[1]])
  plist = do.call(order,plist)
  ret = pts[plist,]
  }else{
    fr = pattern.freqs(data);
    fr = cbind(fr,traitobj$trait);
    print(fr);
    fr = full.pattern.expand(pts = fr,expandcol = ncol(data)+1) 
    ret = fr[,ncol(data)+2]
  }
  ret
}

#' expande los patrones (x2 itemfit)
indexPat = function(data,pats){
  comprimData = apply(data,MARGIN=1,FUN=paste,collapse="/")
  comprimPats = apply(pats[,1:ncol(data)],MARGIN=1,FUN=paste,collapse="/")
  index = lapply(comprimPats,FUN = function(x) {which(x == comprimData)})
  index
}


#' Corrige fr. esperadas muy pequeñas para la estadística de orlando
#' @param E, es una lista (de longitud nitems) con fr. esperadas de respuesta incorrectas y correctas para cada item
#' @param O, es una lista (de longitud nitems) con fr. observadas de respuesta incorrectas y correctas para cada item
#' @param mincell, es el minimo numero esperado de personas con score k que contestan o no al item j
#' @export

collapseCells <- function(O, E, mincell = 1){
  for(i in 1L:length(O)){ #para i entre 1 y nitems
    On <- O[[i]]
    En <- E[[i]]
    
    L <- En < mincell & En != 0
    while(any(L, na.rm = TRUE)){
      if(!is.matrix(L)) break
      whc <- min(which(rowSums(L) > 0L)) #selecciona los scores con problemas
      if(whc == 1L){ #soluciona problemas para el primer score
        En[2L,] <- En[2L, ] + En[1L,]
        On[2L,] <- On[2L, ] + On[1L,]
        En <- En[-1L,]; On <- On[-1L,]
      } else if(whc == nrow(En)){#soluciona problemas para el ultimo score
        En[nrow(En)-1L,] <- En[nrow(En)-1L, ] + En[nrow(En),]
        On[nrow(On)-1L,] <- On[nrow(On)-1L, ] + On[nrow(On),]
        En <- En[-nrow(En),]; On <- On[-nrow(On),]
      } else {#soluciona problemas para scores intermedios
        ss <- c(sum(On[whc-1L,]), sum(On[whc+1L,]))
        up <- (min(ss) == ss)[1L]
        pick <- if(up) whc-1L else whc+1L
        En[pick,] <- En[pick, ] + En[whc,]
        On[pick,] <- On[pick, ] + On[whc,]
        En <- En[-whc,]; On <- On[-whc,]
      }
      L <- En < mincell & En != 0
    }
    En[En == 0] <- NA
    E[[i]] <- En
    O[[i]] <- On
  }###aca termina el for
  return(list("O"=O, "E"=E))
}

#' Verosimilitudes de los scores clasicos para calcular fr. esperadas en Orlando
#' Calcula las veorimilitudes
#' @param pr: matriz de probabilidad
#' @param nitems, cantidad de items considerados en la matriz de probabilidad

s_ss=function(pr,nitems,G){
  sact = matrix(0,ncol = nitems +1,nrow = G)   ###con los cambios en los indices se parece mucho mas a mirt
  for(m in 1:G){                 
    sant = rep(0,(nitems))         
    sant[1] = 1 - pr[m,1]    #(6)
    sant[2] = pr[m,1]       #(7)
    
    for(k in 2:(nitems)){  #item "i" añadido(empieza desde 2 por q ya incluyo el primer item en (6) y (7) )     
      
      sact[m,1] = (1-pr[m,k]) * sant[1]   #(8)
      
      for(kk in 2:k){ #scores 1:(i-1)
        sact[m,kk] = pr[m,k] * sant[kk-1] + (1 - pr[m,k]) * sant[kk]  #(9)
      }
      sact[m,(k+1)] = pr[m,k] * sant[k] #score "i"  #(10)
      sant = sact[m,]
    }
    
    
  }
  return(sact)
}

###Multidimensional simulation of tests

#reference
#Robert, C. P. Simulation of truncated normal variables. Statistics and Computing (1995) 5, 121?125


#' @export
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

#' @export
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


#' @name test.plot
#' @title test.plot
#' @author Juan Liberato
#' @description plots a test or an item.
#' @param z : Item parameter list
#' @param i : Optional, item to plot.
#' 
#' @return Void, draws a plot of the ICC or the test's ICC's
#' @export
test.plot = function(z, i = NULL){
  #### Añadir list = T a model.transform
  est = NULL;
  est$z = z;
  itms = length(est$z$a);
  pts  =60
  pts = ((1:pts)/pts)*10-5
  itemplot = lapply(pts, function(x){probability.3pl(est$z,theta=x)})
  itpl.matrix= matrix(unlist(itemplot),nrow=length(itemplot),ncol=length(itemplot[[1]]),byrow = T)
  st = 1;
  if(!is.null(i)){
    st = i;
    itms = i;
  }
  for (item in st:itms) {
    comb.mat2 = cbind(pts,itpl.matrix[,item])
    comb.mat2 = comb.mat2[order(comb.mat2[,1]),]
    comb.mat2 = rbind(comb.mat2[1,],comb.mat2,comb.mat2[nrow(comb.mat2),])
    comb.mat2[1,2] = 0
    comb.mat2[nrow(comb.mat2),2] = 1
    theta = comb.mat2[,1]
    p = comb.mat2[,2]
    if(item == 1 || !is.null(i)){plot(y=p,x=theta,type='l') 
    }else{lines(y=p,x=theta) }
    abline(v=est$z$b[[item]], col = "red");
    abline(h=est$z$c[[item]], col = "green");
    a=est$z$a[[item]]
    b=est$z$b[[item]]
    c=est$z$c[[item]]
    # abline(a=-b/2,b=a/5, col = "blue");
  }
}

#' @name z.item
#' @title z.item
#' @description Reads of a parameter list by item or items.
#' @author Juan Liberato
#' @param i Item or item vector to get.
#' @return The parameter list with the items specified.
#' 
z.item <- function(z,i){
    lapply(z,function(x){x[i]})  
}



