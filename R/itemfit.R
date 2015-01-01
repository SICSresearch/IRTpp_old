model="3PL"
G = 10 
B=100
FUN = median
library(ltm)
library(IRTpp)
data(LSAT7)
LSAT7 = as.matrix(LSAT7)
datos = expand.table(LSAT7)
modelo=tpm(datos)
modelo$coefficients
nitems=ncol(datos)
zita=matrix(0,nrow=nitems,ncol=3)
zita[,1]=modelo$coefficients[,3]
zita[,2]=modelo$coefficients[,2]
zita[,3]=plogis(modelo$coefficients[, 1])
zita=zita[,-4]
zita=parameter.list(zita,"3PL",dpar=T)
z=zita
patterns=factor.scores(modelo)$score.dat[,-c(7,9)]
p=ncol(patterns[,-c(ncol(patterns)-1,ncol(patterns))])


#' itemfit
#' evalua que tan bien se ajusta el modelo a los datos.
#' @param model el modelo implementado.
#' @param z Parametros de los items, el c esta en (0,1) y recibe en el siguiente orden a,d,c OJO:d
#' @param patterns Patrones, frecuencias y trazos
#' @param pval.sim si es verdadero se simularara un valor p con bootstrap
#' @param G es el numero de grupos
#' @param FUN es la funcion con la que se calcula la probabilidad esperada en cada grupo
#' @param B es el numero de iteraciones bootstrap, por defecto es 100
#' @return est Estadistica X² y P-valor
itemfit<-function(model,z,patterns,pval.sim,G,FUN,B=NULL){
  
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
      X.new=simulateTest(model,ncol(patterns)-2,individuals = sum(patterns[,ncol(patterns)-1]),itempars = z)$test[[1]]
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


#################################################################################


#' x2
#' Calcula la estadistica X².
#' @param model el modelo implementado.
#' @param z Parametros de los items
#' @param patterns Matriz con los patrones, las frecuencias y los trazos.
#' @param G Numero de particiones de los trazos.
#' @param FUN funcion de representante de grupo
#' @return x2 Estadistica X²
x2<-function(model,z,patterns,G,FUN){
  ##Expandir los patrones
  #fp = full.pattern.expand(patterns,ncol(patterns)-1);
  #theta = fp[,ncol(patterns)]
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
  z = model.transform(z,"3PL","c","cp")
  prs = matrix(unlist(lapply(thetaG,function(x){probability.3pl(parameter.list(z,dpar=T,cp = T),theta=x)})),ncol=5,byrow = T)
  z = model.transform(z,"3PL","cp","c")
  
  Njs = as.vector(tapply(frec, groups.Ind, sum))
  Obss2 = rowsum(frec * patterns[,-c(ncol(patterns)-1,ncol(patterns))], groups.Ind, reorder = T)/Njs
  
  chi.square = Njs * (Obss2 - prs)^2/(prs * (1 - prs))   #matriz pre estadistica  
  x2 = colSums(chi.square, na.rm = TRUE)        #estadistica
  return(x2)
}


####################################################

itemfit(model = "3PL",patterns = patterns,z = z,pval.sim = T,FUN = median,G=10)

