#' itemfit
#' evalua que tan bien se ajusta el modelo a los datos.
#' @param model el modelo implementado.
#' @param z Parametros de los items
#' @param patterns Patrones y frecuencias
#' @param theta Trazos por patron
#' @return est Estadistica X² y P-valor
itemfit<-function(model,z,patterns,theta){
  
  
}


#' x2
#' Calcula la estadistica X².
#' @param model el modelo implementado.
#' @param z Parametros de los items
#' @param patterns Matriz con los patrones, las frecuencias y los trazos.
#' @param G Numero de particiones de los trazos.
#' @param FUN funcion de representante de grupo
#' @return x2 Estadistica X²
x2<-function(model,z,patterns,G){
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


model.transform<-function(z,model,src,target){
  z = parameter.matrix(z);
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


pattern.expand = function(pts){
  pts[rep(1:nrow(pts),pts[,ncol(pts)]),-ncol(pts)]
}


full.pattern.expand = function(pts,expandcol){
  pts[rep(1:nrow(pts),pts[,expandcol]),]
}

pattern.freqs = function(data, traitobj){
  d <- data.frame(data)
  ones <- rep(1,nrow(d))
  pts = aggregate(ones,by=as.list(d),FUN=sum)
  plist = lapply(names(pts)[-length(pts)],function(x){c(pts[x])}[[1]])
  plist = do.call(order,plist)
  pts[plist,]
}


probability.3pl<-function (z=NULL, a = z$a, b = z$b, c = z$c , d=z$d , cp = z$cp , theta){
  if(is.null(d)){
    d = -a*b;
  }
  if(is.null(b)){
    b = -d / a;
  }
  if (is.null(cp)) {
    c + ((1 - c)/(1 + exp(-a * (theta - b))))
  }
  else {
    exp(cp)/(1 + exp(cp)) + (1 - (exp(cp)/(1 + exp(cp)))) * 
      (1 + exp(-(a * theta + d)))^(-1)
  }
}


parameter.list<-function (pars, model = "3PL" , dpar = NULL , cp = NULL) 
{
  if(is.null(dpar)){b = "b"}else{b="d"}
  if(is.null(dpar)){c = "c"}else{c="cp"}
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
