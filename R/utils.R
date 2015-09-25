
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
parameter.matrix<-function(pars, model="3PL"){
  cols=0;
  if(model=="1PL"){cols=3}
  if(model=="2PL"){cols=3}
  if(model=="3PL"){cols=3}
  if(model=="Rasch"){cols=3}
  if(cols==0){stop("No valid model provided")}
  matrix(unlist(pars,use.names=F,recursive=T),ncol=cols)
}

#' Returns a parameter list from a parameter matrix.
#' @param pars The parameter matrix
#' @param model The model whose the matrix refers to.
#' @export
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

#' @export
pattern.expand = function(pts){
  pts[rep(1:nrow(pts),pts[,ncol(pts)]),-ncol(pts)]
}

#' @export
full.pattern.expand = function(pts,expandcol){
  pts[rep(1:nrow(pts),pts[,expandcol]),]
}

#' @export
pattern.freqs = function(data, traitobj){
  d <- data.frame(data)
  ones <- rep(1,nrow(d))
  pts = aggregate(ones,by=as.list(d),FUN=sum)
  plist = lapply(names(pts)[-length(pts)],function(x){c(pts[x])}[[1]])
  plist = do.call(order,plist)
  pts[plist,]
}
