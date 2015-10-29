
#' Autoapply a lapply to a function if the parameter is a list instead of a parameter.
#' @param clist Argument that is possibly a list
#' @param fun The function
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
print.sentence<-function(...,sep=" ",verbose=T){
  if(verbose){print(paste(sep=sep,...))}
}

#' Returns a parameter matrix from a parameter list.
#' @param pars The parameter list
#' @param model The model whose the list refers to.
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
parameter.list<-function(pars,model="3PL"){
  names=NULL
  if(model=="1PL"){names=c("a","b","c")}
  if(model=="2PL"){names=c("a","b","c")}
  if(model=="3PL"){names=c("a","b","c")}
  if(model=="Rasch"){names=c("a","b","c")}
  if(is.null(names)){stop("No Valid model provided")}
  pars=list(pars[,1],pars[,2],pars[,3])
  names(pars)<-names
  pars
}




#' Return a valid model string for a representation of a model.
#' @param model Model representation in a integer, string, character, or list with a model element
#' @param asnumber Boolean. Set to true if you need the model as a integer (i.e. To interface with cpp)
#' @return model The valid string model for this model.
#' 
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
#' @param stop Optional, If false, this function wont throw a stop. 
#' @return No return, in case that the model is not valid throws a stop, if error is false, Only prints a message
check.model<-function(model,msg="",error=T){
  checkModel(model,msg,error)
}
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
irtpp.models<-function(){
  c("1PL","2PL","3PL","Rasch")
}


