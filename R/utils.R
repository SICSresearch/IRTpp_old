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

#' Check Model
#' Checks a test according to the model library to validate if it can be estimated or abort the current procedure
#' @param model The model to check
#' @param msg Additional error message
#' @param stop Optional, If false, this function wont throw a stop. 
#' @return No return, in case that the model is not valid throws a stop, if error is false, Only prints a message
checkModel<-function(model,msg="",error=T){
  #check if a model is valid according to the model list.
  if(!(model=="1PL"||model=="2PL"||model=="3PL"||model=="1PLAD"))
    if(error){
      stop(model," Is not a valid model",call.=F,msg)
    }
  else {
    print.sentence(model,"Is not a valid model",msg,sep=" ");
  }
}


#' Apply a function to the indices of a list
#' @param X a vector
#' @param FUN a function
#' @param idx a vector with the indices to consider from the list
#' @examples
#' k=list(3,4,5,6)
#' idx=list(2,4)
#' iapply(function(x) x*x,k,idx)
iapply <- function(FUN,X,idx,...){
  #extract a list with the values to treat and compute the function for each value of the list
  ia.values = lapply(lapply(idx,function(x) X[[x]]),FUN,...)
  #copy the new value to the X list value that is the correspondent
  mapply(function(x,y) {X[[x]]<<-ia.values[[y]]},idx,seq(length(idx)))
  #return the list modified
  X
}
#' Apply a function over indices of multiple lists and store in the original array
#' Warning, some functions need dots to be added to pass the additional arguments
#' @examples
#' k=list(3,4,5,6)
#' y=list(1,0,1,0)
#' idx=list(2,4)
#' p=mapply(function(x,y) x*y,k,y)
#' p
#' y=list(1,1,1,1)
#' imapply(p,function(x,y,...){x*y},idx,k,y)
imapply<- function(X,FUN,idx,...){
  ima.values = mapply(FUN,idx,...)
  mapply(function(x,y) {X[[x]]<<-ima.values[[y]]},idx,seq(length(idx)))
  X
}

#' Return the elements of a list on some indices
#'@examples
#'k=list(3,4,5,6)
#'search.index(k,c(2,3))
search.index<-function(x,idx){
  lapply(idx,function(y) x[[y]])
}


