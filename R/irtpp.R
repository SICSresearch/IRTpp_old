#' Estimate a test item parameters according to Item Response Theory.
#' @param dataset The matrix with the responses from the individuals
#' @param model The model used to calibrate the parameters
#' @return The item parameters in a matrix.
irtpp <- function(dataset,model){
  if(is.list(dataset)){
    print.sentence("irtpp in list mode")
    ret = autoapply(dataset,irtpp,model)
  }
  else{
  model = irtpp.model(model,asnumber=T);
  cuads = as.matrix(read.table(system.file("extdata","Cuads.csv",package="IRTpp"),sep=",",header=T))
  est = irtppinterface(dataset,model,cuads);
  est = unlist(est)
  ret = matrix(est,ncol=3)
  }
  ret
}

#' Estimate the latent traits of the individuals in a test with some given item parameters
#' @param dataset The matrix with the responses from the individuals
#' @param model The model used to calibrate the parameters
#' @param itempars The item parameters for the model.
#' @param method The method to estimate traits
#' @return A list with the patterns and the estimated latent traits
individual.traits<-function(dataset,model,itempars,method){
  model = irtpp.model(model,asnumber=T);
  cuads = as.matrix(read.table(system.file("extdata","Cuads.csv",package="IRTpp"),sep=",",header=T))
  est.method = ifelse(method == "EAP", eapinterface, mapinterface);
  est = est.method(zita_par=itempars,dat=dataset,e_model=model,quads=cuads)
  est = list(matrix(est[[1]],ncol=dim(dataset)[[2]],byrow=T),est[[2]])
  names(est) <- c("patterns","trait")
  est
}