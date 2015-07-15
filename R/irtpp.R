#' Estimate a test item parameters according to Item Response Theory.
#' @param dataset The matrix with the responses from the individuals
#' @param model The model used to calibrate the parameters
#' @return The item parameters in a matrix.
irtpp <- function(dataset,model){
  model = irtpp.model(model,asnumber=T);
  cuads = as.matrix(read.table(system.file("extdata","Cuads.csv",package="IRTpp"),sep=",",header=T))
  est = irtppinterface(dataset,model,cuads);
  est = unlist(est)
  matrix(est,ncol=3)
}

#' Estimate the latent traits of the individuals in a test with some given item parameters
#' @param dataset The matrix with the responses from the individuals
#' @param model The model used to calibrate the parameters
#' @param itempars The item parameters for the model.
#' @return A list with the patterns and the estimated latent traits
individual.traits<-function(dataset,model,itempars){
  model = irtpp.model(model,asnumber=T);
  cuads = as.matrix(read.table(system.file("extdata","Cuads.csv",package="IRTpp"),sep=",",header=T))
  est = eapinterface(zita_par=itempars,dat=dataset,e_model=model,quads=cuads);
  est = list(matrix(est[[1]],ncol=dim(dataset)[[2]],byrow=T),est[[2]])
  names(est) <- c("patterns","trait")
  est
}