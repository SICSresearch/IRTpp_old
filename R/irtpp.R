#' Estimate a test item parameters according to Item Response Theory.
#' 
irtpp <- function(dataset,model){
  cuads= as.matrix(read.table(system.file("extdata","Cuads.csv",package="IRTpp"),sep=",",header=T))
  est = irtppinterface(dataset,model,cuads);
  est = unlist(est)
  matrix(est,ncol=3)
}