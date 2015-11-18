#library(FactoMineR)
#TEST=read.table("~/Desktop/2015-1/TRI/parcial 2/y.txt")
#DIMTEST=c(1,20,21,41,42,65,66,80,81,109)  #5 dimensiones



#' fix.items
#' Establishes tests to fix in the multidimensional model.
#' @param test : A response matrix.
#' @param dimtest: A vector with the ranges of the clusters in the test, the ranges must be specified sequentially \cr
#' For instance in the case fo three clusters \eqn{(a,b,c)} , the vector will be \eqn{(a_l,a_u,b_l,b_u,c_l,c_u)}. Where \eqn{(a_l)} is the lower
#' boundary of the \eqn{a} cluster, and \eqn{a_u} is the upper boundary.
#' @param simTest : a Simulated test object which already specifies the clusters. \cr
#' Warning : this function requires the package Factominer.
#' @export 
fix.items=function(TEST,DIMTEST,simTest){
  if(!is.null(simTest))
  {
    TEST = simTest$test;
    DIMTEST = simTest$clusterlist
  }  
  
  acp=list()
  fijados=NULL
  i=1
  while(i<length(DIMTEST)){
    acp[[i]]=PCA(TEST[,DIMTEST[i]:DIMTEST[i+1]],graph=F)
    cor=acp[[i]]$var$cor[,1]
    fijados[i]=as.numeric(which(cor==max(cor)))+ifelse(i==1,0,DIMTEST[i-1])
    i=i+2
  }
  fijados=na.omit(fijados)
  return(fijados)
}