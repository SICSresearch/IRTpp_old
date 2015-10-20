#library(FactoMineR)
#TEST=read.table("~/Desktop/2015-1/TRI/parcial 2/y.txt")
#DIMTEST=c(1,20,21,41,42,65,66,80,81,109)  #5 dimensiones



#' fix.items
#' establece que items fijar en el modelo multi-unidimensional
#' @param TEST: el test que se quiere modelar
#' @param DIMTEST: un vector que indique que items corresponden a que dimension
#' @export 


fix.items=function(TEST,DIMTEST){

acp=list()
fijados=NULL
i=1

while(i<length(DIMTEST)){
acp[[i]]=PCA(y[,DIMTEST[i]:DIMTEST[i+1]],graph=F)
cor=acp[[i]]$var$cor[,1]
fijados[i]=names(cor[cor==max(cor)])
i=i+2
}
fijados=na.omit(fijados)
return(fijados)

}

#fix.items(TEST,DIMTEST)
