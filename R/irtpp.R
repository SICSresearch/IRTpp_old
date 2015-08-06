#' Estimate a test item parameters according to Item Response Theory.
#' @param dataset The matrix with the responses from the individuals
#' @param model The model used to calibrate the parameters
#' @param initialvalues The matrix with the initial values for the optimization process
#' @return The item parameters in a matrix.
irtpp <- function(dataset=NULL,model, initialvalues = NULL, filename=NULL){
  if(is.null(dataset)){
    if(is.null(filename)){
      stop("Please provide a dataset to irtpp")
    }
  }
  if(is.list(dataset)){
    print.sentence("irtpp in list mode")
    ret = autoapply(dataset,irtpp,model,initialvalues,filename)
  }
  else{
    model = irtpp.model(model,asnumber=T)
    cuads = as.matrix(read.table(system.file("extdata","Cuads.csv",package="IRTpp"),sep=",",header=T))
    if(is.null(filename)){
      
      if(is.null(initialvalues))
        ret = irtppinterface(dataset,model,cuads)
      if(!is.null(initialvalues))
        ret = irtppinterfacevalues(dataset,model,cuads,initialvalues)
    }
    else{
      dataset = filename;
      if(is.null(initialvalues))
        ret = irtppinterfacefile(dataset,model,cuads)
      if(!is.null(initialvalues))
        ret = irtppinterfacefilevalues(dataset,model,cuads,initialvalues)
    }
  }
  ret
}

#' Estimate the latent traits of the individuals in a test with some given item parameters
#' @param dataset The matrix with the responses from the individuals
#' @param model The model used to calibrate the parameters
#' @param itempars The item parameters for the model.
#' @param method The method to estimate traits
#' @return A list with the patterns and the estimated latent traits
individual.traits<-function(dataset=NULL,model,itempars,method, filename=NULL){
  if(is.null(filename)){
    if(is.null(dataset)){
      stop("Please provide a dataset or filename")
    }
    model = irtpp.model(model,asnumber=T)
    cuads = as.matrix(read.table(system.file("extdata","Cuads.csv",package="IRTpp"),sep=",",header=T))
    est.method = ifelse(method == "EAP", eapinterface, mapinterface)
    est = est.method(zita_par=itempars,dat=dataset,e_model=model,quads=cuads)
    est = list(matrix(est[[1]],ncol=dim(dataset)[[2]],byrow=T),est[[2]])
    names(est) <- c("patterns","trait")
  }
  else{
    model = irtpp.model(model,asnumber=T)
    cuads = as.matrix(read.table(system.file("extdata","Cuads.csv",package="IRTpp"),sep=",",header=T))
    est.method = ifelse(method == "EAP", eapinterfacefile, mapinterfacefile)
    est = est.method(zita_par=itempars,dat=filename,e_model=model,quads=cuads)
    est = list(matrix(est[[1]],ncol=2,byrow=T),est[[2]])
    names(est) <- c("patterns","trait")
  }
  est
}
