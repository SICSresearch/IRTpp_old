#' Estimate a test item parameters according to Item Response Theory.
#' @param dataset The matrix with the responses from the individuals
#' @param model The model used to calibrate the parameters
#' @param dims The dimensions to use on the estimation, remember to use the initial parameters if you want highquality estimation
#' @param initialvalues The matrix with the initial values for the optimization process
#' @param filename Optional argument specifying a CSV file to read instead of a dataset in memory
#' @param output Optional. Additonal arguments that need to be documented by cristian
#' @param fixeditems Optional. Items to fix.
#' @return The item parameters in a matrix.
#' @export
irtpp <- function(dataset=NULL,model, dims =1 ,initialvalues = NULL,
                  filename=NULL, output=NULL, restricted.items=c()){
  
  if(dims > 1){
    print("Entering in the multidim.")
    model = irtpp.model(model,asnumber=T)
    cuads = as.matrix(read.table(system.file("extdata","Cuads10.csv",package="IRTpp"),sep=",",header=T))
    ## Initial must be provided at this point.
    print(typeof(dataset));
    print(cuads);
    print(typeof(initialvalues));
    print("Hallooo miss haimering")
    ret = irtppmultidim(dataset,model,cuads,initialvalues, dims , restricted.items)
  }
  else{
  if(is.null(dataset)){
    if(is.null(filename)){
      stop("Please provide a dataset to irtpp")
    }
  }
  if(is.list(dataset)){
    print.sentence("irtpp in list mode")
    ret = autoapply(dataset,irtpp,model,dims,initialvalues,filename,output)
  }
  else{
    model = irtpp.model(model,asnumber=T)
    cuads = as.matrix(read.table(system.file("extdata","Cuads.csv",package="IRTpp"),sep=",",header=T))
    if(is.null(filename)){
      
      if(is.null(initialvalues))
        ret = irtppinterface(dataset,model,cuads,!is.null(output),ifelse(is.null(output), "", output))
      if(!is.null(initialvalues))
        ret = irtppinterfacevalues(dataset,model,cuads,initialvalues,!is.null(output),ifelse(is.null(output), "", output))
    }
    else{
      dataset = filename;
      if(is.null(initialvalues))
        ret = irtppinterfacefile(dataset,model,cuads,!is.null(output),ifelse(is.null(output), "", output))
      if(!is.null(initialvalues))
        ret = irtppinterfacefilevalues(dataset,model,cuads,initialvalues,!is.null(output),ifelse(is.null(output), "", output))
    }
  }
    }
  ret
  
}

#' Estimate the latent traits of the individuals in a test with some given item parameters
#' @param dataset The matrix with the responses from the individuals
#' @param model The model used to calibrate the parameters
#' @param itempars The item parameters for the model.
#' @param method The method to estimate traits
#' @param filename The input filename instead of a in-memory dataset
#' @param probability_matrix The pointer returned in the estimation to the probability matrix in case it does not need to be recalculated
#' @param output Optional. Additonal arguments that need to be documented by cristian
#' @return A list with the patterns and the estimated latent traits
#' @export
individual.traits<-function(model,
                            itempars,
                            method,
                            dataset             = NULL,
                            filename            = NULL,
                            output              = NULL,
                            probability_matrix  = NULL)
{
  
  model = irtpp.model(model,asnumber=T)
  cuads = as.matrix(read.table(system.file("extdata","Cuads.csv",package="IRTpp"),sep=",",header=T))

  if(is.null(filename)){
    if(is.null(dataset)){
      stop("Please provide a dataset or filename")
    }
    est.method = ifelse(method == "EAP", eapinterface, mapinterface)
  }else{
    est.method = ifelse(method == "EAP", eapinterfacefile, mapinterfacefile)
  }

  est = est.method(zita_par     = itempars,
                   dat          = dataset,
                   e_model      = model,
                   quads        = cuads,
                   to_file_flag = !is.null(output),
                   output_path  = ifelse(is.null(output), "", output),
                   matrix_flag  = !is.null(probability_matrix),
                   prob_matrix  = ifelse(is.null(probability_matrix), "", probability_matrix)
                   )

  est = individual.traits.aux(output=output, dataset=dataset, est=est)

  est
}

individual.traits.aux <- function(output, dataset, est){
  if(is.null(output)){
    est = list(matrix(est[[1]],ncol=dim(dataset)[[2]],byrow=T),est[[2]])
    names(est) <- c("patterns","trait")
  }else{
    est = list(est[[3]])
    names(est) <- c("path")
  }
  est
}
