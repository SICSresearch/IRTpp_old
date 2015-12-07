#' SimulateTest
#' Simulates a test according to a irt model
#' 
#' @export
#' @description Simulates a test according to a IRT Model
#' @author Juan Liberato
#' @param model. The IRT model. in Multidimensional is 3PL.
#' @param items. The items to simulate in the test.
#' @param latentTraits. The latentTraits in the UIRT case.
#' @param individuals. The individuals to simulate in the test.
#' @param The boundaries for parameter generation. See the \code{\link{simulateItemParameters}} documentation for help
#' @param dims. The dimensionality of the test.
#' @param itempars. A list with teh parameters of the items if provided.
#' @param verbose. Verbosibity for more information about the process, use when simulating long tests.
#' @param threshold. A boundary for the percent of answered questions to a test. i.e. 0.05 indicates the individuals will answer maximum the 95% of the questions and minimum the 5% of the questions. UIRT only.
#' @param seed. The seed for reproducibility.
#' @param clusters. In the multidimensional case, the clusters to simulate.
#' @param repetition. Seed of a single set , this seed does not affect the itempars generated.
#' @return List with model, seed  , itempars, and the test simulated. check the documentation of \code{\link{simulateTestMD}} for information on the return in the Multidimensional case.
simulateTest = function(model = "3PL" , items = 10 , latentTraits=NULL ,individuals = 1000
                        , boundaries = NULL, dims = 1 , itempars = NULL , verbose = F ,
                        threshold = 0, seed = 1 , clusters = NULL, repetition = 1){
  ret = NULL;
  
  if(dims > 1){
    if(is.null(clusters)){
      clusters = dims;
    }
   ret =  simulateTestMD(items , individuals , dims , clusters , seed , itempars , repetition)
  }else{
  
  ret$model = model;
  seed = ua(seed,floor(runif(1)*10000000))
  set.seed(seed);
  ret$seed = seed;
  check.model(model);
  z = ua(itempars,simulateItemParameters(items,model,dims,boundaries));
  ret$itempars = z;
  #Generate the individual parameters (assume normal for now, change later)
  th=rnorm(individuals*dims,0,1)
  th=(th-mean(th))/sd(th)
  th = matrix(th,ncol=dims);
  ret$latentTraits = ua(latentTraits,th)
  th=ret$latentTraits;
  gc()
  ##Here th must be exactly the latent traits of these individuals in this test.
  ret$prob=do.call(rbind,lapply(th,function(x,z) probability.3pl(theta=x,z=z),z=z))
  gc()
  coins=matrix(runif(items*individuals),ncol=items);
  gc()
  ret$test = ifelse(ret$prob>coins,1,0)
  coins=NULL
  gc()
  if(verbose){print("Simulation finished ... ")}
    }
  ret
}

#cor(rowSums(ret$prob),rowSums(ret$test))

#' SimulateTest.file
#' Simulates a test according to a model and saves it to files, can be slow due to disk usage.
#' 
#' 
#' @export
#' @description This function simulates tests according to a IRT model.
#' @author Juan Liberato
#' @return A List with the model, the seed , itempars the item parameters
#' @param model A string with the model to simulate, please refer to the model documentation in irtpp documentation.
#' @param items the number of items to simulate
#' @param individuals the number of individuals to simulate
#' @param threshold The threshold that indicates the boundaries on the individual scores (to avoid nearly perfect or nearly )
#' @param reps The number of tests to generate with this settings
#' @param filename A name to give the tests.
#' @param directory The directory to output the tests
#' @param latentTraits A set of latent traits to set them for the individuals
#' @param dims Optional. The number of dimensions to simulate in the test if the model is multidimensional TODO (Untested in multidimensional, please do not use this parameter for now)
#' @param boundaries Optional. The kind of boundaries that are specified for the parameters. 
#' @param itempars Optional. Item parameters to be used in the simulation. When the parameters are not generated, the item parameters must be specified.
#' @param seed Optional. Seed to use to generate all the data
#' @param verbose Optional. If true, output is made to know the status of the algorithm
#' @examples
#' k=simulateTest(items=20,individuals=2000,threshold=0.01,dims=1,reps=3,model="3PL")
simulateTest.file<-function(model="2PL",items=10,individuals=1000,reps=1,dims=1,filename="test",directory=NULL,boundaries=NULL,itempars=NULL,latentTraits=NULL,seed=NULL,verbose=F,threshold=0)
{ 
  dirflag=F;
  if(!is.null(directory)){
    print.sentence("Outputting to directory",directory,verbose=verbose)
    dirflag = T;
  }
  
  cells = items*individuals*reps;
  print.sentence("Total cells to simulate : ",cells,verbose=verbose)
  groups=ceiling(cells/10000000)
  gsize = floor(individuals/groups)
  lsize = individuals - (groups*gsize);
  if(lsize==0) lsize=gsize;
  oind = individuals;
  if(groups == 1) lsize=individuals;
  print.sentence("Groups in simulation : ",groups, "size : ", gsize, verbose=verbose , " last group size : ", lsize)
  
  model=irtpp.model(model)
  #TODO Implement multidimensional test simulation
  dims=1
  ret = NULL;
  ret$model = model;
  #set the seed if not set
  seed = ua(seed,floor(runif(1)*10000000))
  set.seed(seed);
  ret$seed = seed;
  check.model(model);
  #Generate the item parameters (or read)
  z = ua(itempars,simulateItemParameters(items,model,dims,boundaries));
  ret$itempars = z;
  #Generate the individual parameters (assume normal for now, change later)
  th=rnorm(individuals*dims,0,1)
  th=(th-mean(th))/sd(th)
  th = matrix(th,ncol=dims);
  ret$latentTraits = ua(latentTraits,th)
  th=NULL;
  gc()
  individuals = gsize;
  ##Break the simulation here in files.
  if(dirflag){setwd(dir=directory)}
  for (j in 1:reps){
  fname = paste0(filename,j,".csv");
  if(file.exists(fname)){
    print.sentence("Deleting file",fname);
    file.remove(fname)}
  }
  
  
  for (i in 1:groups){
    
    if(i == groups){
      individuals = lsize
    }
    ##Select only the thetas in this file reduction.
    #print.sentence("Length of traits",length(ret$latentTraits))
    b1 = ((i-1)*gsize)+1
    b2 = (((i-1)*gsize))+individuals
    if(i == groups) b2 = oind;
    th = ret$latentTraits[b1:b2];
    #print.sentence(b1," : ",b2)
    
    
    if(verbose){print("Starting simulation")}
    ##Here th must be exactly the latent traits of these individuals in this test.
    ret$prob=replicate(reps,do.call(rbind,lapply(th,function(x,z) probability.3pl(theta=x,z=z),z=z)),simplify=F)
    gc()
    coins=replicate(reps,matrix(runif(items*individuals),ncol=items),simplify=F);
    gc()
    ret$test=mapply(function(x,y){ifelse(x>y,1,0)},ret$prob,coins,SIMPLIFY=F);
    coins=NULL
    gc()
    if(verbose){print("Simulation finished ... ")}
    #Impute the test to exclude individuals in the threshold
    
    repeat{
      #scores of individuals and items
      if(verbose){print("")}
      individual.scores = lapply(ret$test,function(x) {
        rowSums(x)/items});
      #outlier scores
      outliers.flags = lapply(individual.scores,function(x) ifelse(x<threshold | x>(1-threshold),T,F))
      outliers.indices = lapply(outliers.flags,function(x) as.list(which(x)))
      outliers.missing = lapply(outliers.indices,length)
      outliers.total = Reduce(sum,outliers.missing)
      if(outliers.total<2){
        print.sentence("No outliers left",verbose=verbose)
        gc()
        break
      }
      else{
        gc()
        print.sentence("Outliers left",outliers.total,verbose=verbose)
      }
      #resimulate the coins
      if(verbose){print("Resimulating coins ...")}
      newcoins = sapply(outliers.missing,function(x){matrix(runif(x*items),ncol=items)},simplify=F)
      probs = mapply(function(x,y){x[as.numeric(y),]},ret$prob,outliers.indices,SIMPLIFY=F)
      if(verbose){print("Calculating new scores ...")}
      newscores=mapply(function(x,y){ifelse(x>y,1,0)},probs,newcoins,SIMPLIFY=F);
      #assign the new scores in the the old test
      if(verbose){print("Re-assigning new scores ...")}
      mapply(function(x,y,z){
        if(outliers.missing[[z]]>0){
          mapply(function(a,b){
            ret$test[[z]][a,]<<-y[b,]
          },x,seq(length(x)),SIMPLIFY=F);
        }
      },outliers.indices,newscores,seq(length(outliers.indices)),SIMPLIFY=F);
    }
    gc()
    ##Handle the return
    ##Output to the file.
    print.sentence("Outputting to files ",verbose=verbose)
    for (j in 1:reps){
      fname = paste0(filename,j,".csv");
      write.table(ret$test[[j]],file=fname,append=T,sep=",",col.names=F,row.names=F)
    }
    print.sentence("... ",verbose=verbose)
    
  }
  ##Delete the files if option is not enabled.
  ##top("hammertime")
  if(dirflag){
    ret$test=directory;
  }
  else {
    for (j in 1:reps){
        fname = paste0(filename,j,".csv");
        path = paste0(getwd(),"/",fname)
        ts=read.table(file=fname,header=F,sep=",")
        ts = data.matrix(ts,rownames.force=F)
        ret$test[[j]]=ts;
        ts=NULL;
        gc()
        file.remove(fname);
    }
  }
  #Generate the tests
  ret$prob = NULL;
  gc();
  ret
}


#' Simulates item parameters depending on a model
#' @name SimulateItemParameters
#' @description Simulates item parameters for UIRT models.
#' @export
#' @param items , Number of items to generate
#' @param model A string with the model to simulate, please refer to the model documentation in irtpp documentation.
#' @param dims Optional. The number of dimensions to simulate in the test if the model is multidimensional
#' @param boundaries Optional. The kind of boundaries that are specified for the parameters.
simulateItemParameters<- function(items, model, dims=1, boundaries=NULL){
  model = irtpp.model(model);
  bd = boundaries;
  bd$b_lower = ua(bd$b_lower,-4); 
  bd$b_upper = ua(bd$b_upper,4); 
  bd$a_upper = ua(bd$a_upper,5); 
  bd$a_lower = ua(bd$a_lower,0.0001); 
  bd$c_upper = ua(bd$c_upper,0.35); 
  bd$c_lower = ua(bd$c_lower,0);
  if(dims == 1){
  
  b = rnorm(items);
  if(model == "3PL"){
    a = rlnorm(items,meanlog=0,sdlog=1/4)
    c = runif(items,min=bd$c_lower,max=bd$c_upper)
  }
  if(model == "2PL"){
    a = rlnorm(items,meanlog=0,sdlog=1/4)
    c = rep(0,items)
  }
  if(model == "Rasch"){
    temp = rlnorm(1,meanlog=0,sdlog=1/4)
    a = rep(temp,items)
    c = rep(0,items) 
  }
  if(model == "1PL"){
    a = rep(1,items)
    c = rep(0,items)
  }
  ret = list(a=a,b=b,c=c);
  }
  if(dims > 1){
    ##Multidimensional
    a=matrix(runif(dims * items, min = 0, max = 7), nrow = items) #a_j
    b=rnorm(items, mean = 0, sd = 0.7)#b
    d=rep(NA, items)
    for (i in 1:items){                                                                                                                                                                                                                   
      d[i]=-b[i]*sqrt(sum(a[i, ]^2))#d
    }
    c <- runif(items, min = 0, max = 0.25) #c
    
    if(model=="2PL"){c=rep(0,items)}
    if(model=="1PL"){c=rep(0,items);a=matrix(1,ncol = dim,nrow = items)}
    if(model=="3PL"){c=c;a=a}
    ret = list(a=a,d=d,c=c);
  }
  ret$model = model;
  ret
}

#FUNCION PARA CALCULAR LA PRBABILIDAD EN EL CASO MULTIDIMENSIONAL

#' prob.

#' @description calcula la probabilidad en el caso multidimensional
#' @author Juan David Cortés
#' @return la probabilidad
#' @param theta: eltrazo latente multidimensional 
#' @param a: parametro de discrimnacion en el caso multidimensional
#' @param d: parametro d en el caso multidimensional
#' @param c: parametro de azar en el caso multidimensional

prob <- function(theta,a,d,c)
{
  prob <- c + (1 - c)/(1 + exp(-(sum(theta*a)+d)))
  return(prob)
}

#FUNCION PARA SIMULAR EL TEST

#' testmulti

#' @description simula un test multidimensional
#' @author Juan David Cortés
#' @return el test y los parámeros poblacionales con los que se simuló
#' @param nitems: numero de items
#' @param ninds: numero de individuos
#' @param dim: dimension del trazo
#' @param model: 1pl, 2pl o 3pl

testmulti=function(nitems,ninds,dim,model){
  
  #GENERACIÓN DE PARÁMETROS
  ret = NULL;
  a=matrix(runif(dim * nitems, min = 0, max = 7), nrow = nitems) #a_j
  b=rnorm(nitems, mean = 0, sd = 0.7)#b
  d=rep(NA, nitems)
  for (i in 1:nitems){
    d[i]=-b[i]*sqrt(sum(a[i, ]^2))#d
  }
  c <- runif(nitems, min = 0, max = 0.25)#c
  
  
  mu=rep(0,dim)
  sigma=matrix(0,dim,dim)
  corr=runif(((dim^2-dim)/2),min = 0,max = .6)
  sigma[lower.tri(sigma)]=corr
  sigma=t(sigma)+sigma
  diag(sigma)=1
  
  
  theta <- mvrnorm(n = ninds, mu = mu, Sigma = sigma)
  
  #PROBABILIDAD
  
  if(model=="2PL"){c=rep(0,nitems)}
  if(model=="1PL"){c=rep(0,nitems);a=matrix(1,ncol = dim,nrow = nitems)}
  if(model=="3PL"){c=c;a=a}
  
  psics <- matrix(NA, nrow = ninds, ncol = nitems)
  for (j in 1:ninds){
    for (i in 1:nitems){
      psics[j, i] <- prob(theta = theta[j, ], a = a[i, ], d = d[i], c = c[i])    
    }
  }
  
  #TEST
  
  test=matrix(NA, nrow = ninds, ncol = nitems)
  for (j in 1:ninds){
    for (i in 1:nitems){
      test[j, i] <- ifelse(runif(1)<psics[j,i],1,0)    
    }
  }
  

  param=list("a"=a,"b"=b,"c"=c,"d"=d)
  lista=list("test"=test,"param"=param)
  return(lista)
}


#' @name Simulate test Multidimensional
#' @author SICS Research team
#' @title Simulate test Multidimensional
#' @param items. The items to simulate 
#' @param individuals. The individuals that answer the test.
#' @param dims. The dimensions of the test.
#' @param seed. A seed for reproducibility.
#' @param z. A parameter list to feed the function.
#' @param repetition. A number of repetition to simulate with the same item parameters different tests.
#' @return 
#' This function returns the following data in a named list :
#' \itemize{
#'   \item test : The simulated test
#'   \item z : A IRT parameter list.
#'   \item clusters : The list of cluster sizes per cluster.
#'   \item direction : The directions in the dimensions for each cluster.
#'   \item clustinit : The principal item of each cluster
#'   \item clusterlist : The range list of the clusters.
#' }
#' @export
simulateTestMD <- function(items = 10, individuals = 1000, dims = 3, clusters = 4 , seed = 10, z = NULL , repetition=1)
{
  ####Start here
  ### Decide number of items per cluster.
  set.seed(seed)
  rem = items%%clusters
  nitems = items - rem
  itemlist = rep(nitems/clusters,clusters)
  itemlist[[clusters]] = itemlist[[clusters]] + rem
  ##split
  #print(itemlist)
  
  ##determinar direcciones principales
  idnoisy = diag(dims)+matrix(rnorm(dims*dims,0.15,0.05),nrow=dims,ncol=dims);
  idnoisy = idnoisy * ifelse(idnoisy < 1 , 1, 0) + diag(dims)
  idnoisy = normalize(idnoisy)
  beta_w =rbind(idnoisy,matrix(rnorm(dims*(clusters-dims),0.3,0.05),nrow = (clusters-dims), dims))
  beta_w
  
  noise <- seq(0.1,by=.25 / clusters, length.out=clusters)
  
  
  # Matrix of the beta directions
  dir_beta = matrix(NA,sum(itemlist),dims)
  
  # seed for reproducible experiments
  set.seed(seed)
  
  #Perturb the dims
  ##List item limits
  ends = cumsum(itemlist)
  inits = (c(0,cumsum(itemlist))+1)[1:length(itemlist)]
  dir_beta[items,]
  i = 1
  for (i in 1:clusters) {
    dir_beta[inits[i]:ends[i],] = (matrix(beta_w[i,], itemlist[i], dims, byrow=TRUE) 
                                   + matrix(runif(itemlist[i]*dims,-noise[i],noise[i]), itemlist[i],dims))
  }
  
  dir_beta = normalize(dir_beta)
  
  ## Quitar negativos
  dir_beta = ifelse(dir_beta<0,0,dir_beta)
  
  ## Definir vectores de identificacion.
  dir_beta[inits[1:dims],] = diag(dims)
  
  
  ##True reference directions
  true_w = matrix(NA,clusters,dims)
  for (i in 1:clusters) {
    true_w[i,] <- abs(eigen(t(dir_beta[inits[i]:ends[i],])%*%dir_beta[inits[i]:ends[i],])$vectors)[,1] 
  }
  
  ### Now simulate itempars
  if(is.null(z)){
    l_a=0.25
    
    u_a = Inf
    Alpha <- rtnorm(items, mean = 0, sd = 1.0, lower = l_a,  upper = u_a)#genera los alphas
    Alpha[inits] = 1
    length(Alpha)
    a = dir_beta * matrix(Alpha, items,dims, byrow=FALSE)
    
    #sqrt(t(dir_beta[2,])%*%(dir_beta[2,]))
    
    # B parametters
    sd.gamma <-1
    Gamma <- rnorm(items,0,sd.gamma)
    Gamma <- Gamma*Alpha
    Gamma[inits[1:dims]] = 0;
    ## C parameters
    guessing=runif(items,0,0.25)
  }
  else {
    a = z$a;
    Gamma = z$d;
    guessing = z$c;
  }
  
  theta <- matrix(rnorm(individuals*dims,0,1),individuals,dims)
  # this line is to garantee the the covariance matrix is the identity
  theta <- theta %*% solve((chol(cov(theta))))
  cov(theta)
  theta.true <- theta
  
  ## Setting prob matrix
  eta  <- theta %*% t(a) -  matrix(Gamma,individuals,items,byrow=TRUE)
  P = guessing + (1-guessing)/(1+exp(-eta))
  ## Coins
  coinseed = seed;
  if(!is.null(repetition)){
    coinseed = repetition*3 + seed;
  }
  set.seed(as.integer(coinseed));
  nna = (repetition*repetition*coinseed)%%100
  dd = runif(nna);
  U <- matrix(runif(items*individuals),individuals,items);
  
  responses <- ifelse(P>U,1,0)
  
  t.score = rowSums(P);
  c.score = rowSums(responses);
  
  cor(t.score,c.score)
  
  restrict = matrix(0,items,dims+2)
  ###restrict a's of the first dims. and the b's
  for (i in 1:ncol(restrict)) {
    for (j in 1:nrow(restrict)) {
      if(i <= dims && j <= dims + 1){
        restrict[i,j]=1;
      }
    }
  }
  ##Return everything
  
  upper = inits;
  lower = inits + itemlist -1; 
  
  ##Estos son los clusters para fix.items
  
  fixedclusts  = c(t(matrix(c(upper,lower),nrow = length(itemlist))))
  
  return (list("test"= responses,"z" = list("a"=a,"d"=Gamma,"c"=guessing),"clusters"=itemlist,"direction" = beta_w,
               "clustinit"=inits, "coinseed"= coinseed, "clusterlist"= fixedclusts
  ))
}