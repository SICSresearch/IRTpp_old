



### Items Reclassify (determine items to drop, function for internal use.)
items.reclassify<-function(n, tmpd , axis){
  POS.PROY<- vector("list")   ## Vector of higher projection in the submatrix
  CAMBIO<- vector()	
  for(k in 1:(n))
  {
    proyec<- matrix(NA, ncol=dim(tmpd[[k]])[2], nrow=(n))
    for(i in 1:dim(tmpd[[k]])[2])
    {
      for(j in 1:(n))
      {
        proyec[j,i]<-t(tmpd[[k]][,i])%*%axis[,j]
      }
    }
    pos.proy<- vector()
    for(i in 1:dim(tmpd[[k]])[2])
    {
      pos.proy[i]<- which(proyec[,i] == max(proyec[,i]))
    }
    POS.PROY[[k]] <- pos.proy
    ant.cluster<- rep(k, dim(tmpd[[k]])[2])
    CAMBIO[k]<- sum(ifelse(pos.proy == ant.cluster, 1, 0))/length(ant.cluster)
  }
  POS.PROY[[k+1]]<-sum(CAMBIO)
  return(POS.PROY)
  
}



## Make clusters of items.
## Iterative PCA.

### y.datF1 is the data 0's and 1's
#library(FactoMineR)
#library(optimx)
#library(numDeriv)

#library(IRTpp)
#t = simulateTestMD(items = 20,individuals = 1000,dims = 3,clusters = 4)
#test = t$test


#acp = ACPI(test,0.8,0.9,6,15)
#y.datF1 = read.table("/home/liberato/Downloads/Aplicacion_UN/respon.original.txt", header=T) 

#names(acp)
#acp[[1]]
#acp[[2]]
#dim(acp[[3]]) ##Latent traits
#dim(y.datF1)
#acp[[4]] ## Items per each dimension
#length(acp[[5]])
#acp[[5]] ## Submatrices for the test.
#lapply(acp[[5]][[2]],dim)
#typeof(acp[[5]][[2]])
#lapply(acp)
#do.call(cbind,acp[[5]])

##Fuse submatrices




ACPI<-function(y.datF1,h7,h2,hn,angu){
  y.datF1<-y.datF1
  originalF1<- y.datF1
  yc.datF1=y.datF1
  tmp.datF1<-vector('list')
  tempF1<- vector()
  axisF1<- vector('list')
  i=1
  vc7<-40 ## Seventh eigenvalue
  vc1<-100 ## First eigenvalue
  while((vc7/vc1)<=h7){
    vp1<- 100
    vp2<- 91 ## Second eigenvalue
    tmp<- angu ## Determinated angle
    cluster<-y.datF1
    print(dim(cluster))
    c=0
    while((vp2/vp1)>=h2)   #### OJO 0.9
    {
      ACP<-PCA(cluster,ncp=3, graph=F)
      hipot<- sqrt(ACP$var$coord[,1]^2 + ACP$var$coord[,2]^2)  
      corr =  ACP$var$coord[,1]/hipot 
      ang = acos(corr)*180/pi	
      pos.ang.1<-as.numeric(which(ang<=tmp))
      while(length(pos.ang.1)<hn)   ### OJO  6
      {
        tmp=tmp+1
        pos.ang.1<- as.numeric(which(ang<=tmp))
      }
      clusterc<-t(t(y.datF1)[pos.ang.1,])
      print(dim(clusterc))
      ACPclust<-PCA(clusterc, graph=F)
      #plot.PCA(ACPclust, cex=0.6, choix="var")
      vp1<-ACPclust$eig[1,1]
      vp2<-ACPclust$eig[2,1] 
      if((vp2/vp1)>=h2){    #### OJO 0.9
        cluster<-cluster
      }else{
        cluster<-clusterc
      }
      tmp=tmp-1
      print(tmp)
      c=c+1
    }
    tempF1[i] = tmp
    axisF1[[i]] <- as.numeric(ACP$ind$coord[,1])
    axisF1[[i]] <- axisF1[[i]]/sd(axisF1[[i]])	
    tmp.datF1[[i]]<-cluster
    print(dim(tmp.datF1[[i]]))
    opos.datF1<- t(t(y.datF1)[-pos.ang.1,])
    print(dim(opos.datF1))
    ACP1<-PCA(opos.datF1,ncp=3, graph=F)
    vc1<-ACP1$eig[1,1]
    vc7<-ACP1$eig[7,1] 
    vc7/vc1
    i<-i+1
    y.datF1<-opos.datF1
    print(i)
  }
  
  ## Number of submatrix
  
  n.clustF1I<-i-1
  n.clustF1<-(i-1)
  
  ## Size of each submatrix
  
  sF1<-vector()
  for (i in 1:n.clustF1){
    sF1[i]<-dim(tmp.datF1[[i]])[2]
  }
  sF1
  
  ## Pure Axis
  
  AXISF1<- matrix(NA, nrow=dim(originalF1)[1], ncol=n.clustF1)
  
  for(i in 1:(n.clustF1))
  {
    AXISF1[,i]<- axisF1[[i]]
  }
  
  ############################################
  #############     RECLASSIFICATION 
  ############################################
  
  POS.PROYF1<-items.reclassify(n.clustF1, tmp.datF1, AXISF1)
  o<-POS.PROYF1[[n.clustF1+1]]
  
  while(o!=(n.clustF1)){
    
    POS.PROY<-vector('list')
    for(i in 1:(n.clustF1)){
      POS.PROY[[i]]=POS.PROYF1[[i]]
    }
    
    ######## Add elements of the other submatrix
    
    tmp.dat1<- tmp.datF1
    POS.PROY1<- POS.PROYF1
    for(i in 1:(n.clustF1))
    {
      for(k in 1:(n.clustF1))
      {
        if(i!=k)
        {
          t<- tmp.dat1[[k]][,POS.PROY1[[k]]==i]
          if(is.vector(t)==TRUE)
          {
            m<-which(POS.PROY1[[k]]==i)
            t<-as.matrix(t)
            colnames(t)<-colnames(head(tmp.dat1[[k]]))[m]
            rownames(t)<-rownames(tmp.dat1[[k]])
            tmp.dat1[[i]]<-cbind(tmp.dat1[[i]], t)
            POS.PROY1[[i]]<- c(POS.PROY1[[i]], i)
          }else
          {
            if(dim(t)[2]!=0)
            {
              tmp.dat1[[i]]<-cbind(tmp.dat1[[i]], tmp.dat1[[k]][,POS.PROY1[[k]]==i])
              POS.PROY1[[i]]<-  c(POS.PROY1[[i]],rep(i,dim(t)[2]))
            }
          }
        }
      }
    }
    
    ######## Remove elements that not belong to the submatrix
    
    for(i in 1:(n.clustF1)){
      tmp.dat1[[i]]<- tmp.dat1[[i]][,POS.PROY1[[i]]==i]
      POS.PROY1[[i]]<- POS.PROY1[[i]][POS.PROY1[[i]]==i] 
    }
    tmp.dat1.fin<-vector("list")
    POS.PROY1.fin<-vector("list")
    k=0
    for (i in 1:(n.clustF1)){
      if(length(POS.PROY1[[i]])==0||length(POS.PROY1[[i]])==1){
      }else{
        k=k+1
        POS.PROY1.fin[[k]]<-rep(k,length(POS.PROY1[[i]]))
        tmp.dat1.fin[[k]]<-tmp.dat1[[i]]
      }
    }
    n.clustF1<-k
    
    AXIS.tmp<- matrix(NA, nrow=dim(originalF1)[1], ncol=(n.clustF1))
    for(i in 1:(n.clustF1)){
      AXIS.tmp[,i]=as.numeric(PCA(tmp.dat1.fin[[i]], graph=F)$ind$coord[,1])
      AXIS.tmp[,i] <- AXIS.tmp[,i]/sd(AXIS.tmp[,i])
    }
    tmp.datF1<-tmp.dat1.fin
    AXISF1<-AXIS.tmp
    POS.PROYF1<-items.reclassify(n.clustF1, tmp.datF1, AXISF1)
    o<-POS.PROYF1[[n.clustF1+1]]
  }
  
  ###### End of the reclassification
  
  ## Number of final submatrix
  
  ##n.clustF1
  
  #### Items in each submatrix
  
  sF1.1<-vector()
  for(i in 1:n.clustF1){
    sF1.1[i]<-dim(tmp.datF1[[i]])[2]
  }
  
  ## Write finally axis of individuals
  
  Nomarti<-rownames(tmp.datF1[[1]])
  colnames(AXISF1)<- paste("v",seq(1,n.clustF1, by=1))
  row.names(AXISF1)<-Nomarti
  
  ##y<-round(runif(6,1,408),0)
  ##xtable(head(AXIS[y,]),digits=4)
  
  
  ## Write the items in each submatrix
  
  items.subm<-matrix(NA, nrow=max(sF1.1), ncol=n.clustF1)
  for(i in 1:n.clustF1){
    items.subm[(1:dim(tmp.datF1[[i]])[2]),i]<-colnames(tmp.datF1[[i]])
  }
  
  #SFin<-matrix(NA,ncol=2,nrow=length(sF1))
  #SFin[,1]<-sF1
  #SFin[,2]<-
  
  tore<-vector('list')
  tore[[1]]=c(n.clustF1I, n.clustF1)
  tore[[2]]=c(sF1, sF1.1)
  tore[["traits"]]<-AXISF1
  tore[["itemdims"]]<-items.subm
  tore[["submatrices"]]<-tmp.datF1
  print(sF1)
  print(sF1.1)
  return(tore)
}



