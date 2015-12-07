#data<-read.table("respon.original.txt", header=T)
#head(data)
#head(tes$test)
#IPCA(data, clas=FALSE, reclas=FALSE)



############################################
#############  FUNCTION   CHANGE OF RECLASSIFICATION 
############################################

change<-function(n, tmpd , axis){
  
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

############################################
#############  FUNCTION   INITIAL
############################################



PCAprin<-function(y.datF1,h7,h2,hn,angu){
  opos.datF1 = NULL;
  y.datF1<-y.datF1
  originalF1<- y.datF1
  yc.datF1=y.datF1
  tmp.datF1<-vector('list')
  tempF1<- vector()
  axisF1<- vector('list')
  i=1
  vc7<-40 ## Seventh eigenvalue
  vc1<-100 ## First eigenvalue
  tmpv=c(angu)
  #repeat{
  while((vc7/vc1)<=h7){
    vp1<- 100
    vp2<- 91 ## Second eigenvalue
    tmp<- angu ## Determinated angle
    cluster<-y.datF1
    c=0
    while((vp2/vp1)>=h2)   
    {
      ACP<-PCA(cluster, ncp=3, graph=F)
      hipot<- sqrt(ACP$var$coord[,1]^2 + ACP$var$coord[,2]^2)  
      corr =  ACP$var$coord[,1]/hipot 
      ang = acos(corr)*180/pi	
      pos.ang.1<-as.numeric(which(ang<=tmp))
      while(length(pos.ang.1)<hn)  
      {
        tmp=tmp+1
        pos.ang.1<- as.numeric(which(ang<=tmp))
      }
      clusterc<-t(t(y.datF1)[pos.ang.1,])
      ACPclust<-PCA(clusterc, graph=F)
      #plot.PCA(ACPclust, cex=0.6, choix="var")
      vp1<-ACPclust$eig[1,1]
      vp2<-ACPclust$eig[2,1] 
      if((vp2/vp1)>=h2){    
        cluster<-cluster
      }else{
        cluster<-clusterc
      }
      tmpv<-c(tmpv,tmp)
      #print(tmp)
      tmp=tmp-1
      c=c+1
      if(length(tmpv)>=5&&sum(diff(tmpv[(length(tmpv)-5):length(tmpv)])^2)==0){
        result1<-vector('list')
        result1[[1]]<-i
        result1[[2]]<-tmp.datF1
        result1[[3]]<-axisF1
        result1[[4]]<-opos.datF1
        return (result1)
        stop("Error Omit")
      }
    }
    tempF1[i] = tmp
    axisF1[[i]] <- as.numeric(ACP$ind$coord[,1])
    axisF1[[i]] <- axisF1[[i]]/sd(axisF1[[i]])	
    tmp.datF1[[i]]<-cluster
    opos.datF1<- t(t(y.datF1)[-pos.ang.1,])
    ACP1<-PCA(opos.datF1,ncp=3, graph=F)
    vc1<-ACP1$eig[1,1]
    vc7<-ACP1$eig[hn,1] 
    vc7/vc1
    i<-i+1
    y.datF1<-opos.datF1
  }
  result1<-vector('list')
  result1[[1]]<-i
  result1[[2]]<-tmp.datF1
  result1[[3]]<-axisF1
  result1[[4]]<-opos.datF1
  return (result1)
}



############################################
#############  FUNCTION   CLUSTERING 
############################################

class<-function(opos.datF1, AXISF1, n.clustF1, tmp.datF1){
  cor.matrix.nois<-cor(opos.datF1, AXISF1)
  axisF1<-matrix(NA, nrow=nrow(opos.datF1), ncol=n.clustF1)
  for(i in 1:ncol(opos.datF1)){
    asig<-sort(cor.matrix.nois[i,], index.return=T)$ix[n.clustF1]
    nam.p<-colnames(tmp.datF1[[asig]])
    tmp.datF1[[asig]]<-cbind(tmp.datF1[[asig]], opos.datF1[,i])
    colnames(tmp.datF1[[asig]])<-c(nam.p,colnames(opos.datF1)[i])
  }
  ACP=PCA(tmp.datF1, graph=F)
  for(i in 1:n.clustF1){
    axisF1[,i]<- as.numeric(ACP$ind$coord[,1])
    axisF1[,i] <- axisF1[,i]/sd(axisF1[,i])	
  }
  result3<-vector('list')
  result3[[1]]<-tmp.datF1
  result3[[2]]<-axisF1
  return(resutl3)
}



############################################
#############  FUNCTION   RECLASSIFICATION 
############################################



reclass<-function(n.clustF1, tmp.datF1, AXISF1, originalF1){
  originalF1=originalF1
  tmp.datF1=tmp.datF1
  AXISF1=AXISF1
  n.clustF1<-n.clustF1
  
  
  POS.PROYF1<-change(n.clustF1, tmp.datF1, AXISF1)
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
        name<-paste("Submatrix", k, sep="")
        POS.PROY1.fin[[k]]<-rep(k,length(POS.PROY1[[i]]))
        tmp.dat1.fin[[name]]<-tmp.dat1[[i]]
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
    POS.PROYF1<-change(n.clustF1, tmp.datF1, AXISF1)
    o<-POS.PROYF1[[n.clustF1+1]]
  }
  return2<-vector('list')
  return2[[1]]=n.clustF1
  return2[[2]]=tmp.datF1
  return2[[3]]=AXISF1
  return(return2) 
}


############################################
#############  FUNCTION IMPRETION VV 
############################################

resultshowVV<-function(tore1, tore2.1, tore2.2, tore2.3, tore3, tore4, tore4C, tore5, tore6){
  result<-list(n.submatrix=tore1, n.ipm.befrec=tore2.1, n.ipm.aftrec=tore2.2, n.ipm.clas=tore2.3,
               p.axes=tore3, ipm=tore4, ipmc=tore4C, DPM=tore5, inc=tore6)
  menu<-matrix(NA, nrow=length(result), ncol=2)
  menu[,1]<-paste("$",names(result), sep="")
  menu[,2]<-c("Total of submatrices Before and After reclassification",
              "Number of items per matrix before reclassification", 
              "Number of items per matrix after reclassification",  
              "Number of items per matrix after classification", "Principal Axis of each submatrix", 
              "Matrix of items on each submatrix before classification",
              "Matrix of items on each submatrix after classification", 
              "Data of each submatrix (List object)", "Items not classified in a principle")
  colnames(menu)<-c("Name", "Description")
  rownames(menu)<-seq(1,length(result))
  
  message<-matrix(NA, ncol=1,  nrow=1)
  message[1,]<-c("**Results for the Dimensionality Reduction algorithm based in PCA method**")
  print(message)
  print(menu)
  return(result)
}


############################################
#############  FUNCTION IMPRETION VF
############################################

resultshowVF<-function(tore1, tore2.1, tore2.2, tore3, tore4, tore4C, tore5, tore6, y.clas){
  result<-list(n.submatrix=tore1, n.ipm.befrec=tore2.1, n.ipm.aftrec=tore2.2,
               p.axes=tore3, ipm=tore4, ipmc=tore4C, DPM=tore5, inc=tore6)
  menu<-matrix(NA, nrow=length(result), ncol=2)
  if(y.clas==1){
    menu[,1]<-paste("$",names(result), sep="")
    menu[,2]<-c("Total of clusters",
                "Number of items per matrix before classification", 
                "Number of items per matrix after classification",  
                "Principal Axis of each cluster", "Matrix of items on each cluster before classification",
                "Matrix of items on each cluster after classification", "Data of each cluster (List object)", 
                "Non classified items in a principle")
  }else{
    menu[,1]<-paste("$",names(result), sep="")
    menu[,2]<-c("Total of submatrix",
                "Number of items per matrix before reclassification", 
                "Number of items per matrix after reclassification",  
                "Principal Axis of each submatrix", "Matrix of items on each cluster before reclassification",
                "Matrix of items on each cluster after reclassification", "Data of each submatrix (List object)", 
                "Non classified items")
  }
  colnames(menu)<-c("Name", "Description")
  rownames(menu)<-seq(1,length(result))
  
  message<-matrix(NA, ncol=1,  nrow=1)
  message[1,]<-c("**Results for the Dimensionality Reduction algorithm based in PCA method**")
  print(message)
  print(menu)
  return(result)
}


############################################
#############  FUNCTION IMPRETION FF
############################################

resultshowFF<-function(tore1, tore2.1, tore3, tore4, tore5, tore6){
  result<-list(n.submatrix=tore1, n.ipm.befrec=tore2.1,
               p.axes=tore3, ipm=tore4, DPM=tore5, inc=tore6)
  menu<-matrix(NA, nrow=length(result), ncol=2)
  menu[,1]<-paste("$",names(result), sep="")
  menu[,2]<-c("Total of submatrices",
              "Number of items per matrix", 
              "Principal Axis of each submatrix",  
              "Matrix of items on each submatrix", "Data of each submatrix (List object)", 
              "Non classified items")
  colnames(menu)<-c("Name", "Description")
  rownames(menu)<-seq(1,length(result))
  
  message<-matrix(NA, ncol=1,  nrow=1)
  message[1,]<-c("**Results for the Dimensionality Reduction algorithm based in PCA method**")
  print(message)
  print(menu)
  return(result)
}


#' @title IPCA
#' @name IPCA
#' @description Iterative Principal Components Analysis (IPCA)  Computes the Iterative components analysis to find dimensionality of a huge data matrix with entries 0 and 1.
#' @param test: Dichotomic matrix data. 
#' @param h7: Threshold to define the items not classificated
#' @param h2: Threshold to define the unidimensionality 
#' @param hn: number of items in each cluster
#' @param angu: Threshol to define the angle of each cluster (it could change)
#' @param clas: If is TRUE then the classification is made, default is TRUE 
#' @param reclas: If is TRUE then the reclassification is made, default is TRUE
#' If clas=TRUE and reclas=TRUE then the results are:
#' @return n.submatrix: Total of submatrices Before and After reclassification 
#' @return n.ipm.befrec: Number of items per matrix before reclassification. If clas=TRUE and reclas=FALSE or clas=FALSE and reclas=FALSE, then this result is NULL.
#' @return n.ipm.aftrec: Number of items per matrix after reclassification. If clas=TRUE and reclas=FALSE or clas=FALSE and reclas=FALSE, then this result is NULL.
#' @return n.ipm.clas: Number of items per matrix after classification. If clas=FALSE and reclas=TRUE or clas=FALSE and reclas=FALSE, then this result is NULL.
#' @return p.axes: Principal Axis of each submatrix
#' @return ipm: Matrix of items on each submatrix before classification. If clas=FALSE and reclas=TRUE or clas=FALSE and reclas=FALSE, then this result only show the first classification.
#' @return ipmc: Matrix of items on each submatrix after classification. If clas=FALSE and reclas=TRUE or clas=FALSE and reclas=FALSE, then this result is NULL.
#' @return DPM: Data of each submatrix (List object)
#' @return inc: Items not classified in a principle
#' @section  For textual data analysis is recomended to make ortogonalization of axes of each dimension to have interpretation of them.
#' @author SICS Research Group, Universidad Nacional de Colombia \email{ammontenegrod@@unal.edu.co}
#' @export 
#' @importFrom FactoMineR PCA
#'
#' @seealso
#' \code{\link{z3_itemf}}, \code{\link{x2_itemf}}
#'
#' @references
#'
#' Montenegro, A.M. & Ordo\~nez, M.F & P\'aez, S (2015)  A Novel Methodology Based on Item Responde Theory to Find Dimensionality and Emergent Information from Textual Data. \emph{Journal of Infometrix(Submited)}
#'
#' @examples
#'
#'irtpp()
IPCA<-function(data,h7=0.9,h2=0.2,hn=10,angu=15,clas=TRUE,reclas=TRUE){
  
  y.datF1<-data
  
  originalF1<-y.datF1
  
  prim<-PCAprin(y.datF1,h7,h2,hn,angu)
  i<-prim[[1]]
  n.clustF1I=i-1
  tmp.datF1I<-prim[[2]]
  axisF1<-prim[[3]]
  opos.datF1<-prim[[4]]
  
  AXISF1<- matrix(NA, nrow=dim(originalF1)[1], ncol=n.clustF1I)
  
  for(i in 1:(n.clustF1I))
  {
    AXISF1[,i]<- axisF1[[i]]
  }
  
  if(n.clustF1I==1){
   if(clas==TRUE){
     tore4<-colnames(tmp.datF1I)
     tore2.1<-dim(tmp.datF1I)
     tore2.2<-NULL
     tmp.datF1C<-originalF1
     tore1<-n.clustF1I
     ACP=PCA(originalF1, graph=F)
     axisF1=as.numeric(ACP$ind$coord[,1])
     tore3=as.matrix(axisF1)
     colnames(tore3)<-c("Axis")
     rownames(tore3)<-seq(1,nrow(opos.datF1))
     tore4C<-colnames(tmp.datF1C)
     tore5=tmp.datF1C
     tore6=opos.datF1
     resultfin<-resultshowVF(tore1, tore2.1, tore2.2, tore3, tore4, tore4C, tore5, tore6, 1)
   }else{
   }
  }else{
  
  
  sF1<-vector()
  for(i in 1:n.clustF1I){
    sF1[i]<-dim(tmp.datF1I[[i]])[2]
  }
  
  tore2.1<-cbind(sF1)
  colnames(tore2.1)<-c("ipm.befrec")
  rownames(tore2.1)<-paste("Sumb", seq(1,length(sF1)))
  
  tore3<-AXISF1
  
  if(clas==TRUE&&reclas==TRUE){
    
    resrecla<-reclass(n.clustF1I,tmp.datF1I, AXISF1, originalF1)
    ###
    n.clustF1=resrecla[[1]]
    tmp.datF1=resrecla[[2]]
    sF1.1<-vector()
    for(i in 1:n.clustF1){
      sF1.1[i]<-dim(tmp.datF1[[i]])[2]
    }
    
    ####
    AXISF1R=resrecla[[3]]
    Nomarti<-rownames(originalF1)
    colnames(AXISF1R)<- paste("Axis",seq(1,n.clustF1, by=1))
    row.names(AXISF1R)<-Nomarti
    
    #####
    items.subm<-matrix(NA, nrow=max(sF1.1), ncol=n.clustF1)
    for(i in 1:n.clustF1){
      items.subm[(1:dim(tmp.datF1[[i]])[2]),i]<-colnames(tmp.datF1[[i]])
    }
    ####
    resclasi<-class(opos.datF1, AXISF1R, n.clustF1, tmp.datF1)
    tmp.datF1C<-resclasi[[1]]
    AXISF1C<-resclasi[[2]]
    
    
    ####
    sF1.1C<-vector()
    for(i in 1:n.clustF1){
      sF1.1C[i]<-dim(tmp.datF1C[[i]])[2]
    }
    
    items.submC<-matrix(NA, nrow=max(sF1.1C), ncol=n.clustF1)
    for(i in 1:n.clustF1){
      items.submC[(1:dim(tmp.datF1C[[i]])[2]),i]<-colnames(tmp.datF1C[[i]])
    }
    ####
    tore1<-cbind(n.clustF1I, n.clustF1)
    tore2.2<-cbind(sF1.1)
    colnames(tore2.2)<-c("ipm.befrec")
    rownames(tore2.2)<-paste("Sumb", seq(1,length(sF1.1)))
    tore2.3<-cbind(sF1.1C)
    colnames(tore2.3)<-c("ipm.befclas")
    rownames(tore2.3)<-paste("Sumb", seq(1,length(sF1.1C)))
    tore3<-AXISF1C
    tore4<-items.subm
    colnames(tore4)<-paste("IPM", seq(1,length(sF1.1)))
    rownames(tore4)<-seq(1,max(sF1.1))
    ####
    tore4C<-items.submC
    colnames(tore4C)<-paste("IPMC", seq(1,length(sF1.1C)))
    rownames(tore4C)<-seq(1,max(sF1.1C))
    tore5<-tmp.datF1C
    tore6<-opos.datF1
    resultfin<-resultshowVV(tore1, tore2.1, tore2.2, tore2.3, tore3, tore4, tore4C, tore5, tore6)
    ##################################################################################
  }else if(clas==TRUE&&reclas==FALSE){
    recsclasi<-class(opos.datF1, AXISF1, n.clustF1I, tmp.datF1I)
    tmp.datF1C<-recsclasi[[1]]
    AXISF1C<-recsclasi[[2]]
    tore3<-AXISF1C
    colnames(tore3)<-paste("IPMC", seq(1,n.clustF1I))
    rownames(tore3)<-seq(1,nrow(originalF1))
    sF1.1C<-vector()
    for(i in 1:n.clustF1I){
      sF1.1C[i]<-dim(tmp.datF1C[[i]])[2]
    }
    tore1<-c(n.clustF1I)
    tore2.2<-cbind(sF1.1C)
    colnames(tore2.2)<-c("ipm.aftercl")
    rownames(tore2.2)<-paste("Sumb", seq(1,length(sF1.1C)))
    items.subm<-matrix(NA, nrow=max(sF1), ncol=n.clustF1I)
    for(i in 1:n.clustF1I){
      items.subm[(1:dim(tmp.datF1I[[i]])[2]),i]<-colnames(tmp.datF1I[[i]])
    }
    tore4<-items.subm
    items.submC<-matrix(NA, nrow=max(sF1.1C), ncol=n.clustF1I)
    for(i in 1:n.clustF1I){
      items.submC[(1:dim(tmp.datF1C[[i]])[2]),i]<-colnames(tmp.datF1C[[i]])
    }
    tore4C<-items.submC
    dim(tore4C)
    colnames(tore4C)<-paste("IPMC", seq(1,length(sF1.1C)))
    rownames(tore4C)<-seq(1,max(sF1.1C))
    tore5<-tmp.datF1C
    tore6<-opos.datF1
    resultfin<-resultshowVF(tore1, tore2.1, tore2.2, tore3, tore4, tore4C, tore5, tore6,y.clas=1)
    #######################################################################
  }else if(clas==FALSE&&reclas==TRUE){
    resrecla<-reclass(n.clustF1I,tmp.datF1I, AXISF1, originalF1)
    ###
    n.clustF1=resrecla[[1]]
    tmp.datF1=resrecla[[2]]
    sF1.1<-vector()
    for(i in 1:n.clustF1){
      sF1.1[i]<-dim(tmp.datF1[[i]])[2]
    }
    tore1<-cbind(n.clustF1I, n.clustF1)
    tore2.2<-cbind(sF1.1)
    colnames(tore2.2)<-c("ipm.aftrec")
    rownames(tore2.2)<-paste("Sumb", seq(1,length(sF1.1)))
    items.subm<-matrix(NA, nrow=max(sF1), ncol=n.clustF1I)
    for(i in 1:n.clustF1I){
      items.subm[(1:dim(tmp.datF1I[[i]])[2]),i]<-colnames(tmp.datF1I[[i]])
    }
    tore4<-items.subm
    items.submR<-matrix(NA, nrow=max(sF1.1), ncol=n.clustF1)
    for(i in 1:n.clustF1){
      items.submR[(1:dim(tmp.datF1[[i]])[2]),i]<-colnames(tmp.datF1[[i]])
    }
    tore4R<-items.submR
    colnames(tore4R)<-paste("IPMC", seq(1,length(sF1.1)))
    rownames(tore4R)<-seq(1,max(sF1.1))
    tore5<-tmp.datF1
    tore6<-opos.datF1
    resultfin<-resultshowVF(tore1, tore2.1, tore2.2, tore3, tore4, tore4R, tore5, tore6,y.clas=0)
    ############################################################################
  }else{
    tore1<-c(n.clustF1I)
    items.subm<-matrix(NA, nrow=max(sF1), ncol=n.clustF1I)
    for(i in 1:n.clustF1I){
      items.subm[(1:dim(tmp.datF1I[[i]])[2]),i]<-colnames(tmp.datF1I[[i]])
    }
    tore4<-items.subm
    tore5<-tmp.datF1I
    tore6<-opos.datF1
    resultfin<-resultshowFF(tore1, tore2.1, tore3, tore4, tore5, tore6)
  }
  }
  return(resultfin)
}










