
#' Orlando
#' Calcula la estadística de Orlando
#' @param patterns: matriz con los patrones de respuesta, las frecuencias y los trazos
#' @param G numero de puntos de cuadratura
#' @param zita: MATRIZ de estimaciones de parametros de los items (a,b,c)
#' @param FUN : funcion de representante de cada grupo
#' @return estadística de Orlando

#####objetos de sics para ir probando:

orlando_itemf=function(patterns,G,zita,FUN){
  
  pats = patterns[,-ncol(patterns)]   #patrones sin frecuencia
  frec = patterns[,"x"]     #fr. de los patrones de respuesta
  patsSinFrec = pats[,-ncol(pats)]   #patrones sin frecuencia
  nitems=ncol(patsSinFrec)
  
  
  seq=seq(-6,6,length=61)
  pesos=dnorm(seq)/sum(dnorm(seq))
  Cuad=matrix(c(seq,pesos),byrow=F,ncol=2)# cuadraturas
  
  theta=Cuad
  w.cuad = theta[,2] #pesos
  thetaG = theta[,1] #nodos
  
  zita=parameter.list(zita,"3PL")  
  pr = lapply(thetaG,FUN=function(x){probability.3pl(theta=x,z=zita)}) #probabilidad para cada punto de cuadratura
  pr=matrix(unlist(pr),ncol=nitems,byrow=T)
  
  ### TOTALES POR SCORE
  
  score = rowSums(patsSinFrec )  #suma de los patrones(scores sin agrupar)
  Nk=NULL
  for(i in 1:nitems - 1){  #recorriendo los scores(fijando un score)
    inds = print(which(score == i)) #para agrupar llos patrones q determinen un mismo score (i)
    patsCoin = pats[inds,]   #patrones q determinan el score "i": agrupados
    if(class(patsCoin) == "matrix"){  
      if(dim(patsCoin)[1] != 0){     #(si hay 1 patron o mas)
        Nk[i] = sum(patsCoin[,ncol(patsCoin)])
      }else{
        Nk[i] = 0
      }        
    }else{
      if(class(patsCoin) == "numeric"){
        Nk[i] = patsCoin[length(patsCoin)]
      }
    }
  }
  
  
  ### FRECUENCIAS OBSERVADAS
  
  
  O=list()
  print("nitems")
  print(nitems)
  Oik = matrix(0,ncol = nitems -1 ,nrow = nitems)  # 4 scores(columnas) y nitems(filas)
  for(i in 1:nitems - 1){  #recorriendo los scores(fijando un score)
    inds = print(which(score == i)); #para agrupar los patrones q determinen un mismo score (i)
    for(j in 1:(nitems)){  #recorriendo los items, para llenar la matriz por items fiijado un score (i)
      patsCoin = pats[inds,]   #patrones q determinan el score "i": agrupados
      if(class(patsCoin) == "matrix"){  
        if(dim(patsCoin)[1] != 0){     #(si hay 1 patron o mas)
          Oik[j,i] = sum(apply(X = patsCoin,MARGIN = 1,FUN = function(x){ifelse(x[j] == 1,yes = x[nitems + 1],0)})) 
        }else{
          Oik[j,i] = 0
        }        
      }else{
        if(class(patsCoin) == "numeric"){
          Oik[j,i] = ifelse(patsCoin[j] == 1,yes = patsCoin[nitems + 1],0)
        }
      }
      O[[j]]=cbind(Nk-Oik[j,],Oik[j,])
    }
  }
  
  ### S SIN MOÑO
  
  sact=s_ss(pr,nitems=nitems)
  Denom = colSums(matrix(rep(w.cuad,nitems -1 ),ncol = nitems - 1) * sact[,-c(1,ncol(sact))])  
  
  
  ### S MOÑO
  
  nitems = nitems - 1
  smoño=list()
  
  for(p in 1:(nitems+1)){
    smoño[[p]]=s_ss(pr[,-p],nitems=nitems) 
  }
  
  ### Eik   
  
  nitems=ncol(patsSinFrec)
  E=list()
  for(i in 1:length(smoño)){
    E[[i]]=cbind(1-colSums(smoño[[i]][,-ncol(smoño[[i]])]*(pr[,i]*w.cuad))/Denom,colSums(smoño[[i]][,-ncol(smoño[[i]])]*(pr[,i]*w.cuad))/Denom)
  }
  
  
  
  #Estadística
  
  S_X2=NULL
  df.S_X2=NULL
  for (i in 1:nitems) E[[i]] <- E[[i]] * Nk
  coll <- collapseCells(O, E, mincell = 1)
  O <- coll$O
  E <- coll$E
  for (i in 1:length(O)) {
    S_X2[i] <- sum((O[[i]] - E[[i]])^2/E[[i]], na.rm = TRUE)
    df.S_X2[i] <- sum(!is.na(E[[i]])) - nrow(E[[i]]) - 3
  }
  
  pval=pchisq(S_X2,df.S_X2,lower.tail = F)
  
  lista=cbind("S_X2"=S_X2,"df.SX2"=df.S_X2,"p.val.S_X2"=pval)
  
  return(lista)
  
}
