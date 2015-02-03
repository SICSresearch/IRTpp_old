#Item.fit de SICS
item.fit.sics = function(pats,zita,theta,G = 10,FUN = median,p.val.sim = FALSE,boot.num = 100){
  frec = pats[,ncol(pats)]
  nitems = nrow(zita)
  patsSinFrec = pats[,-ncol(pats)]
  groups = quantile(theta,seq(0, 1, length.out = G + 1))
  groups[1] = groups[1] - 0.1
  groups[G + 1] = groups[G + 1] + 0.1
  groups.ind = findInterval(theta,groups)
  groups.ind = factor(groups.ind, levels = unique(groups.ind))
  thetaG = tapply(rep(theta, frec), rep(groups.ind, frec), FUN = FUN)
  pr = matrix(NA,ncol = nitems,nrow = G)
  for(i in 1:G){
    for(j in 1:nitems){
      pr[i,j] = gg(a = zita[j,1],d = zita[j,2],cp = qlogis(zita[j,3]),theta = thetaG[i])
    }
  }
  Nj = as.vector(tapply(frec, groups.ind, sum))
  Obs = rowsum(frec * patsSinFrec, groups.ind, reorder = FALSE)/rep(Nj,nitems)
  chi.square = Nj * (Obs - pr)^2/(pr * (1 - pr))
  Tobs = colSums(chi.square, na.rm = TRUE)
  df = G - 3 
  pvals <- pchisq(Tobs, df = df, lower.tail = FALSE)
  salida = matrix(c(Tobs,pvals),ncol= 2)
  salida
}

item.fit.sics(pats,est$zita,thetas,G=5)



#Estadística Orlando SICS
item.fit.sics = function(pats,zita,theta,G = 41,FUN = median,p.val.sim = FALSE,boot.num = 100){
  frec = pats[,ncol(pats)]
  nitems = nrow(zita)
  patsSinFrec = pats[,-ncol(pats)]
  w.cuad = theta[,2]
  thetaG = theta[,1]
  
  pr = matrix(NA,ncol = nitems,nrow = G)
  for(i in 1:G){
    for(j in 1:nitems){
      pr[i,j] = gg(a = zita[j,1],d = zita[j,2],cp = qlogis(zita[j,3]),theta = thetaG[i])
    }
  }
  
  #Matriz de S_k para cada thetaG
  sact = matrix(0,ncol = nitems +1,nrow = G)
  for(m in 1:G){
    sant = rep(0,nitems)
    sant[1] = 1 - pr[m,1]
    sant[2] = pr[m,1]
    
    #sk = (1 - pr[2,1]) * s0m
    for(k in 2:(nitems-1)){
      #s0
      sact[m,1] = (1-pr[m,k]) * sant[1]
      for(kk in 2:k){
        #sk hasta k-1
        sact[m,kk] = pr[m,kk] * sant[kk-1] + (1 - pr[m,kk]) * sant[kk]
      }
      sact[m,k+1] = pr[m,k] * sant[k]
      sant = sact[m,]
    }
    sact[m,nitems + 1 ] = pr[m,nitems] * sant[nitems-1] + (1 - pr[m,nitems]) * sant[nitems]
  }
  
  #sact = sact[,-1]
  
  #Matriz de S_k para cada thetaG con items faltantes
  #dimensión x = s_k
  #dimensión y = grupo G del theta
  #dimensión z = calculo de matriz con item z faltante
  #donde dim = dim(y,x,z)
  nitems2 = nitems - 1
  sactInc = array(0,dim = c(G,nitems,nitems))
  for(j in 1:nitems){
    prInc = pr[,-j]
    for(m in 1:G){
      sant = rep(0,nitems2)
      sant[1] = 1 - prInc[m,1]
      sant[2] = prInc[m,1]
      #sk = (1 - pr[2,1]) * s0m
      for(k in 2:(nitems2-1)){
        #s0
        sactInc[m,1,j] = (1-prInc[m,k]) * sant[1]
        for(kk in 2:k){
          #sk hasta k-1
          sactInc[m,kk,j] = prInc[m,kk] * sant[kk-1] + (1 - prInc[m,kk]) * sant[kk]
        }
        sactInc[m,k+1,j] = prInc[m,k] * sant[k]
        sant = sactInc[m,,j]
      }
      sactInc[m,nitems,j] = prInc[m,nitems2] * sant[nitems2-1] + (1 - prInc[m,nitems2]) * sant[nitems2]
    }
    #sact = sact[,-1]
  }
  #print(sact)
  #print(sactInc)
  
  Denom = colSums(matrix(rep(w.cuad,nitems -1 ),ncol = nitems - 1) * sact[,-c(1,ncol(sact))])
  print(Denom)
  
  #print(Denom)
  #print(pr)
  #print(c(sactInc[,1,1]))
  Eik = matrix(0,ncol = nitems -1,nrow = nitems)
  for(i in 1:nitems){
    for(j in 1:(nitems - 1)){
      Eik[i,j] = sum(pr[,i] * sactInc[,j,i] * w.cuad) / Denom[j]
      print(pr[,i] * sactInc[,j,i] * w.cuad)    
    }
  }
  print(pr > 1)
  print(Eik)
  score = rowSums(patsSinFrec )
  
  Oik = matrix(0,ncol = nitems -1 ,nrow = nitems)
  for(i in 1:nitems - 1){
    inds = print(which(score == i))
    for(j in 1:(nitems)){
      patsCoin = pats[inds,]
      if(class(patsCoin) == "matrix"){
        if(dim(patsCoin)[1] != 0){
          Oik[j,i] = sum(apply(X = patsCoin,MARGIN = 1,FUN = function(x){ifelse(x[j] == 1,yes = x[nitems + 1],0)})) / i
        }else{
          Oik[j,i] = 0
        }        
      }else{
        if(class(patsCoin) == "numeric"){
          Oik[j,i] = ifelse(patsCoin[j] == 1,yes = patsCoin[nitems + 1],0) / i
        }
      }
    }
  }
  
  Oik = Oik / sum(pats[,nitems + 1 ])
  print(Oik)
  SXCuad = rowSums((Oik - Eik)^2 / (Eik * (1-(Eik))))
  SXCuad
}

#est$pats
#scoresSics

pt.cuad = read.table("/home/mirt/Trabajo IRT/Algoritmo SICS/PWcuad.csv",dec=".",sep = " ",header = T)


inicio = Sys.time()
a=item.fit.sics(pats = est$pats,zita = est$zita,theta = pt.cuad,G=41)
Sys.time() - inicio
sum(a)
a
