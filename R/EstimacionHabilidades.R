#--------------------------Segundo método-------------------------------------------------
#Método EAP

#Probabilidad del patrón
ggpat = function(pat,zita,nitems,cuad){
  p = rep(0,nitems)
  for(i in 1:nitems){
    p[i] = gg(a = zita[i,1],d = zita[i,2],cp = qlogis(zita[i,3]),theta = cuad)
    if(pat[i] == 0){
      p[i] = 1 - p[i] 
    }
  }
  prod(p)
}

pats = est$pats
npats = nrow(est$pats)
nitems = ncol(est$pats) -1

#EAP
thetaEst3 = rep(0,npats)
for(j in 1:npats){
  sumNum = 0
  sumDen = 0
  for(k in 1:41){
    pat = pats[j,][1:nitems]
    sumNum = sumNum + (pt.cuad[k] * w.cuad[k] * ggpat(pat = pat,zita = est$zita,nitems = nitems,cuad = pt.cuad[k]))
    sumDen = sumDen + (w.cuad[k] * ggpat(pat = pat,zita = est$zita,nitems = nitems,cuad = pt.cuad[k]))
  }
  thetaEst3[j] = sumNum / sumDen
}


#Comparaciones con MIrt
pats2 = pats
pats2 = cbind(pats2,thetaEst3)
pats2

thetaEstMirt = fscores(fit)
ncol(thetaEstMirt)
nrow(thetaEstMirt)
nrow(pats2)
max(thetaEstMirt[,11] - pats2[,12])
mean(thetaEstMirt[,11] - pats2[,12])

#--------------------------Tercer método-------------------------------------------------
#Método MAP
#A la matriz de patrones le quita la columna frecuencia y añade una columna de ceros a modo de valores iniciales
patsTheta = cbind(est$pats[,-ncol(est$pats)],rep(0,nrow(est$pats)))
patsTheta


logL = function(theta,pat,zita,nitems){
  probPat = log(ggpat(pat,zita,nitems,theta)) - ((theta^2)/2)
  -probPat
}


#Uso de la función nlp para oprimizar funciones de R -> R
#No se puede usar optim debido a que optimiza funciones de R^p -> R
MAP = function(patsTheta,zita,nitems){
  pat = patsTheta[1:nitems]
  theta = patsTheta[nitems+1]
  opt = nlm(f = logL,p=theta,pat=pat,zita=zita,nitems=nitems)
  opt$estimate
}

#Apply para cada uno de los patrones
scoresSics = apply(patsTheta,1,FUN = MAP,zita = est$zita,nitems =10)
scoresMirt = fscores(fit)

#Comparaciones con MIRT
scoresMirt[,12] - scoresSics
mean(scoresMirt[,12] - scoresSics)
max(scoresMirt[,12] - scoresSics)

thetas = pats2[,12]
thetas
