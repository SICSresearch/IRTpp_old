getwd()
setwd("git/IRTpp/")
yy=read.csv("irtppdebug/farr.log", header = F)
yy=read.csv("irtppdebug/args.log",header = F)
yy=yy[,1:(ncol(yy)-1)]
dim(yy)

persp()

sm.ind = NULL;
i=1
for (i in 1:23) {
  sm = matrix(as.numeric(yy[i,]),ncol=10,nrow=10)
  #image(sm , col = topo.colors(100))
  persp(sm, theta = (i))
  sm.ind[[i]] = sm
  print(i)
  print(sum(sm.ind[[i]]))
  Sys.sleep(0.1)
}

sum(sm.ind[[3]])

for (i in 1:((920-80)/40)) {
  i=i*40;
  sm = matrix(as.numeric(yy[i,]),ncol=10,nrow=10)
  sm = matrix(as.numeric(unlist(yy[i:(i+40),])),ncol=4,nrow=41);
  image(sm , col = topo.colors(100))
  sm.ind[[i]] = sm
  print(i)
  Sys.sleep(0.5)
}

rowSums(yy)
plot(colSums(yy))
?image
f.mat = matrix(fff,ncol = 100,byrow = T)
plot(f.mat[,1])
rowSums(f.mat)
