
############################################3 Only one


library(IRTpp);library(FactoMineR)
dims=2;items =40
t = simulateTestMD(items=items,individuals = 1000, dims = dims, clusters = 2)
test = t$test
initialvalues = matrix(3,nrow = items, ncol = dims+2);initialvalues[,1:dims] = 0.9;
initialvalues[,dims+1] = 0;initialvalues[,dims+2] = 0.2
testt = t;
upper = testt$clustinit;lower = testt$clustinit + testt$clusters -1 
##Estos son los clusters para fix.items

fixedclusts  = c(t(matrix(c(upper,lower),nrow = length(testt$clusters))))
fi = fix.items(test,fixedclusts)
##Items reportados por fix.items.

initiala = initialvalues[,1:dims]
for (ii in 1:dims) {
  a.tmp = initiala[upper[[ii]]:lower[[ii]],-ii];
  initiala[upper[[ii]]:lower[[ii]],-ii] = rep(0.3,length(a.tmp))
}
initialvalues[,1:dims] = initiala

for (ii in 1:dims) {
  initialvalues[t$clustinit[ii],1:dims] = 0;
  initialvalues[t$clustinit[ii],ii] = 1
}

#est = irtpp(dataset = test,model = "3PL", dims = 2 , initialvalues = initialvalues, restricted.items = fi )
est = irtpp(dataset = test,model = "3PL", dims = 2 , initialvalues = initialvalues, restricted.items = t$clustinit[1:dims] )
p3 = parameter.matrix(pars = est$zita,model = "3PL",dims = 2)
p3
