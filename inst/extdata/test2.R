library(IRTpp)
t = simulateTest(model="2PL", individuals=1000, items=10, seed=1)
uirtestimate(t$test, 2)
