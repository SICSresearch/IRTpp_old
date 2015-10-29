library(IRTpp)
context("irtpp test")


test_that("irtpp type , dimension and accuracy ",{
  ts = simulateTest(model="2PL",items=5,individuals=100,seed=1)
  est = irtpp(ts$test[[1]],2)
  expect_is(ts$test[[1]],"matrix")
  expect_equal(length(ts$test[[1]]),500)
  expect_equal(est[1,1],1.14853190978)  
})