---
title: "Simulating Tests with IRTpp"
author: "Juan Liberato"
date: "2015-07-15"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Simulating Tests with IRTpp}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

Simulating Item Response Theory tests and item parameters
========================================================

Sometimes it is useful to simulate tests and items in Item Response Theory, this vignette specifies how to simulate tests and use or interpret the output, using the IRTpp package.

To simulate a test use the function `simulateTest` as follows

```r
library(IRTpp)
test <- simulateTest(model="2PL",items=10,individuals=1000)
```

This runs a simulation with a 2 parameter logistic model, using 10 items and 1000 individuals. A list is returned with the item parameters, the test (Correct responses of the individuals are marked as a 1, in dichotomic models).
For the test we just simulated, item parameters are : 


```r
test$itempars
```

```
## $a
##  [1] 0.9146088 1.5047161 0.6867993 1.0090911 1.3239394 1.0760605 1.1743445
##  [8] 1.0996683 1.0078389 1.3468080
## 
## $b
##  [1]  1.5042396  0.2194698 -0.3738645  0.1858426  1.7736034  0.4768768
##  [7] -0.3025159  0.7919127 -1.2915661  0.3137360
## 
## $c
##  [1] 0 0 0 0 0 0 0 0 0 0
```

Where `test$itempars$a` indicates the parameters belonging to the discrimination parameter, and so on. Notice that the c parameter is a vector of 0, since we are simulating a 2 parameter logistic model.

The test element of the list returned by `simulateTest` contain a list of tests generated, notice that if you simulate only one test, the response data will be in :


```r
responses <-test$test[[1]]
responses[1:10,]
```

```
##       [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
##  [1,]    0    0    0    0    0    1    1    1    1     0
##  [2,]    0    0    1    0    1    0    0    1    1     1
##  [3,]    0    1    0    0    0    0    0    0    1     1
##  [4,]    0    0    0    1    0    0    1    0    1     1
##  [5,]    1    0    1    1    0    1    1    1    1     1
##  [6,]    0    0    1    0    1    1    1    0    1     0
##  [7,]    0    0    1    0    0    1    0    0    0     0
##  [8,]    0    1    1    1    1    1    1    1    1     1
##  [9,]    0    1    1    0    0    1    1    0    1     0
## [10,]    0    0    1    0    0    1    1    1    1     1
```

`responses[1:10,]` are the first 10 rows of the simulated test.


## Simulating multiple tests

IRTpp has built-in capabilities for simulating multiple tests at once, use the option `reps` to simulate multiple tests at once. All the tests are simulated with the same individual and item parameters.

For example simulating 3 tests of the one parameter logistic model :

```r
t = simulateTest(model="1PL",items=4,individuals=5,reps=3)
length(t$test)
```

```
## [1] 3
```

## Adding boundaries to the simulation

Sometimes you want the item parameters to fall in a range, use the boundaries option to achieve this.
The boundaries option receives a list of the parameter bounds you want to implement.
Let's look at the parameter boundaries in an ordinary 3PL simulated test.

```r
t3 = simulateTest(model="3PL",items=500,individuals=10);
summary(t3$itempars$c)
```

```
##      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
## 0.0002223 0.0914500 0.1716000 0.1739000 0.2615000 0.3472000
```

As you can see, the item parameters will most likely fall inside 0 and 0.35 which are the package defaults, however we can impose a lower boundary for the c parameter as follows:


```r
bd = list(c_lower=0.2)
t3 = simulateTest(model="3PL",items=500,individuals=10,boundaries=bd);
summary(t3$itempars$c)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##  0.2000  0.2370  0.2770  0.2756  0.3142  0.3500
```

Notice how the lower boundary has been imposed in the c parameter.

## Keeping individual scores in a given threshold

When simulating tests, sometimes, perfect responses by an individual can affect procedures made with the tests themselves, or it is required for the simulation that the individuals do not answer less (or more) than a percent of the items, to impose restrictions in the individual classic scores use the threshold parameter.


```r
t3 = simulateTest(model="3PL",items=10,individuals=100,threshold=0.2);
```

This threshold ensures that the individuals do not answer less than 20% of the answers or more than 80% of the answers.


```r
response <- t3$test[[1]]
summary(rowSums(response))
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##    2.00    5.00    6.00    5.65    7.00    9.00
```
