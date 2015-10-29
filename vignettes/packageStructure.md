---
title: "IRTpp package structure and naming conventions"
author: "Juan Liberato"
date: "2015-07-15"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{IRTpp package structure and naming conventions}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

Package Structure
==================================
Currently the package has the following main functionalities.

1. Estimating Item and Individual parameters according to an IRT model.
2. Exploratory statistics for tests.
3. Test Simulation and Goodness of fit statistics.

Some common objects between the functions are : 

* Item parameters
* Test (responses)
* Response pattern
* Model

### Item parameters
Item parameters can be passed around in two ways : 

*As a list of named parameters


```r
library(IRTpp) 
t = simulateItemParameters(items=10, model="2PL" , dims= 5)
t
```

```
## $a
##  [1] 0.7621606 1.0330272 1.0086665 0.6783732 1.2551742 0.9292022 1.0038345
##  [8] 0.7303077 1.0364596 1.0634311
## 
## $b
##  [1]  0.7835553  0.5206096 -0.7594416 -0.5060052  0.4095894 -0.4387058
##  [7] -1.7435894 -2.1817079 -0.9354447 -0.8256422
## 
## $c
##  [1] 0 0 0 0 0 0 0 0 0 0
```

*As a matrix or dataframe with equal column number.


```r
parameter.matrix(t,model="2PL")
```

```
##            [,1]       [,2] [,3]
##  [1,] 0.7621606  0.7835553    0
##  [2,] 1.0330272  0.5206096    0
##  [3,] 1.0086665 -0.7594416    0
##  [4,] 0.6783732 -0.5060052    0
##  [5,] 1.2551742  0.4095894    0
##  [6,] 0.9292022 -0.4387058    0
##  [7,] 1.0038345 -1.7435894    0
##  [8,] 0.7303077 -2.1817079    0
##  [9,] 1.0364596 -0.9354447    0
## [10,] 1.0634311 -0.8256422    0
```

For casting between these types `parameter.list` and `parameter.matrix` functions are provided

### Test (responses)

Functions can receive either a single test or a test vector for performing the same computation on many tests, in the function please check if test is a vector of tests instead of a single test and either vectorise the computation or provide an error message to the user.


```r
tests = simulateTest(model="3PL",items=10,individuals=3,reps=4)
typeof(tests$test)
```

```
## [1] "list"
```

```r
typeof(tests$test[[1]])
```

```
## [1] "double"
```

```r
is.matrix(tests$test[[1]])
```

```
## [1] TRUE
```

A test must be a R numeric or integer matrix, and typeof double, if checked with `is.matrix` it must return TRUE.

A test can also be expressed as the response patterns of the test and the frequencies of these patterns

### Response patterns

In IRT response patterns are very important because these represent a set of individuals who answered the test in the same manner and must have the same latent trait.

Patterns are usually represented as a ordered matrix, when it is returned from the `individual.traits` function. 


```r
t = simulateTest(model="2PL",items=3,individuals=10)
individual.traits(dataset=t$test[[1]],model="2PL",itempars=parameter.matrix(t$itempars))
```

```
## Error in eval(expr, envir, enclos): could not find function "individual.traits"
```

### Models

Models in IRT can be very diverse, a list with the current valid models must be mantained to instruct users in what models can they use.

A model is represented by a single string with the model standard name in the package, in UIRT the three main models are "1PL" , "2PL" , and "3PL".

In the IRTpp package, multiple functions are given to help receiving models upon user input.

#### Helper functions for model

`check.model` verifies if a model is valid, if not it throws an error to the user.

`irtpp.model` Casts possible unambiguous inputs for a model into a single standardized model string.

If a function receives the model as a parameter this is usually named modispla
del. And the function irtpp.model can cast the user input.
irtpp.models can list all the available models.


```r
irtpp.model(2)
```

```
## [1] "2PL"
```

```r
irtpp.model("2")
```

```
## [1] "2PL"
```

```r
irtpp.model("2pl")
```

```
## [1] "2PL"
```

```r
irtpp.model("2PL")
```

```
## [1] "2PL"
```

```r
irtpp.models()
```

```
## [1] "1PL"   "2PL"   "3PL"   "1PLAD"
```

#### Probability functions within models

The naming of probability functions is usually probability.<model> in lowercase, the function irtpp.p is also provided which gives the model probability function for a specific model.


```r
probability.3pl
```

```
## function (z, a = z$a, b = z$b, c = z$c, theta, d = -a * b, cp = NULL) 
## {
##     if (is.null(cp)) {
##         c + ((1 - c)/(1 + exp(-a * (theta - b))))
##     }
##     else {
##         exp(cp)/(1 + exp(cp)) + (1 - (exp(cp)/(1 + exp(cp)))) * 
##             (1 + exp(-(a * theta + d)))^(-1)
##     }
## }
## <environment: namespace:IRTpp>
```

```r
irtpp.p(model="3PL")
```

```
## NULL
```
