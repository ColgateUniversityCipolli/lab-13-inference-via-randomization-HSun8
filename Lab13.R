# Henry Sun 
# Lab 13 
library(tidyverse)
library(e1071)
#########################################################################
# Task 1
# part a
(zebrafinch.data <- read_csv("zebrafinches.csv"))
further.data <- zebrafinch.data$further
(n <- length(further.data))
mu0 <- 0
# hedge's g + CI
# hedges_g(x = x.further, mu = mu0, alternative = "less")
# interpret_hedges_g(-1.51)

# t.test
further.t.test <- t.test(x=further.data, mu = mu0, alternative = "less")
(t.stat <- further.t.test$statistic[[1]])

(error.num <- skewness(further.data) * (2*t.stat^2 + 1) * dnorm(t.stat))
(error.denom <- 6 * sqrt(n))
(potential.error <- error.num/error.denom)

# part b
for (i in seq(-5,5,length.out=1000)){
  
}
