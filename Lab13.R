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
t.errors <- rep(NA, length.out = 1000)
# NOTE: I don't like the way I'm doing this for loop ...
j <- 1
for (i in seq(-10,10,length.out=1000)){
  num <- skewness(further.data) * (2*i^2 + 1) * dnorm(i)
  denom <- 6 * sqrt(n)
  t.errors[j] <- num/denom
  j <- j + 1
}
ggplot()+
  geom_line(aes(x=seq(-10,10,length.out=1000), y = t.errors))

# part c
alpha0.1 <- error.formula
(min.nsize <- (error.num/(6 * alpha0.1))^2)


