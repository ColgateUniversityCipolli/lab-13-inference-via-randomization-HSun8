# Henry Sun 
# Lab 13 
library(tidyverse)
library(e1071)
#########################################################################
# Task 1
# part a
(zebrafinch.data <- read_csv("zebrafinches.csv"))
mu0 <- 0
further.data <- zebrafinch.data$further
(n <- length(further.data))

# t.test and t.stat
further.t.test <- t.test(x=further.data, mu = mu0, alternative = "less")
(t.further <- further.t.test$statistic[[1]])

# potential error calculation
(error.num <- skewness(further.data) * (2*t.further^2 + 1) * dnorm(t.further))
(error.denom <- 6 * sqrt(n))
(potential.error <- error.num/error.denom)

# part b
gg.errors <- rep(NA, length.out = 1000)
gg.tvals <- seq(-10,10,length.out=1000)

# create data for errors
for (i in 1:length(gg.tvals)){
  num <- skewness(further.data) * (2*gg.tvals[i]^2 + 1) * dnorm(gg.tvals[i])
  denom <- 6 * sqrt(n)
  gg.errors[i] <- num/denom
}

# plot
errors.plot <- ggplot()+
  geom_line(aes(x= gg.tvals, y = gg.errors))+
  theme_bw()+
  ylab("Potential Error")+
  xlab("T")+
  ggtitle("Potential Error for T, from -10 to 10")+
  geom_vline(aes(xintercept = t.further), color = "red")

# part c
alpha <- 0.05
t.alpha <- qnorm(alpha)

(min.nsize <- ((skewness(further.data) * (2*t.alpha^2 + 1) * dnorm(t.alpha))/
                 (6 * 0.1* alpha))^2)

# Task 2
# part a
R <- 1000
resamples <- tibble(xbars =rep(NA, R))
for (i in 1:R){
  curr.resample <- sample(x = further.data,
                          size= length(further.data),
                          replace = T)
  resamples$xbars[i] <- mean(curr.resample)
}
resamples <- resamples |>
  mutate(t = (xbars-0)/(sd(further.data)/sqrt(n)))
(mean(resamples$t))


(delta <- mean(resamples$t) - 0)
(delta <- abs(mean(resamples$t) - mu0))
(delta <- -mean(resamples$t) + mu0)
(delta <- -(mean(resamples$t) - mu0))

#low <- 
resamples <- resamples |>
  mutate(t.shifted = t - delta)

ggplot(resamples)+
  geom_histogram(aes(x=t, y= after_stat(density)))+
  geom_density(aes(x=t))+
  geom_histogram(aes(x=t.shifted, y = after_stat(density)))+
  geom_density(aes(x=t.shifted))

# part b
further.t.test <- t.test(x=further.data, mu = mu0, alternative = "less")
(t.stat <- further.t.test$statistic[[1]])
low <- t.stat

resamples |>
  summarize(mean = mean(t.shifted),
            p.low = mean(t.shifted <= low))
            
