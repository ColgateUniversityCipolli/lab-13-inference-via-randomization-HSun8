# Henry Sun 
# Lab 13 
library(tidyverse)
library(e1071)
library(xtable)
library(boot)
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

################################################################################
# Task 2
# part a
R <- 1000
# data
further.data <- zebrafinch.data$further
closer.data <- zebrafinch.data$closer
diff.data <- zebrafinch.data$diff
# set.seed() to make results reproducible
set.seed(13345)
# resampling
resamples <- tibble(resamps.further =rep(NA, R),
                    resamps.closer =rep(NA, R),
                    resamps.diff =rep(NA, R))
for (i in 1:R){
  # further
  further.resample <- sample(x = further.data,
                          size= length(further.data),
                          replace = T)
  resamples$resamps.further[i] <- (mean(further.resample - mu0))/
                                  (sd(further.data)/sqrt(n))
  
  # closer
  closer.resample <- sample(x = closer.data,
                             size= length(closer.data),
                             replace = T)
  resamples$resamps.closer[i] <- (mean(closer.resample - mu0))/
                                  (sd(closer.data)/sqrt(n))
  
  # diff
  diff.resample <- sample(x = diff.data,
                            size= length(diff.data),
                            replace = T)
  resamples$resamps.diff[i] <- (mean(diff.resample - mu0))/
                               (sd(diff.data)/sqrt(n))
}

# shifting data
resamples <- resamples |>
 mutate(resamps.further.null = resamps.further - mean(resamps.further),
        resamps.closer.null = resamps.closer - mean(resamps.closer),
        resamps.diff.null = resamps.diff - mean(resamps.diff))

ggplot(resamples)+
  geom_density(aes(x=resamps.further))+
  geom_density(aes(x=resamps.further.null))

(delta <- mean(resamples$t) - mu0)
(delta <- abs(mean(resamples$t) - mu0))
(delta <- -mean(resamples$t) + mu0)
(delta <- -(mean(resamples$t) - mu0))

#low <- 
# resamples <- resamples |>
#   mutate(t.shifted = t - delta)
# 
# ggplot(resamples)+
#   geom_histogram(aes(x=t, y= after_stat(density)))+
#   geom_density(aes(x=t))+
#   geom_histogram(aes(x=t.shifted, y = after_stat(density)))+
#   geom_density(aes(x=t.shifted))

# part b
# further.t.test <- t.test(x=further.data, mu = mu0, alternative = "less")
# (t.stat <- further.t.test$statistic[[1]])
# low <- t.stat
# 
# resamples |>
#   summarize(mean = mean(t.shifted),
#             p.low = mean(t.shifted <= low))
#            

# part b 
boot.pvals <- resamples |>
  summarize(p.further = mean(resamps.further.null <= mean(resamps.further)),
            p.closer = mean(resamps.closer.null >= mean(resamps.closer)),
            p.low.diff = mean(resamps.diff.null <= -mean(resamps.diff)),
            p.high.diff = mean(resamps.diff.null >= mean(resamps.diff)),
            pdiff = p.low.diff + p.high.diff)

# p vals from t test
# further
t.pval.further <- t.test(x=further.data, mu = mu0, alternative = "less")$p.value
# closer 
t.pval.closer <- t.test(x=closer.data, mu = mu0, alternative = "greater")$p.value
# diff
t.pval.further <- t.test(x=diff.data, mu = mu0, alternative = "two.sided")$p.value

all.pvals <- tibble(test = c("T-test", "Bootstrapping"),
                    p.further = c(t.pval.further,boot.pvals$p.further),
                    p.closer = c(t.pval.closer, boot.pvals$p.closer),
                    p.diff = c(t.pval.closer, boot.pvals$pdiff))
xtable(all.pvals)

# part c
# 5th percentile for t test (same for all)
t.5th <- qt(p=0.05, df=n-1)

# 5th percentile for bootstrapping 
# further
boot.f5th <- quantile(resamples$resamps.further.null, 0.05)

# closer
boot.c5th <- quantile(resamples$resamps.closer.null, 0.05)

# diff
boot.d5th <- quantile(resamples$resamps.diff.null, 0.05)

# create table
all.5thp <- tibble(Type = c("T-test", "Bootstrapping"),
                   further = c(t.5th, boot.f5th),
                   closer = c(t.5th, boot.c5th),
                   diff = c(t.5th, boot.d5th))
xtable(all.5thp)

# bca helper
# boot.mean <- function(d, i) {
#   mean(d[i])
# }
# closer.helper <- boot(data = closer.xbars,
#                       statistic = boot.mean,
#                       R=R)
# (boot.closer.bca <- boot.ci(closer.helper, type="bca"))

R <- 1000

# part d
# need to keep track of mean (xbar) now
further.xbars <- rep(NA, R)
closer.xbars <- rep(NA, R)
diff.xbars <- rep(NA, R)

for (i in 1:R){
  # further
  further.resample <- sample(x = further.data,
                             size= length(further.data),
                             replace = T)
  further.xbars[i] <- mean(further.resample)
  
  # closer
  closer.resample <- sample(x = closer.data,
                            size= length(closer.data),
                            replace = T)
  closer.xbars[i] <- mean(closer.resample)
  
  # diff
  diff.resample <- sample(x = diff.data,
                          size= length(diff.data),
                          replace = T)
  diff.xbars[i] <- mean(diff.resample)
}

# bca test
# helper function
boot.mean <- function(d, i){
  mean(d[i])
}
R <- 1000
# further
further.helper <- boot(data = further.data,
                       statistic = boot.mean,
                       R = R)

# closer
closer.helper <- boot(data = closer.data,
                      statistic = boot.mean,
                      R = R)

# diff
diff.helper <- boot(data = diff.data,
                    statistic = boot.mean,
                    R = R)

# CI
# further
(boot.further.CI <- quantile(further.xbars, c(0.025, 0.975)))
(t.further.CI <- t.test(x=further.data, mu=mu0, alternative = "two.sided")$conf.int)
(further.bca <- boot.ci(further.helper, type="bca"))

# closer
(boot.closer.CI <- quantile(closer.xbars, c(0.025, 0.975)))
(t.closer.CI <- t.test(x=closer.data, mu=mu0, alternative = "two.sided")$conf.int)

# diff
(boot.diff.CI <- quantile(diff.xbars, c(0.025, 0.975)))
(t.diff.CI <- t.test(x=diff.data, mu=mu0, alternative = "two.sided")$conf.int)
(diff.bca <- boot.ci(diff.helper, type="bca"))


################################################################################
# task 3
# part a
R <- 1000
rand <- tibble(xbars.further = rep(NA, R),
               xbars.closer = rep(NA, R),
               xbars.diff = rep(NA, R))
# since mean 0 (mu0 = 0) under H0 is given, no need to shift

# RANDOMIZE / SHUFFLE
set.seed(13345)
for(i in 1:R){
  # further
  further.rand <- further.data *
    sample(x = c(-1, 1),
           size = length(further.data),
           replace = T)
  
  rand$xbars.further[i] <- mean(further.rand)
  
  # closer
  closer.rand <- closer.data *
    sample(x = c(-1, 1),
           size = length(closer.data),
           replace = T)
  
  rand$xbars.closer[i] <- mean(closer.rand)
  
  # diff
  diff.rand <- diff.data *
    sample(x = c(-1, 1),
           size = length(diff.data),
           replace = T)
  
  rand$xbars.diff[i] <- mean(diff.rand)
}

# no need to shift back either!

# task b
# p-value 
(rand.pvals <- rand |>
  summarize(p.further = mean(xbars.further <= mean(further.data)),
            p.closer = mean(xbars.closer >= mean(closer.data)),
            p.low.diff = mean(xbars.diff <= -mean(diff.data)),
            p.high.diff = mean(xbars.diff >= mean(diff.data)),
            pdiff = p.low.diff + p.high.diff))

# task c
# further 
R <- 10000
mu0.iterate.further <- 0.01
starting.point.further <- mean(further.data)

mu.lower.further <- starting.point.further
repeat{
  rand <- tibble(xbars = rep(NA, R))
  
  # PREPROCESSING: shift the data to be mean 0 under H0
  x.shift <- further.data - mu.lower.further
  # RANDOMIZE / SHUFFLE
  for(i in 1:R){
    curr.rand <- x.shift *
      sample(x = c(-1, 1),
             size = length(x.shift),
             replace = T)
    
    rand$xbars[i] <- mean(curr.rand)
  }
  # Thinking is hard
  rand <- rand |>
    mutate(xbars = xbars + mu.lower.further) # shifting back
  
  # p-value 
  p.val <-mean(rand$xbars >= mean(further.data))
  
  if(p.val < 0.05){
    break
  }else{
    mu.lower.further <- mu.lower.further - mu0.iterate.further
  }
}

mu.upper.further <- starting.point.further
repeat{
  rand <- tibble(xbars = rep(NA, R))
  
  # PREPROCESSING: shift the data to be mean 0 under H0
  x.shift <- further.data - mu.upper.further
  # RANDOMIZE / SHUFFLE
  for(i in 1:R){
    curr.rand <- x.shift *
      sample(x = c(-1, 1),
             size = length(x.shift),
             replace = T)
    
    rand$xbars[i] <- mean(curr.rand)
  }
  # Thinking is hard
  rand <- rand |>
    mutate(xbars = xbars + mu.upper.further) # shifting back
  
  # p-value 
  p.val <-mean(rand$xbars <= mean(further.data))
  
  if(p.val < 0.05){
    break
  }else{
    mu.upper.further <- mu.upper.further + mu0.iterate.further
  }
}

c(mu.lower.further, mu.upper.further)

# closer
R <- 10000
mu0.iterate.closer <- 0.01
starting.point.closer <- mean(closer.data)

mu.lower.closer <- starting.point.closer
repeat{
  rand <- tibble(xbars = rep(NA, R))
  
  # PREPROCESSING: shift the data to be mean 0 under H0
  x.shift <- closer.data - mu.lower.closer
  # RANDOMIZE / SHUFFLE
  for(i in 1:R){
    curr.rand <- x.shift *
      sample(x = c(-1, 1),
             size = length(x.shift),
             replace = T)
    
    rand$xbars[i] <- mean(curr.rand)
  }
  # Thinking is hard
  rand <- rand |>
    mutate(xbars = xbars + mu.lower.closer) # shifting back
  
  # p-value 
  p.val <-mean(rand$xbars >= mean(closer.data))
  
  if(p.val < 0.05){
    break
  }else{
    mu.lower.closer <- mu.lower.closer - mu0.iterate.closer
  }
}

mu.upper.closer <- starting.point.closer
repeat{
  rand <- tibble(xbars = rep(NA, R))
  
  # PREPROCESSING: shift the data to be mean 0 under H0
  x.shift <- closer.data - mu.upper.closer
  # RANDOMIZE / SHUFFLE
  for(i in 1:R){
    curr.rand <- x.shift *
      sample(x = c(-1, 1),
             size = length(x.shift),
             replace = T)
    
    rand$xbars[i] <- mean(curr.rand)
  }
  # Thinking is hard
  rand <- rand |>
    mutate(xbars = xbars + mu.upper.closer) # shifting back
  
  # p-value 
  p.val <-mean(rand$xbars <= mean(closer.data))
  
  if(p.val < 0.05){
    break
  }else{
    mu.upper.closer <- mu.upper.closer + mu0.iterate.closer
  }
}

c(mu.lower.closer, mu.upper.closer)

# diff
R <- 10000
mu0.iterate.diff <- 0.01
starting.point.diff <- mean(diff.data)

mu.lower.diff <- starting.point.diff
repeat{
  rand <- tibble(xbars = rep(NA, R))
  
  # PREPROCESSING: shift the data to be mean 0 under H0
  x.shift <- diff.data - mu.lower.diff
  # RANDOMIZE / SHUFFLE
  for(i in 1:R){
    curr.rand <- x.shift *
      sample(x = c(-1, 1),
             size = length(x.shift),
             replace = T)
    
    rand$xbars[i] <- mean(curr.rand)
  }
  # Thinking is hard
  rand <- rand |>
    mutate(xbars = xbars + mu.lower.diff) # shifting back
  
  # p-value 
  delta <- abs(mean(diff.data) - mu.lower.diff)
  low <- mu.lower.diff - delta # mirror
  high <- mu.lower.diff + delta   # xbar
  p.val <- mean(rand$xbars <= low) +
      mean(rand$xbars >= high)
  
  if(p.val < 0.05){
    break
  }else{
    mu.lower.diff <- mu.lower.diff - mu0.iterate.diff
  }
}

mu.upper.diff <- starting.point.diff
repeat{
  rand <- tibble(xbars = rep(NA, R))
  
  # PREPROCESSING: shift the data to be mean 0 under H0
  x.shift <- diff.data - mu.upper.diff
  # RANDOMIZE / SHUFFLE
  for(i in 1:R){
    curr.rand <- x.shift *
      sample(x = c(-1, 1),
             size = length(x.shift),
             replace = T)
    
    rand$xbars[i] <- mean(curr.rand)
  }
  # Thinking is hard
  rand <- rand |>
    mutate(xbars = xbars + mu.upper.diff) # shifting back
  
  # p-value 
  delta <- abs(mean(diff.data) - mu.upper.diff)
  low <- mu.upper.diff - delta # mirror
  high <- mu.upper.diff + delta   # xbar
  p.val <- mean(rand$xbars <= low) +
    mean(rand$xbars >= high)
  
  if(p.val < 0.05){
    break
  }else{
    mu.upper.diff <- mu.upper.diff + mu0.iterate.diff
  }
}

c(mu.lower.diff, mu.upper.diff)
