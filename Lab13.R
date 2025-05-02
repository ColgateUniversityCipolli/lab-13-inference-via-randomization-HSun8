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
                               (sd(further.data)/sqrt(n))
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

# CI
# further
(boot.further.CI <- quantile(further.xbars, c(0.025, 0.975)))
(t.further.CI <- t.test(x=further.data, mu=mu0, alternative = "two.sided")$conf.int)

# closer
(boot.closer.CI <- quantile(closer.xbars, c(0.025, 0.975)))
(t.closer.CI <- t.test(x=closer.data, mu=mu0, alternative = "two.sided")$conf.int)

# diff
(boot.diff.CI <- quantile(diff.xbars, c(0.025, 0.975)))
(t.diff.CI <- t.test(x=diff.data, mu=mu0, alternative = "two.sided")$conf.int)

