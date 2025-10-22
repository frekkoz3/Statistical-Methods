#################################
# R code for LAB3 - 20/10/2025  #
#                               #
# ! This is the pre-LAB version #
# You must replace XXX properly #
#################################



#   _____________________________________________________________________________________
#   Interval estimation: for $\mu$ under the Gaussian case, when $\sigma^2$ is known ####


B <- 1000     # Number of replications
n <- 10       # sample size
mu <- 5       # True population mean 
sigma <- 2    # True population standard deviation

alpha <- 0.05 # Confidence level: 1- alpha

# CI is matrix where we save the confidence intervals for each replication:
# -) first column: lower bound
# -) second column: upper bound
CI <- matrix(0, B, 2)

# "l" is a vector whose elements assume TRUE (1) or FALSE(0) depending on 
# whether the true parameter value lies within the interval
l <- rep(0, B)


set.seed(1234)
q0975 <- qnorm(1 - alpha/2) # quantile of order 1-alpha/2 of N(0,1)
for(i in 1 : B) {
  x <- rnorm(n, mu, sigma)
  CI[i,] <- mean(x) + c(-q0975, q0975) * sigma/sqrt(n)
  l[i] <- (mu > CI[i,1] & mu < CI[i,2])
}

#Empirical coverage probability: 
mean(l)
par(mfrow= c(1,1))
#Plot the first 100 c.i.:
# black: intervals not including mu
# red: intervals including  mu
plot(1, xlim = c(0, 10), ylim = c(0, 11), type = "n", 
     xlab = expression(mu), ylab = "", yaxt  = "n", 
     main = paste("100 IC for the mean (known variance)"), cex.main = 1.2)
abline(v = mu)

d <- 0
for(i in 1 : 100){
  d <- d + 0.1 
  lines(seq(CI[i, 1], CI[i, 2], length = 100), rep(d, 100), col = (l[i] + 1))
}

# number of intervals (out the 100) including the true parameter value 
sum(l[1 : 100]) 

# Interval estimation: for $\mu$ under the Gaussian case, when $\sigma^2$ is unknown ----
# consider the same setting as above: 
# mu = 5, sigma=2, n=10, B=1000, alpha=0.05
CI <- matrix(0, B, 2)
l <- rep(0, B)

set.seed(1234)
q0975 <- qt(1 - alpha/2, n - 1) # quantile of order 1-alpha/2  of t(n-1)

###############################################################################
# TODO:
# -) Generate samples from the N(mu, sigma^2)
# -) Obtain the confidence interval for each replication 
# -) Use a flag to detect if the confidence interval 
#    includes the true parameter value
for(i in 1 : B) {
  x <- rnorm(n, mu, sigma)
  CI[i,] <- #XXX
  l[i] <- #XXX
}

#Empirical coverage probability:

#XXX
###############################################################################


# Plot of the (first 100) CIs
plot(1, xlim = c(0,10), ylim = c(0,11), type = "n", 
     xlab = expression(mu), ylab = "", yaxt  = "n", 
     main = paste("100 IC for the mean (unknown variance)"), cex.main = 1.2)
abline(v = mu)

d <- 0
for (i in 1 : 100){
  d <- d + 0.1 
  lines(seq(CI[i, 1], CI[i, 2], length = 100), rep(d, 100), col = (l[i] + 1))
}

# number of intervals (out the 100) including the true parameter value 
sum(l[1 : 100]) 





#   ____________________________________________________________________________
#   Interval estimation: difference between two means                       ####

Anor <- read.table("http://stat4ds.rwth-aachen.de/data/Anorexia.dat", 
                   header=TRUE)
# Get difference post-pre treatment for the group cb and c 
cogbehav <- Anor$after[Anor$therapy == "cb"] - Anor$before[Anor$therapy == "cb"]
control <- Anor$after[Anor$therapy == "c"] - Anor$before[Anor$therapy == "c"]


# Get the 95% CI via t.test function   
res <- t.test(cogbehav,control,var.equal=TRUE,conf.level=0.95)
res$conf.int

###############################################################################
# TODO:
# Obtain the result above by hand
# Fill XXX
n1 <- XXX
n2 <- XXX

s2 <- XXX

CI <- XXX 
###############################################################################


 
#   ____________________________________________________________________________
#   Interval estimation: difference between two proportions                 ####

success <- c(315, 304)
total <-  c(604, 597)
res <- prop.test(success, total, conf.level = 0.95, correct = FALSE)
res$conf.int

# By hand
p1 <- success[1]/total[1]
p2 <- success[2]/total[2]

p1 - p2 + c(-1, 1) * qnorm(0.975) * sqrt(p1 * (1 - p1)/total[1] + p2 * (1 - p2)/total[2])
