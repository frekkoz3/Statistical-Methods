#################################
# R code for LAB4 - 27/10/2025  #
#                               #
# ! This is the pre-LAB version #
# You must replace XXX properly #
#################################


Anor <- read.table("http://stat4ds.rwth-aachen.de/data/Anorexia.dat", 
                   header=TRUE)

# Get difference post-pre treatment for the group cb and c 
cogbehav <- Anor$after[Anor$therapy == "cb"] - Anor$before[Anor$therapy == "cb"]
control <- Anor$after[Anor$therapy == "c"] - Anor$before[Anor$therapy == "c"]


# Get the 95% CI via t.test function   
res <- t.test(cogbehav,control,var.equal=TRUE,conf.level=0.95)
res$conf.int


# Obtain the result above by hand: we will see together during next lab
n1 <- length(cogbehav)
n2 <- length(control)

var(cogbehav)
var(control)

s2 <- ((n1- 1) * var(cogbehav)  + (n2 - 1) * var(control))/(n1 + n2 -2)



#   ____________________________________________________________________________
#   Hypothesis testing                                                      ####


##  ............................................................................
##  Test for the mean difference                                            ####


# Two-sided two-sample test

res.two  <- t.test(cogbehav, control, var.equal = TRUE)
res.two

# PErform th test by hand (complete XXX)
testStat <- (mean(cogbehav) - mean(control))/sqrt(s2 * (1/n1 + 1/n2))

# p-value
p.value.two <- 2 * pt(testStat, df = n1 + n2 -2, lower.tail =  FALSE)

2 * (1 - pt(testStat, df = n1 + n2 -2))

# One-sided two-sample test
res.one  <- t.test(cogbehav, control, var.equal = TRUE, alternative = "greater")
res.one

## one-sided
p.value.one <- pt(testStat, df = n1 + n2 -2, lower.tail =  FALSE)


library(RColorBrewer)
plotclr <- brewer.pal(6, "YlOrRd")

curve(dt(x, n1 + n2 - 2), xlim = c(-5, 5), ylim = c(0, 0.4), 
      main = "p-values and rejection region", col = "blue", 
      lwd = 2, xlab = "x-y",  ylab = expression(t[13]),  yaxs="i")
cord.x <- c(qt(0.95, n1 + n2 - 2),seq(qt(0.95, n1 + n2 - 2), 5, 0.01), 5)
cord.y <- c(0, dt(seq(qt(0.95, n1 + n2 - 2), 5, 0.01), 13), 0)
polygon(cord.x, cord.y, col = plotclr[3], border = NA )

abline(v = res.one$statistic, lty = 2, lwd = 2, col = "red")
text(0, 0.2, paste("Accept", expression(H0)))
text(2.7, 0.08, paste("Reject", expression(H0)))
text(as.double(res.one$statistic) - 0.15, 0.02, "t", col = "red", cex = 1.2)
# On the plot above note that the value of observed test statistic can vary as
# function of the sample (only for this particular case is too close to alpha
# that you cannot detect a graphical difference)



##  ............................................................................
##  Test for equality of variance                                           ####

# Using the var.test function
var.test(cogbehav, control, alternative = "two.sided") 

# Test for equality of variance by hand
ratiovar <- var(cogbehav)/var(control) # Test statistic 
pv_bi <- 2 * min(pf(ratiovar, n1 - 1, n2 - 1, lower = FALSE), 
                 1 - pf(ratiovar, n1 - 1, n2 - 1, lower = FALSE))
pv_bi


##  ............................................................................
##  Pearson's chi-squared test                                              ####

## Test for independence
n <- 760
obs_freq <- matrix(c(50, 70, 30, 100,
                     114, 30, 10, 100, 
                     116, 27, 13, 100),3,4,T)
colnames(obs_freq) <- c("A", "B", "AB", "O")
rownames(obs_freq) <- c("South", "Central", "North")
chisq.test(obs_freq)


mx <- colSums(obs_freq)/n; mx
my <- rowSums(obs_freq)/n; my

# !!!!: Note that the result of outer()
# can be obtained by doing 
# t(mx %*% t(my))
# my%*%t(mx)
# tcrossprod(my,mx)

exp_freq <- outer(my, mx) * n; exp_freq
chi2 <- sum((obs_freq - exp_freq)^2/exp_freq); chi2
chi2
pchisq(chi2, (ncol(obs_freq) - 1) * (nrow(obs_freq) - 1), lower.tail = FALSE)

## Goodness of fit test

child <- 0:5
fam <- c(52, 60, 55, 18, 8, 7)
lambda <-  1.5
plot(0:5, fam/200, ylim = c(0,0.35))
segments(0:5,rep(0,5), 0:5, 
         c(dpois(0:4, lambda),ppois(4, lambda, lower = FALSE)))

c(dpois(0:4, lambda), ppois(4, lambda, lower = FALSE))


exp <- c(dpois(0:4, lambda),ppois(4, lambda, lower = FALSE)) * 200 
round(exp, 4)
chisq_el <- (fam-exp)^2/exp
round((fam-exp)^2/exp, 4)


chisq.obs <- sum(chisq_el)
chisq.obs 
pchisq(chisq.obs, df = 5, lower=FALSE)


chisq.test(fam, p = exp/200)



#   ____________________________________________________________________________
#   Likelihood Inference                                                    ####


##  ............................................................................
##  Binomial model                                                          ####

# Let's generate some data (n = 100, p = 0.6)
set.seed(13)
n <- 100 
p <- 0.6
y <- rbinom(n, 1, p)

# (log) likelihood function: equivalent formulations
# Using dbinom
llik_bin <- function(theta, data){
  sum(dbinom(x = data, size = 1, prob = theta, log = TRUE))
}

# Writing the expression by hand 
llik_bin2 <- function(theta, data){
  sum(data) * log(theta) + (length(data) - sum(data)) * log(1-theta)
}

# Check for a fixed value of p 
llik_bin(0.5,y)
llik_bin2(0.5,y)

# ML estimate
MLEp  <- mean(y)
MLEp

# Value of the log likelihood at the MLE
llik_bin(MLEp, y)

# Check
pgrid <- seq(0.01, 0.99, by = 0.01)
pgrid[which(llik_bin2(pgrid, y) == max(llik_bin2(pgrid, y)))]

#  Note that the last check make use of the function we built manually
# This because this function is vectorised, in contrast to  first one.
# Indeed 
llik_bin(c(0.5, 0.75, 0.99), y) # Wrong 

# The results above is due to this operation
sum(dbinom(x = y[seq(1,100,3)], size = 1, prob = c(0.5), log = TRUE)) + 
  sum(dbinom(x = y[seq(2,100,3)], size = 1, prob = c(0.75), log = TRUE)) + 
  sum(dbinom(x = y[seq(3,100,3)], size = 1, prob = c(0.99), log = TRUE))

llik_bin2(c(0.5, 0.75, 0.99), y) # Correct

# However, we can vectorise the function  by means of 
llik_bin_v <- Vectorize(llik_bin, 'theta')

# So we can visualise the log likelihood function (using both the approaches)
par(mfrow = c(1,2))
curve(llik_bin2(theta=x, data = y), 0, 1, 
      xlab= expression(p), ylab = expression(l(p)),
      main = "Log-likelihood function")
abline(v = p, col = "red")
abline(v = MLEp, col = "blue")
legend("topleft", legend = c(expression(p), expression(hat(p))),
       col = c("red", "blue"), lty = c(1, 1))


curve(llik_bin_v(theta=x, data = y), 0, 1, 
      xlab= expression(p), ylab = expression(l(p)),
      main = "Log-likelihood function")
abline(v = p, col = "red")
abline(v = MLEp, col = "blue")
legend("topleft", legend = c(expression(p), expression(hat(p))),
       col = c("red", "blue"), lty = c(1,1 ))

