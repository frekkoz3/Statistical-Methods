#################################
# R code for LAB5 - 06/10/2025  #
#                               #
#################################


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


#   ____________________________________________________________________________
#   Reparametrisation                                                       ####


# logit(theta) 
psi <- function(theta){
  log(theta/(1-theta))
}

# theta  
theta <- function(psi){
  exp(psi)/(1 + exp(psi))
}

# log likelihood under reparametrisation
llik_bin_rep <- function(param, data) llik_bin2(theta(param), data)

par(mfrow = c(1,1))
curve(llik_bin_rep(param = x, data = y), -10, 10, 
      xlab= expression(psi), ylab = expression(l(psi)),
      main = "Log-likelihood function")
abline(v = psi(p), col = "red")
abline(v = psi(MLEp), col = "blue")
legend("topleft", legend = c(expression(psi), expression(hat(psi))),
       col = c("red", "blue"), lty = c(1, 1))


# Monte Carlo simulation: asymptotic normality of MLE
set.seed(1234)
R <- 5000
p <- 0.6
n2 <- 1000 

prop <- rep(0, R)

for(i in 1:R){
  sample <- rbinom(n2, 1, p)
  prop[i] <- mean(sample)
}


par(mfrow=c(1,2))
hist(prop, freq = F, breaks = 30, main = expression(p))
curve(dnorm(x, mean = p, sd = sqrt((p * (1 - p)/n2))), 
      lwd = 2, xlab = "", ylab = "", add = T)

hist(psi(prop), freq = F, breaks = 30, main = expression(psi))
curve(dnorm(x, mean = psi(p), sd = 1/sqrt((p * (1 - p) * n2) )), 
      lwd = 2, xlab = "", ylab = "", add = T)

## Optimisation ----
# Since many optimisers allow only minimisation, 
# we implement the negative lok-likelihood (also considering the reparametrisation) 
nllik_bin <- function(theta, data){
  -sum(data) * log(theta) - (length(data) - sum(data)) * log(1-theta)
}
nllik_bin_rep <- function(param, data) nllik_bin(theta(param), data)

# Score function
score <- function(theta, data){
  sum(data)/theta - (length(data) - sum(data))/(1 - theta)
}

# Observed information
obs_info <- function(theta, data){
  sum(data)/(theta^2) + (length(data) - sum(data))/((1 - theta)^2)
}

# Let's do a check of the score and the obs info by using numerical derivatives
library(numDeriv)
grad(nllik_bin, 0.5, data = y)
-score(0.5, data = y)

hessian(nllik_bin, 0.5, data = y)
obs_info(0.5, data = y)


# Minimise the negative log-likelihood 
bin.nlm_start1 <- nlm(f = nllik_bin, p = 0.5,
                      data = y, hessian = TRUE)
bin.nlm_start1

nllik_bin(bin.nlm_start1$estimate, data = y) 
obs_info(bin.nlm_start1$estimate, data = y)

# Minimise the negative log-likelihood using optim 
bin.optim_start1 <- optim(par = 0.5, fn= nllik_bin, hessian = TRUE,
                          data = y, method = "L-BFGS-B", lower = 1e-7, upper = 1-1e-7)
bin.optim_start1

bin.optim_start2 <- optim(par = 0.001, fn = nllik_bin,
                          data = y, method = "L-BFGS-B", lower = 1e-7, upper =1-1e-7)
bin.optim_start2
optimHess(bin.optim_start2$par, nllik_bin, data = y)

# Optimise under reparametrisation
bin.nlm_start1_rep <- nlm(f = nllik_bin_rep , p = 0,
                          data = y, hessian = TRUE)
bin.nlm_start1_rep

nllik_bin_rep (bin.nlm_start1_rep$estimate, data = y) 

c(bin.nlm_start1_rep$estimate, psi(bin.nlm_start1$estimate), psi(mean(y)))

hessian(nllik_bin_rep, bin.nlm_start1_rep$estimate, data = y)



#   ____________________________________________________________________________
#   Linear regression                                                       ####

# Load the data
library(DAAG)
data(nihills)
n <- nrow(nihills)

# Structure, first and last six rows of the data 
str(nihills)
head(nihills)
tail(nihills)

# Summary of the variables  
summary(nihills)

# Focus on time of males 

### Simple linear regression ----
# Consider time ~ dist +error 

# Scatterplot 
par(mfrow=c(1,1))
plot(nihills$dist, nihills$time, pch = 19, xlab = "dist", ylab = "time")
# Equivalently
# with(nihills, plot(dist, time, pch = 19))

# Fit the model 
lm_dist <- lm(time ~ dist, data = nihills)
lm_dist

# A detailed model output including several information
summary(lm_dist)

# Residual plots
par(mfrow = c(2, 2))
plot(lm_dist)

# Scatterplot + (estimated) regression line
par(mfrow = c(1, 1))
with(nihills, plot(dist, time, pch = 19))
abline(coef(lm_dist), col = "red", lty = "solid")
text(13, 3, expression(time == hat(beta)[0] + hat(beta)[1] * dist), col = "red")

points(nihills$dist, predict(lm_dist), col = "red", pch = 19, cex = 0.8)

nihills[17,]

segments(nihills[17,]$dist, nihills[17,]$time, 
         nihills[17,]$dist, fitted(lm_dist)[17], lty = "dashed")

### The summary output of a lm object ----

### Save the results of the fitted model summary in a object 
sum_lm_dist <- summary(lm_dist)
### Get the needed quantities 
X <- model.matrix(~ dist, data = nihills)
y <- nihills$time
p <- ncol(X)

### Estimated coefficients 
beta.hat <- solve(t(X) %*% X) %*% t(X) %*% y
beta.hat
sum_lm_dist$coefficients[,1] # Check 

### Fitted values 
y.hat <- X %*% beta.hat
max(abs(predict(lm_dist) - y.hat)) # Check

### Residuals 
residuals <- y - y.hat 
max(abs(residuals(lm_dist) - residuals)) # Check
cbind(summary(residuals), round(summary(sum_lm_dist$residuals), 5)) # Check

### Estimate of variance of error distribution using residuals 
s2 <- sum(residuals^2)/(n - p); s2
sum_lm_dist$sigma^2 # Check 

### Estimates of stdev of regression coefficient estimator
stdb <- sqrt(diag(s2 * solve(t(X) %*% X))); stdb
sum_lm_dist$coefficients[, 2] ## Check 

### Observed test statistics 
beta.hat/stdb
as.numeric(sum_lm_dist$coefficients[, 3])  ## Check 

### p-values 
2 * pt(abs(beta.hat/stdb), df = n - p, lower.tail = FALSE)
as.numeric(sum_lm_dist$coefficients[,4]) ## Check

### r-squared 
Tot_SS <- sum((y - mean(y))^2) 
Res_SS <- sum((y.hat - y)^2)
Mod_SS <- sum((y.hat - mean(y))^2)

R_sq <- Mod_SS/Tot_SS; 
R_sq <- 1 - Res_SS/Tot_SS
R_sq
with(nihills, cor(time, dist))^2
sum_lm_dist$r.squared  # Check 

### adjusted r-squared
1 - (Res_SS/(n - p))/(Tot_SS/(n - 1))
R_sq - (1 - R_sq)*(p - 1)/(n - p)     # Check 
sum_lm_dist$adj.r.squared # Check

### Null model 
X0 <- X[, 1]
p0 <- 1
beta.hat0 <- solve(t(X0) %*% X0) %*% t(X0) %*% y
beta.hat0
beta.hat0 <- mean(y)

y.hat0 <- mean(y)
residuals0 <- y - y.hat0 

### F-statistic 
F <- ((sum(residuals0^2) - sum(residuals^2))/(p - p0))/(sum(residuals^2)/(n - p))
F
sum_lm_dist$coefficients[2, 3]^2   # Check

sum_lm_dist$fstatistic # Check

### p-value
pf(F, p0, n - p, lower.tail = FALSE)

