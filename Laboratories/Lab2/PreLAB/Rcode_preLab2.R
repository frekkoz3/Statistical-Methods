#################################
# R code for LAB2 - 13/10/2025  #
#                               #
# ! This is the pre-LAB version #
# You must replace XXX properly #
#################################


#   ____________________________________________________________________________
#   Distribution of the sample mean                                         ####

# Case 1) Generate from N(0,1)
# Case 2) Generate from t-student with 3 degrees of freedom 
# Case 3) Generate from (continuous) uniform in (0,1)


set.seed(1234) # Set a seed 
R <- 1000      # Fix the number of simulations
n <- 30        # Set the sample size 

# Generate the sample of size n (for the 3 distributions) R times 
samples <- array(0, c(3, n, R))
for(i in 1 : R){
  samples[1, ,i] <- rnorm(n, 0, 1)
  samples[2, ,i] <- rt(n, df = 3)
  samples[3, ,i] <- runif(n, 0, 1)
}


# Getting the sample mean
sample_stat <- apply(samples, c(1,3), mean)

# Equivalently (creating a)
sample_stat2 <- matrix(0, nrow = 3, ncol = R)
for(i in 1 : R) {
  sample_stat2[, i] <- apply(samples[, , i], 1, mean)
}
max(abs(sample_stat - sample_stat2)) # Check 


# Visualise the sample distributions by means of hiostograms 
par (mfrow=c(1,3))

hist(sample_stat[1,], nclass = 30, probability = TRUE, 
     xlab="y", main= "N(0,1)", cex.main = 1.5)
hist(sample_stat[2,], nclass = 30, probability = TRUE, 
     xlab = "y", main = expression(t[3]), cex.main = 1.5)
hist(sample_stat[3,], nclass = 30, probability = TRUE, 
     xlab = "y", main = "U(0,1)", cex.main = 1.5)


###############################################################################
# TODO: Overlap the density of the proper normal distribution to the histogram:

# Case 1) N( , )

# Case 2) N( , )

# Case 3) N( , )

###############################################################################


# Other graphichal tools for this assesment 

###############################################################################
# TODO: Overlap the density of the proper normal distribution substituing the XXX:

# ECDF vs CDF 
par(mfrow = c(1, 3))
plot(ecdf(sample_stat[1,]), xlab="y", main= "N(0,1)", cex.main = 1.5)
curve(XXX, add = TRUE, col = "red", lty = 2, lwd = 2)

plot(ecdf(sample_stat[2,]), xlab="y", main= expression(t[3]), cex.main = 1.5)
curve(XXX, add = TRUE, col = "red", lty = 2, lwd = 2)

plot(ecdf(sample_stat[3,]), xlab="y", main= "U(0,1)", cex.main = 1.5)
curve(XXX, add = TRUE, col = "red", lty = 2, lwd = 2)


# QQ-Norm 
par(mfrow = c(1, 3))

help(abline)
qqnorm(sample_stat[1,])
abline(XXX, XXX)

qqnorm(sample_stat[2,])
abline(XXX, XXX)

qqnorm(sample_stat[3,])
abline(XXX, XXX)
###############################################################################

par (mfrow = c(1, 2))
sigma <- 1 

# Recall the the samples are arleady generated and stored in the array samples
sample_var <- apply(samples, c(1,3), var)[1,]

# Histogram 
hist(sample_var, nclass = 30, probability = TRUE, 
     xlab = expression(s^2), main = "Case 1: N(0,1)", cex.main = 1.5)
curve((n-1)/sigma^2 * dchisq((n-1)*x/sigma^2, df = n - 1),  add = TRUE, col = "red", lwd = 2)

# ECDF vs CDF 0
plot(ecdf(sample_var), xlab = expression(s^2), main = "Case 1: N(0,1)", cex.main = 1.5)
curve( pchisq((n-1)*x/sigma^2, df = n - 1),  add = TRUE, col = "red", lwd = 2)




#   ____________________________________________________________________________
#   Point estimation: Comparison of four estimators of the mean parameter   ####

R <- 1000   # Number of replications
n <- 10     # Sample size
mu <- 2     # Population mean
sigma <- 2  # Population standard deviation

est <- matrix(0, R, 4)
label_est <- c("Mean", "Median", "(Min + Max)/2", "10pTrimMean")

set.seed(1234)

###############################################################################
# TODO: Generate values and fill the matrix est and compute the bias, variance and MSE
for (i in 1 : R) {
  # XXX generate values and fill the matrix est XXX
}

par(mfrow = c(1, 1), xaxt = "n")
boxplot(est, main="Comparison between four estimators")
par(xaxt = "s")
axis(1, 1 : 4, label_est)
abline(h = mu, lwd = 2, col = "blue")

# Compute the Bias
bias <- #XXX 
bias

# Compute the Variances
variance <- #XXX
variance

# Compute the MSE 
MSE <- #XXX
MSE 
###############################################################################


# Point estimation: Assess properties increasing n ----
R <- 1000
n <- c(10, 200, 1000)
mu <- 2
sigma <- 2

est <- array(0, dim = c(R, 4, 3))
label_est <- c("Mean", "Median", "(Min + Max)/2", "10pTrimMean")

set.seed(13)


for(j in 1 : length(n)){
  for (i in 1 : R) {
    x <- rnorm(n[j], mu, sigma)
    est[i, 1, j] <- mean(x)
    est[i, 2, j] <- median(x)
    est[i, 3, j] <- (max(x) + min(x))/2
    est[i, 4, j] <- mean(x, trim = 0.1)
  }
}  

# Note that the plots obtained here are related to the last case 
# -) You can use the left arrow to see the other plots 
# -) Produce the plots for each k (do not consider the for loop and fix k)
for(k in 1: 4){
  par(mfrow = c(1, 3))
  for (j in 1 : length(n)) {
    hist(est[,k, j], nclass = 30,  probability = T, xlab = "", xlim = c(-2, 6),
         main = label_est[k], cex.main = 2)
    abline(v = mu, col = "blue", lwd = 2)
  }
}  

# Variances
variance <- matrix(0, 3, 4)
for(j in 1:length(n)){
  variance[j,] <- apply(est[,,j], 2, var)
}
# Bias
bias <- matrix(0, 3, 4)
for(j in 1:length(n)){
  bias[j,] <- apply(est[,,j], 2, mean) - mu 
}

# MSE 
MSE <- variance + bias^2

colnames(bias) <- colnames(variance) <- colnames(MSE) <- label_est
rownames(bias) <- rownames(variance) <- 
  rownames(MSE) <- c("n = 10", "n = 200", "n = 1000" )

variance 
bias
MSE


#   ____________________________________________________________________________
#   Point estimation: Comparison of unbiased and biased sample variance estimators ####

# Check the biased nature of the (biased) sample variance estimator via MC simulation, 
# by generating n = 10 iid values from N(0,1), and compare the results with the unbiased estimator 

#Initial settings
set.seed(2)
R <- 1000
n <- 10
mu <- 0
sigma <- 1

# Save the results in two vectors:
# s2: unbiased sample variance estimates 
# s2_b: biased sample variance estimates
s2 <- rep(0, R)
s2_b <- rep(0, R)

# For each replication we generate 10 samples from
# a normal r.v. with mean mu and variance sigma^2
# and we compute the four sample variance estimates 
for(i in 1 : R) {
  y <- rnorm(n, mu, sigma)
  s2[i] <- var(y) 
  s2_b[i] <- var(y) * (n - 1)/n
}

# Get a point estimates 
s2_mean <- mean(s2)
s2_b_mean <- mean(s2_b)

#plot s2
par(mfrow = c(1, 2), oma = c(0, 0, 0, 0))
hist(s2, breaks = 50, xlab = expression(s^2), probability = TRUE,
     main = expression(s^2), cex.main = 1.5)
#in red the true mean, in blue the estimated mean
abline(v = sigma^2, col = "red", lwd = 2)
abline(v = s2_mean, col = "blue", lwd = 3, lty = 3) 

#plot s2 biased
hist(s2_b, breaks = 50, xlab = expression(s[b]^2), probability = TRUE,
     main = expression(s[b]^2), cex.main = 1.5)
#in red the true mean, in blue the estimated mean
abline(v = sigma^2, col = "red", lwd = 3)
abline(v = s2_b_mean, col = "blue", lwd = 3, lty = 2) 

