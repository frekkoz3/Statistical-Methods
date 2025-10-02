## Exercises on Prob1 ----

### ========== Prob1 ==========
## Random variables

## Build.in functions in R to work with rvs ------------------------------------
#?dnorm #pdf/pmf
#?pnorm #cdf
#?qnorm #quantile
#?rnorm #random numbers generation

#other examples: dgamma, dbinom, dpois, dexp, dchisq, dt,...

## Compute P(X = 0) for X ∼ Pois(3) --------------------------------------------

dpois(0, lambda = 3)

## Compute P(X = x) for X ∼ Pois(3), for x = 0,1,...,5 -------------------------

x <- 0:5
probs <- dpois(x, lambda = 3)
data.frame(x, probs)

## Compute P(X<=4) for X ~ N(3,9) and make a plot of the normal distribution 
## highlighing the probability you have computed ----------

p <- pnorm(4, mean = 3, sd = 3)
p

x <- seq(-5, 12, length=200)
y <- dnorm(x, mean=3, sd=3)

plot(x, y, type="l", lwd=2, main="N(3,9)")
polygon(c(x[x<=4],4), c(y[x<=4],0), col="skyblue")


## Compute P(1<=X<=10) for X ~ Gamma(shpae=3,scale=2) and make a plot of the 
## gamma distribution highlighing the probability you have computed ----------

p <- pgamma(10, shape=3, scale=2) - pgamma(1, shape=3, scale=2)
p

x <- seq(0, 15, length=200)
y <- dgamma(x, shape=3, scale=2)

plot(x, y, type="l", lwd=2, main="Gamma(3,2)")
polygon(c(x[x>=1 & x<=10], 10, 1),
        c(y[x>=1 & x<=10], 0, 0),
        col="lightgreen")

## For X ~ N(3,10), what is the quantile of order 0.9?

qnorm(0.9, mean=3, sd=sqrt(10))

##repeat the points above with another continuous rv of your choice

pexp(4, rate=0.5)
qexp(0.9, rate=0.5)

x<- seq(0, 15, length=200)
y <- dexp(x, rate=0.5)

plot(x, y, type="l", lwd=2, main="Exponential(rate=0.5)")
polygon(c(x[x<=4],4), c(y[x<=4],0), col="orange")

## qqplot normal vs student-t  -------------------------------------------------
set.seed(123)
par(mfrow=c(1,2))
x <- rt(1000,2);
qqnorm(x, pch = 16, main = "")
qqline(x)

y <- rnorm(1000)
qqnorm(y, pch = 16, main = "")
qqline(y)

## make the ecdf plot from scratch following the steps and using the simulated data ----
## 1. sort the data
## 2. compute the empirical cumulative distribution function
## 3. plot the sorted data vs the ecdf
## 4. overlap the theoretical cumulative distribution function

set.seed(123)
par(mfrow=c(1,2))
n <- 11 #sample size
x <- rnorm(n) #simulate from a std normal

## Step 1: sort the data
x_sorted <- sort(x)

## Step 2: compute ECDF values = (1:n)/n
ecdf_vals <- (1:n)/n

## Step 3: plot the ECDF (points + step function style)
plot(x_sorted, ecdf_vals, pch=16, col="blue",
     main="ECDF vs Theoretical CDF",
     xlab="x", ylab="Cumulative Probability",
     ylim=c(0,1))
# Add ECDF as steps
for (i in 1:n) {
  if (i==1) {
    segments(-Inf,0, x_sorted[i], 0, col="blue")
  } else {
    segments(x_sorted[i-1], ecdf_vals[i-1], x_sorted[i], ecdf_vals[i-1], col="blue")
  }
  segments(x_sorted[i], ecdf_vals[i]-1/n, x_sorted[i], ecdf_vals[i], col="blue")
}
segments(x_sorted[n], 1, Inf, 1, col="blue")

## Step 4: overlay the theoretical CDF
curve(pnorm(x), col="red", lwd=2, add=TRUE)
legend("topleft", legend=c("ECDF","Theoretical CDF"),
       col=c("blue","red"), lwd=c(1,2), bty="n")
