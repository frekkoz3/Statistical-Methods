#################################
# R code for LAB6 - 10/11/2025  #
#                               #
#################################



#   ____________________________________________________________________________
#   Multiple linear regression                                              ####

# Load the data
library(DAAG)
data(nihills)
n <- nrow(nihills)

lm_dist <- lm(time ~ dist, data = nihills)

### Polinomial regression ---
# Fit the model: time ~ dist + dist^2 + error 
lm_dist2 <- lm(time ~ dist + I(dist^2), data = nihills)
# Explore the summary, residuals, and visualise the fitted regression curve
summary(lm_dist2) 
par(mfrow = c(2,2))
plot(lm_dist2)

par(mfrow = c(1,1))
with(nihills, plot(dist, time, pch=19))
curve(predict(lm_dist2, data.frame(dist = x)), col = "red", 
      lty = "solid", lwd = 2, add = TRUE)
text(10, 3, expression(time == hat(beta)[0] + hat(beta)[1] * dist + 
                         hat(beta)[1] * dist^2), col = "red")
points(nihills$dist, predict(lm_dist2), col = "red", pch = 19, cex = 0.8)

### Adding other predictors ---

# Scatterplot matrix 
library(PerformanceAnalytics)
chart.Correlation(nihills[, c("dist", "climb", "time")])

# Fit the model: time ~ dist + climb  + error 
lm_distclimb <- lm(time ~ dist + climb, data = nihills)
summary(lm_distclimb)

par(mfrow = c(2, 2))
plot(lm_distclimb)

# Plot the residuals w.r.t. the explanatory variables 
par(mfrow = c(1,2))
plot(nihills$dist, lm_distclimb$residuals)
abline(h = 0, lty = "dashed")
plot(nihills$climb, lm_distclimb$residuals)
abline(h = 0, lty = "dashed")

# Fit the model: time ~ dist + dist^2 +  climb  + error 
lm_distclimb2 <- lm(time ~ dist + I(dist^2) + climb, data = nihills)
summary(lm_distclimb2) 
par(mfrow = c(2,2))
plot(lm_distclimb2)

### Transforming variables ----

# Get back to the first graphical plots
par(mfrow = c(1, 3))
hist(nihills$time, probability = TRUE, breaks = 15)
hist(nihills$climb, probability = TRUE, breaks = 15)
hist(nihills$dist, probability = TRUE, breaks = 15)

# Create and add to the dataset the log variables
nihills$log_time <- log(nihills$time)
nihills$log_climb <- log(nihills$climb)
nihills$log_dist <- log(nihills$dist)
chart.Correlation(nihills[, c("log_climb", "log_dist", "log_time")])


lm_ldist_lclimb <- lm(log_time ~ log_dist + log_climb, data=nihills)
summary(lm_ldist_lclimb)
par(mfrow = c(2, 2))
plot(lm_ldist_lclimb)


### Linear model with categorical variables and interaction
# Starting point: we create a variable for exercise 
qdist <- with(nihills, quantile(dist, prob = c(1/3, 2/3)))

nihills$distdisc <- rep("Low", n)
nihills$distdisc[nihills$dist > qdist[1] & nihills$dist <= qdist[2]] <- "Med"
nihills$distdisc[nihills$dist > qdist[2]] <- "High"
nihills$distdisc <- factor(nihills$distdisc, level = c("Low", "Med", "High"))  

model.matrix(~ distdisc, data = nihills)[1:10,]

# Any difficulty?
lm_interaction <- lm(time ~ distdisc * climb, data = nihills)



#   ____________________________________________________________________________
#   Model selection                                                         ####

# Comparison summary() and anova() using the first model 
summary(lm_dist)
anova(lm_dist)


# Comparison between nested models
anova(lm_dist,lm_distclimb)
summary(lm_distclimb)$coefficients[3,3]^2

anova(lm_dist,lm_distclimb2)


# Information criterion 
AIC(lm_dist)
-2 * logLik(lm_dist)[1]    + 2 * (length(coef(lm_dist)) + 1) 

BIC(lm_dist)
-2 * logLik(lm_dist)[1]  + log(nrow(nihills)) * (length(coef(lm_dist)) + 1) 

# Model comparison using AIC

AIC <- rbind(AIC(lm_dist),  AIC(lm_dist2), AIC(lm_distclimb), AIC(lm_distclimb2))
BIC <- rbind(AIC(lm_dist, k = log(n)), AIC(lm_dist2, k = log(n)),
             AIC(lm_distclimb, k = log(n)), AIC(lm_distclimb2, k = log(n)))
cbind(AIC, BIC)


AIC(lm_ldist_lclimb) + sum(2*log(nihills$time))



# Extra: computational aspects
set.seed(123)
error <- rnorm(1000)
X <- cbind(1, matrix(rnorm(1000 * 9), 1000, 9))
beta <- rnorm(10, 0, 2)
y <- X %*% beta + error

library(microbenchmark)
microbenchmark(t(X) %*% X, crossprod(X))

w1 <- function(){
  XTX <- t(X) %*% X; return(solve(XTX) %*% t(X) %*% y)
}
w2 <- function(){
  XTX <- crossprod(X); return(solve(XTX) %*% t(X) %*% y)
}
w3 <- function(){
  XTX <- crossprod(X); return(solve(XTX) %*% (t(X) %*% y))
}
w4 <- function(){
  XTX <- crossprod(X); return(solve(XTX, t(X) %*% y)) 
}
microbenchmark(w1(), w2(), w3(), w4())
