#################################
# R code for LAB7 - 20/11/2025  #
#                               #
#################################


#   ____________________________________________________________________________
#   Logistic regression                                                     ####

### Wells data ----
data <- read.csv2("Wells.csv", header = TRUE)
str(data)
summary(data)
par(mfrow = c(1, 2))
boxplot(dist ~ switch,
  data = data, horizontal = TRUE,
  xlab = "Distance (in meters)", ylab = "Switch"
)
boxplot(arsenic ~ switch,
  data = data, horizontal = TRUE,
  xlab = "Arsenic (in hundreds of mg per liter)", ylab = "Switch"
)

par(mfrow = c(1, 1))
barplot(prop.table(table(data$switch, data$educ), margin = 1),
  beside = TRUE, legend.text = c("no switch", "switch"), xlab = "Education (years)"
)

### Single predictor (distance): logit(p_i) = beta_1 + beta_2 * dist_i

fit1 <- glm(switch ~ dist, data = data, family = binomial(link = "logit"))
res <- summary(fit1)
res

# Results of model output
X <- model.matrix(~dist, data = data)

# Std.error
res$coefficients[, 2]
sqrt(diag(solve((t(X) %*% diag(fit1$weights) %*% X) / fit1$family$dispersion)))

fit1$weights[1] # First element of the W matrix
predict(fit1)[1]
help(predict.glm)
hat_mu1 <- predict(fit1, type = "response")[1]
hat_mu1 * (1 - hat_mu1)

# z value
z.val <- res$coefficients[, 1] / res$coefficients[, 2]
z.val

# Pr(>|z|)
p.val <- 2 * pnorm(abs(z.val), lower = FALSE)
p.val

# Confidence interval for beta2
coef(fit1)[2] + c(-1, 1) * summary(fit1)$coefficients[2, 2] * qnorm(0.975)
confint.default(fit1, data = data)[2, ]

# Deviance residuals (only the first)
with(data, sign(switch[1
] - hat_mu1) * sqrt(2 * switch[1
] * log(switch[1
] / hat_mu1)))
residuals(fit1, type = "deviance")[1]

# Null deviance
fit0 <- glm(switch ~ 1, data = data, family = binomial(link = "logit"))
sum(residuals(fit0, type = "deviance")^2)
with(data, 2 * (sum(switch[switch == 1
] * log(1 / mean(switch))) +
  sum((1 - switch[switch == 0
  ]) * log(1 / (1 - mean(switch))))))

# Residual deviance
pred1 <- predict(fit1, type = "response")
sum(residuals(fit1, type = "deviance")^2)
with(data, 2 * (sum(switch[switch == 1
] * log(1 / pred1[switch == 1])) +
  sum((1 - switch[switch == 0
  ]) * log(1 / (1 - pred1[switch == 0])))))

# anova() function
anova(fit1)

# AIC
AIC(fit1)
-2 * logLik(fit1) + 2 * length(fit1$coefficients)


# Probability of switching well if you are very close to a safe well (dist = 0)
InvLogit <- function(eta) 1 / (1 + exp(-eta))

as.numeric(InvLogit(c(1, 0) %*% coef(fit1)))
predict(fit1, newdata = data.frame(dist = 0), type = "response")

# OR: the odds when dist = c+ 1 is OR times the odds when dist = c
exp(coef(fit1)[2])

# Difference on the probability of switch well due to an increase
# of one unit on the predictor from the mean
as.numeric(InvLogit(c(1, mean(data$dist) + 1) %*% coef(fit1)) -
  InvLogit(c(1, mean(data$dist)) %*% coef(fit1)))

# Equivalently we can use the predict function
predict(fit1, newdata = data.frame(dist = mean(data$dist) + 1), type = "response") -
  predict(fit1, newdata = data.frame(dist = mean(data$dist)), type = "response")

# Difference on the probability of switch well between 300 m and 0 m
as.numeric(InvLogit(c(1, 300) %*% coef(fit1)) - InvLogit(c(1, 0) %*% coef(fit1)))


# Difference on the probability of switch well due to an increase
# of one unit on the predictor from the 99th percentile
as.numeric(InvLogit(c(1, quantile(data$dist, 0.99) + 1) %*% coef(fit1)) -
  InvLogit(c(1, quantile(data$dist, 0.99)) %*% coef(fit1)))
res
# Scaling the dist variable and fit again the model: better interpretation
data$dist100 <- data$dist / 100
fit2 <- glm(switch ~ dist100, data = data, family = binomial(link = "logit"))
summary(fit2)

# Difference on probability scale due to an increase of 1 unit
# on the predictor  dist100 from the mean
as.numeric(InvLogit(c(1, mean(data$dist100) + 1) %*% coef(fit2)) -
  InvLogit(c(1, mean(data$dist100)) %*% coef(fit2)))

# Note that is the same difference due to an increase of 100 unit on the predictor
# dist from the mean (i.e. considering fit1)
as.numeric(InvLogit(c(1, mean(data$dist) + 100) %*% coef(fit1)) -
  InvLogit(c(1, mean(data$dist)) %*% coef(fit1)))

# Plot of the fitted model
switch.jitter <- jitter(data$switch, factor = 0.10)

par(mfrow = c(1, 1))
with(data, plot(dist, switch.jitter,
  xlab = "Distance (in meters) to nearest safe well",
  ylab = "Pr(Swiching)", type = "n"
))
curve(InvLogit(coef(fit1)[1] + coef(fit1)[2] * x), lwd = 1, add = TRUE)
points(data$dist, switch.jitter, pch = 20, cex = .1)

### Adding predictors (arsenic) ----
par(mfrow = c(1, 1))
hist(data$arsenic,
  freq = TRUE, xlab = "Arsenic concentration",
  breaks = seq(0, 0.25 + max(data$arsenic), 0.25), main = ""
)

fit3 <- glm(switch ~ dist100 + arsenic, data = data, family = binomial("logit"))
summary(fit3)

par(mfrow = c(1, 2))
# Pr(Switching) vs distance (consider Arsenic = 0.5 and Arsenic = 1.5)
with(data, plot(dist100, switch.jitter, xlab = "Distance", ylab = "Pr(Switching)", type = "n"))
curve(InvLogit(coef(fit3)[1] + coef(fit3)[2] * x + coef(fit3)[3] * 1.5), add = TRUE, col = "blue")
curve(InvLogit(coef(fit3)[1] + coef(fit3)[2] * x + coef(fit3)[3] * 0.5), add = TRUE, col = "red")
points(data$dist100, switch.jitter, pch = 20, cex = .1)
text(0.50, .27, "if As = 0.5", adj = 0, cex = .8, col = "red")
text(0.75, .50, "if As = 1.5", adj = 0, cex = .8, col = "blue")
# Pr(Switching) vs Arsenic (Consider distance = 0 and distance = 50)
with(data, plot(arsenic, switch.jitter, xlab = "Arsenic", ylab = "Pr(Switching)", type = "n"))
curve(InvLogit(coef(fit3)[1] + coef(fit3)[3] * x), add = TRUE, col = "red")
text(1.5, 0.8, "if dist=0", col = "red")
curve(InvLogit(coef(fit3)[1] + coef(fit3)[2] * 0.5 + coef(fit3)[3] * x), add = TRUE, col = "blue")
text(4, 0.6, "if dist=50", col = "blue")
points(data$arsenic, switch.jitter, cex = 0.01, pch = 20)

### Including interaction terms ----
fit4 <- glm(switch ~ dist100 * arsenic,
  data = data, family = binomial(link = "logit")
)

# Equivalently
# fit4 <- glm(switch ~ dist100 + arsenic + dist100:arsenic,
#            data = data, family = binomial(link = "logit"))
summary(fit4)

# Probability of switching for an average distance and
# an average concentration of arsenic
as.numeric(InvLogit(c(
  1, mean(data$dist100), mean(data$arsenic),
  mean(data$dist100) * mean(data$arsenic)
) %*% coef(fit4)))
# Equivalently
# predict(fit4, newdata=data.frame(dist100 = mean(data$dist100),
#                                 arsenic = mean(data$arsenic)),
#        type = "response")


# Difference in probability due to an increase of 100 meters
# in the distance (from the mean), keeping arsenic at its mean
mean_d100 <- mean(data$dist100)
mean_ars <- mean(data$arsenic)
as.numeric(InvLogit(c(1, mean_d100 + 1, mean_ars, (mean_d100 + 1) * mean_ars) %*% coef(fit4)) -
  InvLogit(c(1, mean_d100, mean_ars, mean_d100 * mean_ars) %*% coef(fit4)))


# Difference in probability due to an increase of 1 unit
# in the arsenic (from the mean), keeping distance at its mean
as.numeric(InvLogit(c(1, mean_d100, mean_ars + 1, mean_d100 * (mean_ars + 1)) %*% coef(fit4)) -
  InvLogit(c(1, mean_d100, mean_ars, mean_d100 * mean_ars) %*% coef(fit4)))

# Difference in logit probability scale when arsenic increases by 1 unit,
# considering dist at its mean level
coef(fit4)[3] + coef(fit4)[4] * mean(data$dist100)

# Difference in logit probability scale when dist increases by 100 meters,
# considering arsenic at its mean level
coef(fit4)[2] + coef(fit4)[4] * mean(data$arsenic)

par(mfrow = c(1, 2))
# Pr(Switching) vs Distance (Consider Arsenic=0.5 and Arsenic = 1)
with(data, plot(data$dist, switch.jitter,
  type = "n",
  xlab = "Distance", ylab = "Pr(Switching)"
))
curve(InvLogit(cbind(1, x / 100, 1.5, 1.5 * x / 100) %*% coef(fit4)), add = TRUE, col = "red")
curve(InvLogit(cbind(1, x / 100, 0.5, 0.5 * x / 100) %*% coef(fit4)), add = TRUE, col = "blue")
text(100, 0.6, "if As = 1.5", col = "red")
text(35, 0.4, "if As = 0.5", col = "blue")
points(data$dist, switch.jitter, pch = 20, cex = 0.01)

# Pr(Switching) vs Arsenic (Consider Distance = 0 and Distance = 50)
with(data, plot(data$arsenic, switch.jitter,
  type = "n",
  xlab = "Arsenic concentration", ylab = "Pr(Switching)"
))
curve(InvLogit(cbind(1, 0, x, 0) %*% coef(fit4)), add = TRUE, col = "red")
curve(InvLogit(cbind(1, 0.5, x, 0.5 * x) %*% coef(fit4)), add = TRUE, col = "blue")
text(2, 0.9, "if Dist = 0", col = "red")
text(3, 0.6, "if Dist = 50", col = "blue")
points(data$arsenic, switch.jitter, pch = 20, cex = 0.01)



#   ____________________________________________________________________________
#   Model Comparison and Model Selection                                    ####

# Analysis of deviance

anova(fit4, test = "LRT") # equivalently test = "Chisq"

# AIC
c(AIC(fit0), AIC(fit1), AIC(fit3), AIC(fit4))


### Extra: Centering the variables and adding social predictors ----
# Cetering the variables
data$c.dist100 <- scale(data$dist100, scale = FALSE)
data$c.arsenic <- scale(data$arsenic, scale = FALSE)

fit5 <- glm(switch ~ c.dist100 + c.arsenic + c.dist100:c.arsenic,
  family = binomial(link = "logit"), data = data
)
summary(fit5)

# Adding social education (centered) predictor
data$c.educ4 <- scale(data$educ / 4, scale = FALSE)
fit6 <- glm(switch ~ c.dist100 * c.arsenic + c.educ4,
  data = data,
  family = binomial(link = "logit")
)
summary(fit6)

fit7 <- glm(switch ~ c.dist100 * c.arsenic + c.dist100 * c.educ4 + c.arsenic * c.educ4,
  data = data, family = binomial(link = "logit")
)
summary(fit7)

### Binned residuals ----
pred4 <- predict(fit4, type = "response")
res4 <- residuals(fit4, type = "response")

plot(c(0, 1), c(-1, 1),
  xlab = "Estimated Pr(Switching)", type = "n",
  ylab = "Observed - Estimated", main = "Residual Plot"
)
abline(h = 0, col = "gray", lwd = 0.5)
points(pred4, res4, pch = 20, cex = 0.2)
abline(c(0, -1), col = "red", lwd = 0.25)
abline(c(1, -1), col = "red", lwd = 0.25)


binned.resids <- function(x, y, nclass = sqrt(length(x))) {
  breaks.index <- floor(length(x) * (1:(nclass)) / nclass)
  breaks <- c(-Inf, sort(x)[breaks.index], Inf)
  output <- NULL
  xbreaks <- NULL
  x.binned <- as.numeric(cut(x, breaks))
  for (i in 1:nclass) {
    items <- (1:length(x))[x.binned == i]
    x.range <- range(x[items])
    xbar <- mean(x[items])
    ybar <- mean(y[items])
    n <- length(items)
    sdev <- sd(y[items])
    output <- rbind(output, c(xbar, ybar, n, x.range, 2 * sdev / sqrt(n)))
  }
  colnames(output) <- c("xbar", "ybar", "n", "x.lo", "x.hi", "2se")
  return(list(binned = output, xbreaks = xbreaks))
}


br4 <- binned.resids(pred4, res4, nclass = 50)$binned
plot(range(br4[, 1]), range(br4[, 2], br4[, 6], -br4[, 6]),
  type = "n",
  xlab = "Estimated Pr(Switching)", ylab = "Average residual", main = "Binned Residual Plot"
)
abline(h = 0, col = "gray", lwd = 0.5)
points(br4[, 1], br4[, 2], pch = 19, cex = 0.5)
lines(br4[, 1], br4[, 6], col = "gray", lwd = 0.5)
lines(br4[, 1], -br4[, 6], col = "gray", lwd = 0.5)

par(mfrow = c(1, 2))
br.dist <- binned.resids(data$dist, res4,
  nclass = 40
)$binned
plot(range(br.dist[, 1]), range(br.dist[, 2], br.dist[, 6], -br.dist[, 6]),
  type = "n",
  xlab = "Distance", ylab = "Average residual", main = "Binned residual plot"
)
abline(h = 0, col = "gray", lwd = 0.5)
lines(br.dist[, 1], br.dist[, 6], col = "gray", lwd = 0.5)
lines(br.dist[, 1], -br.dist[, 6], col = "gray", lwd = 0.5)
points(br.dist[, 1], br.dist[, 2], pch = 19, cex = 0.5)

br.arsenic <- binned.resids(data$arsenic, res4,
  nclass = 40
)$binned
plot(range(br.arsenic[, 1]), range(br.arsenic[, 2], br.arsenic[, 6], -br.arsenic[, 6]),
  type = "n",
  xlab = "Arsenic concentration", ylab = "Average residual", main = "Binned residual plot"
)
abline(h = 0, col = "gray", lwd = 0.5)
lines(br.arsenic[, 1], br.arsenic[, 6], col = "gray", lwd = 0.5)
lines(br.arsenic[, 1], -br.arsenic[, 6], col = "gray", lwd = 0.5)
points(br.arsenic[, 1], br.arsenic[, 2], pch = 19, cex = 0.5)


#   ____________________________________________________________________________
#   Error rate                                                              ####

# Error rate of the null model

fit0 <- glm(switch ~ 1, family = binomial(link = "logit"), data = data)
m0 <- mean(data$switch)
err.rate0 <- min(m0, 1 - m0)
err.rate0

# Error rate of model 4
err.rate4 <- mean((pred4 > 0.5 & data$switch == 0) | ((pred4 < 0.5 & data$switch == 1)))

# Confusion matrix
conf_mat <- table(pred4 > 0.5, data$switch)
conf_mat
sum(diag(conf_mat)) / nrow(data)
1 - sum(diag(conf_mat)) / nrow(data)
err.rate4


#   ____________________________________________________________________________
#   Poisson regression                                                      ####

accidents <- read.table("SingaporeAuto.txt", header = T)
accidents[1:5, ]

accidents$NCD <- as.factor(accidents$NCD)
accidents$AgeCat <- as.factor(accidents$AgeCat)
accidents$VAgeCat <- as.factor(accidents$VAgeCat)


table(accidents$Clm_Count)
barplot(table(accidents$Clm_Count))

table(accidents$Clm_Count) / nrow(accidents)
round(dpois(0:3, mean(accidents$Clm_Count)), 6)


# Fit Poisson regression including all the variables
model1 <- glm(Clm_Count ~ Female + AgeCat + VAgeCat + NCD,
  family = poisson, data = accidents
)
summary(model1)
str(model1)
# Interpretation
exp(coef(model1)[2])

exp(coef(model1)[2] + c(-1, 1) * summary(model1)$coefficients[2, 2] * qnorm(0.975))

exp(coef(model1)[3])

anova(model1, test = "LRT")

model0 <- glm(Clm_Count ~ 1, family = poisson, data = accidents)
model2 <- glm(Clm_Count ~ VAgeCat + NCD,
  family = poisson, data = accidents
)
model3 <- glm(Clm_Count ~ NCD,
  family = poisson, data = accidents
)
AIC(model0, model1, model2, model3)
