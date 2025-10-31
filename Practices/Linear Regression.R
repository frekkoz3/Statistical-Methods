load("linear regression simple examples/fev.RData")

summary(fev) # no missing data
dim(fev)
str(fev)

hist(fev$fev)
boxplot(fev$fev)

library(PerformanceAnalytics)

chart.Correlation(fev)

# using only the height
fit <- lm(fev~height, data = fev)
summary(fit)
par(mfrow = c(2, 2))
plot(fit)
# using the height plus the smoke as addictive variable
fit <- lm(fev~height+smoke, data = fev)
summary(fit)
par(mfrow = c(2, 2))
plot(fit)
#using the height plus the male as addictive variable
fit <- lm(fev~height+male, data = fev)
summary(fit)
par(mfrow = c(2, 2))
plot(fit)
#using the height plus height^2 plus the male as addictive variable
fit <- lm(fev~height+I(height^2)+male, data = fev)
summary(fit)
par(mfrow = c(2, 2))
plot(fit)
