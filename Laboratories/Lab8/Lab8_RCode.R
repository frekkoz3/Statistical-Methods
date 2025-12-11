#################################
# R code for LAB8 - 01/12/2025  #
#                               #
#################################


#   ____________________________________________________________________________
#   Regularisation techniques                                               ####

### Prostate data example ----
load("Prostate.rda") 
str(Prostate)
summary(Prostate)
table(Prostate$gleason)
Prostate$gleason <- as.numeric(factor(Prostate$gleason, 
                              labels = c("Low", "High","High","High")))

pairs(Prostate[,-c(5,7)], col = 1 + Prostate$svi, 
      pch = Prostate$gleason,
      main = paste("Prostate data, n = ", nrow(Prostate)))
par(mfrow = c(1, 2))
plot(lpsa ~ as.factor(gleason), data = Prostate)
plot(lpsa ~ as.factor(svi), data = Prostate)

# Linear model 
fit.lin <- lm(lpsa ~ ., data = Prostate) 
par(mfrow = c(1,2))
plot(fit.lin, which = c(1, 2))
fit1s <- summary(fit.lin) 
fit1s

# VIF index (multicollinearity)
library(car) 
vif(fit.lin)

### Ridge regression ----
library(glmnet) 

# Data preparation: glmnet requires the data to be passed in a specific way 
n <- dim(Prostate)[1] 
p <- dim(Prostate)[2]
X <- as.matrix(Prostate[, 1 : 8]) 
y <- Prostate[, p] 

# Fit ridge
fit_ridge <- glmnet(X, y, alpha = 0)

# Explore the call to the procedure 
print(fit_ridge)

# See the regression coefficients when $lambda=0.1$
round(t(coef(fit_ridge, s = 0.1)), 4)

par(mfrow=c(1,1))
# Plot the coefficients against the L1-norm
plot(fit_ridge, cex.lab = 1.5)
#abline(v=sum(abs(t(coef(fit_ridge, s = 0.1)))[-1]))
#abline(h = 0.4781)
# Cross-validation procedure to select lambda
cvfit_ridge <- cv.glmnet(X, y, alpha = 0)
cvfit_ridge

# See the MSE vs log(lambda)
plot(cvfit_ridge)

cvfit_ridge$lambda.min
cvfit_ridge$lambda.1se

# Let's compare the coefficients 
cbind(coef(fit.lin), coef(cvfit_ridge, s = "lambda.min"), coef(cvfit_ridge, s = "lambda.1se"))

### LASSO regression ----
fit_lasso <- glmnet(X, y, alpha = 1)
print(fit_lasso)
plot(fit_lasso, cex.lab = 1.4)
cvfit_lasso <- cv.glmnet(X,y, alpha = 1)
cvfit_lasso

plot(cvfit_lasso)

round(t(coef(cvfit_lasso, s = "lambda.min")), 4)
round(t(coef(cvfit_lasso, s = "lambda.1se")), 4)

cvfit_lasso <- cv.glmnet(X,y, alpha = 1)
cvfit_lasso

plot(cvfit_lasso)

round(t(coef(cvfit_lasso, s = "lambda.min")), 4)
round(t(coef(cvfit_lasso, s = "lambda.1se")), 4)

res <- cbind(coef(fit.lin), coef(cvfit_ridge, s = "lambda.min"),
             coef(cvfit_lasso, s = "lambda.min"), coef(cvfit_lasso, s = "lambda.1se"))
colnames(res) <- c("LS", "RIDGE", "LASSOmin", "LASSO1se")
res

# Model performance on a test set
set.seed(23)
idx_train <- sample(1 : nrow(Prostate), 0.90 * nrow(Prostate), replace = FALSE)
y_train <- y[idx_train]
y_test <- y[-idx_train]
X_train <- X[idx_train,]
X_test <- X[-idx_train,]
fit.lin <- lm(lpsa ~ ., data = Prostate[idx_train,])
pred_lm <- as.numeric(predict(fit.lin, newdata = Prostate[-idx_train,])) 
cvfit_ridge <- cv.glmnet(X_train, y_train, alpha = 0)
cvfit_ridge
pred_ridge1se <- as.numeric(predict(cvfit_ridge, newx = X_test, s = "lambda.1se"))
pred_ridgemin <- as.numeric(predict(cvfit_ridge, newx = X_test, s = "lambda.min"))
cvfit_lasso <- cv.glmnet(X_train, y_train, alpha = 1)
cvfit_lasso
pred_lasso1se <- as.numeric(predict(cvfit_lasso, newx = X_test, s = "lambda.1se"))
pred_lassomin <- as.numeric(predict(cvfit_lasso, newx = X_test, s = "lambda.min"))
MSE <- function(y, pred_y){
  return(mean((y-pred_y)^2))
}
MSE_lm <- MSE(y_test, pred_lm)
MSE_ridge1se <- MSE(y_test, pred_ridge1se)
MSE_lasso1se <- MSE(y_test, pred_lasso1se)
MSE_ridgemin <- MSE(y_test, pred_ridgemin)
MSE_lassomin <- MSE(y_test, pred_lassomin)
c(MSE_lm, MSE_ridge1se, MSE_lasso1se)
c(MSE_lm, MSE_ridgemin, MSE_lassomin)



#   ____________________________________________________________________________
#   EXTRA                                                                   ####

### Diabetes example ----
library(mlbench); data("PimaIndiansDiabetes2", package = "mlbench")
str(PimaIndiansDiabetes2)
apply(is.na(PimaIndiansDiabetes2), 2, sum)
PimaIndiansDiabetes2 <- na.omit(PimaIndiansDiabetes2)

# Fit the GLM
mod_diabetes <- glm(diabetes ~ ., family = binomial, data = PimaIndiansDiabetes2)
summary(mod_diabetes)

# Predict the probability of having diabetes for a certain woman 
predict(mod_diabetes, 
        newdata = data.frame(pregnant= 1, 
                             glucose = 102, pressure = 80, triceps = 68, 
                             insulin = 156, mass = 33, pedigree = 0, age = 60), type = "response")
# Ridge and lasso regression
n <- dim(PimaIndiansDiabetes2)[1]
p <- dim(PimaIndiansDiabetes2)[2]
X <- as.matrix(PimaIndiansDiabetes2[, 1 : 8])
y <- PimaIndiansDiabetes2[, p]
cvfit_ridge <- cv.glmnet(X, y, alpha = 0, family = "binomial")
cvfit_ridge
plot(cvfit_ridge)

# Lasso regression 
cvfit_lasso <- cv.glmnet(X,y, alpha = 1, family = "binomial")
cvfit_lasso

plot(cvfit_lasso)

res <- cbind(coef(mod_diabetes), coef(cvfit_ridge, s = "lambda.min"), 
             coef(cvfit_ridge, s = "lambda.1se"),
             coef(cvfit_lasso, s = "lambda.min"), 
             coef(cvfit_lasso, s = "lambda.1se"))
colnames(res) <- c("LS", "RIDGEmin", "RIDGE1se", "LASSOmin", "LASSO1se")
res


#   ____________________________________________________________________________
#   Spline functions                                                        ####

# Boston housing  
library(splines)
library(ggplot2)
library(MASS)
data(Boston) # From MASS package
str(Boston)
summary(Boston)

# Scatterplot 
ggplot(Boston, aes(lstat, medv)) + geom_point() + theme_bw()

# Simple linear model 
fit <- lm(medv ~ lstat, data = Boston)

# Residuals 
par(mfrow = c(1, 2))
plot(fit, which = c(1, 2))
hist(fit$residuals, breaks = 50)


# Degree-2 polynomial regression
fit.poly2 <- lm(medv ~ lstat + I(lstat^2), data = Boston)
fit.poly2 <- lm(medv ~ poly(lstat, 2, raw = TRUE), data = Boston) 
par(mfrow = c(1,2))
plot(fit.poly2,  which = c(1,2))

# Degree-5 polynomial regression
fit.poly5 <- lm(medv ~ poly(lstat, 5, raw = TRUE), data = Boston)
par(mfrow=c(1,2))
plot(fit.poly5, which = c(1,2))

# fitted curve
ggplot(Boston, aes(lstat, medv)) + 
  geom_point() + theme_bw() + 
  stat_smooth(method = lm, formula = y ~ x, col = "red") +
  stat_smooth(method = lm, formula = y ~ poly(x, 2, raw = TRUE), col = "green") +
  stat_smooth(method = lm, formula = y ~ poly(x, 5, raw = TRUE))

# cubic regression splines (3 knots - equivalent formulations)
knots <- quantile(Boston$lstat)[2 : 4]
fit.spline <- lm(medv ~ bs(lstat, knots = knots), data = Boston)
summary(fit.spline)
fit.spline <- lm(medv ~ bs(lstat, df = 6), data = Boston)
summary(fit.spline)

# Compare the fitted curve
ggplot(Boston, aes(lstat, medv)) + 
  geom_point() + theme_bw() +
  stat_smooth(method = lm, formula = y ~ poly(x, 2, raw = TRUE) , se = FALSE)+
  stat_smooth(method = lm, formula = y ~ bs(x, df = 6), col = "red", se = FALSE) +
  stat_smooth(method = lm, formula = y ~ bs(x, df = 100), col = "green", se = FALSE ) 

# linear, quadratic and cubic regression splines 
knot <- knots[2]
# knot <- median(Boston$lstat)
par(mfrow=c(1,2))
ggplot(Boston, aes(lstat, medv)) + 
  geom_point() + theme_bw() + 
  stat_smooth(method = lm, formula = y ~ bs(x, knots = knot, degree = 1), 
              col = "green", se = FALSE ) +  
  stat_smooth(method = lm, formula = y ~ bs(x, knots = knot, degree = 2), 
              col = "red", se = FALSE) + 
  stat_smooth(method = lm, formula = y ~ bs(x, knots = knot), se = FALSE) 

c(AIC(lm(lstat ~ bs(medv, knots = knot, degree = 1), data = Boston)), 
  AIC(lm(lstat ~ bs(medv, knots = knot, degree = 2), data = Boston)),
  AIC(lm(lstat ~ bs(medv, knots = knot, degree = 3), data = Boston)))



#   ____________________________________________________________________________
#   GAMs                                                                    ####

### Ozone Example ----
library(mgcv)
library(PerformanceAnalytics)
library(faraway)
data(ozone)
ozone$log_O3 <- log(ozone$O3) 
summary(ozone)

# Restrict the analysis to temperature, visibility, wind speed and day of the year
chart.Correlation(ozone[, c("vis", "doy", "vh", "wind", "log_O3")])

# Model 1
gamfit <- gam(log_O3 ~ s(temp) + s(vis) + s(doy) + s(wind), data = ozone)
summary(gamfit)
length(gamfit$coef)
gamfit$sp

par(mfrow=c(1,2))
plot(gamfit, select = 1, scale = FALSE)
plot(gamfit, select = 2, scale = FALSE)

par(mfrow=c(1,2))
plot(gamfit, select = 3, scale = FALSE)
plot(gamfit, select = 4, scale = FALSE)

plot(gamfit, residuals = TRUE, pch = 19, pages = 1)

# Part of the plot obtained via gam routines 
pr_data <- cbind(ozone$temp, residuals(gamfit, type="pearson") +
                   gam(log_O3 ~ s(temp) + s(vis) + s(doy) + s(wind), data = ozone,
                       fit=FALSE)$X[,2:10] %*%coef(gamfit)[2:10]) 
par(mfrow = c(1,2))
plot(pr_data[,1],pr_data[,2], pch = 19)
rug(pr_data[,1])
plot(gamfit, residuals = TRUE, pch = 19, select = 1, se = FALSE)


# Model 2
gamfit2 <- gam(log(O3)~s(temp) + s(vis) + s(doy) + wind, data = ozone)
summary(gamfit2)

# Model 3
glmfit_viagam <- gam(log(O3) ~ temp + vis + doy + wind, data = ozone)
summary(glmfit_viagam)

# Model 4
glmfit <- glm(log(O3) ~ temp + vis + doy + wind, data = ozone)
summary(glmfit)

# Consider all the covariates 
gamfit_ext <- gam(log(O3) ~ s(vh) + s(wind) + s(humidity) + s(temp) + s(ibh) + 
                    s(dpg) + s(ibt) + s(vis) + s(doy), data = ozone)
summary(gamfit_ext)

round(gamfit_ext$sp,4)
par(mfrow = c(1, 3))
par(mfrow=c(3,3))
for(j in 1 : 9) plot(gamfit_ext, select = j, scale = FALSE)

# Residuals
plot(gamfit_ext, residuals = TRUE, pch = 19)

# Extended model 2
gamfit_ext2 <- gam(log(O3) ~ vh + wind + s(humidity) + s(temp) + s(ibh) + 
                     s(dpg) + ibt + s(vis) + s(doy), data=ozone)
summary(gamfit_ext2)

# Extended model 3 (GLM)
gamfit_ext3 <- gam(log(O3) ~ vh + wind + humidity + temp + ibh + 
                     dpg + ibt + vis + doy, data=ozone)
summary(gamfit_ext3)

cbind(AIC(gamfit_ext,gamfit_ext2,gamfit_ext3))


# Smoothing parameters selection
sp <- gamfit_ext$sp # smoothing parameter obtained for the extended model 1

# Let's define a grid of value multiplying the sp 
tuning.scale<-c(1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3, 1e4, 1e5)
scale.exponent <- log10(tuning.scale)

n.tuning <- length(tuning.scale)
min2ll <- aic <- bic <- gcv <- rep(NA, n.tuning)

# By means of sp argument of the gam function we are fixing the sp to a specific value 
# For each value of the grid we optimise the penalized log-likelihood by fixing the sp
# In the end we extract the values of -2nll, AIC, BIC and GCV
for (i in 1:n.tuning) {
  gamfit_opt<-gam(log(O3) ~ s(vh) + s(wind) + s(humidity) + s(temp) + s(ibh)+
                    s(dpg) + s(ibt) + s(vis) + s(doy), data = ozone, sp = tuning.scale[i] * sp)
  min2ll[i]<- -2 * logLik(gamfit_opt)
  aic[i] <- AIC(gamfit_opt)
  bic[i] <- BIC(gamfit_opt)
  gcv[i] <- gamfit_opt$gcv.ubre
}

par(mfrow = c(2, 2)) 
plot(scale.exponent, min2ll, type = "b",main = "-2 log likelihood")
plot(scale.exponent, aic, type = "b", main = "AIC")
plot(scale.exponent, bic, type = "b", main = "BIC")
plot(scale.exponent, gcv, type = "b", main = "GCV")

# A comparison between the model selected via BIC e the default (GCV)
min.bic <- 1e100
opt.tuning.scale <- NULL
for (i in 1 :n.tuning) {
  if (bic[i] < min.bic) {
    min.bic <- bic[i]
    opt.tuning.scale <- tuning.scale[i]
  }
}
opt.sp <- opt.tuning.scale * sp  # Smoothing parameter selected via BIC

gamobj <- gam(log(O3) ~ s(vh) + s(wind) + s(humidity) + s(temp) + s(ibh) + 
                s(dpg) + s(ibt) + s(vis) + s(doy),
              data = ozone, sp = opt.sp)
summary(gamobj)

summary(gamfit_ext)

AIC(gamobj, gamfit_ext)

# Splitting into train and test 
ozoneTrain <- ozone[ozone$doy <= 360,] 
ozoneTest <- ozone[ozone$doy > 360,]
y_train <- ozoneTrain$log_O3
y_test <- ozoneTest$log_O3
X_train <- as.matrix(ozoneTrain[,-c(1,11)])
X_test <- as.matrix(ozoneTest[,-c(1,11)])

# Extended model 1
gamfit_ext <- gam(log(O3) ~ s(vh) + s(wind) + s(humidity) + s(temp) + s(ibh) + 
                    s(dpg) + s(ibt) + s(vis) + s(doy), data = ozoneTrain)
# Extended model 2
gamfit_ext2 <- gam(log(O3) ~ vh + wind + s(humidity) + s(temp) + s(ibh) + 
                     s(dpg) + ibt + s(vis) + s(doy), data=ozoneTrain)

# Extended model 3 (GLM)
gamfit_ext3 <- gam(log(O3) ~ vh + wind + humidity + temp + ibh + 
                     dpg + ibt + vis + doy, data=ozoneTrain)
cvfit_ridge <- cv.glmnet(X_train, y_train, alpha = 0) #Ridge
cvfit_lasso <- cv.glmnet(X_train, y_train, alpha = 1) #LASSO

# Prredicted values 
pred_GAM1 <- predict(gamfit_ext, newdata = ozoneTest)
pred_GAM2 <- predict(gamfit_ext2, newdata = ozoneTest)
pred_GAM3 <- predict(gamfit_ext3, newdata = ozoneTest)
pred_ridgemin <- as.numeric(predict(cvfit_ridge, newx = X_test, s = "lambda.min"))
pred_lassomin <- as.numeric(predict(cvfit_lasso, newx = X_test, s = "lambda.min"))

# MSE 
MSE_GAM1 <- MSE(ozoneTest$log_O3, pred_GAM1)
MSE_GAM2 <- MSE(ozoneTest$log_O3, pred_GAM2)
MSE_GAM3 <- MSE(ozoneTest$log_O3, pred_GAM3)
MSE_ridgemin <- MSE(y_test, pred_ridgemin)
MSE_lassomin <- MSE(y_test, pred_lassomin)

c(MSE_GAM1, MSE_GAM2, MSE_GAM3, MSE_ridgemin,  MSE_lassomin)

par(mfrow= c(1,1))
plot(ozoneTest$doy, ozoneTest$log_O3, ylim = c(0,2.5))
points(ozoneTest$doy, pred_GAM2, col = "blue", pch = 20)
points(ozoneTest$doy, pred_GAM3, col = "red", pch = 21)

### Electricity Demand Example ----
load("GEF14_d4.Rdata")
str(GEF14_d4[, c("dow", "doy", "temp95_h17", "load24_h17", "load_h17")])
summary(GEF14_d4[, c("temp95_h17", "load24_h17", "load_h17")])
chart.Correlation(GEF14_d4[, c("doy", "temp95_h17", "load24_h17", "load_h17")])
plot(load_h17 ~ factor(dow, labels= c("Sun", "Mon", "Thu", "Wed", "Thur", "Fri", "Sat")), 
     xlab = "Dow", data = GEF14_d4)


gam1_gef14 <- gam(load_h17 ~ dow + load24_h17 + s(temp95_h17) + s(doy), data = GEF14_d4)
summary(gam1_gef14)
summary(gam1_gef14)$p.table
summary(gam1_gef14)$s.table

# A model with a higher number of basis functions for doy and using cubic ciclyc-splines
gam2_gef14 <- gam(load_h17 ~ dow + load24_h17 + 
                    s(temp95_h17) + s(doy, k = 20), data = GEF14_d4)
summary(gam2_gef14)

# As above but considering cubic ciclyc splines
gam3_gef14 <- gam(load_h17 ~ dow + load24_h17 + 
                    s(temp95_h17) + s(doy, k = 20, bs = "cc"), data = GEF14_d4)
summary(gam3_gef14)

# Just a linear model 
glm_gef14 <- gam(load_h17 ~ dow + load24_h17 +temp95_h17 + doy, data = GEF14_d4)
summary(glm_gef14)

# Comparison 
AIC(gam1_gef14, gam2_gef14, gam3_gef14, glm_gef14)

# Splitting into train and test 
GEF14_d4_Train <- GEF14_d4[GEF14_d4$year != 2011 | GEF14_d4$doy<=300,] 
GEF14_d4_Test <- GEF14_d4[GEF14_d4$year == 2011 & GEF14_d4$doy > 300,]


# Extended model 1
gam1_gef14 <- gam(load_h17 ~ dow + load24_h17 + s(temp95_h17) + s(doy), data = GEF14_d4_Train)

# Extended model 2
gam2_gef14 <- gam(load_h17 ~ dow + load24_h17 + 
                    s(temp95_h17) + s(doy, k = 20), data = GEF14_d4_Train)
# Extended model 3 
gam3_gef14 <- gam(load_h17 ~ dow + load24_h17 + 
                    s(temp95_h17) + s(doy, k = 20, bs = "cc"), data = GEF14_d4_Train)


# Just a linear model 
glm_gef14 <- gam(load_h17 ~ dow + load24_h17 +temp95_h17 + doy, data = GEF14_d4_Train)


pred_GAM1 <- predict(gam1_gef14, newdata = GEF14_d4_Test)
pred_GAM2 <- predict(gam2_gef14, newdata = GEF14_d4_Test)
pred_GAM3 <- predict(gam3_gef14, newdata = GEF14_d4_Test)
pred_GAM4 <- predict(glm_gef14, newdata = GEF14_d4_Test)

MSE_GAM1 <- MSE(GEF14_d4_Test$load_h17, pred_GAM1)
MSE_GAM2 <- MSE(GEF14_d4_Test$load_h17, pred_GAM2)
MSE_GAM3 <- MSE(GEF14_d4_Test$load_h17, pred_GAM3)
MSE_GAM4 <- MSE(GEF14_d4_Test$load_h17, pred_GAM4)

c(MSE_GAM1,MSE_GAM2,MSE_GAM3,MSE_GAM4)
#[1] 120.1685 118.3201 115.8620 470.0182


#   ____________________________________________________________________________
#   EXTRA                                                                   ####

### Example Retynopathy ----
library(gss)
data(wesdr)
str(wesdr)
wesdr$ret <- as.factor(wesdr$ret)
par(mfrow = c(1,3))
plot(dur ~ ret, data = wesdr)
plot(dur ~ ret, data = wesdr)
plot(dur ~ ret, data = wesdr)

# GAMs
gam_log <- gam(ret ~ s(dur) + s(gly) + s(bmi), family = binomial(), data = wesdr)
summary(gam_log)

# Visualise the effects
par(mfrow = c(1,3))
for(j in 1 : 3) {
  plot(gam_log, select = j)
}

# Removing the spline for gly
gam_log2 <- gam(ret ~ s(dur) + gly + s(bmi), family = binomial(), data = wesdr)
summary(gam_log2)

# GLM
gam_log3 <- gam(ret ~ dur + gly + bmi, family = binomial(), data = wesdr)
summary(gam_log3)

# AIC
AIC(gam_log, gam_log2, gam_log3)

# Correct classification on all the data 
sum(diag(table(wesdr$ret,predict(gam_log, type = "response", newdata = wesdr)>0.5)))/nrow(wesdr)
sum(diag(table(wesdr$ret,predict(gam_log2, type = "response", newdata = wesdr)>0.5)))/nrow(wesdr)
sum(diag(table(wesdr$ret,predict(gam_log3, type = "response", newdata = wesdr)>0.5)))/nrow(wesdr)

# Assess model performance on a test set 
set.seed(1)
idx_train <- sample(1:nrow(wesdr), nrow(wesdr)*0.75, replace = FALSE)
wesdr_train <- wesdr[idx_train,]
wesdr_test <- wesdr[-idx_train,]

gam_log <- gam(ret ~ s(dur) + s(gly) + s(bmi), family = binomial(), data = wesdr_train)
gam_log2 <- gam(ret ~ dur + gly + bmi, family = binomial(), data = wesdr_train)
gam_log3 <- gam(ret ~ s(dur) + s(bmi), family = binomial(), data = wesdr_train)


table(wesdr_test$ret,predict(gam_log, type = "response", newdata = wesdr_test)>0.5)
table(wesdr_test$ret,predict(gam_log2, type = "response", newdata = wesdr_test)>0.5)
table(wesdr_test$ret,predict(gam_log3, type = "response", newdata = wesdr_test)>0.5)

library(ROCR)
pred1 <- prediction(as.numeric(predict(gam_log, type = "response", newdata = wesdr_test)), 
                    wesdr_test$ret)
pred2 <- prediction(as.numeric(predict(gam_log2, type = "response", newdata = wesdr_test)), 
                    wesdr_test$ret)
pred3 <- prediction(as.numeric(predict(gam_log3, type = "response", newdata = wesdr_test)), 
                    wesdr_test$ret)

# TPR and FPR
perf1 <- performance(pred1, "tpr", "fpr")
perf1@x.values[[1]][55] # FPR
perf1@y.values[[1]][55] # TPR
perf1@alpha.values[[1]][55] 
11/(11+76) # FPR
43/(38 + 43) #TPR

perf2 <- performance(pred2, "tpr", "fpr")
perf3 <- performance(pred3, "tpr", "fpr")
# ROC curve
plot(perf1, col = "darkblue", lwd = 2, main = "ROC Curve with ROCR")
plot(perf2, col = "darkred", lwd = 2, add = TRUE)
plot(perf3, col = "darkgreen", lwd = 2, add = TRUE)
abline(a = 0, b = 1, col = "red", lty = 2)

