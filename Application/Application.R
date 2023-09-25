# ==============================================================================
#      Application: Dimension Adaptive Estimation
# ==============================================================================

# Author: Chencheng Fang, ccfang@uni-bonn.de
# Supervisor: Prof. Dr. Joachim Freyberger

# Note: I have designed an R package 'dimada', which helps to implement dimension 
# adaptive estimation easily, and I have made it available on GitHub at ccfang2.

# -----------libraries ---------------------------------------------------------
options(warn=-1)
rm(list=ls())

# install package 'dimada' from my GitHub
# devtools::install_github("ccfang2/dimada")
library(dimada)

library(reticulate) 
library(stargazer)

# ----------source -------------------------------------------------------------
setwd("~/Desktop/bonn master/thesis/Application/data")

# ----------data ---------------------------------------------------------------
# dataset is sourced from Kaggle, and the original data is in npy format.
# https://www.kaggle.com/datasets/madhavmalhotra/car-price-regression-preprocessed
# X and y in R format is also available in 'ApplicationDataset.RData'.
np <- import("numpy")
X <- np$load("X.npy")
y<-np$load("y.npy")

names.X <- c("year", "mileage", "tax", "mpg", "engine","diesel","electricity", "hybrid", "petrol")
colnames(X) <- names.X

#dataset <- cbind(price=y,X)
#stargazer(dataset, summary=TRUE, type="latex", median=TRUE)

X<-X[,1:5]

# rescale dataset by z-score standardization
y.scale<- scale(y)
X.scale <-X;  X.scale[,1] <- X.scale[,1]-min(X.scale[,1])
X.scale[,1:5]<- apply(X.scale[,1:5], 2, scale)

# split dataset
set.seed(200)
train.index <- sample(1:nrow(X.scale), size=floor(nrow(X.scale)*0.6), replace = FALSE)
X.scale.train <- X.scale[train.index,]
X.scale.test <- X.scale[-train.index,]

y.scale.train <- y.scale[train.index]
y.scale.test <- y.scale[-train.index]

# -------------------ols estimator ---------------------------------------------
ols <- glm(y.scale.train~X.scale.train, family="gaussian")
(ols.mse <- mean((as.matrix(cbind(1,X.scale.test))%*%ols$coefficients-as.vector(y.scale.test))^2, na.rm=TRUE))

# --------------------additive estimator ---------------------------------------
# power series as basis
add.power.sieve <- poly.gen(data=X.scale.train, test.data = X.scale.test, max.interaction = 1, legendre=FALSE, n.basis = 20)

set.seed(200)
add.dimada.power <- dimada(y=as.vector(y.scale.train), x=NULL, x.sieve = add.power.sieve$train,
                       methods = c("Lasso","adaLasso"), family = "gaussian")

(add.power.Lasso.mse <- mean((as.matrix(cbind(1,add.power.sieve$test[,add.dimada.power$Lasso$coefs.final$terms, drop=FALSE]))%*%ifelse(is.na(add.dimada.power$Lasso$post.lm$coefficients),0,add.dimada.power$Lasso$post.lm$coefficients)-as.vector(y.scale.test))^2))
(add.power.adaLasso.mse <- mean((as.matrix(cbind(1,add.power.sieve$test[,add.dimada.power$adaLasso$coefs.final$terms, drop=FALSE]))%*%ifelse(is.na(add.dimada.power$adaLasso$post.lm$coefficients),0,add.dimada.power$adaLasso$post.lm$coefficients)-as.vector(y.scale.test))^2))

# legendre polynomials as basis
add.legendre.sieve <- poly.gen(data=X.scale.train, test.data = X.scale.test, max.interaction = 1, legendre = TRUE,n.basis=20)

set.seed(200)
add.dimada.legendre <- dimada(y=as.vector(y.scale.train), x=NULL, x.sieve = add.legendre.sieve$train,
                          methods = c("Lasso","adaLasso"), family = "gaussian")
                     
(add.legendre.Lasso.mse <- mean((as.matrix(cbind(1,add.legendre.sieve$test[,add.dimada.legendre$Lasso$coefs.final$terms, drop=FALSE]))%*%ifelse(is.na(add.dimada.legendre$Lasso$post.lm$coefficients),0,add.dimada.legendre$Lasso$post.lm$coefficients)-as.vector(y.scale.test))^2))
(add.legendre.adaLasso.mse <- mean((as.matrix(cbind(1,add.legendre.sieve$test[,add.dimada.legendre$adaLasso$coefs.final$terms, drop=FALSE]))%*%ifelse(is.na(add.dimada.legendre$adaLasso$post.lm$coefficients),0,add.dimada.legendre$adaLasso$post.lm$coefficients)-as.vector(y.scale.test))^2))

# B-Splines as basis
add.bspline.sieve <- bspline.gen(data=X.scale.train, test.data = X.scale.test, max.interaction = 1, n.basis=10)

set.seed(200)
add.dimada.bspline <- dimada(y=as.vector(y.scale.train), x=NULL, x.sieve = add.bspline.sieve$train,
                         methods = c("Lasso","adaLasso"), family = "gaussian")

(add.bspline.Lasso.mse <- mean((as.matrix(cbind(1,add.bspline.sieve$test[,add.dimada.bspline$Lasso$coefs.final$terms, drop=FALSE]))%*%ifelse(is.na(add.dimada.bspline$Lasso$post.lm$coefficients),0,add.dimada.bspline$Lasso$post.lm$coefficients)-as.vector(y.scale.test))^2))
(add.bspline.adaLasso.mse <- mean((as.matrix(cbind(1,add.bspline.sieve$test[,add.dimada.bspline$adaLasso$coefs.final$terms, drop=FALSE]))%*%ifelse(is.na(add.dimada.bspline$adaLasso$post.lm$coefficients),0,add.dimada.bspline$adaLasso$post.lm$coefficients)-as.vector(y.scale.test))^2))

# trigonometric polynomials as basis
add.trig.sieve <- trig.gen(data=X.scale.train, test.data = X.scale.test, max.interaction = 1,n.basis = 20)

set.seed(200)
add.dimada.trig <- dimada(y=as.vector(y.scale.train), x=NULL, x.sieve = add.trig.sieve$train,
                      methods = c("Lasso","adaLasso"), family = "gaussian")

(add.trig.Lasso.mse <- mean((as.matrix(cbind(1,add.trig.sieve$test[,add.dimada.trig$Lasso$coefs.final$terms, drop=FALSE]))%*%ifelse(is.na(add.dimada.trig$Lasso$post.lm$coefficients),0,add.dimada.trig$Lasso$post.lm$coefficients)-as.vector(y.scale.test))^2))
(add.trig.adaLasso.mse <- mean((as.matrix(cbind(1,add.trig.sieve$test[,add.dimada.trig$adaLasso$coefs.final$terms, drop=FALSE]))%*%ifelse(is.na(add.dimada.trig$adaLasso$post.lm$coefficients),0,add.dimada.trig$adaLasso$post.lm$coefficients)-as.vector(y.scale.test))^2))

# ----------dimension adaptive estimator ---------------------------------------
# power series as basis
power.sieve <- poly.gen(data=X.scale.train, test.data = X.scale.test, max.interaction = 5, legendre=FALSE, n.basis = 20)

set.seed(200)
dimada.power <- dimada(y=as.vector(y.scale.train), x=NULL, x.sieve = power.sieve$train,
                       methods = c("Lasso","adaLasso"), family = "gaussian")

(power.Lasso.mse <- mean((as.matrix(cbind(1,power.sieve$test[,dimada.power$Lasso$coefs.final$terms, drop=FALSE]))%*%ifelse(is.na(dimada.power$Lasso$post.lm$coefficients),0,dimada.power$Lasso$post.lm$coefficients)-as.vector(y.scale.test))^2))
(power.adaLasso.mse <- mean((as.matrix(cbind(1,power.sieve$test[,dimada.power$adaLasso$coefs.final$terms, drop=FALSE]))%*%ifelse(is.na(dimada.power$adaLasso$post.lm$coefficients),0,dimada.power$adaLasso$post.lm$coefficients)-as.vector(y.scale.test))^2))

# legendre polynomials as basis
legendre.sieve <- poly.gen(data=X.scale.train, test.data = X.scale.test, max.interaction = 5, legendre = TRUE,n.basis=20)

set.seed(200)
dimada.legendre <- dimada(y=as.vector(y.scale.train), x=NULL, x.sieve = legendre.sieve$train,
                          methods = c("Lasso","adaLasso"), family = "gaussian")

(legendre.Lasso.mse <- mean((as.matrix(cbind(1,legendre.sieve$test[,dimada.legendre$Lasso$coefs.final$terms, drop=FALSE]))%*%ifelse(is.na(dimada.legendre$Lasso$post.lm$coefficients),0,dimada.legendre$Lasso$post.lm$coefficients)-as.vector(y.scale.test))^2))
(legendre.adaLasso.mse <- mean((as.matrix(cbind(1,legendre.sieve$test[,dimada.legendre$adaLasso$coefs.final$terms, drop=FALSE]))%*%ifelse(is.na(dimada.legendre$adaLasso$post.lm$coefficients),0,dimada.legendre$adaLasso$post.lm$coefficients)-as.vector(y.scale.test))^2))

# B-Splines as basis
bspline.sieve <- bspline.gen(data=X.scale.train, test.data = X.scale.test, max.interaction = 5, n.basis=10)

set.seed(200)
dimada.bspline <- dimada(y=as.vector(y.scale.train), x=NULL, x.sieve = bspline.sieve$train,
                         methods = c("Lasso","adaLasso"), family = "gaussian")

(bspline.Lasso.mse <- mean((as.matrix(cbind(1,bspline.sieve$test[,dimada.bspline$Lasso$coefs.final$terms, drop=FALSE]))%*%ifelse(is.na(dimada.bspline$Lasso$post.lm$coefficients),0,dimada.bspline$Lasso$post.lm$coefficients)-as.vector(y.scale.test))^2))
(bspline.adaLasso.mse <- mean((as.matrix(cbind(1,bspline.sieve$test[,dimada.bspline$adaLasso$coefs.final$terms, drop=FALSE]))%*%ifelse(is.na(dimada.bspline$adaLasso$post.lm$coefficients),0,dimada.bspline$adaLasso$post.lm$coefficients)-as.vector(y.scale.test))^2))

# trigonometric polynomials as basis
trig.sieve <- trig.gen(data=X.scale.train, test.data = X.scale.test, max.interaction = 5,n.basis = 20)

set.seed(200)
dimada.trig <- dimada(y=as.vector(y.scale.train), x=NULL, x.sieve = trig.sieve$train,
                      methods = c("Lasso","adaLasso"), family = "gaussian")

(trig.Lasso.mse <- mean((as.matrix(cbind(1,trig.sieve$test[,dimada.trig$Lasso$coefs.final$terms, drop=FALSE]))%*%ifelse(is.na(dimada.trig$Lasso$post.lm$coefficients),0,dimada.trig$Lasso$post.lm$coefficients)-as.vector(y.scale.test))^2))
(trig.adaLasso.mse <- mean((as.matrix(cbind(1,trig.sieve$test[,dimada.trig$adaLasso$coefs.final$terms, drop=FALSE]))%*%ifelse(is.na(dimada.trig$adaLasso$post.lm$coefficients),0,dimada.trig$adaLasso$post.lm$coefficients)-as.vector(y.scale.test))^2))

save(dimada.power, dimada.legendre, dimada.bspline, dimada.trig, file="Application.RData")

# -----------------results summary ---------------------------------------------

methods<- c("Post LASSO", "Post Adaptive LASSO")
basis <- c("Power Series", "Legendre Polynomials", "B-Splines", "Trigonometric Polynomials")
estimator <- c("addt","dimada")
mse <- c(add.power.Lasso.mse, add.power.adaLasso.mse, add.legendre.Lasso.mse, add.legendre.adaLasso.mse,
         add.bspline.Lasso.mse, add.bspline.adaLasso.mse, add.trig.Lasso.mse, add.trig.adaLasso.mse,
         power.Lasso.mse, power.adaLasso.mse, legendre.Lasso.mse, legendre.adaLasso.mse,
         bspline.Lasso.mse, bspline.adaLasso.mse, trig.Lasso.mse, trig.adaLasso.mse
         )
result <- cbind(expand.grid(method=methods, basis=basis, estimator=estimator), mse)
result.final <- rbind(result, data.frame(method=NA, basis=NA, estimator="ols", mse=ols.mse)) %>% mutate(ratio=mse/ols.mse)
               
tabular.result <- tabular((Factor(estimator, "Estimator")*Factor(basis, "Basis")*Factor(method,"Method")) ~ identity*(mse+ratio), data=cbind(result,ratio=result.final$ratio[-17]))
# toLatex(tabular.result)

options(warn=0)

