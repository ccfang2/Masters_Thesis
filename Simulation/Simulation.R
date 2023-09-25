# ==============================================================================
#      Simulation: Dimension Adaptive Estimation
# ==============================================================================

# Author: Chencheng Fang, ccfang@uni-bonn.de
# Supervisor: Prof. Dr. Joachim Freyberger

# Note: I have designed an R package 'dimada', which helps to implement dimension 
# adaptive estimation easily, and I have made it available on GitHub at ccfang2.

# -----------libraries ---------------------------------------------------------
rm(list=ls())

# install package 'dimada' from my GitHub
# devtools::install_github("ccfang2/dimada")
library(dimada)

library(ggplot2)
library(ggthemes)
library(ggpubr)
library(dplyr)
library(tables)
library(DT)

# ----------source -------------------------------------------------------------
setwd("~/Desktop/bonn master/thesis")
source("Simulation/dimada.sim.R")

# ==============================================================================
# Part 1: Comparing the Performance of Estimators (n=400)
# ==============================================================================

# ++++++++++++++ Model 1: Parametric +++++++++++++++++++++++++++++++++++++++++++

# -------------- basis: power series -------------------------------------------
# the same seed to ensure the same train data set to be used for different basis functions
set.seed(200)
data.param.power.repl.n400 <- replicate(500,dimada.sim(n=400, n.test=1000, d=5,
                                                      formula="3*X[,1]+1.8*X[,2]+X[,3]+2.5*X[,4]+X[,5]",
                                                      basis="power",est=c("unrestricted","additive","ols")),
                                        simplify = FALSE)
                                                            
# -------------- basis: legendre polynomials -----------------------------------
set.seed(200)
data.param.legendre.repl.n400 <- replicate(500,dimada.sim(n=400, n.test=1000, d=5,
                                                         formula="3*X[,1]+1.8*X[,2]+X[,3]+2.5*X[,4]+X[,5]",
                                                         basis="legendre",est=c("unrestricted","additive","ols")),
                                           simplify = FALSE)
                                                                 
# -------------- basis: B-splines-----------------------------------------------
set.seed(200)
data.param.bspline.repl.n400 <- replicate(500,dimada.sim(n=400, n.test=1000, d=5,
                                                        formula="3*X[,1]+1.8*X[,2]+X[,3]+2.5*X[,4]+X[,5]",
                                                        basis="bspline",est=c("unrestricted","additive","ols")),
                                          simplify = FALSE)
                                                                
# -------------- basis: trigonometric polynomials-------------------------------
set.seed(200)
data.param.trig.repl.n400 <- replicate(500,dimada.sim(n=400, n.test=1000, d=5,
                                                     formula="3*X[,1]+1.8*X[,2]+X[,3]+2.5*X[,4]+X[,5]",
                                                     basis="trig",est=c("unrestricted","additive","ols")),
                                       simplify = FALSE)
                                                        
save(data.param.power.repl.n400, data.param.legendre.repl.n400,
     data.param.bspline.repl.n400,data.param.trig.repl.n400,
     file="param.n400.RData")        

# --------------Table: Comparisons (mean) --------------------------------------
table.param.power.test.n400 <-
  Reduce(function(x, y) inner_join(x, y, by = c("Method","Estimator","Sigma")),
         plyr::llply(data.param.power.repl.n400, function(x) x)) %>% 
  rowwise() %>%
  mutate(Mean.NonZeros = mean(c_across(starts_with("NonZeros")), na.rm = TRUE),
         Mean.MSE = mean(c_across(starts_with("MSE")), na.rm = TRUE)) %>%
  select(Mean.NonZeros, Mean.MSE, Method, Estimator, Sigma) %>%
  mutate(Basis= dplyr::if_else(is.na(Method), NA, "Power Series"))

table.param.legendre.test.n400 <-
  Reduce(function(x, y) inner_join(x, y, by = c("Method","Estimator","Sigma")),
         plyr::llply(data.param.legendre.repl.n400, function(x) x)) %>% 
  rowwise() %>%
  mutate(Mean.NonZeros = mean(c_across(starts_with("NonZeros")), na.rm = TRUE),
         Mean.MSE = mean(c_across(starts_with("MSE")), na.rm = TRUE)) %>%
  select(Mean.NonZeros, Mean.MSE, Method, Estimator, Sigma) %>%
  mutate(Basis= dplyr::if_else(is.na(Method), NA, "Legendre Polynomials"))

table.param.bspline.test.n400 <-
  Reduce(function(x, y) inner_join(x, y, by = c("Method","Estimator","Sigma")),
         plyr::llply(data.param.bspline.repl.n400, function(x) x)) %>% 
  rowwise() %>%
  mutate(Mean.NonZeros = mean(c_across(starts_with("NonZeros")), na.rm = TRUE),
         Mean.MSE = mean(c_across(starts_with("MSE")), na.rm = TRUE)) %>%
  select(Mean.NonZeros, Mean.MSE, Method, Estimator, Sigma) %>%
  mutate(Basis=dplyr::if_else(is.na(Method), NA, "B-Splines"))

table.param.trig.test.n400 <-
  Reduce(function(x, y) inner_join(x, y, by = c("Method","Estimator","Sigma")),
         plyr::llply(data.param.trig.repl.n400, function(x) x)) %>% 
  rowwise() %>%
  mutate(Mean.NonZeros = mean(c_across(starts_with("NonZeros")), na.rm = TRUE),
         Mean.MSE = mean(c_across(starts_with("MSE")), na.rm = TRUE)) %>%
  select(Mean.NonZeros, Mean.MSE, Method, Estimator, Sigma) %>%
  mutate(Basis=dplyr::if_else(is.na(Method), NA, "Trigonometric Polynomials"))

table.param.test.n400 <-
  do.call(rbind, list(table.param.power.test.n400,
                      table.param.legendre.test.n400,
                      table.param.bspline.test.n400,
                      table.param.trig.test.n400)) %>%
  distinct() %>%
  mutate(Method=ordered(Method, levels=c("Post LASSO", "Post Adaptive LASSO")),
         Estimator=ordered(Estimator, levels=c("Unrestricted","Additive","OLS")), 
         Sigma=ordered(Sigma, levels=c("0.05","0.2")),
         Basis=factor(Basis, ordered=TRUE, levels=c("Power Series","Legendre Polynomials","B-Splines","Trigonometric Polynomials"))) 
  # %>%
  # datatable(filter="top", caption = "True Model: Parametric") %>%
  # DT::formatRound(columns=1:2, digits=c(2,7))

tabular.param.test.n400 <- tabular((Factor(Basis, "Basis"))*(Factor(Estimator, "Estimator")) ~ (Factor(Sigma, "Sigma"))*(Factor(Method, "Method"))*identity*(MSE + Terms), 
                                   data=table.param.test.n400 %>% rename(MSE=Mean.MSE, Terms=Mean.NonZeros))
# toLatex(tabular.param.test.n400)

# --------------Table: Comparisons (median) --------------------------------------
table.param.power.test.n400.md <-
  Reduce(function(x, y) inner_join(x, y, by = c("Method","Estimator","Sigma")),
         plyr::llply(data.param.power.repl.n400, function(x) x)) %>% 
  rowwise() %>%
  mutate(Median.NonZeros = median(c_across(starts_with("NonZeros")), na.rm = TRUE),
         Median.MSE = median(c_across(starts_with("MSE")), na.rm = TRUE)) %>%
  select(Median.NonZeros, Median.MSE, Method, Estimator, Sigma) %>%
  mutate(Basis= dplyr::if_else(is.na(Method), NA, "Power Series"))

table.param.legendre.test.n400.md <-
  Reduce(function(x, y) inner_join(x, y, by = c("Method","Estimator","Sigma")),
         plyr::llply(data.param.legendre.repl.n400, function(x) x)) %>% 
  rowwise() %>%
  mutate(Median.NonZeros = median(c_across(starts_with("NonZeros")), na.rm = TRUE),
         Median.MSE = median(c_across(starts_with("MSE")), na.rm = TRUE)) %>%
  select(Median.NonZeros, Median.MSE, Method, Estimator, Sigma) %>%
  mutate(Basis= dplyr::if_else(is.na(Method), NA, "Legendre Polynomials"))

table.param.bspline.test.n400.md <-
  Reduce(function(x, y) inner_join(x, y, by = c("Method","Estimator","Sigma")),
         plyr::llply(data.param.bspline.repl.n400, function(x) x)) %>% 
  rowwise() %>%
  mutate(Median.NonZeros = median(c_across(starts_with("NonZeros")), na.rm = TRUE),
         Median.MSE = median(c_across(starts_with("MSE")), na.rm = TRUE)) %>%
  select(Median.NonZeros, Median.MSE, Method, Estimator, Sigma) %>%
  mutate(Basis=dplyr::if_else(is.na(Method), NA, "B-Splines"))

table.param.trig.test.n400.md <-
  Reduce(function(x, y) inner_join(x, y, by = c("Method","Estimator","Sigma")),
         plyr::llply(data.param.trig.repl.n400, function(x) x)) %>% 
  rowwise() %>%
  mutate(Median.NonZeros = median(c_across(starts_with("NonZeros")), na.rm = TRUE),
         Median.MSE = median(c_across(starts_with("MSE")), na.rm = TRUE)) %>%
  select(Median.NonZeros, Median.MSE, Method, Estimator, Sigma) %>%
  mutate(Basis=dplyr::if_else(is.na(Method), NA, "Trigonometric Polynomials"))

table.param.test.n400.md <-
  do.call(rbind, list(table.param.power.test.n400.md,
                      table.param.legendre.test.n400.md,
                      table.param.bspline.test.n400.md,
                      table.param.trig.test.n400.md)) %>%
  distinct() %>%
  mutate(Method=ordered(Method, levels=c("Post LASSO", "Post Adaptive LASSO")),
         Estimator=ordered(Estimator, levels=c("Unrestricted","Additive","OLS")), 
         Sigma=ordered(Sigma, levels=c("0.05","0.2")),
         Basis=factor(Basis, ordered=TRUE, levels=c("Power Series","Legendre Polynomials","B-Splines","Trigonometric Polynomials"))) 
  # %>%
  # datatable(filter="top", caption = "True Model: Parametric") %>%
  # DT::formatRound(columns=1:2, digits=c(2,7))

tabular.param.test.n400.md <- tabular((Factor(Basis, "Basis"))*(Factor(Estimator, "Estimator")) ~ (Factor(Sigma, "Sigma"))*(Factor(Method, "Method"))*identity*(MSE + Terms), 
                                   data=table.param.test.n400.md %>% rename(MSE=Median.MSE, Terms=Median.NonZeros))
# toLatex(tabular.param.test.n400.md)

# +++++++++++++++ Model 2: Non-Parametric with Additivity ++++++++++++++++++++++

# --------------- basis: power series ------------------------------------------
set.seed(200)
data.add.power.repl.n400 <- replicate(500,dimada.sim(n=400, n.test=1000, d=5,
                                                    formula="sin(4*X[,1])+1.5*log(X[,2])+1/cos(X[,3])+sin(sqrt(X[,4]))+sin(X[,5]^2)",
                                                    basis="power",est=c("unrestricted","additive","ols")),
                                      simplify = FALSE)
                                                      
# --------------- basis: legendre polynomials-----------------------------------
set.seed(200)
data.add.legendre.repl.n400 <- replicate(500,dimada.sim(n=400, n.test=1000, d=5,
                                                       formula="sin(4*X[,1])+1.5*log(X[,2])+1/cos(X[,3])+sin(sqrt(X[,4]))+sin(X[,5]^2)",
                                                       basis="legendre",est=c("unrestricted","additive","ols")),
                                         simplify = FALSE)
                                                               
# --------------- basis: B-splines----------------------------------------------
set.seed(200)
data.add.bspline.repl.n400 <- replicate(500,dimada.sim(n=400, n.test=1000, d=5,
                                                      formula="sin(4*X[,1])+1.5*log(X[,2])+1/cos(X[,3])+sin(sqrt(X[,4]))+sin(X[,5]^2)",
                                                      basis="bspline",est=c("unrestricted","additive","ols")),
                                        simplify = FALSE)
                                                              
# --------------- basis: trigonometric polynomials------------------------------
set.seed(200)
data.add.trig.repl.n400 <- replicate(500,dimada.sim(n=400, n.test=1000, d=5,
                                                   formula="sin(4*X[,1])+1.5*log(X[,2])+1/cos(X[,3])+sin(sqrt(X[,4]))+sin(X[,5]^2)",
                                                   basis="trig",est=c("unrestricted","additive","ols")),
                                     simplify = FALSE)
                                
save(data.add.power.repl.n400, data.add.legendre.repl.n400,
     data.add.bspline.repl.n400,data.add.trig.repl.n400,
     file="add.n400.RData")      

# --------------- Table: Comparisons (mean)-------------------------------------
table.add.power.test.n400 <-
  Reduce(function(x, y) inner_join(x, y, by = c("Method","Estimator","Sigma")),
         plyr::llply(data.add.power.repl.n400, function(x) x)) %>%
  rowwise() %>%
  mutate(Mean.NonZeros = mean(c_across(starts_with("NonZeros")), na.rm = TRUE),
         Mean.MSE = mean(c_across(starts_with("MSE")), na.rm = TRUE)) %>%
  select(Mean.NonZeros, Mean.MSE, Method, Estimator, Sigma) %>%
  mutate(Basis= dplyr::if_else(is.na(Method), NA, "Power Series"))

table.add.legendre.test.n400 <-
  Reduce(function(x, y) inner_join(x, y, by = c("Method","Estimator","Sigma")),
         plyr::llply(data.add.legendre.repl.n400, function(x) x)) %>%
  rowwise() %>%
  mutate(Mean.NonZeros = mean(c_across(starts_with("NonZeros")), na.rm = TRUE),
         Mean.MSE = mean(c_across(starts_with("MSE")), na.rm = TRUE)) %>%
  select(Mean.NonZeros, Mean.MSE, Method, Estimator, Sigma) %>%
  mutate(Basis= dplyr::if_else(is.na(Method), NA, "Legendre Polynomials"))

table.add.bspline.test.n400 <-
  Reduce(function(x, y) inner_join(x, y, by = c("Method","Estimator","Sigma")),
         plyr::llply(data.add.bspline.repl.n400, function(x) x)) %>% 
  rowwise() %>%
  mutate(Mean.NonZeros = mean(c_across(starts_with("NonZeros")), na.rm = TRUE),
         Mean.MSE = mean(c_across(starts_with("MSE")), na.rm = TRUE)) %>%
  select(Mean.NonZeros, Mean.MSE, Method, Estimator, Sigma) %>%
  mutate(Basis= dplyr::if_else(is.na(Method), NA, "B-Splines"))

table.add.trig.test.n400 <-
  Reduce(function(x, y) inner_join(x, y, by = c("Method","Estimator","Sigma")),
         plyr::llply(data.add.trig.repl.n400, function(x) x)) %>% 
  rowwise() %>%
  mutate(Mean.NonZeros = mean(c_across(starts_with("NonZeros")), na.rm = TRUE),
         Mean.MSE = mean(c_across(starts_with("MSE")), na.rm = TRUE)) %>%
  select(Mean.NonZeros, Mean.MSE, Method, Estimator, Sigma) %>%
  mutate(Basis= dplyr::if_else(is.na(Method), NA, "Trigonometric Polynomials"))

table.add.test.n400 <-
  do.call(rbind, list(table.add.power.test.n400,
                      table.add.legendre.test.n400,
                      table.add.bspline.test.n400,
                      table.add.trig.test.n400)) %>%
  distinct() %>%
  mutate(Method=ordered(Method, levels=c("Post LASSO", "Post Adaptive LASSO")),
         Estimator=ordered(Estimator, levels=c("Unrestricted","Additive","OLS")), 
         Sigma=ordered(Sigma, levels=c("0.05","0.2")),
         Basis=factor(Basis, ordered=TRUE, levels=c("Power Series","Legendre Polynomials","B-Splines","Trigonometric Polynomials")))
  # %>%
  # datatable(filter="top", caption = "True Model: Non-parametric with Additivity") %>%
  # DT::formatRound(columns=1:2, digits=c(2,7))

tabular.add.test.n400 <- tabular((Factor(Basis, "Basis"))*(Factor(Estimator, "Estimator")) ~ (Factor(Sigma, "Sigma"))*(Factor(Method, "Method"))*identity*(MSE + Terms), 
                                   data=table.add.test.n400 %>% rename(MSE=Mean.MSE, Terms=Mean.NonZeros))
# toLatex(tabular.add.test.n400)

# --------------- Table: Comparisons (median)-------------------------------------
table.add.power.test.n400.md <-
  Reduce(function(x, y) inner_join(x, y, by = c("Method","Estimator","Sigma")),
         plyr::llply(data.add.power.repl.n400, function(x) x)) %>%
  rowwise() %>%
  mutate(Median.NonZeros = median(c_across(starts_with("NonZeros")), na.rm = TRUE),
         Median.MSE = median(c_across(starts_with("MSE")), na.rm = TRUE)) %>%
  select(Median.NonZeros, Median.MSE, Method, Estimator, Sigma) %>%
  mutate(Basis= dplyr::if_else(is.na(Method), NA, "Power Series"))

table.add.legendre.test.n400.md <-
  Reduce(function(x, y) inner_join(x, y, by = c("Method","Estimator","Sigma")),
         plyr::llply(data.add.legendre.repl.n400, function(x) x)) %>%
  rowwise() %>%
  mutate(Median.NonZeros = median(c_across(starts_with("NonZeros")), na.rm = TRUE),
         Median.MSE = median(c_across(starts_with("MSE")), na.rm = TRUE)) %>%
  select(Median.NonZeros, Median.MSE, Method, Estimator, Sigma) %>%
  mutate(Basis= dplyr::if_else(is.na(Method), NA, "Legendre Polynomials"))

table.add.bspline.test.n400.md <-
  Reduce(function(x, y) inner_join(x, y, by = c("Method","Estimator","Sigma")),
         plyr::llply(data.add.bspline.repl.n400, function(x) x)) %>% 
  rowwise() %>%
  mutate(Median.NonZeros = median(c_across(starts_with("NonZeros")), na.rm = TRUE),
         Median.MSE = median(c_across(starts_with("MSE")), na.rm = TRUE)) %>%
  select(Median.NonZeros, Median.MSE, Method, Estimator, Sigma) %>%
  mutate(Basis= dplyr::if_else(is.na(Method), NA, "B-Splines"))

table.add.trig.test.n400.md <-
  Reduce(function(x, y) inner_join(x, y, by = c("Method","Estimator","Sigma")),
         plyr::llply(data.add.trig.repl.n400, function(x) x)) %>% 
  rowwise() %>%
  mutate(Median.NonZeros = median(c_across(starts_with("NonZeros")), na.rm = TRUE),
         Median.MSE = median(c_across(starts_with("MSE")), na.rm = TRUE)) %>%
  select(Median.NonZeros, Median.MSE, Method, Estimator, Sigma) %>%
  mutate(Basis= dplyr::if_else(is.na(Method), NA, "Trigonometric Polynomials"))

table.add.test.n400.md <-
  do.call(rbind, list(table.add.power.test.n400.md,
                      table.add.legendre.test.n400.md,
                      table.add.bspline.test.n400.md,
                      table.add.trig.test.n400.md)) %>%
  distinct() %>%
  mutate(Method=ordered(Method, levels=c("Post LASSO", "Post Adaptive LASSO")),
         Estimator=ordered(Estimator, levels=c("Unrestricted","Additive","OLS")), 
         Sigma=ordered(Sigma, levels=c("0.05","0.2")),
         Basis=factor(Basis, ordered=TRUE, levels=c("Power Series","Legendre Polynomials","B-Splines","Trigonometric Polynomials")))
  # %>%
  # datatable(filter="top", caption = "True Model: Non-parametric with Additivity") %>%
  # DT::formatRound(columns=1:2, digits=c(2,7))

tabular.add.test.n400.md <- tabular((Factor(Basis, "Basis"))*(Factor(Estimator, "Estimator")) ~ (Factor(Sigma, "Sigma"))*(Factor(Method, "Method"))*identity*(MSE + Terms), 
                                 data=table.add.test.n400.md %>% rename(MSE=Median.MSE, Terms=Median.NonZeros))
# toLatex(tabular.add.test.n400.md)


# +++++++++++++++ Model 3: Non-Parametric ++++++++++++++++++++++++++++++++++++++

# --------------- basis: power series ------------------------------------------
set.seed(200)
data.np.power.repl.n400 <- replicate(500,dimada.sim(n=400, n.test=1000, d=5,
                                                   formula="3*(X[,1]+4*X[,2]+X[,3]*X[,4]*X[,5])^(1/4)+2*sin(X[,4]+X[,5]^2+X[,1]*X[,2]*X[,3])+3*log(X[,3]^2+X[,4]+2*X[,5])",
                                                   basis="power",est=c("unrestricted","additive","ols")),
                                     simplify = FALSE)
                                                      
# --------------- basis: legendre polynomials ----------------------------------
set.seed(200)
data.np.legendre.repl.n400 <- replicate(500,dimada.sim(n=400, n.test=1000, d=5,
                                                      formula="3*(X[,1]+4*X[,2]+X[,3]*X[,4]*X[,5])^(1/4)+2*sin(X[,4]+X[,5]^2+X[,1]*X[,2]*X[,3])+3*log(X[,3]^2+X[,4]+2*X[,5])", 
                                                      basis="legendre",est=c("unrestricted","additive","ols")),
                                        simplify = FALSE)
                                                              
# --------------- basis: B-splines ---------------------------------------------
set.seed(200)
data.np.bspline.repl.n400 <- replicate(500,dimada.sim(n=400, n.test=1000, d=5,
                                                     formula="3*(X[,1]+4*X[,2]+X[,3]*X[,4]*X[,5])^(1/4)+2*sin(X[,4]+X[,5]^2+X[,1]*X[,2]*X[,3])+3*log(X[,3]^2+X[,4]+2*X[,5])",
                                                     basis="bspline",est=c("unrestricted","additive","ols")),
                                       simplify = FALSE)
                                                             
# --------------- basis: trigonometric polynomials -----------------------------
set.seed(200)
data.np.trig.repl.n400 <- replicate(500,dimada.sim(n=400, n.test=1000, d=5,
                                                  formula="3*(X[,1]+4*X[,2]+X[,3]*X[,4]*X[,5])^(1/4)+2*sin(X[,4]+X[,5]^2+X[,1]*X[,2]*X[,3])+3*log(X[,3]^2+X[,4]+2*X[,5])",
                                                  basis="trig",est=c("unrestricted","additive","ols")),
                                    simplify = FALSE)
                                                   
save(data.np.power.repl.n400, data.np.legendre.repl.n400,
     data.np.bspline.repl.n400,data.np.trig.repl.n400,
     file="np.n400.RData")   

# --------------- Table: Comparisons (mean) ------------------------------------
table.np.power.test.n400 <-
  Reduce(function(x, y) inner_join(x, y, by = c("Method","Estimator","Sigma")),
         plyr::llply(data.np.power.repl.n400, function(x) x)) %>% #x$result.test
  rowwise() %>%
  mutate(Mean.NonZeros = mean(c_across(starts_with("NonZeros")), na.rm = TRUE),
         Mean.MSE = mean(c_across(starts_with("MSE")), na.rm = TRUE)) %>%
  select(Mean.NonZeros, Mean.MSE, Method, Estimator, Sigma) %>%
  mutate(Basis= dplyr::if_else(is.na(Method), NA, "Power Series"))

table.np.legendre.test.n400 <-
  Reduce(function(x, y) inner_join(x, y, by = c("Method","Estimator","Sigma")),
         plyr::llply(data.np.legendre.repl.n400, function(x) x)) %>% #x$result.test
  rowwise() %>%
  mutate(Mean.NonZeros = mean(c_across(starts_with("NonZeros")), na.rm = TRUE),
         Mean.MSE = mean(c_across(starts_with("MSE")), na.rm = TRUE)) %>%
  select(Mean.NonZeros, Mean.MSE, Method, Estimator, Sigma) %>%
  mutate(Basis= dplyr::if_else(is.na(Method), NA, "Legendre Polynomials"))

table.np.bspline.test.n400 <-
  Reduce(function(x, y) inner_join(x, y, by = c("Method","Estimator","Sigma")),
         plyr::llply(data.np.bspline.repl.n400, function(x) x)) %>% #x$result.test
  rowwise() %>%
  mutate(Mean.NonZeros = mean(c_across(starts_with("NonZeros")), na.rm = TRUE),
         Mean.MSE = mean(c_across(starts_with("MSE")), na.rm = TRUE)) %>%
  select(Mean.NonZeros, Mean.MSE, Method, Estimator, Sigma) %>%
  mutate(Basis= dplyr::if_else(is.na(Method), NA, "B-Splines"))

table.np.trig.test.n400 <-
  Reduce(function(x, y) inner_join(x, y, by = c("Method","Estimator","Sigma")),
         plyr::llply(data.np.trig.repl.n400, function(x) x)) %>% 
  rowwise() %>%
  mutate(Mean.NonZeros = mean(c_across(starts_with("NonZeros")), na.rm = TRUE),
         Mean.MSE = mean(c_across(starts_with("MSE")), na.rm = TRUE)) %>%
  select(Mean.NonZeros, Mean.MSE, Method, Estimator, Sigma) %>%
  mutate(Basis= dplyr::if_else(is.na(Method), NA, "Trigonometric Polynomials"))

table.np.test.n400 <-
  do.call(rbind, list(table.np.power.test.n400,
                      table.np.legendre.test.n400,
                      table.np.bspline.test.n400,
                      table.np.trig.test.n400)) %>%
  distinct() %>%
  mutate(Method=ordered(Method, levels=c("Post LASSO", "Post Adaptive LASSO")),
         Estimator=ordered(Estimator, levels=c("Unrestricted","Additive","OLS")), 
         Sigma=ordered(Sigma, levels=c("0.05","0.2")),
         Basis=factor(Basis, ordered=TRUE, levels=c("Power Series","Legendre Polynomials","B-Splines","Trigonometric Polynomials")))
  # %>%
  #datatable(filter="top", caption = "True Model: Non-Parametric") %>%
  #DT::formatRound(columns=1:2, digits=c(2,7))

tabular.np.test.n400 <- tabular((Factor(Basis, "Basis"))*(Factor(Estimator, "Estimator")) ~ (Factor(Sigma, "Sigma"))*(Factor(Method, "Method"))*identity*(MSE + Terms), 
                                 data=table.np.test.n400 %>% rename(MSE=Mean.MSE, Terms=Mean.NonZeros))
# toLatex(tabular.np.test.n400)

# --------------- Table: Comparisons (median) ------------------------------------
table.np.power.test.n400.md <-
  Reduce(function(x, y) inner_join(x, y, by = c("Method","Estimator","Sigma")),
         plyr::llply(data.np.power.repl.n400, function(x) x)) %>% #x$result.test
  rowwise() %>%
  mutate(Median.NonZeros = median(c_across(starts_with("NonZeros")), na.rm = TRUE),
         Median.MSE = median(c_across(starts_with("MSE")), na.rm = TRUE)) %>%
  select(Median.NonZeros, Median.MSE, Method, Estimator, Sigma) %>%
  mutate(Basis= dplyr::if_else(is.na(Method), NA, "Power Series"))

table.np.legendre.test.n400.md <-
  Reduce(function(x, y) inner_join(x, y, by = c("Method","Estimator","Sigma")),
         plyr::llply(data.np.legendre.repl.n400, function(x) x)) %>% #x$result.test
  rowwise() %>%
  mutate(Median.NonZeros = median(c_across(starts_with("NonZeros")), na.rm = TRUE),
         Median.MSE = median(c_across(starts_with("MSE")), na.rm = TRUE)) %>%
  select(Median.NonZeros, Median.MSE, Method, Estimator, Sigma) %>%
  mutate(Basis= dplyr::if_else(is.na(Method), NA, "Legendre Polynomials"))

table.np.bspline.test.n400.md <-
  Reduce(function(x, y) inner_join(x, y, by = c("Method","Estimator","Sigma")),
         plyr::llply(data.np.bspline.repl.n400, function(x) x)) %>% #x$result.test
  rowwise() %>%
  mutate(Median.NonZeros = median(c_across(starts_with("NonZeros")), na.rm = TRUE),
         Median.MSE = median(c_across(starts_with("MSE")), na.rm = TRUE)) %>%
  select(Median.NonZeros, Median.MSE, Method, Estimator, Sigma) %>%
  mutate(Basis= dplyr::if_else(is.na(Method), NA, "B-Splines"))

table.np.trig.test.n400.md <-
  Reduce(function(x, y) inner_join(x, y, by = c("Method","Estimator","Sigma")),
         plyr::llply(data.np.trig.repl.n400, function(x) x)) %>% 
  rowwise() %>%
  mutate(Median.NonZeros = median(c_across(starts_with("NonZeros")), na.rm = TRUE),
         Median.MSE = median(c_across(starts_with("MSE")), na.rm = TRUE)) %>%
  select(Median.NonZeros, Median.MSE, Method, Estimator, Sigma) %>%
  mutate(Basis= dplyr::if_else(is.na(Method), NA, "Trigonometric Polynomials"))

table.np.test.n400.md <-
  do.call(rbind, list(table.np.power.test.n400.md,
                      table.np.legendre.test.n400.md,
                      table.np.bspline.test.n400.md,
                      table.np.trig.test.n400.md)) %>%
  distinct() %>%
  mutate(Method=ordered(Method, levels=c("Post LASSO", "Post Adaptive LASSO")),
         Estimator=ordered(Estimator, levels=c("Unrestricted","Additive","OLS")), 
         Sigma=ordered(Sigma, levels=c("0.05","0.2")),
         Basis=factor(Basis, ordered=TRUE, levels=c("Power Series","Legendre Polynomials","B-Splines","Trigonometric Polynomials"))) 
  # %>%
  # datatable(filter="top", caption = "True Model: Non-Parametric") %>%
  # DT::formatRound(columns=1:2, digits=c(2,7))

tabular.np.test.n400.md <- tabular((Factor(Basis, "Basis"))*(Factor(Estimator, "Estimator")) ~ (Factor(Sigma, "Sigma"))*(Factor(Method, "Method"))*identity*(MSE + Terms), 
                                data=table.np.test.n400.md %>% rename(MSE=Median.MSE, Terms=Median.NonZeros))
# toLatex(tabular.np.test.n400.md)

# ==============================================================================
# Part 2: Convergence Rates of Dimension Adaptive Estimator in Different Underlying Models
# n=50, 200, 400, 600, 800, 1000; 
# Legendre Polynomials as Basis
# ==============================================================================

# +++++++++++++++ Model 1: parametric ++++++++++++++++++++++++++++++++++++++++++

# --------------- basis: Legendre polynomials ----------------------------------
set.seed(200)
data.param.legendre.repl.n50.rate <- replicate(500, dimada.sim(n=50, n.test=1000, d=5,
                                                              formula="3*X[,1]+1.8*X[,2]+X[,3]+2.5*X[,4]+X[,5]",
                                                              empirical.mse = TRUE,
                                                              basis="legendre",est="unrestricted"),
                                               simplify = FALSE)

set.seed(200)
data.param.legendre.repl.n100.rate <- replicate(500, dimada.sim(n=100, n.test=1000, d=5,
                                                              formula="3*X[,1]+1.8*X[,2]+X[,3]+2.5*X[,4]+X[,5]",
                                                              empirical.mse = TRUE,
                                                              basis="legendre",est="unrestricted"),
                                               simplify = FALSE)

set.seed(200)
data.param.legendre.repl.n150.rate <- replicate(500, dimada.sim(n=150, n.test=1000, d=5,
                                                               formula="3*X[,1]+1.8*X[,2]+X[,3]+2.5*X[,4]+X[,5]",
                                                               empirical.mse = TRUE,
                                                               basis="legendre",est="unrestricted"),
                                                simplify = FALSE)

set.seed(200)
data.param.legendre.repl.n200.rate <- replicate(500,dimada.sim(n=200, n.test=1000, d=5,
                                                              formula="3*X[,1]+1.8*X[,2]+X[,3]+2.5*X[,4]+X[,5]",
                                                              empirical.mse = TRUE,
                                                              basis="legendre",est="unrestricted"),
                                                simplify = FALSE)

set.seed(200)
data.param.legendre.repl.n250.rate <- replicate(500, dimada.sim(n=250, n.test=1000, d=5,
                                                               formula="3*X[,1]+1.8*X[,2]+X[,3]+2.5*X[,4]+X[,5]",
                                                               empirical.mse = TRUE,
                                                               basis="legendre",est="unrestricted"),
                                                simplify = FALSE)

set.seed(200)
data.param.legendre.repl.n300.rate <- replicate(500, dimada.sim(n=300, n.test=1000, d=5,
                                                              formula="3*X[,1]+1.8*X[,2]+X[,3]+2.5*X[,4]+X[,5]",
                                                              empirical.mse = TRUE,
                                                              basis="legendre",est="unrestricted"),
                                               simplify = FALSE)

set.seed(200)
data.param.legendre.repl.n350.rate <- replicate(500, dimada.sim(n=350, n.test=1000, d=5,
                                                               formula="3*X[,1]+1.8*X[,2]+X[,3]+2.5*X[,4]+X[,5]",
                                                               empirical.mse = TRUE,
                                                               basis="legendre",est="unrestricted"),
                                                simplify = FALSE)
                                                         
set.seed(200)
data.param.legendre.repl.n400.rate <- replicate(500,dimada.sim(n=400, n.test=1000, d=5,
                                                              formula="3*X[,1]+1.8*X[,2]+X[,3]+2.5*X[,4]+X[,5]",
                                                              empirical.mse = TRUE,
                                                              basis="legendre",est="unrestricted"),
                                                simplify = FALSE)

set.seed(200)
data.param.legendre.repl.n450.rate <- replicate(500, dimada.sim(n=450, n.test=1000, d=5,
                                                                formula="3*X[,1]+1.8*X[,2]+X[,3]+2.5*X[,4]+X[,5]",
                                                                empirical.mse = TRUE,
                                                                basis="legendre",est="unrestricted"),
                                                simplify = FALSE)


save(data.param.legendre.repl.n50.rate, data.param.legendre.repl.n100.rate, data.param.legendre.repl.n150.rate, 
     data.param.legendre.repl.n200.rate, data.param.legendre.repl.n250.rate, data.param.legendre.repl.n300.rate,
     data.param.legendre.repl.n350.rate, data.param.legendre.repl.n400.rate, data.param.legendre.repl.n450.rate, 
     file="param.rate.RData")     


# +++++++++++++++ Model 2: non parametric with additivity+++++++++++++++++++++++

# --------------- basis: Legendre polynomials ----------------------------------
set.seed(200)
data.add.legendre.repl.n50.rate <- replicate(500,dimada.sim(n=50, n.test=1000, d=5,
                                                            formula="sin(4*X[,1])+1.5*log(X[,2])+1/cos(X[,3])+sin(sqrt(X[,4]))+sin(X[,5]^2)",
                                                            empirical.mse = TRUE, 
                                                            basis="legendre", est="unrestricted"),
                                              simplify = FALSE)

set.seed(200)
data.add.legendre.repl.n100.rate <- replicate(500,dimada.sim(n=100, n.test=1000, d=5,
                                                           formula="sin(4*X[,1])+1.5*log(X[,2])+1/cos(X[,3])+sin(sqrt(X[,4]))+sin(X[,5]^2)",
                                                           empirical.mse = TRUE, 
                                                           basis="legendre", est="unrestricted"),
                                             simplify = FALSE)

set.seed(200)
data.add.legendre.repl.n150.rate <- replicate(500,dimada.sim(n=150, n.test=1000, d=5,
                                                            formula="sin(4*X[,1])+1.5*log(X[,2])+1/cos(X[,3])+sin(sqrt(X[,4]))+sin(X[,5]^2)",
                                                            empirical.mse = TRUE, 
                                                            basis="legendre", est="unrestricted"),
                                              simplify = FALSE)

set.seed(200)
data.add.legendre.repl.n200.rate <- replicate(500,dimada.sim(n=200, n.test=1000, d=5,
                                                            formula="sin(4*X[,1])+1.5*log(X[,2])+1/cos(X[,3])+sin(sqrt(X[,4]))+sin(X[,5]^2)",
                                                            empirical.mse = TRUE, 
                                                            basis="legendre", est="unrestricted"),
                                              simplify = FALSE)

set.seed(200)
data.add.legendre.repl.n250.rate <- replicate(500,dimada.sim(n=250, n.test=1000, d=5,
                                                            formula="sin(4*X[,1])+1.5*log(X[,2])+1/cos(X[,3])+sin(sqrt(X[,4]))+sin(X[,5]^2)",
                                                            empirical.mse = TRUE, 
                                                            basis="legendre", est="unrestricted"),
                                              simplify = FALSE)

set.seed(200)
data.add.legendre.repl.n300.rate <- replicate(500,dimada.sim(n=300, n.test=1000, d=5,
                                                            formula="sin(4*X[,1])+1.5*log(X[,2])+1/cos(X[,3])+sin(sqrt(X[,4]))+sin(X[,5]^2)",
                                                            empirical.mse = TRUE, 
                                                            basis="legendre", est="unrestricted"),
                                              simplify = FALSE)

set.seed(200)
data.add.legendre.repl.n350.rate <- replicate(500,dimada.sim(n=350, n.test=1000, d=5,
                                                            formula="sin(4*X[,1])+1.5*log(X[,2])+1/cos(X[,3])+sin(sqrt(X[,4]))+sin(X[,5]^2)",
                                                            empirical.mse = TRUE, 
                                                            basis="legendre", est="unrestricted"),
                                              simplify = FALSE)

set.seed(200)
data.add.legendre.repl.n400.rate <- replicate(500,dimada.sim(n=400, n.test=1000, d=5,
                                                            formula="sin(4*X[,1])+1.5*log(X[,2])+1/cos(X[,3])+sin(sqrt(X[,4]))+sin(X[,5]^2)",
                                                            empirical.mse = TRUE, 
                                                            basis="legendre", est="unrestricted"),
                                              simplify = FALSE)


set.seed(200)
data.add.legendre.repl.n450.rate <- replicate(500,dimada.sim(n=450, n.test=1000, d=5,
                                                             formula="sin(4*X[,1])+1.5*log(X[,2])+1/cos(X[,3])+sin(sqrt(X[,4]))+sin(X[,5]^2)",
                                                             empirical.mse = TRUE, 
                                                             basis="legendre", est="unrestricted"),
                                              simplify = FALSE)

save(data.add.legendre.repl.n50.rate, data.add.legendre.repl.n100.rate, data.add.legendre.repl.n150.rate, 
     data.add.legendre.repl.n200.rate, data.add.legendre.repl.n250.rate, data.add.legendre.repl.n300.rate,
     data.add.legendre.repl.n350.rate, data.add.legendre.repl.n400.rate, data.add.legendre.repl.n450.rate, 
     file="add.rate.RData") 


# +++++++++++++++ Model 3: non parametric ++++++++++++++++++++++++++++++++++++++

# --------------- basis: Legendre polynomials ----------------------------------
set.seed(200)
data.np.legendre.repl.n50.rate <- replicate(500,dimada.sim(n=50, n.test=1000, d=5,
                                                       formula="3*(X[,1]+4*X[,2]+X[,3]*X[,4]*X[,5])^(1/4)+2*sin(X[,4]+X[,5]^2+X[,1]*X[,2]*X[,3])+3*log(X[,3]^2+X[,4]+2*X[,5])",
                                                       empirical.mse = TRUE,
                                                       basis="legendre", est="unrestricted"),
                                         simplify = FALSE)

set.seed(200)
data.np.legendre.repl.n100.rate <- replicate(500,dimada.sim(n=100, n.test=1000, d=5,
                                                          formula="3*(X[,1]+4*X[,2]+X[,3]*X[,4]*X[,5])^(1/4)+2*sin(X[,4]+X[,5]^2+X[,1]*X[,2]*X[,3])+3*log(X[,3]^2+X[,4]+2*X[,5])",
                                                          empirical.mse = TRUE,
                                                          basis="legendre", est="unrestricted"),
                                            simplify = FALSE)

set.seed(200)
data.np.legendre.repl.n150.rate <- replicate(500,dimada.sim(n=150, n.test=1000, d=5,
                                                          formula="3*(X[,1]+4*X[,2]+X[,3]*X[,4]*X[,5])^(1/4)+2*sin(X[,4]+X[,5]^2+X[,1]*X[,2]*X[,3])+3*log(X[,3]^2+X[,4]+2*X[,5])",
                                                          empirical.mse = TRUE,
                                                          basis="legendre", est="unrestricted"),
                                            simplify = FALSE)

set.seed(200)
data.np.legendre.repl.n200.rate <- replicate(500,dimada.sim(n=200, n.test=1000, d=5,
                                                           formula="3*(X[,1]+4*X[,2]+X[,3]*X[,4]*X[,5])^(1/4)+2*sin(X[,4]+X[,5]^2+X[,1]*X[,2]*X[,3])+3*log(X[,3]^2+X[,4]+2*X[,5])",
                                                           empirical.mse = TRUE,
                                                           basis="legendre", est="unrestricted"),
                                             simplify = FALSE)

set.seed(200)
data.np.legendre.repl.n250.rate <- replicate(500,dimada.sim(n=250, n.test=1000, d=5,
                                                          formula="3*(X[,1]+4*X[,2]+X[,3]*X[,4]*X[,5])^(1/4)+2*sin(X[,4]+X[,5]^2+X[,1]*X[,2]*X[,3])+3*log(X[,3]^2+X[,4]+2*X[,5])",
                                                          empirical.mse = TRUE,
                                                          basis="legendre", est="unrestricted"),
                                            simplify = FALSE)

set.seed(200)
data.np.legendre.repl.n300.rate <- replicate(500,dimada.sim(n=300, n.test=1000, d=5,
                                                          formula="3*(X[,1]+4*X[,2]+X[,3]*X[,4]*X[,5])^(1/4)+2*sin(X[,4]+X[,5]^2+X[,1]*X[,2]*X[,3])+3*log(X[,3]^2+X[,4]+2*X[,5])",
                                                          empirical.mse = TRUE,
                                                          basis="legendre", est="unrestricted"),
                                            simplify = FALSE)

set.seed(200)
data.np.legendre.repl.n350.rate <- replicate(500,dimada.sim(n=350, n.test=1000, d=5,
                                                          formula="3*(X[,1]+4*X[,2]+X[,3]*X[,4]*X[,5])^(1/4)+2*sin(X[,4]+X[,5]^2+X[,1]*X[,2]*X[,3])+3*log(X[,3]^2+X[,4]+2*X[,5])",
                                                          empirical.mse = TRUE,
                                                          basis="legendre", est="unrestricted"),
                                            simplify = FALSE)
                                                       
set.seed(200)
data.np.legendre.repl.n400.rate <- replicate(500,dimada.sim(n=400, n.test=1000, d=5,
                                                           formula="3*(X[,1]+4*X[,2]+X[,3]*X[,4]*X[,5])^(1/4)+2*sin(X[,4]+X[,5]^2+X[,1]*X[,2]*X[,3])+3*log(X[,3]^2+X[,4]+2*X[,5])",
                                                           empirical.mse = TRUE,
                                                           basis="legendre", est="unrestricted"),
                                             simplify = FALSE)



set.seed(200)
data.np.legendre.repl.n450.rate <- replicate(500,dimada.sim(n=450, n.test=1000, d=5,
                                                            formula="3*(X[,1]+4*X[,2]+X[,3]*X[,4]*X[,5])^(1/4)+2*sin(X[,4]+X[,5]^2+X[,1]*X[,2]*X[,3])+3*log(X[,3]^2+X[,4]+2*X[,5])",
                                                            empirical.mse = TRUE,
                                                            basis="legendre", est="unrestricted"),
                                             simplify = FALSE)
                                                   
save(data.np.legendre.repl.n50.rate, data.np.legendre.repl.n100.rate, data.np.legendre.repl.n150.rate,
     data.np.legendre.repl.n200.rate, data.np.legendre.repl.n250.rate, data.np.legendre.repl.n300.rate,
     data.np.legendre.repl.n350.rate, data.np.legendre.repl.n400.rate, data.np.legendre.repl.n450.rate, 
     file="np.rate.RData")     

# +++++++++++++++++++++ Comparison Plots +++++++++++++++++++++++++++++++++++++++

sigma <- c(0.05,0.2)
methods <- c("Post LASSO", "Post Adaptive LASSO")
models <- c("m1","m2", "m3")

param.rate.all <- rbind(do.call(rbind, data.param.legendre.repl.n50.rate) %>% mutate(N=50),
                        do.call(rbind, data.param.legendre.repl.n100.rate) %>% mutate(N=100),
                        do.call(rbind, data.param.legendre.repl.n150.rate) %>% mutate(N=150),
                        do.call(rbind, data.param.legendre.repl.n200.rate) %>% mutate(N=200),
                        do.call(rbind, data.param.legendre.repl.n250.rate) %>% mutate(N=250),
                        do.call(rbind, data.param.legendre.repl.n300.rate) %>% mutate(N=300),
                        do.call(rbind, data.param.legendre.repl.n350.rate) %>% mutate(N=350),
                        do.call(rbind, data.param.legendre.repl.n400.rate) %>% mutate(N=400),
                        do.call(rbind, data.param.legendre.repl.n450.rate) %>% mutate(N=450)) %>% mutate(Model="m1")

add.rate.all <- rbind(do.call(rbind, data.add.legendre.repl.n50.rate) %>% mutate(N=50),
                      do.call(rbind, data.add.legendre.repl.n100.rate) %>% mutate(N=100),
                      do.call(rbind, data.add.legendre.repl.n150.rate) %>% mutate(N=150),
                      do.call(rbind, data.add.legendre.repl.n200.rate) %>% mutate(N=200),
                      do.call(rbind, data.add.legendre.repl.n250.rate) %>% mutate(N=250),
                      do.call(rbind, data.add.legendre.repl.n300.rate) %>% mutate(N=300),
                      do.call(rbind, data.add.legendre.repl.n350.rate) %>% mutate(N=350),
                      do.call(rbind, data.add.legendre.repl.n400.rate) %>% mutate(N=400),
                      do.call(rbind, data.add.legendre.repl.n450.rate) %>% mutate(N=450)) %>% mutate(Model="m2")

np.rate.all <- rbind(do.call(rbind, data.np.legendre.repl.n50.rate) %>% mutate(N=50),
                     do.call(rbind, data.np.legendre.repl.n100.rate) %>% mutate(N=100),
                     do.call(rbind, data.np.legendre.repl.n150.rate) %>% mutate(N=150),
                     do.call(rbind, data.np.legendre.repl.n200.rate) %>% mutate(N=200),
                     do.call(rbind, data.np.legendre.repl.n250.rate) %>% mutate(N=250),
                     do.call(rbind, data.np.legendre.repl.n300.rate) %>% mutate(N=300),
                     do.call(rbind, data.np.legendre.repl.n350.rate) %>% mutate(N=350),
                     do.call(rbind, data.np.legendre.repl.n400.rate) %>% mutate(N=400),
                     do.call(rbind, data.np.legendre.repl.n450.rate) %>% mutate(N=450)) %>% mutate(Model="m3")

options(scipen = 999)

# ---------------- comparison (mean) -------------------------------------------
rate.all <- rbind(param.rate.all, add.rate.all, np.rate.all) %>%
  mutate(Model=factor(Model, ordered=TRUE, levels=models),
         Method=factor(Method, ordered=TRUE, levels=methods)) %>%
  group_by(Method, Model, Sigma, N) %>%
  summarise(mean.mse=mean(MSE), .groups="drop") %>%
  mutate(mean.true.mse=round(ifelse(Sigma==sigma[1],mean.mse-sigma[1]^2, mean.mse-sigma[2]^2), digits=8)) %>%
  group_by(Method, Model, Sigma) %>%
  filter(N!=50) %>%
  mutate(mean.true.mse.n100=mean.true.mse[N==100]) %>%
  mutate(ratio.true.mse=mean.true.mse/mean.true.mse.n100) %>%
  select(-mean.true.mse.n100)

plot.sigma <- rate.all %>% 
  ggplot(data=., aes(N, ratio.true.mse, lty=Model))+
  geom_point(shape=1)+
  geom_line()+
  facet_wrap(vars(Sigma, Method), scales = "fixed")+
  labs(x="Sample Size", y="Ratio of True MSE")+
  #lims(y=c(0,1))+
  theme_bw()+
  theme(strip.background =element_rect(fill="white"))

tabular.rate <- tabular((Factor(N, "N")) ~ (Factor(Sigma, "Sigma"))*(Factor(Method, "Method"))*(Factor(Model, "Model"))*identity*(mean.true.mse+ratio.true.mse), data=rate.all)
# toLatex(tabular.rate)

# ---------------- comparison (median) -------------------------------------------
rate.all.md <- rbind(param.rate.all, add.rate.all, np.rate.all) %>%
  mutate(Model=factor(Model, ordered=TRUE, levels=models),
         Method=factor(Method, ordered=TRUE, levels=methods)) %>%
  group_by(Method, Model, Sigma, N) %>%
  summarise(median.mse=median(MSE), .groups="drop") %>%
  mutate(median.true.mse=round(ifelse(Sigma==sigma[1],median.mse-sigma[1]^2, median.mse-sigma[2]^2), digits=8)) %>%
  group_by(Method, Model, Sigma) %>%
  filter(N!=50) %>% 
  mutate(median.true.mse.n100=median.true.mse[N==100]) %>%
  mutate(ratio.true.mse=median.true.mse/median.true.mse.n100) %>%
  select(-median.true.mse.n100)

plot.sigma.md <- rate.all.md %>% 
  ggplot(data=., aes(N, ratio.true.mse, lty=Model))+
  geom_point(shape=1)+
  geom_line()+
  facet_wrap(vars(Sigma, Method), scales = "fixed")+
  labs(x="Sample Size", y="Ratio of True MSE")+
  #lims(y=c(0,1))+
  theme_bw()+
  theme(strip.background =element_rect(fill="white"))

tabular.rate.md <- tabular((Factor(N, "N")) ~ (Factor(Sigma, "Sigma"))*(Factor(Method, "Method"))*(Factor(Model, "Model"))*identity*(median.true.mse+ratio.true.mse), data=rate.all.md)
# toLatex(tabular.rate.md)

options(scipen = 0)
