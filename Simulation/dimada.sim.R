

dimada.sim <- function(n=200,
                       n.test=1000,
                       d=5,
                       formula,
                       empirical.mse=TRUE,
                       sigma=c(0.05,0.2),
                       basis="bspline",
                       est=c("unrestricted","additive","ols"),
                       family="gaussian",
                       nfolds=10,
                       parallel=FALSE,
                       lambda.min.ratio=1e-10)
  {

  # ---------------------------
  # Check Arguments
  # ---------------------------

  # check if 'formula' is a single character
  if (!is.character(formula) | length(formula)!=1) stop("The argument 'formula' has to be a single character.")

  # check if 'empirical.mse' is a single logical value
  if (!is.logical(empirical.mse) | length(empirical.mse)!=1) stop("The argument 'empirical.mse' has to be a single logical value.")

  # check if 'sigma' is a vector of non-negative numbers
  if (!is.numeric(sigma) | !is.vector(sigma) | any(sigma <0)) stop("The argument 'sigma' should be a vector of non-negative numbers")

  # check if 'basis' is defined correctly
  if (any(length(basis)!=1, !(basis %in% c("power","legendre","bspline","cosine","sine","trig","haar","daubechies")))) stop("The argument 'basis' should be 'power', 'legendre', 'bspline', 'cosine', 'sine', 'trig', 'haar' or 'daubechies'.")

  # check if 'est' is defined correctly
  if (!("unrestricted" %in% est)) stop("The baseline estimator 'unrestricted' has to be in the argument 'est'.")
  if (any(!(est %in% c("unrestricted","additive","ols")))) stop("Elements in the argument 'est' have to be among 'unrestricted','additive' and 'ols'.")

  # check if 'family' is defined correctly
  if (any(length(family)!=1, !(family %in% c("gaussian", "binomial", "poisson", "multinomial", "cox", "mgaussian")))) stop("Elements in the argument 'family' have to be 'gaussian', 'binomial', 'poisson', 'multinomial', 'cox' or 'mgaussian'.")

  # check if 'nfolds' is a single integer between 3 and the number of observations
  if (any(nfolds%%1!=0, nfolds>n, nfolds<3, length(nfolds)!=1)) stop("The argument 'nfolds' has to be a single integer between 3 and the number of observations.")

  # check if 'parallel' is a single logical value
  if (!is.logical(parallel) | length(parallel)!=1) stop("The argument 'parallel' has to be a single logical value.")

  # check if 'lambda.min.ratio' is defined correctly
  if (any(lambda.min.ratio>1, lambda.min.ratio<0, length(lambda.min.ratio)!=1)) stop("The argument 'lambda.min.ratio' has to be a single ratio number.")

  # ---------------------------
  # Generate original data set
  # ---------------------------

  X <- matrix(runif(n*d, 0, 1), n, d)
  test.data <- matrix(runif(n.test*d, 0, 1), n.test, d)
  epsilon <- rnorm(n, 0, 1)
  epsilon.test <- rnorm(n.test,0,1)

  Y.true <- eval(parse(text=formula))
  Y <- plyr::llply(1:length(sigma), function(x) Y.true + epsilon*sigma[x])
  Y.names <- paste("Y.sigma", 1:length(sigma), sep = "")
  names(Y) <- Y.names

  X.temp <- X; X <- test.data

  Y.true.test <- eval(parse(text=formula))
  X<- X.temp

  Y.test <- plyr::llply(1:length(sigma), function(x) Y.true.test + epsilon.test*sigma[x])
  Y.test.names <- paste("Y.test.sigma", 1:length(sigma), sep = "")
  names(Y.test) <- Y.test.names

  # ---------------------------
  # Generate sieves
  # ---------------------------

  if (basis=="power") {
    # unrestricted underlying model
    sieve.X.free <- poly.gen(data=X, test.data=test.data, n.basis=NULL, max.interaction=d, legendre=FALSE)
    # additive underlying model
    if ("additive" %in% est) sieve.X.add <- poly.gen(data=X, test.data=test.data, n.basis=NULL, max.interaction=1, legendre=FALSE)
    # parametric underlying model
    # sieve.X.param <- poly.gen(data=X, test.data=test.data, n.basis=4, max.interaction=d, legendre=FALSE)
  }

  if (basis=="legendre") {
    # unrestricted underlying model
    sieve.X.free <- poly.gen(data=X, test.data=test.data, n.basis=NULL, max.interaction=d, legendre=TRUE)
    # additive underlying model
    if ("additive" %in% est)  sieve.X.add <- poly.gen(data=X, test.data=test.data, n.basis=NULL, max.interaction=1, legendre=TRUE)
    # parametric underlying model
    # sieve.X.param <- poly.gen(data=X, test.data=test.data, n.basis=4, max.interaction=d, legendre=TRUE)
  }

  if (basis=="bspline") {
    # unrestricted underlying model
    sieve.X.free <- bspline.gen(data=X, test.data=test.data, n.basis=NULL, max.interaction=d, spline.degree=2)
    # additive underlying model
    if ("additive" %in% est)  sieve.X.add <- bspline.gen(data=X, test.data=test.data, n.basis=NULL, max.interaction=1, spline.degree=2)
    # parametric underlying model
    # sieve.X.param <- bspline.gen(data=X, test.data=test.data, n.basis=4, max.interaction=d, spline.degree=2)
  }

  if (basis=="cosine") {
    # unrestricted underlying model
    sieve.X.free <- cosine.gen(data=X, test.data=test.data, n.basis=NULL, max.interaction=d)
    # additive underlying model
    if ("additive" %in% est)  sieve.X.add <- cosine.gen(data=X, test.data=test.data, n.basis=NULL, max.interaction=1)
    # parametric underlying model
    # sieve.X.param <- cosine.gen(data=X, test.data=test.data, n.basis=4, max.interaction=d)
  }

  if (basis=="sine") {
    # unrestricted underlying model
    sieve.X.free <- sine.gen(data=X, test.data=test.data, n.basis=NULL, max.interaction=d)
    # additive underlying model
    if ("additive" %in% est)  sieve.X.add <- sine.gen(data=X, test.data=test.data, n.basis=NULL, max.interaction=1)
    # parametric underlying model
    # sieve.X.param <- sine.gen(data=X, test.data=test.data, n.basis=4, max.interaction=d)
  }

  if (basis=="trig") {
    # unrestricted underlying model
    sieve.X.free <- trig.gen(data=X, test.data=test.data, n.basis=NULL, max.interaction=d)
    # additive underlying model
    if ("additive" %in% est)  sieve.X.add <- trig.gen(data=X, test.data=test.data, n.basis=NULL, max.interaction=1)
    # parametric underlying model
    # sieve.X.param <- trig.gen(data=X, test.data=test.data, n.basis=4, max.interaction=d)
  }

  if (basis=="haar") {
    # unrestricted underlying model
    sieve.X.free <- haar.gen(data=X, test.data=test.data, n.basis=NULL, max.interaction=d)
    # additive underlying model
    if ("additive" %in% est)  sieve.X.add <- haar.gen(data=X, test.data=test.data, n.basis=NULL, max.interaction=1)
    # parametric underlying model
    # sieve.X.param <- haar.gen(data=X, test.data=test.data, n.basis=4, max.interaction=d)
  }

  if (basis=="daubechies") {
    # unrestricted underlying model
    sieve.X.free <- daubechies.gen(data=X, test.data=test.data, n.basis=NULL, max.interaction=d)
    # additive underlying model
    if ("additive" %in% est)  sieve.X.add <- daubechies.gen(data=X, test.data=test.data, n.basis=NULL, max.interaction=1)
    # parametric underlying model
    # sieve.X.param <- daubechies.gen(data=X, test.data=test.data, n.basis=4, max.interaction=d)
  }

  if ("additive" %in% est) {sieve.X.list <- list(sieve.X.free$train, sieve.X.add$train)} else {sieve.X.list <- list(sieve.X.free$train)}
  if ("additive" %in% est) {sieve.X.list.test <- list(sieve.X.free$test, sieve.X.add$test)} else {sieve.X.list.test <- list(sieve.X.free$test)}
  if ("additive" %in% est) {estimators<- c("Unrestricted","Additive")} else {estimators<- c("Unrestricted")}

  # ---------------------------
  # OLS
  # ---------------------------

  # as a comparison, estimate ols with intercept
  if ("ols" %in% est) {
    ols <-plyr::llply(Y, function(x) lm(x~X))
    if (empirical.mse) {
      ols.test.mse <- plyr::ldply(1:length(sigma), function(x) data.frame(mse.test=mean(( cbind(1,test.data)%*%ols[[x]]$coefficients-Y.test[[x]])^2, na.rm=TRUE),
                                                                          estimator="OLS",
                                                                          sigma=sigma[x]))
    } else {
      ols.test.mse <- plyr::ldply(1:length(sigma), function(x) data.frame(mse.test=mean(( cbind(1,test.data)%*%ols[[x]]$coefficients-Y.true.test)^2, na.rm=TRUE),
                                                                          estimator="OLS",
                                                                          sigma=sigma[x]))
    }
  }

  # -------------------------------
  # Dimension Adaptive Estimator
  # -------------------------------

  data <- plyr::llply(Y,function(y) plyr::llply(sieve.X.list,
                                                function(x) dimada(y=as.vector(y), x=NULL, basis=basis, x.sieve=x, methods=c("Lasso","adaLasso"), family=family,
                                                                   nfolds=nfolds, parallel=parallel, lambda.min.ratio=lambda.min.ratio, s="lambda.min")))

  #--------Lasso: test data based-----------

  if(empirical.mse) {
    post.lasso.test.mse <- plyr::ldply(1:length(sigma), function(x) data.frame(plyr::ldply(1:length(estimators),
                                                                                           function(y) data.frame(nzeros=data[[x]][[y]]$Lasso$parameters.final$nzeros,
                                                                                                                  lambdas=data[[x]][[y]]$Lasso$parameters.final$lambdas,
                                                                                                                  mse.test=mean((as.matrix(cbind(1,sieve.X.list.test[[y]][,data[[x]][[y]]$Lasso$coefs.final$terms, drop=FALSE]))%*%ifelse(is.na(data[[x]][[y]]$Lasso$post.lm$coefficients),0,data[[x]][[y]]$Lasso$post.lm$coefficients)-Y.test[[x]])^2, na.rm=TRUE),
                                                                                                                  estimator=estimators[y])),
                                                                               sigma=sigma[x]))
  } else {
    post.lasso.test.mse <- plyr::ldply(1:length(sigma), function(x) data.frame(plyr::ldply(1:length(estimators),
                                                                                           function(y) data.frame(nzeros=data[[x]][[y]]$Lasso$parameters.final$nzeros,
                                                                                                                  lambdas=data[[x]][[y]]$Lasso$parameters.final$lambdas,
                                                                                                                  mse.test=mean((as.matrix(cbind(1,sieve.X.list.test[[y]][,data[[x]][[y]]$Lasso$coefs.final$terms, drop=FALSE]))%*%ifelse(is.na(data[[x]][[y]]$Lasso$post.lm$coefficients),0,data[[x]][[y]]$Lasso$post.lm$coefficients)-Y.true.test)^2, na.rm=TRUE),
                                                                                                                  estimator=estimators[y])),
                                                                               sigma=sigma[x]))
  }

  #--------Adaptive Lasso: test data based-----------


  if(empirical.mse){
    post.adaLasso.test.mse <- plyr::ldply(1:length(sigma), function(x) data.frame(plyr::ldply(1:length(estimators),
                                                                                              function(y) data.frame(nzeros=data[[x]][[y]]$adaLasso$parameters.final$nzeros,
                                                                                                                     lambdas=data[[x]][[y]]$adaLasso$parameters.final$lambdas,
                                                                                                                     mse.test=mean((as.matrix(cbind(1,sieve.X.list.test[[y]][,data[[x]][[y]]$adaLasso$coefs.final$terms, drop=FALSE]))%*%ifelse(is.na(data[[x]][[y]]$adaLasso$post.lm$coefficients),0,data[[x]][[y]]$adaLasso$post.lm$coefficients)-Y.test[[x]])^2, na.rm=TRUE),
                                                                                                                     estimator=estimators[y])),
                                                                                  sigma=sigma[x]))
  } else {
    post.adaLasso.test.mse <- plyr::ldply(1:length(sigma), function(x) data.frame(plyr::ldply(1:length(estimators),
                                                                                              function(y) data.frame(nzeros=data[[x]][[y]]$adaLasso$parameters.final$nzeros,
                                                                                                                     lambdas=data[[x]][[y]]$adaLasso$parameters.final$lambdas,
                                                                                                                     mse.test=mean((as.matrix(cbind(1,sieve.X.list.test[[y]][,data[[x]][[y]]$adaLasso$coefs.final$terms, drop=FALSE]))%*%ifelse(is.na(data[[x]][[y]]$adaLasso$post.lm$coefficients),0,data[[x]][[y]]$adaLasso$post.lm$coefficients)-Y.true.test)^2, na.rm=TRUE),
                                                                                                                     estimator=estimators[y])),
                                                                                  sigma=sigma[x]))
  }

  # -------------------------------
  # Results
  # -------------------------------

  if ("ols" %in% est) {
    result.test <- dplyr::bind_rows(post.lasso.test.mse %>% mutate(method="Post LASSO"),
                                    post.adaLasso.test.mse %>% mutate(method="Post Adaptive LASSO"),
                                    ols.test.mse %>% mutate(method=NA)) %>%
      mutate(estimator=as.factor(estimator), sigma=as.factor(sigma), method=as.factor(method)) %>%
      dplyr::rename(NonZeros=nzeros, MSE=mse.test,Method=method,Estimator=estimator,Sigma=sigma)
  } else {
    result.test <- dplyr::bind_rows(post.lasso.test.mse %>% mutate(method="Post LASSO"),
                                    post.adaLasso.test.mse %>% mutate(method="Post Adaptive LASSO")) %>%
      mutate(estimator=as.factor(estimator), sigma=as.factor(sigma), method=as.factor(method)) %>%
      dplyr::rename(NonZeros=nzeros, MSE=mse.test,Method=method,Estimator=estimator,Sigma=sigma)
  }

  print("*")
  return(result.test)

}
