log_profile_restricted_likelihood_IEB <- function(x.tilde,kappa,y.tilde,d){
  # This function computes the log profile restricted likelihood for the case
  # when there are multiple regressors.
  # This is for GWAS structure with kinship random effects.
  # Arguments to the function:
  #    kappa: "variance" parameter for the kinship random effects.
  #    y.tilde: dependent variable in the spectral domain.
  #    x.tilde: regressor in the spectral domain.
  #    d: vector with the eigenvalues of the kinship matrix.
  n <- nrow(x.tilde)
  p <- ncol(x.tilde)
  D.star.inv <- 1/(1 + kappa*d)
  xty <- t(x.tilde) %*% (D.star.inv*y.tilde)
  xtx <- t(x.tilde) %*% (D.star.inv*x.tilde)
  aux <- eigen(xtx,symmetric = TRUE)
  P <- aux$vectors
  eigH <- aux$values
  eigH[abs(eigH) < 10^(-10)] <- 0
  eigH <- 1/eigH
  eigH[!is.finite(eigH)] <- 0
  beta.reml <- P%*%(eigH*t(P)) %*% xty
  residual <- y.tilde - x.tilde %*% beta.reml
  sum.squares <- sum(D.star.inv * residual^2)
  p <- sum(eigH > 0)
  sigma.reml <- as.numeric(sum.squares/(n - p))
  as.numeric(-0.5*(n - p)*log(2*pi*sigma.reml)  + 0.5 * sum(log(D.star.inv)) - (n-p)/2 + 0.5*sum(log(eigH[eigH != 0])))
}

log_profile_restricted_likelihood_IEB_1Dregressor <- function(x.tilde,kappa,y.tilde,d){
  # This function computes the log profile restricted likelihood for the case
  # when the regressor is one-dimensional.
  # This is for GWAS structure with kinship random effects.
  # Arguments to the function:
  #    kappa: "variance" parameter for the kinship random effects.
  #    y.tilde: dependent variable in the spectral domain.
  #    x.tilde: regressor in the spectral domain.
  #    d: vector with the eigenvalues of the kinship matrix.
  n <- length(y.tilde)     # sample size
  p <- 1
  D.star.inv <- 1/(1 + kappa*d)
  xty <- sum(D.star.inv * x.tilde * y.tilde)
  xtx <- sum(D.star.inv * x.tilde^2)
  beta.reml <- xty / xtx
  residual <- y.tilde - x.tilde * beta.reml
  sum.squares <- sum(D.star.inv * residual^2)
  sigma.reml <- as.numeric(sum.squares/(n - p))
  as.numeric(-.5*(n - p)*log(2*pi*sigma.reml)  + .5 * sum(log(D.star.inv)) - (n-p)/2 - .5*log(xtx))
}

REML_estimation_IEB <- function(x.tilde,y.tilde,d){
  # This function computes REML estimates for the case
  # when there are multiple regressors.
  # This is for GWAS structure with kinship random effects.
  # Arguments to the function:
  #    y.tilde: dependent variable in the spectral domain.
  #    x.tilde: regressor in the spectral domain.
  #    d: vector with the eigenvalues of the kinship matrix.
  n <- length(y.tilde)     # sample size
  p <- ncol(x.tilde)
  kappa <- stats::optimize(log_profile_restricted_likelihood_IEB, interval = c(0,100), x.tilde = x.tilde, y.tilde = y.tilde, d = d, maximum = TRUE)$maximum
  D.star.inv <- 1/(1 + kappa*d)
  xty <- t(x.tilde) %*% (D.star.inv*y.tilde)
  xtx <- t(x.tilde) %*% (D.star.inv*x.tilde)
  aux <- eigen(xtx,symmetric = TRUE)
  P <- aux$vectors
  eigH <- aux$values
  eigH[abs(eigH) < 10^(-10)] <- 0
  eigH <- 1/eigH
  eigH[!is.finite(eigH)] <- 0
  beta.reml <- P%*%(eigH*t(P)) %*% xty
  residual <- y.tilde - x.tilde %*% beta.reml
  sum.squares <- sum(D.star.inv * residual^2)
  p <- sum(eigH > 0)
  sigma.reml <- as.numeric(sum.squares/(n - p))

  return(list(beta=beta.reml, sigma2=sigma.reml, kappa=kappa))
}

REML_estimation_IEB_1Dregressor <- function(x.tilde,y.tilde,d){
  # This function computes REML estimates for the case
  # when the regressor is one-dimensional.
  # This is for GWAS structure with kinship random effects.
  # Arguments to the function:
  #    y.tilde: dependent variable in the spectral domain.
  #    x.tilde: regressor in the spectral domain.
  #    d: vector with the eigenvalues of the kinship matrix.
  n <- length(y.tilde)     # sample size
  p <- 1
  kappa <- stats::optimize(log_profile_restricted_likelihood_IEB_1Dregressor,interval = c(0,100),x.tilde = x.tilde, y.tilde = y.tilde, d = d, maximum = TRUE)$maximum
  D.star.inv <- 1/(1 + kappa*d)
  xty <- sum(D.star.inv * x.tilde * y.tilde)
  xtx <- sum(D.star.inv * x.tilde^2)
  beta.reml <- xty / xtx
  residual <- y.tilde - x.tilde * beta.reml
  sum.squares <- sum(D.star.inv * residual^2)
  sigma.reml <- as.numeric(sum.squares/(n - p))

  return(list(beta=beta.reml, sigma2=sigma.reml, kappa=kappa))
}

log_marginal_likelihood_P3D_null_IEB <- function(one.tilde,y.tilde,D.star.inv){
  # This function computes the log marginal likelihood for the case
  # when the regressor is a vector of ones. This assumes population
  # pre-determined parameters (P3D).
  # This is for GWAS structure with kinship random effects.
  # Arguments to the function:
  #    sigma2: variance of the error term. This is estimated somewhere else.
  #    one.tilde: vector of ones in the spectral domain.
  #    y.tilde: dependent variable in the spectral domain.
  #    D.star.inv: vector with the elements of a diagonal matrix that has the eigenvalues of the precision matrix of the observations.
  n <- length(y.tilde)     # sample size
  p <- 1
  xty <- sum(D.star.inv * one.tilde * y.tilde)
  xtx <- sum(D.star.inv * one.tilde^2)
  beta.reml <- xty / xtx
  residual <- y.tilde - one.tilde * beta.reml
  sum.squares <- sum(D.star.inv * residual^2)
  sigma2 <- sum.squares/(n - p)
  as.numeric(-.5*(n - p)*log(2*pi*sigma2)  + .5 * sum(log(D.star.inv)) - 0.5*sum.squares/sigma2 - .5*log(xtx))
}


log.marginal.likelihood.unit <- function(k,sigma2,xtx,xty,y.tilde,D.star.inv){
  n <- length(y.tilde)
  return(-0.5*n*log(2*pi*sigma2) - 0.5*k*log(n + 1) - 0.5*sum(log(D.star.inv[D.star.inv != 0])) - (0.5/sigma2)*(t(y.tilde)%*%(D.star.inv * y.tilde) - (n/(n + 1))*t(xty)%*%solve(xtx)%*%xty))
}

unit_oneiter <- function(y,X,P,D,FDR.threshold = 0.95,maxiterations = 1000, runs_til_stop = 100, v1=1, v2=1, Xc,indices_previous,initial_pi0,iter_count,kinship){
  # prior of pi0: beta(v1=2, v2=2)
  # fix pi0 in model selection
  y.tilde <- t(P)%*%y
  one <- matrix(1,ncol = 1,nrow = length(y.tilde))
  one.tilde <- t(P)%*%one
  y_MS <- y - one%*%REML_estimation_IEB_1Dregressor(one.tilde,y.tilde, D)$beta
  D_MS <- D
  n <- length(y)

  if(is.null(Xc)){
    Xc <- one
    Xc.tilde <- one.tilde
    p_fixed <- ncol(Xc.tilde)
  }else{
    Xc <- cbind(one, Xc)
    Xc.tilde <- t(P)%*%Xc
    p_fixed <- ncol(Xc.tilde)
  }

  if(Matrix::rankMatrix(Xc.tilde)[1] < ncol(Xc.tilde)){
    dropped_cols <- caret::findLinearCombos(Xc.tilde)$remove
    REML.estimates <- REML_estimation_IEB(x.tilde = Xc.tilde[,-dropped_cols,drop = FALSE],y.tilde = y.tilde,d = D)
    alpha <- REML.estimates$beta
    kappa <- REML.estimates$kappa

    y_SMA <- y - Xc[,-dropped_cols,drop = FALSE]%*%alpha

    V <- diag(n) + kappa*kinship
    aux <- eigen(V - Xc[,-dropped_cols,drop = FALSE] %*% solve(t(Xc[,-dropped_cols,drop = FALSE])%*%solve(V)%*%Xc[,-dropped_cols,drop = FALSE])%*%t(Xc[,-dropped_cols,drop = FALSE]),symmetric = TRUE)
    p_fixed <- ncol(Xc[,-dropped_cols,drop = FALSE])
  }else{
    REML.estimates <- REML_estimation_IEB(x.tilde = Xc.tilde,y.tilde = y.tilde,d = D)
    alpha <- REML.estimates$beta
    kappa <- REML.estimates$kappa

    y_SMA <- y - Xc%*%alpha

    V <- diag(n) + kappa*kinship
    aux <- eigen(V - Xc %*% solve(t(Xc)%*%solve(V)%*%Xc)%*%t(Xc),symmetric = TRUE)
  }

  P <- aux$vectors
  eigH <- aux$values
  eigH[abs(eigH) < 10^(-10)] <- 0
  D <- eigH
  rm(aux);rm(eigH)

  y.tilde_SMA <- t(P)%*%y_SMA
  X.tilde <- t(P)%*%X

  D.star.inv <- 1/D
  D.star.inv[!is.finite(D.star.inv)] <- 0

  xj.t.xj <- apply(X.tilde*(D.star.inv*X.tilde),2,sum)
  xj.t.y <- t(X.tilde) %*% (D.star.inv*y.tilde_SMA)
  beta.hat <- xj.t.y / xj.t.xj
  sigma.star2 <- as.numeric(((beta.hat^2)*xj.t.xj - 2*beta.hat*c(xj.t.y) + as.numeric(t(y.tilde_SMA)%*%(D.star.inv*y.tilde_SMA)))/(n - p_fixed - 1))
  var.beta.hat <- sigma.star2*(1/xj.t.xj)

  return_dat <- cbind(beta.hat,var.beta.hat)
  return_dat <- as.data.frame(return_dat)
  colnames(return_dat) <- c("Beta_Hat","Var_Beta_Hat")

  indices_bad <- abs(beta.hat) < 10e-12 | (xj.t.xj < 10e-8 & abs(beta.hat) > 1000)

  if(!is.null(indices_previous)){
    indices_bad[indices_previous] <- TRUE
    beta.hat<- beta.hat[!indices_bad]
    var.beta.hat <- var.beta.hat[!indices_bad]
  }else{
    beta.hat<- beta.hat[!indices_bad]
    var.beta.hat <- var.beta.hat[!indices_bad]
  }

  log_predictive_density = function(param) # function(y, m)
  {
    # Computes log marginal likelihood for 1D regressors as a function of tau and pi0.
    # This assumes that the prior for the regression coefficient is a mixture of a point of
    # mass at zero and a nonlocal prior. This also assumes that the estimates of the different
    # regression coefficients are (approximately) independent. Thus, this is kind of a
    # pseudo likelihood.
    # Maximization of this function leads to empirical Bayes estimates of tau and pi0.
    # This is much a faster (and maybe more precise) alternative to the Langaas' method.
    pi0 <- param[1]
    return(sum(log(pi0*stats::dnorm(beta.hat,mean= 0,sd=sqrt(var.beta.hat)) +
                     (1-pi0)*(2*pi*(n + 1)*var.beta.hat)^(-0.5)*exp(-0.5*(beta.hat^2)*(1/((n + 1)*var.beta.hat)))
    ))+log(stats::dbeta(x = pi0, shape1 = v1, shape2 = v2)))
  }

  pi0.hat <- stats::optimize(log_predictive_density,lower = 0.5,upper = 1,maximum = TRUE)$maximum

  # Compute posterior probability of beta_j different than 0:
  numerator <- (1-pi0.hat)*(2*pi*(n + 1)*var.beta.hat)^(-0.5)*exp(-0.5*(beta.hat^2)*(1/((n + 1)*var.beta.hat)))
  denominator <- pi0.hat*stats::dnorm(beta.hat,mean= 0,sd=sqrt(var.beta.hat)) +
    (1-pi0.hat)*(2*pi*(n + 1)*var.beta.hat)^(-0.5)*exp(-0.5*(beta.hat^2)*(1/((n + 1)*var.beta.hat)))

  postprob <- vector(mode = "numeric",length = ncol(X.tilde))
  postprob[!indices_bad] <- numerator / denominator
  postprob[indices_bad] <- 0

  ###############  Bayesian FDR  #################

  order.postprob <- order(postprob, decreasing=TRUE)
  postprob.ordered <- postprob[order.postprob]

  FDR.Bayes <- cumsum(postprob.ordered) / 1:ncol(X.tilde)
  if(sum(FDR.Bayes > FDR.threshold) == 0){
    return_dat <- cbind(return_dat,postprob,FALSE)
    return_dat <- as.data.frame(return_dat)
    colnames(return_dat) <- c("Beta_Hat","Var_Beta_Hat","PostProb","Significant")
  }else{
    return_dat <- cbind(return_dat,postprob,postprob >= postprob.ordered[max(which(FDR.Bayes > FDR.threshold))])
    return_dat <- as.data.frame(return_dat)
    colnames(return_dat) <- c("Beta_Hat","Var_Beta_Hat","PostProb","Significant")
  }

  if(iter_count != 1){
    pi0.hat <- initial_pi0
  }

  if(sum(return_dat$Significant) > 0){
    indices.significant <- which(return_dat$Significant)
    X.tilde.significant <- X.tilde[, indices.significant,drop = FALSE]

    if(!is.null(Xc)){
      X.tilde.significant <- cbind(Xc.tilde[,-1],X.tilde.significant)
      indices.significant <- c(indices_previous,indices.significant)
    }

    if(Matrix::rankMatrix(cbind(one.tilde,X.tilde.significant))[1] < ncol(cbind(one.tilde,X.tilde.significant))){
      dropped_cols <- caret::findLinearCombos(cbind(one.tilde,X.tilde.significant))$remove
      REML.estimates.ms <- REML_estimation_IEB(cbind(one.tilde,X.tilde.significant)[,-dropped_cols,drop = FALSE], y.tilde, D_MS)
    }else{
      REML.estimates.ms <- REML_estimation_IEB(cbind(one.tilde,X.tilde.significant), y.tilde, D_MS)
    }
    kappa.hat <- REML.estimates.ms$kappa

    V <- diag(n) + kappa.hat*kinship
    aux <- eigen(V - one %*% solve(t(one)%*%solve(V)%*%one)%*%t(one),symmetric = TRUE)
    P <- aux$vectors
    eigH <- aux$values
    D <- rep(0,length(eigH))
    D[1:(n - p_fixed)] <- eigH[1:(n - p_fixed)]
    rm(aux);rm(eigH)

    X.tilde.significant <- t(P)%*% X[, indices.significant,drop = FALSE]
    y.tilde_MS <- t(P)%*%y_MS
    one.tilde <- t(P)%*%one

    D.star.inv <- rep(0,length(D))
    D.star.inv[1:(n - p_fixed)] <- (1/D)[1:(n - p_fixed)]
    X.significant.ty <- t(X.tilde.significant) %*% (D.star.inv*y.tilde_MS)
    X.significant.t.X.significant <- t(X.tilde.significant) %*% (D.star.inv*X.tilde.significant)

    total.p <- ncol(X.tilde.significant)

    if(total.p < 16){
      # Do full model search
      total.models <- 2^total.p
      log.unnormalized.posterior.probability <- rep(NA, total.models)
      log.unnormalized.posterior.probability[1] <- total.p * log(pi0.hat) + log_marginal_likelihood_P3D_null_IEB(one.tilde,y.tilde_MS,D.star.inv)
      dat <- rep(list(0:1), total.p)
      dat <- as.matrix(expand.grid(dat))
      for (i in 1:(total.models-1)){
        model <- unname(which(dat[i + 1,] == 1))
        Xsub <- X.tilde.significant[,model,drop = FALSE]
        k <- length(model)
        if(Matrix::rankMatrix(Xsub)[1] < ncol(Xsub)){
          dropped_cols <- caret::findLinearCombos(Xsub)$remove
          model <- model[-dropped_cols]
        }
        sub_xtx <- matrix(X.significant.t.X.significant[model, model],ncol = length(model),nrow = length(model))
        sub_xty <- matrix(X.significant.ty[model,],ncol = 1)
        sub_x <- matrix(X.tilde.significant[, model],ncol = length(model))
        beta.reml <- solve(sub_xtx) %*% sub_xty
        residual <- y.tilde_MS - sub_x %*% beta.reml
        sum.squares <- sum(D.star.inv * residual^2)
        sigma2.hat <- as.numeric(sum.squares/(n - length(model) - 1))
        log.unnormalized.posterior.probability[i+1] <-
          k*log(1-pi0.hat) + (total.p-k)*log(pi0.hat) +
          log.marginal.likelihood.unit(k, sigma2.hat,
                                       sub_xtx,
                                       sub_xty,
                                       y.tilde_MS,
                                       D.star.inv)
      }
      log.unnormalized.posterior.probability <- log.unnormalized.posterior.probability - max(log.unnormalized.posterior.probability)
      unnormalized.posterior.probability <- exp(log.unnormalized.posterior.probability)
      posterior.probability <- unnormalized.posterior.probability/sum(unnormalized.posterior.probability)

    } else {
      # Do model search with genetic algorithm
      fitness_ftn <- function(string){
        if(sum(string) == 0){
          return(total.p * log(pi0.hat) + log_marginal_likelihood_P3D_null_IEB(one.tilde,y.tilde_MS,D.star.inv))
        }
        model <- which(string==1)
        Xsub <- X.tilde.significant[,model,drop = FALSE]
        k <- length(model)
        if(Matrix::rankMatrix(Xsub)[1] < ncol(Xsub)){
          dropped_cols <- caret::findLinearCombos(Xsub)$remove
          model <- model[-dropped_cols]
        }
        sub_xtx <- matrix(X.significant.t.X.significant[model, model],ncol = length(model),nrow = length(model))
        sub_xty <- matrix(X.significant.ty[model,],ncol = 1)
        sub_x <- matrix(X.tilde.significant[, model],ncol = length(model))
        beta.reml <- solve(sub_xtx) %*% sub_xty
        residual <- y.tilde_MS - sub_x %*% beta.reml
        sum.squares <- sum(D.star.inv * residual^2)
        sigma2.hat <- as.numeric(sum.squares/(n - length(model) - 1))

        return(k*log(1-pi0.hat) + (total.p-k)*log(pi0.hat) +
                 log.marginal.likelihood.unit(k, sigma2.hat,
                                              sub_xtx,
                                              sub_xty,
                                              y.tilde_MS,
                                              D.star.inv))
      }
      if(total.p > 99){
        suggestedsol <- diag(total.p)
        tmp_log.unnormalized.posterior.probability <- vector()
        for(i in 1:total.p){
          model <- which(suggestedsol[i,]==1)
          sub_xtx <- matrix(X.significant.t.X.significant[model, model],ncol = length(model),nrow = length(model))
          sub_xty <- matrix(X.significant.ty[model,],ncol = 1)
          sub_x <- matrix(X.tilde.significant[, model],ncol = length(model))
          beta.reml <- solve(sub_xtx) %*% sub_xty
          residual <- y.tilde_MS - sub_x %*% beta.reml
          sum.squares <- sum(D.star.inv * residual^2)
          sigma2.hat <- as.numeric(sum.squares/(n - length(model) - 1))
          k <- length(model)
          tmp_log.unnormalized.posterior.probability[i] <- k*log(1-pi0.hat) + (total.p-k)*log(pi0.hat) +
            log.marginal.likelihood.unit(k, sigma2.hat,
                                         sub_xtx,
                                         sub_xty,
                                         y.tilde_MS,
                                         D.star.inv)
        }
        suggestedsol <- rbind(0,suggestedsol[order(tmp_log.unnormalized.posterior.probability,decreasing = TRUE)[1:99],])
      }else{
        suggestedsol <- rbind(0,diag(total.p))
      }

      fitness_ftn <- memoise::memoise(fitness_ftn)
      ans <- GA::ga("binary", fitness = fitness_ftn, nBits = total.p,maxiter = maxiterations,popSize = 100,elitism = min(c(10,2^total.p)),run = runs_til_stop,suggestions = suggestedsol,monitor = FALSE)
      memoise::forget(fitness_ftn)
      dat <- ans@population
      dupes <- duplicated(dat)
      dat <- dat[!dupes,]
      ans@fitness <- ans@fitness[!dupes]
      log.unnormalized.posterior.probability <- ans@fitness - max(ans@fitness)
      unnormalized.posterior.probability <- exp(log.unnormalized.posterior.probability)
      posterior.probability <- unnormalized.posterior.probability/sum(unnormalized.posterior.probability)
    }

    inclusion_prb <- unname((t(dat)%*%posterior.probability)/sum(posterior.probability))

    model <- dat[which.max(posterior.probability),]
    model_dat <- cbind(indices.significant,model,inclusion_prb)
    model_dat <- as.data.frame(model_dat)
    colnames(model_dat) <- c("SNPs","BestModel","Inclusion_Prob")

    tmp <- stats::cor(cbind(t(P)%*%X[,indices_bad,drop = FALSE],X.tilde.significant))
    indices.correlated <- which(tmp == 1 & lower.tri(tmp), arr.ind = TRUE, useNames = FALSE)
    indices.correlated <- unique(c(indices.correlated[indices.correlated[,1] %in% which(model == 1),2],indices.correlated[indices.correlated[,2] %in% which(model == 1),1]))

    return(list(prescreen = return_dat,modelselection = model_dat,correlated_SNPs = c(which(indices_bad),indices.significant)[indices.correlated],pi_0_hat = pi0.hat))

  }else{
    return(list(prescreen = return_dat,modelselection = "No significant in prescreen1",correlated_SNPs = "None",pi_0_hat = pi0.hat))
  }
}

#' @title  IEB for Gaussian Phenotypes
#' @description  Performs IEB
#' @param Y The observed numeric phenotypes
#' @param SNPs The SNP matrix, where each column represents a single SNP encoded as the numeric coding 0, 1, 2. This is entered as a matrix object.
#' @param FDR_Nominal The nominal false discovery rate for which SNPs are selected from in the screening step.
#' @param kinship The observed kinship matrix, has to be a square positive semidefinite matrix. Defaulted as the identity matrix. The function used
#' to create the kinship matrix used in the BICOSS paper is A.mat() from package rrBLUP.
#' @param maxiterations The maximum iterations the genetic algorithm in the model selection step iterates for.
#' Defaulted at 400 which is the value used in the BICOSS paper simulation studies.
#' @param runs_til_stop The number of iterations at the same best model before the genetic algorithm in the model selection step converges.
#' Defaulted at 40 which is the value used in the BICOSS paper simulation studies.
#' @return The column indices of SNPs that were in the best model identified by BICOSS.
#' @examples
#' library(GWAS.BAYES)
#' IEB(Y = Y, SNPs = SNPs, kinship = kinship,
#'     FDR_Nominal = 0.05,
#'     maxiterations = 400,runs_til_stop = 40)
#' @export
IEB <- function(Y,SNPs,kinship,FDR_Nominal = 0.05,maxiterations = 400,runs_til_stop = 40){

  if(FDR_Nominal > 1 | FDR_Nominal < 0){
    stop("FDR_Nominal has to be between 0 and 1")
  }
  if(!is.numeric(Y)){
    stop("Y has to be numeric")
  }
  if(!is.matrix(SNPs)){
    stop("SNPs has to be a matrix object")
  }
  if(!is.numeric(SNPs)){
    stop("SNPs has to contain numeric values")
  }
  if(nrow(kinship) != ncol(kinship)){
    stop("kinship has to be a square matrix")
  }
  if(nrow(kinship) != nrow(SNPs)){
    stop("kinship has to have the same number of rows as SNPs")
  }
  if(!is.numeric(kinship)){
    stop("kinship has to contain numeric values")
  }
  if(maxiterations-floor(maxiterations)!=0){
    stop("maxiterations has to be a integer")
  }
  if(runs_til_stop-floor(runs_til_stop)!=0){
    stop("runs_til_stop has to be a integer")
  }
  if(maxiterations < runs_til_stop){
    stop("maxiterations has to be larger than runs_til_stop")
  }

  y <- Y
  rm(Y)

  X <- scale(SNPs)
  rm(SNPs)

  FDR.threshold <- 1 - FDR_Nominal
  rm(FDR_Nominal)

  v1 = 1;v2 = 1

  aux <- eigen(kinship,symmetric = TRUE)
  P <- aux$vectors
  eigH <- aux$values
  D <- c(eigH[-length(eigH)],0)
  rm(aux);rm(eigH)

  count <- 1
  p <- 1
  indices_previous <- NULL
  indices_sig_list <- list()
  indices_cor_list <- list()
  Xc <- NULL
  initial_pi0 <- NULL

  while(p > 0 & count < 10){
    tmp <- unit_oneiter(y,X,P,D,FDR.threshold = FDR.threshold,maxiterations = maxiterations, runs_til_stop = runs_til_stop, v1=v1, v2=v2,Xc = Xc,indices_previous = indices_previous,initial_pi0 = initial_pi0,iter_count = count,kinship = kinship)

    if(is.character(tmp$modelselection)){
      break
    }

    indices_sig_list[[count]] <- tmp$modelselection$SNPs[tmp$modelselection$BestModel == 1]
    indices_previous <- tmp$modelselection$SNPs[tmp$modelselection$BestModel == 1]

    if(count == 1){
      p <- length(tmp$modelselection$SNPs[tmp$modelselection$BestModel == 1])
      initial_pi0 <- tmp$pi_0_hat
    }else{
      p <- length(tmp$modelselection$SNPs[tmp$modelselection$BestModel == 1][!(tmp$modelselection$SNPs[tmp$modelselection$BestModel == 1]%in%indices_sig_list[[count - 1]])])
    }

    Xc <- X[,tmp$modelselection$SNPs[tmp$modelselection$BestModel == 1]]

    indices_cor_list[[count]] <- tmp$correlated_SNPs

    count <- count + 1
  }

  if(count == 1 & is.character(tmp$modelselection)){
    return(list(best_model = "No SNPs found in Screening"))
  }else{

    tmp <- stats::cor(cbind(X[,unlist(indices_cor_list)],X[,indices_sig_list[[count - 1]]]))
    indices.correlated <- which(tmp == 1 & lower.tri(tmp), arr.ind = TRUE, useNames = FALSE)
    indices.correlated <- unique(c(indices.correlated[,2],indices.correlated[,1]))
    indices.correlated <- c(unlist(indices_cor_list),indices_sig_list[[count - 1]])[indices.correlated]
    indices.correlated <- indices.correlated[!(indices.correlated %in% indices_sig_list[[count - 1]])]

    return(list(best_model = indices_sig_list[[count - 1]]))
  }
}
