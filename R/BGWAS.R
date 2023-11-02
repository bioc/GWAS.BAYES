svd1 <- function(x){
  p <- ncol(x)
  d <- svd(x,nu = p,nv = p)$d
  return(list(max(d)/min(d),p - length(d)))
}

log_profile_restricted_likelihood_1Dregressor_BGWAS <- function(x.tilde,kappa,y.tilde,d){
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

REML_estimation_1Dregressor_BGWAS <- function(x.tilde,y.tilde,d){
  # This function computes REML estimates for the case
  # when the regressor is one-dimensional.
  # This is for GWAS structure with kinship random effects.
  # Arguments to the function:
  #    y.tilde: dependent variable in the spectral domain.
  #    x.tilde: regressor in the spectral domain.
  #    d: vector with the eigenvalues of the kinship matrix.
  n <- length(y.tilde)     # sample size
  p <- 1
  kappa <- stats::optimize(log_profile_restricted_likelihood_1Dregressor_BGWAS,interval = c(0,100),x.tilde = x.tilde, y.tilde = y.tilde, d = d, maximum = TRUE)$maximum
  D.star.inv <- 1/(1 + kappa*d)
  xty <- sum(D.star.inv * x.tilde * y.tilde)
  xtx <- sum(D.star.inv * x.tilde^2)
  beta.reml <- xty / xtx
  residual <- y.tilde - x.tilde * beta.reml
  sum.squares <- sum(D.star.inv * residual^2)
  sigma.reml <- as.numeric(sum.squares/(n - p))

  return(list(beta=beta.reml, sigma2=sigma.reml, kappa=kappa))
}

log_profile_restricted_likelihood <- function(x,t,y,d){
  n <- nrow(x)
  p <- ncol(x)
  estar <- 1/(1 + t*d)
  xty <- t(x) %*% (estar*y)
  xtx <- t(x) %*% (estar*x)
  aux <- eigen(xtx,symmetric = TRUE)
  P <- aux$vectors
  eigH <- aux$values
  eigH[abs(eigH) < 10^(-10)] <- 0
  eigH <- 1/eigH
  eigH[!is.finite(eigH)] <- 0
  beta.reml <- P%*%(eigH*t(P)) %*% xty
  residual <- y - x %*% beta.reml
  sum.squares <- sum(estar * residual^2)
  p <- sum(eigH > 0)
  sigma.reml <- as.numeric(sum.squares/(n - p))
  as.numeric(-0.5*(n - p)*log(2*pi*sigma.reml)  + 0.5 * sum(log(estar)) - (n-p)/2 - 0.5*determinant(xtx,logarithm=TRUE)$modulus[1])
}

REML_estimation <- function(x,y,d){
  n <- length(y)
  p <- ncol(x)
  t <- stats::optimize(log_profile_restricted_likelihood, interval = c(0,100), x = x, y = y, d = d, maximum = TRUE)$maximum
  estar <- 1/(1 + t*d)
  xty <- t(x) %*% (estar*y)
  xtx <- t(x) %*% (estar*x)
  aux <- eigen(xtx,symmetric = TRUE)
  P <- aux$vectors
  eigH <- aux$values
  eigH[abs(eigH) < 10^(-10)] <- 0
  eigH <- 1/eigH
  eigH[!is.finite(eigH)] <- 0
  beta.reml <- P%*%(eigH*t(P)) %*% xty
  residual <- y - x %*% beta.reml
  sum.squares <- sum(estar * residual^2)
  p <- sum(eigH > 0)
  sigma.reml <- as.numeric(sum.squares/(n - p))

  return(list(beta=beta.reml, sigma2=sigma.reml, t=t))
}

log_marginal_likelihood_P3D_null <- function(one.tilde,y.tilde,D.star.inv){
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

log.marginal.likelihood.nondiagonal_1sim = function(k,tau,sigma2,xtx,xty,x.tilde,y.tilde,D.star.inv,nboot)
{
  # Computes log marginal likelihood for a model based on a pMOM nonlocal prior.
  # Considers a random effect model for GWAS data with kinship random effects
  # and multiple regressors.
  # Uses the approximation: E(\prod{\beta_j^2}) \approx \prod E({\beta_j})^2

  n <- length(y.tilde)     # sample size
  #  k <- ncol(x.tilde)
  nt_1 <- n*tau + 1
  C.k <- (nt_1/(n*tau))*xtx
  C.k_1 <- solve(C.k)

  beta.k.hat <- C.k_1 %*% xty
  R.k <- t(y.tilde)%*%(D.star.inv * y.tilde) - t(beta.k.hat)%*%C.k%*%beta.k.hat

  E2 <- MASS::mvrnorm(nboot,mu = beta.k.hat,Sigma = sigma2*C.k_1)
  E1 <- sqrt(nt_1)*sweep(E2,2,c(beta.k.hat))

  E2 <- mean(apply(E2^2, 1, prod))
  E1 <- mean(apply(E1^2, 1, prod))

  return(-0.5*n*log(2*pi*sigma2) - 0.5*sum(log(D.star.inv)) - (0.5*k)*log(nt_1) - 0.5*R.k/sigma2 + log(E2) - log(E1))
}

all_Marco_nondiagonal_pi0_known_expected <- function(y,X,P,D,FDR.threshold = 0.95,maxiterations = 1000, runs_til_stop = 100,nboot = 100, TAU, v1=2, v2=2){
  # optimize tau and pi0 unconstrained
  # prior of pi0: beta(v1=2, v2=2)
  # prior dist of pi0 in model selection

  n <- length(y)

  y.tilde <- t(P) %*% y
  X.tilde <- t(P) %*% X

  one.tilde <- t(P) %*% rep(1,length(y))

  REML.estimates <- REML_estimation_1Dregressor_BGWAS(x.tilde = one.tilde,y.tilde = y.tilde,d = D)
  alpha <- REML.estimates$beta
  kappa <- REML.estimates$kappa

  y.tilde <- y.tilde - one.tilde*alpha

  D.star.inv <- 1/(1 + kappa*D)

  xj.t.xj <- apply(X.tilde*(D.star.inv*X.tilde),2,sum)
  xj.t.y <- t(X.tilde) %*% (D.star.inv*y.tilde)
  beta.hat <- xj.t.y / xj.t.xj
  sigma.star2 <- as.numeric(((beta.hat^2)*xj.t.xj - 2*beta.hat*c(xj.t.y) + as.numeric(t(y.tilde)%*%(D.star.inv*y.tilde)))/(length(y) - 2))
  var.beta.hat <- sigma.star2 / xj.t.xj

  return_dat <- cbind(beta.hat,var.beta.hat)
  return_dat <- as.data.frame(return_dat)
  colnames(return_dat) <- c("Beta_Hat","Var_Beta_Hat")

  tau <- TAU

  log.marginal.likelihood_tau = function(param) # function(y, m)
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
                     (1-pi0)*(2*pi*var.beta.hat)^(-0.5) * (n*tau+1)^(-1.5) *
                     (1+n*tau*beta.hat^2/((n*tau+1)*var.beta.hat)) *
                     exp(-beta.hat^2 / (2*var.beta.hat*(n*tau+1)))
    ))+log(stats::dbeta(x = pi0, shape1 = v1, shape2 = v2)))
  }

  pi0.hat <- stats::optimize(log.marginal.likelihood_tau,lower = 0.5,upper = 1,maximum = TRUE)$maximum
  tau.hat <- tau

  # Compute posterior probability of beta_j different than 0:
  numerator <- (1-pi0.hat)*(2*pi*var.beta.hat)^(-0.5) * (n*tau.hat+1)^(-1.5) *
    (1+n*tau.hat*beta.hat^2/((n*tau.hat+1)*var.beta.hat)) *
    exp(-beta.hat^2 / (2*var.beta.hat*(n*tau.hat+1)))
  denominator <- pi0.hat*stats::dnorm(beta.hat,mean= 0,sd=sqrt(var.beta.hat)) +
    (1-pi0.hat)*(2*pi*var.beta.hat)^(-0.5) * (n*tau.hat+1)^(-1.5) *
    (1+n*tau.hat*beta.hat^2/((n*tau.hat+1)*var.beta.hat)) *
    exp(-beta.hat^2 / (2*var.beta.hat*(n*tau.hat+1)))

  postprob <- numerator / denominator

  ###############  Bayesian FDR  #################

  order.postprob <- order(postprob, decreasing=TRUE)
  postprob.ordered <- postprob[order.postprob]

  FDR.Bayes <- cumsum(postprob.ordered) / 1:ncol(X)
  if(sum(FDR.Bayes > FDR.threshold) == 0){
    return_dat <- cbind(return_dat,postprob,FALSE)
    return_dat <- as.data.frame(return_dat)
    colnames(return_dat) <- c("Beta_Hat","Var_Beta_Hat","PostProb","Significant")
  }else{
    return_dat <- cbind(return_dat,postprob,postprob >= postprob.ordered[max(which(FDR.Bayes > FDR.threshold))])
    return_dat <- as.data.frame(return_dat)
    colnames(return_dat) <- c("Beta_Hat","Var_Beta_Hat","PostProb","Significant")
  }

  if(sum(return_dat$Significant) > 0){
    indices.significant <- which(return_dat$Significant)
    X.tilde.significant <- X.tilde[, indices.significant,drop = FALSE]

    if(Matrix::rankMatrix(X.tilde.significant)[1] < ncol(X.tilde.significant)){
      dropped_cols <- caret::findLinearCombos(X.tilde.significant)$remove
      REML.estimates.ms <- REML_estimation(X.tilde.significant[,-dropped_cols,drop = FALSE], y.tilde, D)
    }else{
      REML.estimates.ms <- REML_estimation(X.tilde.significant, y.tilde, D)
    }
    kappa.hat <- REML.estimates.ms$t

    D.star.inv <- 1/(1 + kappa.hat*D)
    X.significant.ty <- t(X.tilde.significant) %*% (D.star.inv*y.tilde)
    X.significant.t.X.significant <- t(X.tilde.significant) %*% (D.star.inv*X.tilde.significant)

    total.p <- ncol(X.tilde.significant)

    if(total.p < 16){
      # Do full model search
      total.models <- 2^total.p
      log.unnormalized.posterior.probability <- rep(NA, total.models)
      log.unnormalized.posterior.probability[1] <- total.p * log(pi0.hat) + log_marginal_likelihood_P3D_null(one.tilde,y.tilde,D.star.inv)
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
        residual <- y.tilde - sub_x %*% beta.reml
        sum.squares <- sum(D.star.inv * residual^2)
        sigma2.hat <- as.numeric(sum.squares/(n - length(model) - 1))
        log.unnormalized.posterior.probability[i+1] <-
          k*log(1-pi0.hat) + (total.p-k)*log(pi0.hat) +
          log.marginal.likelihood.nondiagonal_1sim(k, tau.hat, sigma2.hat,
                                                   sub_xtx,
                                                   sub_xty,
                                                   sub_x,
                                                   y.tilde,
                                                   D.star.inv,
                                                   nboot = nboot)
      }
      log.unnormalized.posterior.probability <- log.unnormalized.posterior.probability - max(log.unnormalized.posterior.probability)
      unnormalized.posterior.probability <- exp(log.unnormalized.posterior.probability)
      posterior.probability <- unnormalized.posterior.probability/sum(unnormalized.posterior.probability)

    } else {
      # Do model search with genetic algorithm
      fitness_ftn <- function(string){
        if(sum(string) == 0){
          return(total.p * log(pi0.hat) + log_marginal_likelihood_P3D_null(one.tilde,y.tilde,D.star.inv))
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
        residual <- y.tilde - sub_x %*% beta.reml
        sum.squares <- sum(D.star.inv * residual^2)
        sigma2.hat <- as.numeric(sum.squares/(n - length(model) - 1))
        return(k*log(1-pi0.hat) + (total.p-k)*log(pi0.hat) +
                 log.marginal.likelihood.nondiagonal_1sim(k, tau.hat, sigma2.hat,
                                                          sub_xtx,
                                                          sub_xty,
                                                          sub_x,
                                                          y.tilde,
                                                          D.star.inv,
                                                          nboot = nboot))
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
          residual <- y.tilde - sub_x %*% beta.reml
          sum.squares <- sum(D.star.inv * residual^2)
          sigma2.hat <- as.numeric(sum.squares/(n - length(model) - 1))
          k <- length(model)
          tmp_log.unnormalized.posterior.probability[i] <- k*log(1-pi0.hat) + (total.p-k)*log(pi0.hat) +
            log.marginal.likelihood.nondiagonal_1sim(k, tau.hat, sigma2.hat,
                                                     sub_xtx,
                                                     sub_xty,
                                                     sub_x,
                                                     y.tilde,
                                                     D.star.inv,
                                                     nboot = nboot)
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

    tmp <- stats::cor(X.tilde.significant) - diag(ncol(X.tilde.significant))
    tmp <- tmp[,model == 1,drop = FALSE]
    indices.correlated <- vector()
    for(i in 1:ncol(tmp)){
      indices.correlated <- c(indices.correlated,which(tmp[,i] == 1))
    }
    model_dat <- cbind(model_dat,0)
    colnames(model_dat) <- c("SNPs","BestModel","Inclusion_Prob","CorrelatedSNPs")
    model_dat$CorrelatedSNPs[indices.correlated] <- 1

    return(list(prescreen = return_dat,modelselection = model_dat,pi_0_hat = pi0.hat))

  }else{
    return(list(prescreen = return_dat,modelselection = "No significant in prescreen1",pi_0_hat = pi0.hat))
  }
}

all_Marco_nondiagonal_optim_tau_pi0_priors <- function(y,X,P,D,FDR.threshold = 0.95,maxiterations = 1000, runs_til_stop = 100,nboot = 100, v1=2, v2=2){
  # optimize tau and pi0 unconstrained
  # prior of pi0: beta(v1=2, v2=2)
  # fix pi0 in model selection

  n <- length(y)

  y.tilde <- t(P) %*% y
  X.tilde <- t(P) %*% X

  one.tilde <- t(P) %*% rep(1,length(y))

  REML.estimates <- REML_estimation_1Dregressor_BGWAS(x.tilde = one.tilde,y.tilde = y.tilde,d = D)
  alpha <- REML.estimates$beta
  kappa <- REML.estimates$kappa

  y.tilde <- y.tilde - one.tilde*alpha

  D.star.inv <- 1/(1 + kappa*D)

  xj.t.xj <- apply(X.tilde*(D.star.inv*X.tilde),2,sum)
  xj.t.y <- t(X.tilde) %*% (D.star.inv*y.tilde)
  beta.hat <- xj.t.y / xj.t.xj
  sigma.star2 <- as.numeric(((beta.hat^2)*xj.t.xj - 2*beta.hat*c(xj.t.y) + as.numeric(t(y.tilde)%*%(D.star.inv*y.tilde)))/(length(y) - 2))
  var.beta.hat <- sigma.star2 / xj.t.xj

  return_dat <- cbind(beta.hat,var.beta.hat)
  return_dat <- as.data.frame(return_dat)
  colnames(return_dat) <- c("Beta_Hat","Var_Beta_Hat")

  log.marginal.likelihood_tau = function(param) # function(y, m)
  {
    # Computes log marginal likelihood for 1D regressors as a function of tau and pi0.
    # This assumes that the prior for the regression coefficient is a mixture of a point of
    # mass at zero and a nonlocal prior. This also assumes that the estimates of the different
    # regression coefficients are (approximately) independent. Thus, this is kind of a
    # pseudo likelihood.
    # Maximization of this function leads to empirical Bayes estimates of tau and pi0.
    # This is much a faster (and maybe more precise) alternative to the Langaas' method.
    pi0 = exp(param[1]) / (1+exp(param[1]))
    tau = exp(param[2])
    return(sum(log(pi0*stats::dnorm(beta.hat,mean= 0,sd=sqrt(var.beta.hat)) +
                     (1-pi0)*(2*pi*var.beta.hat)^(-0.5) * (n*tau+1)^(-1.5) *
                     (1+n*tau*beta.hat^2/((n*tau+1)*var.beta.hat)) *
                     exp(-beta.hat^2 / (2*var.beta.hat*(n*tau+1)))
    ))+log(stats::dbeta(x = pi0, shape1 = v1, shape2 = v2)) + stats::dgamma(1/tau,shape = 0.55/0.022 + 1,rate = 0.55,log = TRUE))
  }

  param.start <- c(2.944439,0)
  result <- stats::optim(param.start, fn=log.marginal.likelihood_tau, method = "L-BFGS-B", hessian=TRUE, control = list(fnscale=-1))

  pi0.hat <- exp(result$par[1]) / (1+exp(result$par[1]))
  tau.hat <- exp(result$par[2])

  # Compute posterior probability of beta_j different than 0:
  numerator <- (1-pi0.hat)*(2*pi*var.beta.hat)^(-0.5) * (n*tau.hat+1)^(-1.5) *
    (1+n*tau.hat*beta.hat^2/((n*tau.hat+1)*var.beta.hat)) *
    exp(-beta.hat^2 / (2*var.beta.hat*(n*tau.hat+1)))
  denominator <- pi0.hat*stats::dnorm(beta.hat,mean= 0,sd=sqrt(var.beta.hat)) +
    (1-pi0.hat)*(2*pi*var.beta.hat)^(-0.5) * (n*tau.hat+1)^(-1.5) *
    (1+n*tau.hat*beta.hat^2/((n*tau.hat+1)*var.beta.hat)) *
    exp(-beta.hat^2 / (2*var.beta.hat*(n*tau.hat+1)))

  postprob <- numerator / denominator

  ###############  Bayesian FDR  #################

  order.postprob <- order(postprob, decreasing=TRUE)
  postprob.ordered <- postprob[order.postprob]

  FDR.Bayes <- cumsum(postprob.ordered) / 1:ncol(X)
  if(sum(FDR.Bayes > FDR.threshold) == 0){
    return_dat <- cbind(return_dat,postprob,FALSE)
    return_dat <- as.data.frame(return_dat)
    colnames(return_dat) <- c("Beta_Hat","Var_Beta_Hat","PostProb","Significant")
  }else{
    return_dat <- cbind(return_dat,postprob,postprob >= postprob.ordered[max(which(FDR.Bayes > FDR.threshold))])
    return_dat <- as.data.frame(return_dat)
    colnames(return_dat) <- c("Beta_Hat","Var_Beta_Hat","PostProb","Significant")
  }

  if(sum(return_dat$Significant) > 0){
    indices.significant <- which(return_dat$Significant)
    X.tilde.significant <- X.tilde[, indices.significant,drop = FALSE]

    if(Matrix::rankMatrix(X.tilde.significant)[1] < ncol(X.tilde.significant)){
      dropped_cols <- caret::findLinearCombos(X.tilde.significant)$remove
      REML.estimates.ms <- REML_estimation(X.tilde.significant[,-dropped_cols,drop = FALSE], y.tilde, D)
    }else{
      REML.estimates.ms <- REML_estimation(X.tilde.significant, y.tilde, D)
    }
    kappa.hat <- REML.estimates.ms$t

    D.star.inv <- 1/(1 + kappa.hat*D)
    X.significant.ty <- t(X.tilde.significant) %*% (D.star.inv*y.tilde)
    X.significant.t.X.significant <- t(X.tilde.significant) %*% (D.star.inv*X.tilde.significant)

    total.p <- ncol(X.tilde.significant)

    if(total.p < 16){
      # Do full model search
      total.models <- 2^total.p
      log.unnormalized.posterior.probability <- rep(NA, total.models)
      log.unnormalized.posterior.probability[1] <- total.p * log(pi0.hat) + log_marginal_likelihood_P3D_null(one.tilde,y.tilde,D.star.inv)
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
        residual <- y.tilde - sub_x %*% beta.reml
        sum.squares <- sum(D.star.inv * residual^2)
        sigma2.hat <- as.numeric(sum.squares/(n - length(model) - 1))
        log.unnormalized.posterior.probability[i+1] <-
          k*log(1-pi0.hat) + (total.p-k)*log(pi0.hat) +
          log.marginal.likelihood.nondiagonal_1sim(k, tau.hat, sigma2.hat,
                                                   sub_xtx,
                                                   sub_xty,
                                                   sub_x,
                                                   y.tilde,
                                                   D.star.inv,
                                                   nboot = nboot)
      }
      log.unnormalized.posterior.probability <- log.unnormalized.posterior.probability - max(log.unnormalized.posterior.probability)
      unnormalized.posterior.probability <- exp(log.unnormalized.posterior.probability)
      posterior.probability <- unnormalized.posterior.probability/sum(unnormalized.posterior.probability)

    } else {
      # Do model search with genetic algorithm
      fitness_ftn <- function(string){
        if(sum(string) == 0){
          return(total.p * log(pi0.hat) + log_marginal_likelihood_P3D_null(one.tilde,y.tilde,D.star.inv))
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
        residual <- y.tilde - sub_x %*% beta.reml
        sum.squares <- sum(D.star.inv * residual^2)
        sigma2.hat <- as.numeric(sum.squares/(n - length(model) - 1))
        return(k*log(1-pi0.hat) + (total.p-k)*log(pi0.hat) +
                 log.marginal.likelihood.nondiagonal_1sim(k, tau.hat, sigma2.hat,
                                                          sub_xtx,
                                                          sub_xty,
                                                          sub_x,
                                                          y.tilde,
                                                          D.star.inv,
                                                          nboot = nboot))
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
          residual <- y.tilde - sub_x %*% beta.reml
          sum.squares <- sum(D.star.inv * residual^2)
          sigma2.hat <- as.numeric(sum.squares/(n - length(model) - 1))
          k <- length(model)
          tmp_log.unnormalized.posterior.probability[i] <- k*log(1-pi0.hat) + (total.p-k)*log(pi0.hat) +
            log.marginal.likelihood.nondiagonal_1sim(k, tau.hat, sigma2.hat,
                                                     sub_xtx,
                                                     sub_xty,
                                                     sub_x,
                                                     y.tilde,
                                                     D.star.inv,
                                                     nboot = nboot)
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

    tmp <- stats::cor(X.tilde.significant) - diag(ncol(X.tilde.significant))
    tmp <- tmp[,model == 1,drop = FALSE]
    indices.correlated <- vector()
    for(i in 1:ncol(tmp)){
      indices.correlated <- c(indices.correlated,which(tmp[,i] == 1))
    }
    model_dat <- cbind(model_dat,0)
    colnames(model_dat) <- c("SNPs","BestModel","Inclusion_Prob","CorrelatedSNPs")
    model_dat$CorrelatedSNPs[indices.correlated] <- 1

    return(list(prescreen = return_dat,modelselection = model_dat,pi_0_hat = pi0.hat,tau = tau.hat))

  }else{
    return(list(prescreen = return_dat,modelselection = "No significant in prescreen1",pi_0_hat = pi0.hat,tau = tau.hat))
  }
}

#' Performs BGWAS analysis as in the BGWAS manuscript (Williams, J., Xu, S. and Ferreira, M.A., 2023. BGWAS: Bayesian variable selection in linear mixed models with nonlocal priors for genome-wide association studies. BMC Bioinformatics, 24(1), pp.1-20.).
#'
#'
#'
#' @param Y The observed numeric phenotypes
#' @param SNPs The SNP matrix, where each column represents a single SNP encoded as the numeric coding 0, 1, 2. This is entered as a matrix object.
#' @param FDR_Nominal The nominal false discovery rate for which SNPs are selected from in the screening step. Defaulted at 0.05.
#' @param kinship The observed kinship matrix, has to be a square positive semidefinite matrix. Defaulted as the identity matrix. The function used
#' to create the kinship matrix used in the BICOSS paper is A.mat() from package rrBLUP.
#' @param tau This specifies the prior or the fixed value for the nonlocal prior on the SNP effects. To specify the Inverse Gamma prior set tau = "IG". To specify fixed values of tau set tau = 0.022. Fixed values explored in the BBGWAS manuscript were 0.022 and 0.348.
#' @param maxiterations The maximum iterations the genetic algorithm in the model selection step iterates for.
#' Defaulted at 400 which is the value used in the BICOSS paper simulation studies.
#' @param runs_til_stop The number of iterations at the same best model before the genetic algorithm in the model selection step converges.
#' Defaulted at 40 which is the value used in the BICOSS paper simulation studies.
#' @return A named list that includes the output from the screening step and the output from the model selection step. The model proposed by the BGWAS method is the output of the model selection step.
#' @examples
#' library(GWAS.BAYES)
#' BGWAS(Y = Y, SNPs = SNPs, kinship = kinship,
#'     FDR_Nominal = 0.05, tau = "IG",
#'     maxiterations = 400,runs_til_stop = 40)
#' @export
BGWAS <- function(Y,SNPs,FDR_Nominal = 0.05,kinship = diag(nrow(SNPs)),tau = "IG",maxiterations = 4000,runs_til_stop = 400){

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
  if(!(tau == "IG" | is.numeric(tau))){
    stop("The only non-numeric input for tau is IG")
  }
  if(is.numeric(tau)){
    if(tau > 1 | tau < 0.01){
      stop("Recommended values of tau should be between 1 and 0.01. We specifically recommend tau = 0.022 or tau = 0.348")
    }
  }

  FDR.threshold <- 1 - FDR_Nominal
  aux <- eigen(kinship,symmetric = TRUE)
  P <- aux$vectors
  eigH <- aux$values
  D <- c(eigH[-length(eigH)],0)
  rm(aux);rm(eigH)

  X <- scale(SNPs)
  rm(SNPs)

  if(tau == "IG"){
    output <- all_Marco_nondiagonal_optim_tau_pi0_priors(y = Y,X = X,P = P,D = D,FDR.threshold = FDR.threshold,
                                                         maxiterations = maxiterations, runs_til_stop = runs_til_stop,
                                                         nboot = 1000, v1=1, v2=1)
    return(list(Screening = output$prescreen,ModelSelection = output$modelselection))
  }else{
    output <- all_Marco_nondiagonal_pi0_known_expected(y = Y,X = X,P = P,D = D,FDR.threshold = FDR.threshold
                                                       ,maxiterations = maxiterations, runs_til_stop = runs_til_stop,
                                                       nboot = 1000, TAU = tau, v1=1, v2=1)
    return(list(Screening = output$prescreen,ModelSelection = output$modelselection))
  }
}
