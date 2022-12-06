svd1 <- function(x){
  p <- ncol(x)
  d <- svd(x,nu = p,nv = p)$d
  return(list(max(d)/min(d),p - length(d)))
}

svd1_modified <- function(x,Xc,one.tilde){
  x <- cbind(one.tilde,x,Xc)
  p <- ncol(x)
  d <- svd(x,nu = p,nv = p)$d
  return(list(max(d)/min(d),p - length(d)))
}

svd1_modified_noXc <- function(x,one.tilde){
  x <- cbind(one.tilde,x)
  p <- ncol(x)
  d <- svd(x,nu = p,nv = p)$d
  return(list(max(d)/min(d),p - length(d)))
}

log_profile_likelihood_MLE <- function(x,t,y,d){
  n <- nrow(x)
  p <- ncol(x)
  estar <- 1/(1 + t*d)
  xtey <- t(x)%*%(estar*y)
  xtex <- t(x)%*%(estar*x)
  aux <- eigen(xtex,symmetric = TRUE)
  P <- aux$vectors
  eigH <- aux$values
  eigH[abs(eigH) < 10^(-10)] <- 0
  eigH <- 1/eigH
  eigH[!is.finite(eigH)] <- 0
  beta.mle <- P%*%(eigH*t(P)) %*% xtey
  sum.squares <- sum(estar * (y - x %*% beta.mle)^2)
  sigma.mle <- as.numeric(sum.squares/(n))
  as.numeric(-.5*(n)*log(2*pi*sigma.mle)  + .5 * sum(log(estar))- n/2)
}

log_profile_likelihood_MLE_1Dregressor <- function(x,t,y,d){
  n <- length(y)
  estar <- 1/(1 + t*d)
  xty <- sum(estar * x * y)
  xtx <- sum(estar * x^2)
  beta.mle <- xty / xtx
  sum.squares <- sum(estar * (y - x * beta.mle)^2)
  sigma.mle <- as.numeric(sum.squares/(n))
  as.numeric(-.5*(n)*log(2*pi*sigma.mle)  + .5 * sum(log(estar)) - n/2)
}

RE_BIC_1Dregressor <- function(x,y,d){
  n <- length(y)
  p <- ncol(x)
  z <- stats::optimize(log_profile_likelihood_MLE_1Dregressor,interval = c(0,100),x = x, y = y, d = d, maximum = TRUE)$objective
  return(-2 * z + p*log(n))
}

RE_BIC_P3D_modelselection <- function(xtx,xty,x,y,d_1,pstar){
  n <- length(d_1)
  aux <- eigen(xtx,symmetric = TRUE)
  P <- aux$vectors
  eigH <- aux$values
  eigH[abs(eigH) < 10^(-10)] <- 0
  eigH <- 1/eigH
  eigH[!is.finite(eigH)] <- 0
  beta.mle <- P%*%(eigH*t(P)) %*% xty
  sum.squares <- sum(d_1 * (y - x %*% beta.mle)^2)
  sigma.mle <- as.numeric(sum.squares/(n))
  z <- as.numeric(-.5*(n)*log(2*pi*sigma.mle)  + .5 * sum(log(d_1[d_1 != 0]))- n/2)
  BIC.cor <- -2 * z + pstar*log(n)
  return(BIC.cor)
}

RE_BIC_modelselection <- function(x,y,d,pstar){
  n <- length(d)
  z <- as.numeric(stats::optimize(log_profile_likelihood_MLE,interval = c(0,100),x = x, y = y,d = d,maximum = TRUE)$objective)
  BIC.cor <- -2*z + pstar*log(n)
  return(BIC.cor)
}

log_profile_restricted_likelihood_NullModel <- function(x,t,y,d){
  n <- length(y)
  p <- 1
  estar <- 1/(1 + t*d)
  xty <- sum(estar * x * y)
  xtx <- sum(estar * x^2)
  beta.reml <- xty / xtx
  residual <- y - x * beta.reml
  sum.squares <- sum(estar * residual^2)
  sigma.reml <- as.numeric(sum.squares/(n - p))
  as.numeric(-.5*(n - p)*log(2*pi*sigma.reml)  + .5 * sum(log(estar)) - (n-p)/2 - .5*log(xtx))
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

RE_p_modified <- function(x,y,d,Xc,one.tilde){
  x <- cbind(one.tilde,x,Xc)
  P <- ncol(Xc)
  n <- nrow(x)
  p <- ncol(x)
  L <- matrix(diag(p)[-c(1,(p - P + 1):p),],nrow = p - length(c(1,(p - P + 1):p)),ncol = p)
  z <- stats::optimize(log_profile_restricted_likelihood,interval = c(0,100),x = x, y = y,d = d,maximum = TRUE)$maximum
  estar <- 1/(1 + z*d)
  aux <- eigen(t(x)%*%(estar*x),symmetric = TRUE)
  P <- aux$vectors
  eigH <- aux$values
  eigH[abs(eigH) < 10^(-10)] <- 0
  eigH <- 1/eigH
  eigH[!is.finite(eigH)] <- 0
  beta2 <- P%*%(eigH*t(P)) %*% t(x) %*% (estar*y)
  residual <- y - x %*% beta2
  sum.squares <- sum(estar * residual^2)
  p <- sum(eigH > 0)
  sigma.reml <- as.numeric(sum.squares/(n - p))
  Fstat <- ((1/sigma.reml)*as.numeric(t(L%*%beta2)%*%solve(L%*%P%*%(eigH*t(P))%*%t(L))%*%L%*%beta2))/(nrow(L))
  return(p_value = stats::pf(Fstat,nrow(L),n - p,lower.tail = FALSE))
}

RE_p_modified_noXc <- function(x,y,d,one.tilde){
  x <- cbind(one.tilde,x)
  n <- nrow(x)
  p <- ncol(x)
  L <- matrix(diag(p)[-1,],nrow = p - 1,ncol = p)
  z <- stats::optimize(log_profile_restricted_likelihood,interval = c(0,100),x = x, y = y,d = d,maximum = TRUE)$maximum
  estar <- 1/(1 + z*d)
  aux <- eigen(t(x)%*%(estar*x),symmetric = TRUE)
  P <- aux$vectors
  eigH <- aux$values
  eigH[abs(eigH) < 10^(-10)] <- 0
  eigH <- 1/eigH
  eigH[!is.finite(eigH)] <- 0
  beta2 <- P%*%(eigH*t(P)) %*% t(x) %*% (estar*y)
  residual <- y - x %*% beta2
  sum.squares <- sum(estar * residual^2)
  p <- sum(eigH > 0)
  sigma.reml <- as.numeric(sum.squares/(n - p))
  Fstat <- ((1/sigma.reml)*as.numeric(t(L%*%beta2)%*%solve(L%*%P%*%(eigH*t(P))%*%t(L))%*%L%*%beta2))/(nrow(L))
  return(stats::pf(Fstat,nrow(L),n - p,lower.tail = FALSE))
}

RE_BIC_modified <- function(x,y,d,one.tilde, Xc){
  x <- cbind(one.tilde,x,Xc)
  n <- length(y)
  p <- ncol(x)
  z <- stats::optimize(log_profile_likelihood_MLE,interval = c(0,100),x = x, y = y, d = d, maximum = TRUE)$objective
  return(-2 * z + p*log(n))
}

RE_BIC_modified_noXc <- function(x,y,d,one.tilde){
  x <- cbind(one.tilde,x)
  n <- length(y)
  p <- ncol(x)
  z <- stats::optimize(log_profile_likelihood_MLE,interval = c(0,100),x = x, y = y, d = d, maximum = TRUE)$objective
  return(-2 * z + p*log(n))
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

REML_NullModel <- function(x,y,d){
  n <- length(y)
  p <- 1
  t <- stats::optimize(log_profile_restricted_likelihood_NullModel,interval = c(0,100),x = x, y = y, d = d, maximum = TRUE)$maximum
  estar <- 1/(1 + t*d)
  xty <- sum(estar * x * y)
  xtx <- sum(estar * x^2)
  beta.reml <- xty / xtx
  residual <- y - x * beta.reml
  sum.squares <- sum(estar * residual^2)
  sigma.reml <- as.numeric(sum.squares/(n - p))

  return(list(beta=beta.reml, sigma2=sigma.reml, t=t))
}

MLE_estimation <- function(x,y,d){
  n <- length(y)
  t <- stats::optimize(log_profile_likelihood_MLE, interval = c(0,100), x = x, y = y, d = d, maximum = TRUE)$maximum
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
  sigma.mle <- as.numeric(sum.squares/(n))

  return(list(beta=beta.reml, sigma2=sigma.mle, t=t))
}

log_mle_P3D_null <- function(one.tilde,y,estar){
  n <- length(y)     # sample size
  p <- 1
  xty <- sum(estar * one.tilde * y)
  xtx <- sum(estar * one.tilde^2)
  beta.reml <- xty / xtx
  residual <- y - one.tilde * beta.reml
  sum.squares <- sum(estar * residual^2)
  sigma2 <- sum.squares/n
  as.numeric(-.5*(n)*log(2*pi*sigma2)  + .5 * sum(log(estar[estar != 0])) - n/2)
}

BICOSS_1iter <- function(y,X,P,D,FDR.threshold = 0.95,maxiterations = 1000, runs_til_stop = 100,P3D,indices_previous,iter_num,initial_pi0,Xc,kinship){
  if(P3D){
    y.tilde <- t(P)%*%y
    X.tilde <- t(P)%*%X
    one <- matrix(1,ncol = 1,nrow = length(y.tilde))
    one.tilde <- t(P)%*%one
    y_MS <- y - one%*%MLE_estimation(one.tilde,y.tilde, D)$beta
    D_MS <- D

    n <- length(y.tilde)

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
      REML_estimates <- REML_estimation(x = Xc.tilde[,-dropped_cols,drop = FALSE],y = y.tilde, D)
      MLE_estimates <- MLE_estimation(x = Xc.tilde[,-dropped_cols,drop = FALSE],y = y.tilde,d = D)

      t_mle <- MLE_estimates$t
      alpha_MLE <- MLE_estimates$beta
      y_SMA_MLE <- y - Xc[,-dropped_cols,drop = FALSE]%*%alpha_MLE

      t_reml <- REML_estimates$t
      alpha_REML <- REML_estimates$beta
      y_SMA_REML <- y - Xc[,-dropped_cols,drop = FALSE]%*%alpha_REML

      V <- diag(n) + t_reml*kinship
      aux_REML <- eigen(V - Xc[,-dropped_cols,drop = FALSE] %*% solve(t(Xc[,-dropped_cols,drop = FALSE])%*%solve(V)%*%Xc[,-dropped_cols,drop = FALSE])%*%t(Xc[,-dropped_cols,drop = FALSE]),symmetric = TRUE)
      p_fixed <- ncol(Xc[,-dropped_cols,drop = FALSE])

      V <- diag(n) + t_mle*kinship
      aux_MLE <- eigen(V - Xc[,-dropped_cols,drop = FALSE] %*% solve(t(Xc[,-dropped_cols,drop = FALSE])%*%solve(V)%*%Xc[,-dropped_cols,drop = FALSE])%*%t(Xc[,-dropped_cols,drop = FALSE]),symmetric = TRUE)
    }else{
      REML_estimates <- REML_estimation(x = Xc.tilde,y = y.tilde, D)
      MLE_estimates <- MLE_estimation(x = Xc.tilde,y = y.tilde,d = D)

      t_mle <- MLE_estimates$t
      alpha_MLE <- MLE_estimates$beta
      y_SMA_MLE <- y - Xc%*%alpha_MLE

      t_reml <- REML_estimates$t
      alpha_REML <- REML_estimates$beta
      y_SMA_REML <- y - Xc%*%alpha_REML

      V <- diag(n) + t_reml*kinship
      aux_REML <- eigen(V - Xc %*% solve(t(Xc)%*%solve(V)%*%Xc)%*%t(Xc),symmetric = TRUE)

      V <- diag(n) + t_mle*kinship
      aux_MLE <- eigen(V - Xc %*% solve(t(Xc)%*%solve(V)%*%Xc)%*%t(Xc),symmetric = TRUE)
    }

    P <- aux_REML$vectors
    eigH <- aux_REML$values
    D <- rep(0,length(eigH))
    D[1:(n - p_fixed)] <- eigH[1:(n - p_fixed)]
    rm(aux_REML);rm(eigH)

    y.tilde_SMA_REML <- t(P)%*%y_SMA_REML
    X.tilde <- t(P)%*%X
    # Xc.tilde <- t(P)%*%Xc

    D.star.inv <- rep(0,length(D))
    D.star.inv[1:(n - p_fixed)] <- (1/D)[1:(n - p_fixed)]
    xj.t.xj <- apply(X.tilde*(D.star.inv*X.tilde),2,sum)
    xj.t.y <- t(X.tilde) %*% (D.star.inv*y.tilde_SMA_REML)
    beta.hat <- xj.t.y / xj.t.xj
    sigma.reml <- as.numeric(((beta.hat^2)*xj.t.xj - 2*beta.hat*c(xj.t.y) + as.numeric(t(y.tilde_SMA_REML)%*%(D.star.inv*y.tilde_SMA_REML)))/(n - p_fixed - 1))
    #var.beta.hat <- sigma.reml*((xj.t.xj)^(-2))*(xj.t.xj - apply(D.star.inv*X.tilde,2,function(x){t(x)%*%Xc.tilde%*%solve(t(Xc.tilde)%*%(D.star.inv*Xc.tilde))%*%t(Xc.tilde)%*%x}))
    var.beta.hat <- sigma.reml*(1/xj.t.xj)

    indices_bad <- abs(beta.hat) < 10e-12

    beta.hat<- beta.hat[!indices_bad]
    var.beta.hat <- var.beta.hat[!indices_bad]
    t.statistic <- beta.hat / sqrt(var.beta.hat)
    pvalues <- vector(mode = "numeric",length = ncol(X.tilde))
    pvalues[!indices_bad] <- 2*stats::pt(abs(t.statistic),n - p_fixed - 1,lower.tail = FALSE)
    pvalues[indices_bad] <- 1
    pvalues[is.na(pvalues)] <- 1

    nullprob <- limma::convest(pvalues[!indices_bad])
    if(nullprob == 1 & iter_num == 1){
      nullprob <- 1 - 100/ncol(X)
    }
    if(iter_num == 1){
      pi0 <- nullprob
    }else{
      pi0 <- initial_pi0
    }
    if(nullprob == 1 & iter_num != 1){
      return(list(best_model = "No Significant SNPs in Screening",pi0 = pi0))
    }

    P <- aux_MLE$vectors
    eigH <- aux_MLE$values
    D <- rep(0,length(eigH))
    D[1:(n - p_fixed)] <- eigH[1:(n - p_fixed)]
    rm(aux_MLE);rm(eigH)

    y.tilde_SMA_MLE <- t(P)%*%y_SMA_MLE
    X.tilde <- t(P)%*%X
    # Xc.tilde <- t(P)%*%Xc

    D.star.inv <- rep(0,length(D))
    D.star.inv[1:(n - p_fixed)] <- (1/D)[1:(n - p_fixed)]
    xj.t.xj <- apply(X.tilde*(D.star.inv*X.tilde),2,sum)
    xj.t.y <- t(X.tilde) %*% (D.star.inv*y.tilde_SMA_MLE)
    beta.hat <- xj.t.y / xj.t.xj
    sigma.mle <- as.numeric(((beta.hat^2)*xj.t.xj - 2*beta.hat*c(xj.t.y) + as.numeric(t(y.tilde_SMA_MLE)%*%(D.star.inv*y.tilde_SMA_MLE)))/n)
    BIC <- vector(mode = "numeric",length = ncol(X.tilde))
    BIC[!indices_bad] <- (-2*as.numeric(-.5*(n)*log(2*pi*sigma.mle[!indices_bad])  + .5 * sum(log(D.star.inv[D.star.inv != 0]))- n/2)) + (p_fixed + 1)*log(n) - 2*1*log(1 - nullprob)
    BIC.null <- -2*(log_mle_P3D_null(t(P)%*%one,y = y.tilde_SMA_MLE,estar = D.star.inv)) + p_fixed*log(n) - 2*log(nullprob)
    BIC[indices_bad] <- BIC.null
    BIC[is.na(BIC)] <- BIC.null

    postprob <- (exp(-.5 * (BIC - apply(cbind(BIC,as.numeric(BIC.null)),1,max))))/(exp(-.5 * (BIC - apply(cbind(BIC,as.numeric(BIC.null)),1,max))) + exp(-.5 * (as.numeric(BIC.null) - apply(cbind(BIC,as.numeric(BIC.null)),1,max))))
    order.postprob <- order(postprob, decreasing=TRUE)
    postprob.ordered <- postprob[order.postprob]
    FDR.Bayes <- cumsum(postprob.ordered) / 1:ncol(X)

    if(sum(FDR.Bayes > FDR.threshold) == 0){
      return_dat <- cbind(postprob,FALSE)
      return_dat <- as.data.frame(return_dat)
      colnames(return_dat) <- c("PostProb","Significant")
    }else{
      return_dat <- cbind(postprob,postprob >= postprob.ordered[max(which(FDR.Bayes > FDR.threshold))])
      return_dat <- as.data.frame(return_dat)
      colnames(return_dat) <- c("PostProb","Significant")
    }

    if(sum(return_dat$Significant) > 0){
      indices.significant <- which(return_dat$Significant == 1)
      X.tilde.significant <- X.tilde[,indices.significant,drop = FALSE]

      if(!is.null(Xc)){
        X.tilde.significant <- cbind(Xc.tilde[,-1],X.tilde.significant)
        indices.significant <- c(indices_previous,indices.significant)
      }

      if(Matrix::rankMatrix(cbind(one.tilde,X.tilde.significant))[1] < ncol(cbind(one.tilde,X.tilde.significant))){
        dropped_cols <- caret::findLinearCombos(cbind(one.tilde,X.tilde.significant))$remove
        REML.estimates.ms <- REML_estimation(cbind(one.tilde,X.tilde.significant)[,-dropped_cols,drop = FALSE], y.tilde, D_MS)
      }else{
        REML.estimates.ms <- REML_estimation(cbind(one.tilde,X.tilde.significant), y.tilde, D_MS)
      }
      kappa.hat <- REML.estimates.ms$t

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
        BICs <- rep(NA, total.models)
        BICs[1] <- -2*(log_mle_P3D_null(one.tilde,y = y.tilde_MS,estar = D.star.inv)) + 1*log(n) - 2*(log(1 - pi0) + (total.p)*log(pi0))
        dat <- rep(list(0:1), total.p)
        dat <- as.matrix(expand.grid(dat))
        for (i in 1:(total.models-1)){
          model <- unname(which(dat[i + 1,] == 1))
          Xsub <- X.tilde.significant[,model,drop = FALSE]
          k <- length(model) + 1
          if(Matrix::rankMatrix(Xsub)[1] < ncol(Xsub)){
            dropped_cols <- caret::findLinearCombos(Xsub)$remove
            model <- model[-dropped_cols]
          }
          sub_xtx <- matrix(X.significant.t.X.significant[model, model],ncol = length(model),nrow = length(model))
          sub_xty <- matrix(X.significant.ty[model,],ncol = 1)
          sub_x <- matrix(X.tilde.significant[, model],ncol = length(model))
          BICs[i+1] <- RE_BIC_P3D_modelselection(xtx = sub_xtx,xty = sub_xty,x = sub_x,y = y.tilde_MS,d_1 = D.star.inv, pstar = k) - 2*(k*log(1 - pi0) + (total.p - k)*log(pi0))
        }

        best_model <- indices.significant[unname(which(dat[which.min(BICs),] == 1))]

      } else {
        # Do model search with genetic algorithm
        fitness_ftn <- function(string){
          if(sum(string) == 0){
            return((-1)*(-2*(log_mle_P3D_null(one.tilde,y = y.tilde_MS,estar = D.star.inv)) + 1*log(n) - 2*(log(1 - pi0) + (total.p)*log(pi0))))
          }
          model <- which(string==1)
          Xsub <- X.tilde.significant[,model,drop = FALSE]
          k <- length(model) + 1
          if(Matrix::rankMatrix(Xsub)[1] < ncol(Xsub)){
            dropped_cols <- caret::findLinearCombos(Xsub)$remove
            model <- model[-dropped_cols]
          }
          sub_xtx <- matrix(X.significant.t.X.significant[model, model],ncol = length(model),nrow = length(model))
          sub_xty <- matrix(X.significant.ty[model,],ncol = 1)
          sub_x <- matrix(X.tilde.significant[, model],ncol = length(model))
          return((-1)*(RE_BIC_P3D_modelselection(xtx = sub_xtx,xty = sub_xty,x = sub_x,y = y.tilde_MS,d_1 = D.star.inv, pstar = k) - 2*(k*log(1 - pi0) + (total.p - k)*log(pi0))))
        }
        if(total.p > 99){
          suggestedsol <- diag(total.p)
          BICs <- vector()
          for(i in 1:total.p){
            model <- which(suggestedsol[i,]==1)
            sub_xtx <- matrix(X.significant.t.X.significant[model, model],ncol = length(model),nrow = length(model))
            sub_xty <- matrix(X.significant.ty[model,],ncol = 1)
            sub_x <- matrix(X.tilde.significant[, model],ncol = length(model))
            k <- length(model) + 1
            BICs[i] <- RE_BIC_P3D_modelselection(xtx = sub_xtx,xty = sub_xty,x = sub_x,y = y.tilde_MS,d_1 = D.star.inv, pstar = k) - 2*(k*log(1 - pi0) + (total.p - k)*log(pi0))
          }
          suggestedsol <- rbind(0,suggestedsol[order(BICs,decreasing = FALSE)[1:99],])
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
        if(length(which.max(ans@fitness))> 1){
          tmp <- dat[which.max(ans@fitness),]
          tmp <- tmp[which.min(rowSums(tmp)),,drop = FALSE]
          best_model <- indices.significant[unname(which(tmp == 1))]
        }else{
          best_model <- indices.significant[unname(which(dat[which.max(ans@fitness),] == 1))]
        }
      }
      return(list(best_model = best_model,pi0 = pi0))
    }else{
      return(list(best_model = "No Significant SNPs in Screening",pi0 = pi0))
    }

  }else{
    y.tilde <- t(P)%*%y
    X.tilde <- t(P)%*%X
    one <- matrix(1,ncol = 1,nrow = length(y.tilde))
    one.tilde <- t(P)%*%one
    n <- length(y.tilde)

    if(is.null(Xc)){
      BIC.null <- RE_BIC_1Dregressor(one.tilde,y = y.tilde,d = D)
      l1 <- unlist(apply(X.tilde,2, svd1_modified_noXc,one.tilde = one.tilde))
      rank_defient <- unname(which(l1[1:length(l1) %% 2 == 0]>0))
      condition_num <- unname(which(l1[1:length(l1) %% 2 == 1]>1e6))
      indices_previous <- unique(c(rank_defient,condition_num))
    }else{
      Xc.tilde <- t(P)%*%Xc
      BIC.null <- RE_BIC_modelselection(cbind(one.tilde,Xc.tilde),y = y.tilde,d = D,pstar = ncol(cbind(one.tilde,Xc.tilde)))
      l1 <- unlist(apply(X.tilde,2, svd1_modified,one.tilde = one.tilde,Xc = Xc.tilde))
      rank_defient <- unname(which(l1[1:length(l1) %% 2 == 0]>0))
      condition_num <- unname(which(l1[1:length(l1) %% 2 == 1]>1e6))
      indices_previous <- unique(c(indices_previous,rank_defient,condition_num))
    }

    if(is.null(indices_previous) | length(indices_previous) == 0){
      pvalues <- apply(X.tilde,2,RE_p_modified_noXc,y = y.tilde, d = D,one.tilde = one.tilde)
      nullprob <- limma::convest(pvalues)
    }else{
      pvalues <- vector(length = ncol(X.tilde))
      pvalues[-indices_previous] <- apply(X.tilde[,-indices_previous],2,RE_p_modified,y = y.tilde,Xc = Xc.tilde, d = D,one.tilde = one.tilde)
      pvalues[indices_previous] <- 1
      nullprob <- limma::convest(pvalues[-indices_previous])
    }

    if(nullprob == 1 & iter_num == 1){
      nullprob <- 1 - 100/ncol(X)
    }
    if(iter_num == 1){
      pi0 <- nullprob
    }else{
      pi0 <- initial_pi0
    }
    if(nullprob == 1 & iter_num != 1){
      return(list(best_model = "No Significant SNPs in Screening",pi0 = pi0))
    }

    if(is.null(indices_previous) | length(indices_previous) == 0){
      BIC <- apply(X.tilde,2,RE_BIC_modified_noXc,y = y.tilde, d = D,one.tilde = one.tilde) - 2*1*log(1 - nullprob)
      BIC.null <- BIC.null - 2*log(nullprob)
    }else{
      BIC <- vector(length = ncol(X.tilde))
      BIC[-indices_previous] <- apply(X.tilde[,-indices_previous],2,RE_BIC_modified,y = y.tilde, d = D,one.tilde = one.tilde,Xc = Xc.tilde) - 2*1*log(1 - nullprob)
      BIC.null <- BIC.null - 2*log(nullprob)
      BIC[indices_previous] <- BIC.null
    }

    postprob <- (exp(-.5 * (BIC - apply(cbind(BIC,as.numeric(BIC.null)),1,max))))/(exp(-.5 * (BIC - apply(cbind(BIC,as.numeric(BIC.null)),1,max))) + exp(-.5 * (as.numeric(BIC.null) - apply(cbind(BIC,as.numeric(BIC.null)),1,max))))
    order.postprob <- order(postprob, decreasing=TRUE)
    postprob.ordered <- postprob[order.postprob]
    FDR.Bayes <- cumsum(postprob.ordered) / 1:ncol(X)

    if(sum(FDR.Bayes > FDR.threshold) == 0){
      return_dat <- cbind(postprob,FALSE)
      return_dat <- as.data.frame(return_dat)
      colnames(return_dat) <- c("PostProb","Significant")
    }else{
      return_dat <- cbind(postprob,postprob >= postprob.ordered[max(which(FDR.Bayes > FDR.threshold))])
      return_dat <- as.data.frame(return_dat)
      colnames(return_dat) <- c("PostProb","Significant")
    }

    if(sum(return_dat$Significant) > 0){
      indices.significant <- which(return_dat$Significant == 1)
      X.tilde.significant <- X.tilde[,indices.significant,drop = FALSE]

      if(!is.null(Xc)){
        X.tilde.significant <- cbind(Xc.tilde,X.tilde.significant)
        indices.significant <- c(indices_previous,indices.significant)
      }

      total.p <- ncol(X.tilde.significant)

      if(total.p < 16){
        # Do full model search
        total.models <- 2^total.p
        BICs <- rep(NA, total.models)
        BICs[1] <- RE_BIC_modelselection(x = one.tilde,y = y.tilde,pstar = 1,d = D) - 2*(log(1 - pi0) + (total.p - 1)*log(pi0))
        dat <- rep(list(0:1), total.p)
        dat <- as.matrix(expand.grid(dat))
        for (i in 1:(total.models-1)){
          model <- unname(which(dat[i + 1,] == 1))
          Xsub <- cbind(one.tilde,X.tilde.significant[,model,drop = FALSE])
          k <- length(model) + 1
          if(Matrix::rankMatrix(Xsub)[1] < ncol(Xsub)){
            dropped_cols <- caret::findLinearCombos(Xsub)$remove
            Xsub <- Xsub[,-dropped_cols,drop = FALSE]
          }
          BICs[i+1] <- RE_BIC_modelselection(x = Xsub,y = y.tilde,d = D, pstar = k) - 2*(k*log(1 - pi0) + (total.p - k)*log(pi0))
        }

        best_model <- indices.significant[unname(which(dat[which.min(BICs),] == 1))]

      } else {
        # Do model search with genetic algorithm
        fitness_ftn <- function(string){
          if(sum(string) == 0){
            return((-1)*(RE_BIC_modelselection(x = one.tilde,y = y.tilde,pstar = 1,d = D) - 2*(log(1 - pi0) + (total.p - 1)*log(pi0))))
          }
          model <- which(string==1)
          Xsub <- cbind(one.tilde,X.tilde.significant[,model,drop = FALSE])
          k <- length(model) + 1
          if(Matrix::rankMatrix(Xsub)[1] < ncol(Xsub)){
            dropped_cols <- caret::findLinearCombos(Xsub)$remove
            Xsub <- Xsub[,-dropped_cols,drop = FALSE]
          }
          return((-1)*(RE_BIC_modelselection(x = Xsub,y = y.tilde,d = D, pstar = k) - 2*(k*log(1 - pi0) + (total.p - k)*log(pi0))))
        }
        if(total.p > 99){
          suggestedsol <- diag(total.p)
          BICs <- vector()
          for(i in 1:total.p){
            model <- which(suggestedsol[i,]==1)
            Xsub <- cbind(one.tilde,X.tilde.significant[,model,drop = FALSE])
            k <- length(model) + 1
            BICs[i] <- RE_BIC_modelselection(x = Xsub,y = y.tilde,d = D, pstar = k) - 2*(k*log(1 - pi0) + (total.p - k)*log(pi0))
          }
          suggestedsol <- rbind(0,suggestedsol[order(BICs,decreasing = FALSE)[1:99],])
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
        if(length(which.max(ans@fitness))> 1){
          tmp <- dat[which.max(ans@fitness),]
          tmp <- tmp[which.min(rowSums(tmp)),,drop = FALSE]
          best_model <- indices.significant[unname(which(tmp == 1))]
        }else{
          best_model <- indices.significant[unname(which(dat[which.max(ans@fitness),] == 1))]
        }
      }
      return(list(best_model = best_model,pi0 = pi0))
    }else{
      return(list(best_model = "No Significant SNPs in Screening",pi0 = pi0))
    }
  }
}



#' @title  BICOSS for Gaussian Phenotypes
#' @description  Performs BICOSS analysis as described in Williams, J., Ferreira, M.A.R. & Ji, T. BICOSS: Bayesian iterative conditional stochastic search for GWAS. BMC Bioinformatics 23, 475 (2022). https://doi.org/10.1186/s12859-022-05030-0.
#' @param Y The observed numeric phenotypes
#' @param SNPs The SNP matrix, where each column represents a single SNP encoded as the numeric coding 0, 1, 2. This is entered as a matrix object.
#' @param FDR_Nominal The nominal false discovery rate for which SNPs are selected from in the screening step.
#' @param kinship The observed kinship matrix, has to be a square positive semidefinite matrix. Defaulted as the identity matrix. The function used
#' to create the kinship matrix used in the BICOSS paper is A.mat() from package rrBLUP.
#' @param maxiterations The maximum iterations the genetic algorithm in the model selection step iterates for.
#' Defaulted at 400 which is the value used in the BICOSS paper simulation studies.
#' @param runs_til_stop The number of iterations at the same best model before the genetic algorithm in the model selection step converges.
#' Defaulted at 40 which is the value used in the BICOSS paper simulation studies.
#' @param P3D Population previous determined, if TRUE BICOSS uses approximated variance parameters estimated from the baseline model when conducting
#' both the screening and the model selection steps. Setting P3D = TRUE is significantly faster. If FALSE, uses exact estimates of the variance
#' parameters all models in both the screening and model selection step.
#' @return The column indices of SNPs that were in the best model identified by BICOSS.
#' @examples
#' library(GWAS.BAYES)
#' BICOSS(Y = Y, SNPs = SNPs, kinship = kinship,
#'     FDR_Nominal = 0.05,P3D = TRUE,
#'     maxiterations = 400,runs_til_stop = 40)
#' @export
BICOSS <- function(Y,SNPs,FDR_Nominal = 0.05,kinship = diag(nrow(SNPs)),maxiterations = 400,runs_til_stop = 40,P3D = TRUE){

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
  if(!is.logical(P3D)){
    stop("P3D has to be logical")
  }


  threshold <- FDR_Nominal
  aux <- eigen(kinship,symmetric = TRUE)
  P <- aux$vectors
  eigH <- aux$values
  D <- c(eigH[-length(eigH)],0)
  rm(aux);rm(eigH)

  if(P3D){
    SNPs <- scale(SNPs)
  }

  count <- 1
  p <- 1
  indices_previous <- NULL
  initial_pi0 <- NULL
  Xc <- NULL
  indices_sig_list <- list()

  while(p > 0 & count < 10){
    tmp <- BICOSS_1iter(y = Y,X = SNPs,P = P,D = D,FDR.threshold = 1 - threshold,maxiterations = maxiterations, runs_til_stop = runs_til_stop,
                        P3D = P3D,indices_previous = indices_previous,iter_num = count,initial_pi0 = initial_pi0,Xc = Xc,kinship = kinship)
    indices_sig_list[[count]] <- tmp$best_model
    indices_previous <- tmp$best_model

    if(is.character(tmp$best_model)){
      break
    }
    Xc <-  SNPs[,tmp$best_model,drop = FALSE]
    if(count == 1){
      p <- length(tmp$best_model)
      initial_pi0 <- tmp$pi0
    }else{
      p <- length(tmp$best_model[!(tmp$best_model%in%indices_sig_list[[count - 1]])])
    }

    count <- count + 1
  }

  if(count == 1 & is.character(tmp$best_model)){
    return(list(best_model = "No SNPs found in Screening"))
  }else{
    return(list(best_model = indices_sig_list[[count - 1]]))
  }
}

#' A. Thaliana Kinship matrix
#'
#' This is a kinship matrix from the TAIR9 genotype information for 328 A. Thaliana Ecotypes from the paper
#' Components of Root Architecture Remodeling in Response to Salt Stress. The kinship matrix was computed using all SNPs with minor allele frequency
#' greater than 0.01.
#'
#' @format ## `kinship`
#' A matrix with 328 rows and 328 columns corresponding to the 328 ecotypes.
"kinship"

#' A. Thaliana Genotype matrix
#'
#' This is a matrix with 328 observations and 9,000 SNPs. Each row is contains 9,000 SNPs from a single A. Thaliana ecotype in the paper
#' Components of Root Architecture Remodeling in Response to Salt Stress.
#'
#' @format ## `SNPs`
#' A matrix with 328 observations and 9,000 SNPs.
"SNPs"

#' A. Thaliana Simulated Phenotype matrix
#'
#' This is a phenotype matrix simulated from the 9,000 SNPs. SNPs at positions 450, 1350, 2250, 3150, 4050,
#' 4950, 5850, 6750, 7650, and 8550 have nonzero coefficients. Further, the data was simulated under the linear mixed model
#' specified in the vignette and the BICOSS manuscript using the kinship matrix (kinship).
#'
#' @format ## `Y`
#' A matrix with 328 rows corresponding to the 328 ecotypes.
"Y"

P3D_pvalues <- function(y,X,kinship){
  X <- scale(X)

  aux <- eigen(kinship,symmetric = TRUE)
  P <- aux$vectors
  eigH <- aux$values
  D <- c(eigH[-length(eigH)],0)
  rm(aux);rm(eigH)
  n <- length(y)

  y.tilde <- t(P) %*% y
  X.tilde <- t(P) %*% X

  one.tilde <- t(P) %*% rep(1,length(y))

  REML.estimates <- REML_NullModel(x = one.tilde,y = y.tilde,d = D)
  alpha <- REML.estimates$beta
  kappa <- REML.estimates$t

  y.tilde <- y.tilde - one.tilde*alpha

  D.star.inv <- 1/(1 + kappa*D)

  xj.t.xj <- apply(X.tilde*(D.star.inv*X.tilde),2,sum)
  xj.t.y <- t(X.tilde) %*% (D.star.inv*y.tilde)
  beta.hat <- xj.t.y / xj.t.xj
  sigma.star2 <- as.numeric(((beta.hat^2)*xj.t.xj - 2*beta.hat*c(xj.t.y) + as.numeric(t(y.tilde)%*%(D.star.inv*y.tilde)))/(length(y) - 2))
  var.beta.hat <- sigma.star2 / xj.t.xj

  t.statistic <- beta.hat / sqrt(var.beta.hat)
  pvalues <- 2*stats::pnorm(abs(t.statistic), mean=0, sd=1, lower.tail=FALSE)
  return(pvalues)
}

NotP3D_pvalues <- function(y,X,kinship){
  X <- scale(X)

  aux <- eigen(kinship,symmetric = TRUE)
  P <- aux$vectors
  eigH <- aux$values
  D <- c(eigH[-length(eigH)],0)
  rm(aux);rm(eigH)
  n <- length(y)

  y.tilde <- t(P) %*% y
  X.tilde <- t(P) %*% X

  one.tilde <- t(P) %*% rep(1,length(y))

  REML.estimates <- REML_NullModel(x = one.tilde,y = y.tilde,d = D)
  alpha <- REML.estimates$beta

  y.tilde <- y.tilde - one.tilde*alpha

  REML_estimation_1Dregressor_p_value <- function(x.tilde,y.tilde,d){
    # This function computes REML estimates for the case
    # when the regressor is one-dimensional.
    # This is for GWAS structure with kinship random effects.
    # Arguments to the function:
    #    y.tilde: dependent variable in the spectral domain.
    #    x.tilde: regressor in the spectral domain.
    #    d: vector with the eigenvalues of the kinship matrix.
    n <- length(y.tilde)     # sample size
    p <- 2
    log_profile_restricted_likelihood <- function(x,t,y,d){
      n <- length(y)
      p <- 2
      estar <- 1/(1 + t*d)
      xty <- sum(estar * x * y)
      xtx <- sum(estar * x^2)
      beta.reml <- xty / xtx
      residual <- y - x * beta.reml
      sum.squares <- sum(estar * residual^2)
      sigma.reml <- as.numeric(sum.squares/(n - p))
      as.numeric(-.5*(n - p)*log(2*pi*sigma.reml)  + .5 * sum(log(estar)) - (n-p)/2 - .5*log(xtx))
    }
    kappa <- stats::optimize(log_profile_restricted_likelihood,interval = c(0,100),x = x.tilde, y = y.tilde, d = d, maximum = TRUE)$maximum
    D.star.inv <- 1/(1 + kappa*d)
    xty <- sum(D.star.inv * x.tilde * y.tilde)
    xtx <- sum(D.star.inv * x.tilde^2)
    beta.reml <- xty / xtx
    residual <- y.tilde - x.tilde * beta.reml
    sum.squares <- sum(D.star.inv * residual^2)
    sigma.reml <- as.numeric(sum.squares/(n - p))

    var.beta.hat <- sigma.reml / xtx
    t.statistic <- beta.reml / sqrt(var.beta.hat)
    pvalues <- 2*stats::pnorm(abs(t.statistic), mean=0, sd=1, lower.tail=FALSE)

    return(pvalues)
  }

  p_values <- apply(X.tilde,2,REML_estimation_1Dregressor_p_value,y.tilde = y.tilde,d = D)
  return(p_values)
}

Pvalues_SLR <- function(y,X){
  n <- length(y)
  p_value_oneSNP <- function(x,y){
    x <- cbind(1,x)
    n <- length(y)
    xtx_1 <- solve(t(x)%*%x)
    beta.hat <- xtx_1%*%t(x)%*%y
    sigma2 <- sum((y - x%*%beta.hat)^2)/(n - 2)
    t <- beta.hat[2]/sqrt(sigma2*xtx_1[2,2])
    return(2*stats::pnorm(abs(t), mean=0, sd=1, lower.tail=FALSE))
  }
  pvalues <- unname(apply(X,2,p_value_oneSNP,y = y))
  return(pvalues)
}

#' Performs Single Marker Association tests for both Linear Mixed Models and Linear models.
#'
#'
#'
#' @param Y The observed numeric phenotypes
#' @param SNPs The SNP matrix, where each column represents a single SNP encoded as the numeric coding 0, 1, 2. This is entered as a matrix object.
#' @param kinship The observed kinship matrix, has to be a square positive semidefinite matrix. Defaulted as the identity matrix. The function used
#' to create the kinship matrix used in the BICOSS paper is A.mat() from package rrBLUP.
#' @param P3D Population previous determined, if TRUE BICOSS uses approximated variance parameters estimated from the baseline model when conducting
#' both the screening and the model selection steps. Setting P3D = TRUE is significantly faster. If FALSE, uses exact estimates of the variance
#' parameters all models in both the screening and model selection step.
#' @return The p-values corresponding to every column provided in SNPs. These p-values can be used with any threshold of your choosing or with
#' p.adjust().
SMA <- function(Y,SNPs,kinship = FALSE,P3D = FALSE){
  if(!is.numeric(Y)){
    stop("Y has to be numeric")
  }
  if(!is.matrix(SNPs)){
    stop("SNPs has to be a matrix object")
  }
  if(!is.numeric(SNPs)){
    stop("SNPs has to contain numeric values")
  }
  if(!is.logical(kinship)){
    if(nrow(kinship) != ncol(kinship)){
      stop("kinship has to be a square matrix")
    }
    if(nrow(kinship) != nrow(SNPs)){
      stop("kinship has to have the same number of rows as SNPs")
    }
    if(!is.numeric(kinship)){
      stop("kinship has to contain numeric values")
    }
  }

  if(is.logical(kinship)){
    return(Pvalues_SLR(y = Y, X = SNPs))
  }else{
    if(P3D){
      return(unname(P3D_pvalues(y = Y,X = SNPs,kinship = kinship)))
    }else{
      return(unname(NotP3D_pvalues(y = Y,X = SNPs,kinship = kinship)))
    }
  }
}
