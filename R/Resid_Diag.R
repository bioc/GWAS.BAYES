resids_diag <- function(Y,SNPs,significant,kinship = FALSE,principal_components = FALSE,plot_it = TRUE){
  requireNamespace("Matrix")
  requireNamespace("caret")

  if(sum(apply(SNPs,2,is.numeric)) != ncol(SNPs)){
    stop("Not every column of SNPs is numeric")
  }
  if(sum(apply(SNPs,2,sum) > nrow(SNPs))){
    stop("Some values of the SNP matrix are not 0 or 1")
  }
  if(!(is.numeric(Y))){
    stop("Y is not a numeric vector")
  }
  if(!(is.numeric(significant))){
    stop("significant is not a numeric vector")
  }
  if(!is.logical(kinship)){
    if(!((nrow(kinship) == nrow(SNPs)) & (ncol(kinship) == nrow(SNPs)))){
      stop("kinship does not have the dimensions nrow(SNPs) x nrow(SNPs)")
    }
  }else{
    if(kinship == TRUE){
      stop("kinship can only take on the values of FALSE or a matrix nrow(SNPs) x nrow(SNPs)")
    }
  }
  if(!is.logical(principal_components)){
    if(!(nrow(principal_components) == nrow(SNPs))){
      stop("principal_components does not have the the same number of rows as SNPs")
    }
    if(sum(apply(principal_components,2,is.numeric)) != ncol(principal_components)){
      stop("Not every column of principal_components is numeric")
    }
  }else{
    if(principal_components == TRUE){
      stop("principal_components can only take on the values of FALSE or a matrix with the same number of rows as SNPs")
    }
  }
  if(!is.logical(plot_it)){
    stop("plot_it is not a logical value")
  }

  Good <- SNPs[,which(significant == 1)]
  oG <- Good

  nX <- nrow(Good)

  oY <- Y

  optim_llik_RE_BETA <- function(x,y,d){
    #This is used for RE model and Kinship SLR model
    n <- nrow(x)
    p <- ncol(x)
    z <- optimize(log_profile_likelihood_REML,interval = c(0,100),x = x, y = y,d = d,maximum = TRUE)$maximum
    estar <- diag(1/(1 + z*diag(d)),ncol = n,nrow = n)
    xtey <- eigenMapMatMult3(t(x),estar,y)
    xtex <- eigenMapMatMult3(t(x),estar,x)
    xtex_1 <- solve(xtex)
    beta.mle <- eigenMapMatMult2(xtex_1,xtey)
    top <- eigenMapMatMult3(t(y - eigenMapMatMult2(x,beta.mle)),estar,y - eigenMapMatMult2(x,beta.mle))
    sigma.mle <- as.numeric(top/(n - p))
    return(list(beta.mle,z))
  }

  optim_llik_SLR_BETA <- function(x,y){
    n <- nrow(x)
    p <- ncol(x)
    xtx <- eigenMapMatMult2(t(x),x)
    xty <- eigenMapMatMult2(t(x),y)
    xtx_1 <- solve(xtx)
    beta.mle <- eigenMapMatMult2(xtx_1,xty)
    return(beta.mle)
  }

  if(is.logical(kinship)){
    if(is.logical(principal_components)){
      #No principal or kinship
      intercept <- matrix(1,nrow = nX,ncol = 1)
      Good <- cbind(intercept,Good)

      if(Matrix::rankMatrix(Good)[1] < ncol(Good)){
        dropped_cols <- caret::findLinearCombos(Good)$remove
        Good <- Good[,-dropped_cols]
      }
      beta.good <- optim_llik_SLR_BETA(x = Good, y = Y)
      yhat <- Good %*% beta.good
      resids <- Y - yhat

      if(plot_it == TRUE){
        hist(resids,main = "Histogram of Residuals",xlab = "Residuals")
        return(shapiro.test(resids))
      }else{
        return(shapiro.test(resids))
      }
    }else{
      #Principal and no kinship
      intercept <- matrix(1,nrow = nX,ncol = 1)
      Good <- cbind(intercept,Good,principal_components)

      if(Matrix::rankMatrix(Good)[1] < ncol(Good)){
        dropped_cols <- caret::findLinearCombos(Good)$remove
        Good <- Good[,-dropped_cols]
      }
      beta.good <- optim_llik_SLR_BETA(x = Good, y = Y)
      yhat <- Good %*% beta.good
      resids <- Y - yhat

      if(plot_it == TRUE){
        hist(resids,main = "Histogram of Residuals",xlab = "Residuals")
        return(shapiro.test(resids))
      }else{
        return(shapiro.test(resids))
      }
    }
  }else{
    if(is.logical(principal_components)){
      #Kinship no principal
      spec.decomp <- eigen(kinship,symmetric = TRUE)
      Q <- spec.decomp$vectors
      Qt <- t(Q)
      D <- diag(spec.decomp$values,nrow = nX,ncol = nX)
      rm(spec.decomp)

      intercept <- matrix(1,nrow = nX,ncol = 1)
      intercept <- do.call(eigenMapMatMult2,list(Qt,intercept))
      Y <- do.call(eigenMapMatMult2,list(Qt,Y)); Good <- do.call(eigenMapMatMult2,list(Qt,Good))

      Good <- cbind(intercept,Good)
      oG <- cbind(1,oG)

      if(Matrix::rankMatrix(Good)[1] < ncol(Good)){
        dropped_cols <- caret::findLinearCombos(Good)$remove
        Good <- Good[,-dropped_cols]
        oG <- oG[,-dropped_cols]
      }

      beta.good <- optim_llik_RE_BETA(x = Good, y = Y, d = D)

      Q <- Q[,-ncol(Q)]
      Qt <- t(Q)
      D <- D[-ncol(D),-ncol(D)]
      uhat <- eigenMapMatMult3(eigenMapMatMult2(Q,1/(diag(1,nrow = nrow(D),ncol = ncol(D)) + (1/beta.good[[2]])*(1/D))),Qt,oY - oG %*% beta.good[[1]])
      resids <- oY - oG %*% beta.good[[1]] - uhat

      if(plot_it == TRUE){
        hist(resids,main = "Histogram of Residuals",xlab = "Residuals")
        return(shapiro.test(resids))
      }else{
        return(shapiro.test(resids))
      }
    }else{
      #Kinship and principal
      ogprinc <- principal_components

      spec.decomp <- eigen(kinship,symmetric = TRUE)
      Q <- spec.decomp$vectors
      Qt <- t(Q)
      D <- diag(spec.decomp$values,nrow = nX,ncol = nX)
      rm(spec.decomp)

      intercept <- matrix(1,nrow = nX,ncol = 1)
      intercept <- do.call(eigenMapMatMult2,list(Qt,intercept))
      Y <- do.call(eigenMapMatMult2,list(Qt,Y)); Good <- do.call(eigenMapMatMult2,list(Qt,Good))
      principal_components <- do.call(eigenMapMatMult2,list(Qt,principal_components))

      Good <- cbind(intercept,Good,principal_components)
      oG <- cbind(1,oG,ogprinc)

      if(Matrix::rankMatrix(Good)[1] < ncol(Good)){
        dropped_cols <- caret::findLinearCombos(Good)$remove
        Good <- Good[,-dropped_cols]
        oG <- oG[,-dropped_cols]
      }

      beta.good <- optim_llik_RE_BETA(x = Good, y = Y, d = D)

      Q <- Q[,-ncol(Q)]
      Qt <- t(Q)
      D <- D[-ncol(D),-ncol(D)]
      uhat <- eigenMapMatMult3(eigenMapMatMult2(Q,1/(diag(1,nrow = nrow(D),ncol = ncol(D)) + (1/beta.good[[2]])*(1/D))),Qt,oY - oG %*% beta.good[[1]])
      resids <- oY - oG %*% beta.good[[1]] - uhat

      if(plot_it == TRUE){
        hist(resids,main = "Histogram of Residuals",xlab = "Residuals")
        return(shapiro.test(resids))
      }else{
        return(shapiro.test(resids))
      }
    }
  }
}
