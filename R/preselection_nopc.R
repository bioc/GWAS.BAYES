preselection_nopc <- function(Y,X,number_cores,frequentist,controlrate,threshold,nullprob,alterprob,kinship = FALSE){
  # These lines need to be run, this essentially cleans the data set
  requireNamespace("parallel")
  ####################

  original_n <- ncol(X)
  Xy <- cbind(Y,X)
  y1 <- c("y",paste0("SNP",1:ncol(X)))
  colnames(Xy) <- y1
  y1 <- y1[y1 %in% paste0("SNP",1:ncol(X))]

  X <- as.matrix(Xy[,y1])
  Y <- matrix(Xy[,1],ncol = 1)
  rm(Xy)

  P <- 2

  ##################

  if(is.logical(kinship)){
    #Simple Linear Regression
    nX <- nrow(X)
    SNPdata_list <- split(t(X),1:ncol(X))
    rm(X)
    if(.Platform$OS.type == "unix"){
      SNPdata_list <- mclapply(SNPdata_list,SNP_data_function_nopcp,mc.cores = number_cores,int = 1)
    } else{
      cl <- makeCluster(number_cores)
      clusterExport(cl,c("SNPdata_list","SNP_data_function_nopcp"),envir=environment())
      SNPdata_list <- parLapply(cl,SNPdata_list,SNP_data_function_nopcp,int = 1)
      stopCluster(cl)
    }

    if(frequentist){
      #Frequentist SLR
      if(.Platform$OS.type == "unix"){
        SNPdata_list <- unlist(mclapply(SNPdata_list,optim_llik_SLR_p,mc.cores = number_cores, y = Y))
      } else{
        cl <- makePSOCKcluster(number_cores)
        clusterExport(cl,c("optimize","solve","eigenMapMatMult3","eigenMapMatMult2"),envir=environment())
        SNPdata_list <- unlist(parLapply(cl,SNPdata_list,optim_llik_SLR_p, y = Y))
        stopCluster(cl)
      }
      Pval_function(p_vals = SNPdata_list,n = original_n,thresh = threshold,control = controlrate)
    }else{
      #Bayesian SLR
      w <- matrix(1,nrow = nX,ncol = 1)
      BIC.null <- optim_llik_SLR_BIC(w,y = Y)

      if(.Platform$OS.type == "unix"){
        SNPdata_list <- unlist(mclapply(SNPdata_list,optim_llik_SLR_BIC,mc.cores = number_cores, y = Y))
      } else{
        cl <- makeCluster(number_cores)
        clusterExport(cl,c("solve","eigenMapMatMult3","eigenMapMatMult2"),envir=environment())
        SNPdata_list <- unlist(parLapply(cl,SNPdata_list,optim_llik_SLR_BIC, y = Y))
        stopCluster(cl)
      }
      p_vec <- (alterprob*exp(-.5 * (SNPdata_list - apply(cbind(SNPdata_list,as.numeric(BIC.null)),1,max))))/(alterprob*exp(-.5 * (SNPdata_list - apply(cbind(SNPdata_list,as.numeric(BIC.null)),1,max))) + nullprob*exp(-.5 * (as.numeric(BIC.null) - apply(cbind(SNPdata_list,as.numeric(BIC.null)),1,max))))
      order_vec <- (1:length(p_vec))[order(p_vec,decreasing = TRUE)]
      p_vec <- p_vec[order(p_vec,decreasing = TRUE)]
      FDR <- vector()
      for(d in 1:length(p_vec)){
        FDR[d] <- sum(1 - p_vec[1:d])/d
      }
      FDR <- FDR[order(order_vec,decreasing = FALSE)]
      p_vec <- p_vec[order(order_vec,decreasing = FALSE)]
      tf_mat <- as.data.frame(cbind(FDR < threshold,p_vec))
      colnames(tf_mat) <- c("Significant","ApprPosteriorProbs")
      return(tf_mat)
    }
  }else{
    #Kinship
    #SLR with kinship component
    nX <- nrow(X)

    spec.decomp <- eigen(kinship,symmetric = TRUE)
    Q <- spec.decomp$vectors
    Qt <- t(Q)
    D <- diag(spec.decomp$values,nrow = nX,ncol = nX)
    rm(spec.decomp)

    intercept <- matrix(1,nrow = nX,ncol = 1)
    intercept <- do.call(eigenMapMatMult2,list(Qt,intercept))
    Y <- do.call(eigenMapMatMult2,list(Qt,Y)); X <- do.call(eigenMapMatMult2,list(Qt,X))

    SNPdata_list <- split(t(X),1:ncol(X))
    rm(X);rm(Qt);rm(Q)

    if(.Platform$OS.type == "unix"){
      SNPdata_list <- mclapply(SNPdata_list,SNP_data_function_nopcp,mc.cores = number_cores,int = intercept)
    } else{
      cl <- makeCluster(number_cores)
      clusterExport(cl,c("SNP_data_function_nopcp"),envir=environment())
      SNPdata_list <- parLapply(cl,SNPdata_list,SNP_data_function_nopcp,int = intercept)
      stopCluster(cl)
    }
    #Frequentist RE model
    if(frequentist){
      if(.Platform$OS.type == "unix"){
        SNPdata_list <- unlist(mclapply(SNPdata_list,optim_llik_RE_p,mc.cores = number_cores, y = Y, d = D))
      } else{
        cl <- makeCluster(number_cores)
        clusterExport(cl,c("optimize","diag","solve","log_profile_likelihood_REML","optim_llik_RE_p","eigenMapMatMult2","eigenMapMatMult3"),envir=environment())
        SNPdata_list <- unlist(parLapply(cl,SNPdata_list,optim_llik_RE_p, y = Y, d = D))
        stopCluster(cl)
      }

      Pval_function(p_vals = SNPdata_list,n = original_n,thresh = threshold,control = controlrate)
    }else{
      #Bayesian RE model
      BIC.null <- optim_llik_RE_BIC(x = intercept,y = Y,d = D)

      if(.Platform$OS.type == "unix"){
        SNPdata_list <- unlist(mclapply(SNPdata_list,optim_llik_RE_BIC,mc.cores = number_cores, y = Y, d = D))
      } else{
        cl <- makeCluster(number_cores)
        clusterExport(cl,c("solve","log_profile_likelihood_REML","optimize","eigenMapMatMult2","eigenMapMatMult3"),envir=environment())
        SNPdata_list <- unlist(parLapply(cl,SNPdata_list,optim_llik_RE_BIC, y = Y, d = D))
        stopCluster(cl)
      }

      p_vec <- (alterprob*exp(-.5 * (SNPdata_list - apply(cbind(SNPdata_list,as.numeric(BIC.null)),1,max))))/(alterprob*exp(-.5 * (SNPdata_list - apply(cbind(SNPdata_list,as.numeric(BIC.null)),1,max))) + nullprob*exp(-.5 * (as.numeric(BIC.null) - apply(cbind(SNPdata_list,as.numeric(BIC.null)),1,max))))
      order_vec <- (1:length(p_vec))[order(p_vec,decreasing = TRUE)]
      p_vec <- p_vec[order(p_vec,decreasing = TRUE)]
      FDR <- vector()
      for(d in 1:length(p_vec)){
        FDR[d] <- sum(1 - p_vec[1:d])/d
      }
      FDR <- FDR[order(order_vec,decreasing = FALSE)]
      p_vec <- p_vec[order(order_vec,decreasing = FALSE)]
      tf_mat <- as.data.frame(cbind(FDR < threshold,p_vec))
      colnames(tf_mat) <- c("Significant","ApprPosteriorProbs")
      return(tf_mat)
    }

  }

}
