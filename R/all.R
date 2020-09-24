standardize <- function(SNPs,method=c("major-minor","alphabetical"),number_cores = 1){
  requireNamespace("parallel")
  match.arg(method)

  if(sum(apply(SNPs, 2, is.character)) == ncol(SNPs)){
    for(j in 1:ncol(SNPs)){
      SNPs[,j] <- as.factor(SNPs[,j])
    }
  }else if(sum(apply(SNPs,2,is.factor)) != ncol(SNPs)){
    stop("Not all columns in SNPs are characters or Factors")
  }

  if(!is.numeric(number_cores)){
    stop("number_cores is not numeric")
  }

  if(method == "major-minor"){
    std_fun <- function(x,snps){
      ret_1 <- integer(nrow(snps))
      ret_1[snps[,x] == names(which.max(table(snps[,x])))] <- 1
      return(ret_1)
    }
  }else{
    std_fun <- function(x,snps){
      ret_1 <- integer(nrow(snps))
      ret_1[snps[,x] == levels(snps[,x])[1]] <- 1
      return(ret_1)
    }
  }

  if(.Platform$OS.type == "unix"){
    z <- mclapply(1:ncol(SNPs),std_fun,mc.cores = number_cores,snps = SNPs)
  } else{
    cl <- makeCluster(number_cores)
    clusterExport(cl,c("std_fun"),envir=environment())
    z <- parLapply(cl,1:ncol(SNPs),std_fun,snps = SNPs)
    stopCluster(cl)
  }
  z <- do.call(cbind, z)
  return(z)
}

pca_function <- function(SNPs,number_components,plot_it = TRUE){
  if(sum(apply(SNPs,2,is.numeric)) != ncol(SNPs)){
    stop("Not every column of SNPs is numeric")
  }
  if(sum(apply(SNPs,2,sum) > nrow(SNPs))){
    stop("Some values of the SNP matrix are not 0 or 1")
  }
  if(!is.logical(plot_it)){
    stop("plot_it is not a logical value")
  }
  if(!is.numeric(number_components)){
    stop("number_components is not numeric")
  }
  SNPs <- as.matrix(SNPs)
  ec2 <- svd(SNPs,number_components,0)
  pca_1 <- as.matrix(ec2$u)
  colnames(pca_1) <- paste0("PC_",1:number_components)
  if(plot_it){
    plot(ec2$d^2/sum(ec2$d^2),xlim = c(0,15),ylab = "Percent Variation Explained",xlab = "Number of Components",type = "b", pch = 16)
    return(pca_1)
  }else{
    return(pca_1)
  }
}

level_function <- function(SNPs,MAF = 0.01){
  if(sum(apply(SNPs,2,is.numeric)) != ncol(SNPs)){
    stop("Not every column of SNPs is numeric")
  }
  if(sum(apply(SNPs,2,sum) >nrow(SNPs))){
    stop("Some values of the SNP matrix are not 0 or 1")
  }
  if(!is.numeric(MAF)){
    stop("MAF is not numeric")
  }
  z <- colMeans(SNPs)
  SNPs <- SNPs[,z <= (1 - MAF)]
  if(sum(!(z <= (1 - MAF))) == 0){
    level_dropped <- "None"
  }else{
    level_dropped <- which(!(z <= (1 - MAF)))
  }
  return(list(SNPs = SNPs,SNPs_Dropped = level_dropped))
}

SNP_data_function_pcp <- function(x,pcp,int){
  x <- as.data.frame(cbind(int,matrix(x,ncol = 1),pcp))
  colnames(x) <- c("Intercept","SNP",paste0("PC",1:ncol(pcp)))
  return(as.matrix(x))
}

SNP_data_function_nopcp <- function(x,int){
  x <- as.data.frame(cbind(int,matrix(x,ncol = 1)))
  colnames(x) <- c("Intercept","SNP")
  return(as.matrix(x))
}

log_profile_likelihood_REML <- function(x,t,y,d){
  n <- nrow(x)
  p <- ncol(x)
  estar <- diag(1/(1 + t*diag(d)),ncol = n,nrow = n)
  xtey <- eigenMapMatMult3(t(x),estar,y)
  xtex <- eigenMapMatMult3(t(x),estar,x)
  beta.mle <- eigenMapMatMult2(solve(xtex),xtey)
  top <- eigenMapMatMult3(t(y - eigenMapMatMult2(x,beta.mle)),estar,y - eigenMapMatMult2(x,beta.mle))
  sigma.mle <- as.numeric(top/(n - p))
  as.numeric(-.5*n*log(2*pi*sigma.mle)  - .5 * sum(log(diag(diag(1,ncol = n,nrow = n) + t*d)))- (1/(2*sigma.mle)) * top - .5*log(det(xtex)) + (p/2)*log(sigma.mle))
}

optim_llik_RE_p <- function(x,y,d){
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
  standard.error <- sigma.mle * xtex_1
  p_value <- 2* pt(abs(as.numeric(beta.mle[2,1]/sqrt(standard.error[2,2]))),n - p,lower.tail = FALSE)
  return(p_value)
}

optim_llik_RE_BIC <- function(x,y,d){
  #This is used for RE model and Kinship SLR model
  n <- ncol(d)
  z <- as.numeric(optimize(log_profile_likelihood_REML,interval = c(0,1000),x = x, y = y,d = d,maximum = TRUE)$objective)
  BIC.cor <- -2 * z + ncol(x)*log(n)
  return(BIC.cor)
}

optim_llik_SLR_BIC <- function(x,y){
  n <- nrow(x)
  p <- ncol(x)
  xtx <- eigenMapMatMult2(t(x),x)
  xty <- eigenMapMatMult2(t(x),y)
  xtx_1 <- solve(xtx)
  beta.mle <- eigenMapMatMult2(xtx_1,xty)
  txy <- y - eigenMapMatMult2(x,beta.mle)
  top <- eigenMapMatMult2(t(txy),txy)
  sigma.mle <- as.numeric(top/(n-p))
  BIC <- -2 *((-n/2)*log(2*pi*sigma.mle) - .5*(1/sigma.mle)*top - .5*log(det(xtx)) + (p/2)*log(sigma.mle)) + (p-1)*log(n)
  return(BIC)
}

optim_llik_SLR_p <- function(x,y){
  n <- nrow(x)
  p <- ncol(x)
  xtx <- eigenMapMatMult2(t(x),x)
  xty <- eigenMapMatMult2(t(x),y)
  xtx_1 <- solve(xtx)
  beta.mle <- eigenMapMatMult2(xtx_1,xty)
  txy <- y - eigenMapMatMult2(x,beta.mle)
  top <- eigenMapMatMult2(t(txy),txy)
  sigma.mle <- as.numeric(top/(n-p))
  standard.error <- sigma.mle * xtx_1
  p_value <- 2* pt(abs(as.numeric(beta.mle[2,1]/sqrt(standard.error[2,2]))),n - p,lower.tail = FALSE)
  return(p_value)
}

Pval_function <- function(p_vals,n,thresh,control){
  p_vals <- as.numeric(p_vals)
  tf_mat <- matrix(NA,ncol = 2,nrow = n)
  #Seeing if the SNP is significant given a multiple comparison error correction and the threshold
  tf_mat[,1] <- p.adjust(p_vals,control) < thresh
  tf_mat[,2] <- p_vals
  tf_mat <- as.data.frame(tf_mat)
  colnames(tf_mat) <- c("Significant","P_values")
  #Returning the p-values as well as if they are significant or not
  return(tf_mat)
}

aggregate_SNPs <- function(SNPs,Y,na.rm = TRUE){

  if(sum(apply(SNPs,2,is.numeric)) != ncol(SNPs)){
    stop("Not every column of SNPs is numeric")
  }
  if(sum(apply(SNPs,2,sum) > nrow(SNPs))){
    stop("Some values of the SNP matrix are not 0 or 1")
  }
  if(!(is.numeric(Y))){
    stop("Y is not a numeric vector")
  }
  if(!is.logical(na.rm)){
    stop("na.rm is not numeric")
  }
  original_n <- ncol(SNPs)
  Xy <- cbind(Y,SNPs)
  y1 <- c("y",paste0("SNP",1:ncol(SNPs)))
  colnames(Xy) <- y1
  y1 <- y1[y1 %in% paste0("SNP",1:ncol(SNPs))]

  if(na.rm){
    Xy <- Xy[is.finite(Xy[,1]),]
  }

  SNPs <- as.matrix(Xy[,y1])
  Y <- matrix(Xy[,1],ncol = 1)
  rm(Xy)

  G <- SNPs
  G[SNPs == 0] <- -1
  G <- round(eigenMapMatMult2(G,t(G)))
  q <- which(G == G[1,1])
  G[q] <- 1
  G[(1:prod(dim(G)))[!(1:prod(dim(G)) %in% q)]] <- 0
  G <- unique(G)
  G <- G / rowSums(G)
  X_small <- round(do.call(eigenMapMatMult2,list(G,SNPs)))
  Y_small <- do.call(eigenMapMatMult2,list(G,Y))
  return(list(SNPs = X_small,Y = Y_small))
}

cor_plot <- function(SNPs,significant,info = FALSE){
  if(sum(apply(SNPs,2,is.numeric)) != ncol(SNPs)){
    stop("Not every column of SNPs is numeric")
  }
  if(sum(apply(SNPs,2,sum) >nrow(SNPs))){
    stop("Some values of the SNP matrix are not 0 or 1")
  }
  if(!(is.numeric(significant))){
    stop("significant is not a numeric vector")
  }
  if(!is.logical(info)){
    if(!((nrow(info) == 2) & (ncol(info) == ncol(SNPs)))){
      stop("info does not have the dimensions 2 x ncol(SNPs)")
    }
  }else{
    if(info == TRUE){
      stop("info can only take on the values of FALSE or a matrix 2 x ncol(SNPs)")
    }
  }
  requireNamespace("reshape2")
  requireNamespace("ggplot2")

  colnames(SNPs) <- paste0("SNP",1:ncol(SNPs))
  SNPs <- SNPs[,significant == 1]
  save1 <- colnames(SNPs)
  if(ncol(SNPs) < 1){
    print("No significant SNPs")
  } else if(ncol(SNPs) == 1){
    print("One significant SNP no plot generated")
  }else if(2 <= ncol(SNPs) & ncol(SNPs) < 12){
    if(!is.logical(info)){
      info <- info[,significant == 1]
      info <- paste0(info[1,],"-",info[2,])
      colnames(SNPs) <- info
    }else{
      colnames(SNPs) <- names(which(significant == 1))
    }
    if(is.null(colnames(SNPs))){
      colnames(SNPs) <- save1
    }
    melted_cormat <- reshape2::melt(round(cor(SNPs),2))
    melted_cormat$Position1 <- as.character(melted_cormat$Var1)
    melted_cormat$Position2 <- as.character(melted_cormat$Var2)

    ggheatmap <- ggplot2::ggplot(melted_cormat, ggplot2::aes(Position1, Position2, fill = value))+ggplot2::geom_tile(color = "white")+ggplot2::scale_fill_gradient2(low = "blue", high = "red", mid = "white",midpoint = 0, limit = c(-1,1), space = "Lab", name="Pearson\nCorrelation")+ggplot2::xlab("Position")+ggplot2::ylab("Position") + ggplot2::theme_minimal() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, size = 12, hjust = 1)) + ggplot2::coord_fixed()
    ggheatmap + ggplot2::geom_text(ggplot2::aes(Position2, Position1, label = value), color = "black", size = 4)
  } else if(11 < ncol(SNPs) & ncol(SNPs) < 25){
    if(!is.logical(info)){
      info <- info[,significant == 1]
      info <- paste0(info[1,],"-",info[2,])
      colnames(SNPs) <- info
    }else{
      colnames(SNPs) <- names(which(significant == 1))
    }
    if(is.null(colnames(SNPs))){
      colnames(SNPs) <- save1
    }
    melted_cormat <- reshape2::melt(round(cor(SNPs),2))
    melted_cormat$Position1 <- as.character(melted_cormat$Var1)
    melted_cormat$Position2 <- as.character(melted_cormat$Var2)

    ggheatmap <- ggplot2::ggplot(melted_cormat, ggplot2::aes(Position1, Position2, fill = value)) + ggplot2::geom_tile(color = "white") + ggplot2::scale_fill_gradient2(low = "blue", high = "red", mid = "white",midpoint = 0, limit = c(-1,1), space = "Lab", name="Pearson\nCorrelation") + ggplot2::xlab("Position")+ggplot2::ylab("Position") + ggplot2::theme_minimal() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, size = 12, hjust = 1)) + ggplot2::coord_fixed()
    ggheatmap
  }else{
    print("Too many SNPs to plot")
  }
}
