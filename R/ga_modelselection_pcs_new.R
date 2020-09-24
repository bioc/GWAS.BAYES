ga_modelselection_pcs_new <- function(Y,X,regions,regionsnames,significant,number_cores,principal_components,maxiterations,runs_til_stop,kinship = FALSE){

  requireNamespace("GA")
  requireNamespace("parallel")
  requireNamespace("doParallel")
  requireNamespace("Matrix")
  requireNamespace("caret")
  requireNamespace("memoise")

  original_n <- ncol(X)
  Xy <- cbind(Y,X,regions,principal_components)
  y1 <- c("y",paste0("SNP",1:ncol(X)),paste0("Regions",1:ncol(regions)),paste0("PC",1:ncol(principal_components)))
  colnames(Xy) <- y1
  y1 <- y1[y1 %in% c(paste0("SNP",1:ncol(X)),paste0("Region",1:ncol(regions)))]
  X <- as.matrix(Xy[,y1])
  Y <- matrix(Xy[,1],ncol = 1)
  principal_components <- as.matrix(Xy[,colnames(Xy)%in%paste0("PC",1:ncol(principal_components))])
  regions <- as.matrix(Xy[,colnames(Xy)%in%paste0("Regions",1:ncol(Regions))])
  rm(Xy)

  P <- 2 + ncol(principal_components)

  if(is.logical(kinship)){
    #SLR
    nX <- nrow(X)
    X <- X[,significant == 1]
    X <- cbind(X,regions)

    fitness_sub <- function(string) {
      inc <- which(string == 1)
      SNPdata_list_sub <- cbind(1,X[,inc],principal_components)

      if(Matrix::rankMatrix(SNPdata_list_sub)[1] < ncol(SNPdata_list_sub)){
        dropped_cols <- caret::findLinearCombos(SNPdata_list_sub)$remove
        SNPdata_list_sub <- SNPdata_list_sub[,-dropped_cols]
        print(paste0("Dropped Columns: ",dropped_cols))
      }
      return((-1)*optim_llik_SLR_BIC(SNPdata_list_sub,Y))
    }
    fitness_sub <- memoise::memoise(fitness_sub)
    ans <- GA::ga("binary", fitness = fitness_sub, nBits = ncol(X),names = y1[c(significant == 1,rep(TRUE,ncol(regions)))],popSize = 100,elitism = 10,maxiter = maxiterations,parallel = number_cores,run = runs_til_stop)
    memoise::forget(fitness_sub)
    dat <- cbind(ans@population,(-1)*ans@fitness)
    dat <- unique(dat)
    dat <- cbind(dat,(exp(-.5 * (dat[,ncol(dat)] - max(dat[,ncol(dat)]))))/(sum(exp(-.5 * (dat[,ncol(dat)] - max(dat[,ncol(dat)]))))))
    if(nrow(ans@solution) <= 10){
      dat <- dat[order(dat[,ncol(dat)],decreasing = TRUE),]
      Solution <- dat[1:10,-c(ncol(dat)-1,ncol(dat))]
      colnames(Solution) <- colnames(ans@solution)
    }else{
      Solution <- ans@solution
      colnames(Solution) <- colnames(ans@solution)
    }
    vec1 <- matrix(0,nrow = nrow(dat),ncol = nrow(Solution))
    for(i in 1:nrow(dat)){
      for(j in 1:nrow(Solution)){
        vec1[i,j] <- sum(dat[i,1:(ncol(dat) - 2)] == Solution[j,])
      }
    }
    vec <- vector()
    for(i in 1:nrow(Solution)){
      vec[i] <- paste0("Model ",i,": y ~ ",paste(names(which(Solution[i,] == 1)),collapse = " + ")," Posterior Prob. of ",round(dat[which(vec1[,i] == (ncol(dat) - 2)),ncol(dat)],digits = 4))
    }
    return(list(Models = vec,Solution = Solution,Regions = regionsnames))
  }else{
    #SLR w Kinship
    nX <- nrow(X)
    X <- X[,significant == 1]
    X <- cbind(X,regions)

    spec.decomp <- eigen(kinship,symmetric = TRUE)
    Q <- spec.decomp$vectors
    Qt <- t(Q)
    D <- diag(spec.decomp$values,nrow = nX,ncol = nX)
    rm(spec.decomp)

    intercept <- matrix(1,nrow = nX,ncol = 1)
    intercept <- do.call(eigenMapMatMult2,list(Qt,intercept))
    Y <- do.call(eigenMapMatMult2,list(Qt,Y)); X <- do.call(eigenMapMatMult2,list(Qt,X)); principal_components <- do.call(eigenMapMatMult2,list(Qt,principal_components))
    rm(Qt);rm(Q)

    fitness_sub <- function(string){
      X_sub <- cbind(intercept,X[,string == 1],principal_components)
      if(Matrix::rankMatrix(X_sub)[1] < ncol(X_sub)){
        dropped_cols <- caret::findLinearCombos(X_sub)$remove
        X_sub <- X_sub[,-dropped_cols]
        print(paste0("Dropped Columns: ",dropped_cols))
      }
      return((-1)*optim_llik_RE_BIC(X_sub,Y,D))
    }
    fitness_sub <- memoise::memoise(fitness_sub)
    ans <- GA::ga("binary", fitness = fitness_sub, nBits = ncol(X),names = y1[c(significant == 1,rep(TRUE,ncol(regions)))],popSize = 100,elitism = 10,maxiter = maxiterations,parallel = number_cores,run = runs_til_stop)
    memoise::forget(fitness_sub)
    dat <- cbind(ans@population,(-1)*ans@fitness)
    dat <- unique(dat)
    dat <- cbind(dat,(exp(-.5 * (dat[,ncol(dat)] - max(dat[,ncol(dat)]))))/(sum(exp(-.5 * (dat[,ncol(dat)] - max(dat[,ncol(dat)]))))))
    if(nrow(ans@solution) <= 10){
      dat <- dat[order(dat[,ncol(dat)],decreasing = TRUE),]
      Solution <- dat[1:10,-c(ncol(dat)-1,ncol(dat))]
      colnames(Solution) <- colnames(ans@solution)
    }else{
      Solution <- ans@solution
      colnames(Solution) <- colnames(ans@solution)
    }
    vec1 <- matrix(0,nrow = nrow(dat),ncol = nrow(Solution))
    for(i in 1:nrow(dat)){
      for(j in 1:nrow(Solution)){
        vec1[i,j] <- sum(dat[i,1:(ncol(dat) - 2)] == Solution[j,])
      }
    }
    vec <- vector()
    for(i in 1:nrow(Solution)){
      vec[i] <- paste0("Model ",i,": y ~ ",paste(names(which(Solution[i,] == 1)),collapse = " + ")," Posterior Prob. of ",round(dat[which(vec1[,i] == (ncol(dat) - 2)),ncol(dat)],digits = 4))
    }
    return(list(Models = vec,Solution = Solution,Regions = regionsnames))
  }
}
