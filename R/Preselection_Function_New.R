preselection <- function(Y,SNPs,number_cores = 1,principal_components = FALSE,frequentist = TRUE,controlrate = "bonferroni",threshold = 0.05,nullprob = NULL,alterprob = NULL,kinship = FALSE,info = FALSE){

  if(sum(apply(SNPs,2,is.numeric)) != ncol(SNPs)){
    stop("Not every column of SNPs is numeric")
  }
  if(sum(apply(SNPs,2,sum) > nrow(SNPs))){
    stop("Some values of the SNP matrix are not 0 or 1")
  }
  if(!(is.numeric(Y))){
    stop("Y is not a numeric vector")
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
  if(!is.logical(info)){
    if(!((nrow(info) == 2) & (ncol(info) == ncol(SNPs)))){
      stop("info does not have the dimensions 2 x ncol(SNPs)")
    }
  }else{
    if(info == TRUE){
      stop("info can only take on the values of FALSE or a matrix 2 x ncol(SNPs)")
    }
  }
  if(!is.numeric(threshold)){
    stop("threshold is not numeric")
  }
  if(threshold > 1 | threshold < 0){
    stop("threshold needs to be between 0 and 1")
  }
  if(!is.logical(frequentist)){
    stop("frequentist needs to be a logical value")
  }
  if(!is.numeric(number_cores)){
    stop("number_cores is not numeric")
  }
  if(!(controlrate %in% p.adjust.methods)){
    stop("control rate needs to be one of p.adjust.methods")
  }
  if(!is.null(nullprob)){
    if(!is.numeric(nullprob)){
      stop("nullprob needs to be numeric")
    }
    if(!is.numeric(alterprob)){
      stop("alterprob needs to be numeric")
    }
    if(alterprob > 1 | alterprob < 0){
      stop("nullprob needs to be between 0 and 1")
    }
    if(alterprob > 1 | alterprob < 0){
      stop("alterprob needs to be between 0 and 1")
    }
    if((nullprob + alterprob) != 1){
      stop("nullprob and alterprob need to sum to 1")
    }
  }

  if(is.logical(info)){
    if(is.logical(principal_components)){
      preselection_nopc(Y = Y,X = SNPs,number_cores = number_cores,frequentist = frequentist,controlrate = controlrate,threshold = threshold,nullprob = nullprob,alterprob = alterprob,kinship = kinship)
    }else{
      preselection_pc(Y = Y,X = SNPs,number_cores = number_cores,principal_components = principal_components,frequentist = frequentist,controlrate = controlrate,threshold = threshold,nullprob = nullprob,alterprob = alterprob,kinship = kinship)
    }
  }else{
    if(is.logical(principal_components)){
      results <- preselection_nopc(Y = Y,X = SNPs,number_cores = number_cores,frequentist = frequentist,controlrate = controlrate,threshold = threshold,nullprob = nullprob,alterprob = alterprob,kinship = kinship)
    }else{
      results <- preselection_pc(Y = Y,X = SNPs,number_cores = number_cores,principal_components = principal_components,frequentist = frequentist,controlrate = controlrate,threshold = threshold,nullprob = nullprob,alterprob = alterprob,kinship = kinship)
    }
    results <- cbind(t(info),results)
    results[,1] <- as.numeric(as.character(results[,1]))
    results[,2] <- as.numeric(as.character(results[,2]))
    colnames(results) <- c("Chromosomes","Positions","Significant","P_values")
    return(results)
  }
}
