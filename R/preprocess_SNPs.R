preprocess_SNPs <- function(SNPs,Y,MAF = 0.01,number_cores,na.rm = FALSE){

  if(!is.logical(na.rm)){
    stop("na.rm is not logical")
  }

  if(na.rm){
    SNPs <- SNPs[!is.na(Y),]
    Y <- Y[!is.na(Y)]
  }

  SNPs <- standardize(SNPs = SNPs,method = "major-minor",number_cores = number_cores)
  list1 <- aggregate_SNPs(SNPs = SNPs, Y = Y)
  SNPs <- list1$SNPs
  Y <- list1$Y
  level_list <- level_function(SNPs = SNPs,MAF = MAF)
  SNPs <- level_list$SNPs
  levels_dropped <- level_list$SNPs_Dropped
  return(list(SNPs = SNPs,Y = Y,SNPs_Dropped = levels_dropped))
}
