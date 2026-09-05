#' @export
efa_preliminary <- function(data=NULL, 
                            formula = ~.,
                            cormatrix=NULL,
                            nobs=NULL,
                            nfactors = 1:3
                            ){
  ## check input
  data_input <- is.data.frame(data)
  cormatrix_input <- is.matrix(cormatrix)
  
  if (!(data_input || cormatrix_input))
    stop("need to specify data or cormatrix")
  
  if (data_input){
    cormat <- cor(data)
  } else if (cormatrix_input){
    cormat <- cormatrix
  }
  
  ## correlation matrix
  cormat[abs(cormat) < 0.3] <- NA
  cormat[lower.tri(cormat)] <- NA

  ## factor models
  names(nfactors) <- nfactors
  if (data_input){
    fmodels <- lapply(nfactors, function(n) factanal(formula, data=data, factors=n))
  } else {
    fmodels <- lapply(nfactors, function(n) factanal(formula, covmat=cormatrix,
                                                     n.obs=nobs, factors=n))
  }
  
  ## chisq, df, p
  test <- do.call(rbind, lapply(fmodels, function(x)
    data.frame("chi"=  x$STATISTIC,
               "df"= x$dof,
               "p" = x$PVAL)))
  test$rounded_p <- round(test$p, 4)
  test$comment <- ifelse(test$p < 0.01, "bad/reject",
                         ifelse(test$p < 0.05, "meh", "ok"))
  
  list(
    ## "models" = fmodels,
    ## "loadings" = loads,
    "readable_cormatrix" = cormat,
    "chisq" = test
  )
}


#' @export
efa_estimate <- function(data=NULL,
                         formula = ~.,
                         cormatrix=NULL,
                         nobs=NULL,
                         nfactors = 1){
  ## check input
  data_input <- is.data.frame(data)
  cormatrix_input <- is.matrix(cormatrix)
  
  if (!(data_input || cormatrix_input))
    stop("need to specify data or cormatrix")

  ## factor models
  if (data_input){
    fmodel <- factanal(formula, data=data, factors=nfactors, rotation="none")
    cor_matrix <- cor(data)
  } else {
    fmodel <- factanal(covmat=cormatrix, n.obs=nobs, factors=nfactors, rotation="none")
    cor_matrix <- cormatrix
  }

  cor_matrix[lower.tri(cor_matrix)] <- NA
  
  ## loadings
  raw_loadings <- loadings(fmodel)

  readable_loadings <- function(loads){
    readable_loadings <- unclass(loads)
    readable_loadings[abs(readable_loadings) < 0.3] <- NA
    readable_loadings
  }
  
  ## communalities and uniqueness
  communalities <- rowSums(raw_loadings^2)
  uniqueness <- 1 - communalities

  ## reproduced_correlation matrix \Lambda \Lambda'
  cor_reproduced <- raw_loadings %*% t(raw_loadings)
  cor_reproduced[lower.tri(cor_reproduced)] <- NA
  ## residual matrix
  cor_residual <-  round(cor_matrix - cor_reproduced, 3)

  ## rotations
  varimax_loadings <- GPArotation::Varimax(raw_loadings)$loadings
  quartimax_loadings <- GPArotation::quartimax(raw_loadings)$loadings
  oblimin_loadings <- GPArotation::oblimin(raw_loadings)$loadings
  
  list(
    "model" = fmodel,
    "communalities" = communalities,
    "uniqueness" = uniqueness,
    "cor_original" = round(cor_matrix, 2),
    "cor_reproduced" = round(cor_reproduced, 2),
    "cor_residual" = round(cor_residual, 2),
    "loadings_raw" = readable_loadings(raw_loadings),
    "loadings_varimax" = readable_loadings(varimax_loadings),
    "loadings_quartimax" = readable_loadings(quartimax_loadings),
    "loadings_oblimin" = readable_loadings(oblimin_loadings)
  )
}
