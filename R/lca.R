lca_preliminary <- function(formula = cbind(A, B, C, D, E, F, G) ~ 1,
                            data = carcinoma,
                            nclass=2:4)
{

  mods <- lapply(nclass, function(n) poLCA::poLCA(
    formula=formula, data=data, nclass=n, nrep=10, verbose = FALSE))

  ## observed/expected frequencies
  gof_OE_freqs <- {
    structure <- mods[[1]]$predcell
    structure <- structure[1:(ncol(structure) - 1)]
    expected <- lapply(mods, function(m) m$predcell[, "expected"])
    names(expected) <- paste0("expected", nclass)
    do.call(cbind,  c(list(structure), expected))
  }

  ## gof criteria table
  llik <- sapply(mods, function(x) x$llik) # loglikelihoods
  npar <- sapply(mods, function(x) x$npar) # ? number of param? (CHECK)
  lrt <- sapply(mods, function(x) x$Gsq)
  chisq <- sapply(mods, function(x) x$Chisq)
  df <- sapply(mods, function(x) x$resid.df)
  chisq_pvalue <- round(1 - pchisq(chisq, df),4)
  lrt_pvalue <- round(1 - pchisq(lrt, df), 4)
  AIC <- sapply(mods, function(x) x$aic)
  BIC <- sapply(mods, function(x) x$bic)
  gof_table <- data.frame("K" = nclass, llik, npar, lrt, chisq, df,
                          chisq_pvalue, lrt_pvalue, AIC, BIC)
  
  list(
    "gof_OE_freqs" = gof_OE_freqs,
    "gof_table" = gof_table
  )
}


