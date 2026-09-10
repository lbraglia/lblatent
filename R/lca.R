#' @export
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
    sparse <- lapply(expected, function(x) x < 5)
    names(expected) <- paste0("expected_", nclass)
    names(sparse) <- paste0("sparse_", nclass)
    do.call(cbind,  c(list(structure), expected, sparse))
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


#' @export
lca_estimate <- function(formula = cbind(A, B, C, D, E, F, G) ~ 1,
                         data = carcinoma,
                         nclass = 3
                         )
{

  ## unordered solution
  unord_mod <- poLCA::poLCA(formula=formula, data=data, nclass=nclass, nrep=10,
                            verbose = FALSE)
  unordered_prior_probs <- round(unord_mod$P,4)
  unordered_conditional_probs <- lapply(unord_mod$probs, round, 2)

  ## ordered probs
  starting_probs <- unord_mod$probs.start
  new_probs <- poLCA::poLCA.reorder(starting_probs, order(unord_mod$P)) # order in the same way the starting values
  mod_ord <- poLCA::poLCA(formula = formula, # re-estimate
                          data = data,
                          nclass = nclass,
                          probs.start = new_probs,
                          verbose = FALSE)
  ordered_prior_probs <- round(mod_ord$P, 4)
  ordered_conditional_probs <- lapply(mod_ord$probs, round, 2)
  posterior_prob <- {
    the_probs <- setNames(data.frame(round(mod_ord$posterior, 3)), nm = seq_len(nclass))
    the_probs$assigned_class <- apply(the_probs, 1, which.max)
    rval <- unique(cbind(data, the_probs))
    rownames(rval) <- NULL
    rval
  }
    
  ## probability of agreement
  agreement_mod <- mod_ord
  agree_class <- 1:2
  names(agree_class) <- agree_class
  ## extract probabilities and add priors
  probs <- lapply(agree_class, function(class) data.frame(do.call(cbind, lapply(agreement_mod$probs, function(x) x[, class]))))
  probs <- lapply(probs, function(x) {x$prior <- agreement_mod$P; x})
  prob_agreement <- sapply(probs, function(x) sum(apply(x, 1, prod)))
  names(prob_agreement) <- {
    tmp <- as.matrix(data)
    dim(tmp) <- NULL
    as.character(sort(unique(tmp)))
  }
  
  pred_class_split <- {
    tmp <- unique(cbind(data, "predclass" = mod_ord$predclass))
    lapply(split(tmp, tmp$predclass), function(x){
      rownames(x) <- NULL
      x
    })
  }

  ricordella_pattern_predizione <- rep(1, ncol(data))
  ricordella_prob_of_observing_ptrn <- poLCA.predcell(
    mod_ord,
    ricordella_pattern_predizione
  )
  ricordella_posterior_prob_of_pattern <- poLCA.posterior(
    mod_ord,
    y = ricordella_pattern_predizione
  )
  
  list(
    ## "unord_mod" = mod,
    ## "unordered" = list("prior_probs" = unordered_prior_probs,
    ##                    "conditional_probs" = unordered_conditional_probs)
    "mod" = mod_ord,
    "probs" = list("prior" = ordered_prior_probs,
                   "conditional" = ordered_conditional_probs,
                   "posterior" = posterior_prob
                   ),
    "agreement" = list(
      "probability" = prob_agreement,
      "expected_count" = prob_agreement * agreement_mod$N
    ),
    "predicted_class" = list(
      "values" = mod_ord$predclass,
      "freqs" = table(mod_ord$predclass),
      "response_patterns" = pred_class_split
      )
  )
}
