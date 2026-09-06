#' @export

ltm_estimate <- function(data, plot_characteristic = FALSE){
  p <- ncol(data)
  des <- ltm::descript(data)
  mod <- ltm::ltm(data  ~ z1, IRT.param=FALSE)
  smod <- summary(mod)
  alpha <- mod$coeff[, 2] # just extract the slopes from coef matrix
  stalpha <- alpha / sqrt(1 + alpha^2)
  
  prob_median_individual <- coef(mod, prob=TRUE, order=TRUE)

  # GOF
  E <- fitted(mod)[,"Exp"] # expected: for each response patterns
  O <- mod$patterns$obs # observed freqs per response pattern
  gof_table <- cbind(setNames(data.frame(mod$pattern$X), nm=names(data)),
                     "O"=O,
                     "E"=E)
  chisq <- sum((E-O)^2/E)
  DOF <- nrow(gof_table) - nrow(coef(smod)) - 1  
  chisq_p <- 1 - pchisq(chisq, DOF)
  lrt <- 2*sum(O*log(O/E))
  lrt_p <- 1 - pchisq(lrt, DOF)
  gof_chisq <- data.frame("chisq"=chisq, "df"=DOF, "p"=chisq_p)
  gof_lrt <- data.frame("LRT"=lrt, "df"=DOF, "p"=lrt_p)
  gof_margins2 <- margins(mod) 
  gof_margins3 <- margins(mod, type="three-way", nprint=2)

  ## factor scores tab
  fs <- factor.scores(mod, method = "EAP")
  comp <- factor.scores(m1, method="Component")
  factor_scores <- cbind(
    fs$score.dat,  # factor scores from EAP
    "component" = comp$score.dat[, "z1"], # Components
    "tot_score" = apply(fs$score.dat[, seq_len(p)], 1, sum)) # total score (sum of agree)
  # order the table by total score: here differently from rasch model, total
  # score is not a sufficient statistics for underlying latent varialbe
  factor_scores <- round(factor_scores[order(factor_scores$tot_score),], 1)

  ## predicted out of sample (ricordella)
  ricordella_predict_outofsample <- factor.scores(mod, resp.pattern = rbind(c(1,0,0,1), c(1,0,1,0)))
  
  ## characteristic curves
  if (plot_characteristic){
    plot(mod, legend = TRUE, cx="bottomright", xlab="latent variable",
         lwd = 3, cex.main = 1.5, cex.lab = 1.3, cex = 1.1)
  }
  
  list(
    "des_perc" = des$perc,
    "des_items" = des$items,
    "des_assoc" = des$pw.ass,
    "mod" = mod,
    "smod" = smod,
    "stalpha" = stalpha,
    "prob_median_individual" = prob_median_individual,
    "gof_table" = gof_table,
    "gof_chisq" = gof_chisq,
    "gof_lrt" = gof_lrt,
    "gof_margins2" = gof_margins2,
    "gof_margins3" = gof_margins3,
    "factor_scores" = factor_scores
    )
}
