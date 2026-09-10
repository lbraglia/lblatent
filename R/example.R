#' @export
example <- function(){
  

  gives_names_to_columns <- TRUE


  library(ltm)
  data(Abortion)
  res_ltm <- ltm_estimate(Abortion)

  library(poLCA)
  data(carcinoma)
  lca_preliminary()
  res_lca <- lca_estimate()

  list(
    "ltm"=res_ltm,
    "lca" = res_lca
  )
}
  
