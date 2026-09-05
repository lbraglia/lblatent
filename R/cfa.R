#' @export
cfa_estimate <- function(model,
                         data=NULL,
                         cormatrix =NULL,
                         nobs = NULL,
                         std.lv = TRUE,
                         orthogonal = FALSE)
{
  data_input <- is.data.frame(data)
  cormatrix_input <- is.matrix(cormatrix)
  
  if (!(data_input || cormatrix_input))
    stop("need to specify data or cormatrix")
  
  if (data_input) {
    fit <- lavaan::cfa(model, data = data, std.lv = std.lv, orthogonal = orthogonal)
  } else {
    fit <- lavaan::cfa(model, sample.cov = cormatrix, sample.nobs = nobs, std.lv = std.lv, orthogonal = orthogonal)
  }
  # browser()
  summary(fit, fit.measures = TRUE)
}
