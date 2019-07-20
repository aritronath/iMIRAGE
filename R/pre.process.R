#' @title Pre-processing gene expression matrix
#'
#' @description Use this function to perform a number of pre-processing steps, including
#' log transformation, normalization and standardization of gene expression matrix.
#' NOTE: If multiple options are set true, the sequence of pre-processing steps
#' performed are variance filtering, log2 transformation, upper quantile normalization
#' and standardization.
#'
#' @param gex a gene expression matrix with genes in columns and samples in rows
#' @param var.filter logical. specify whether genes are filtered based on variance. Default = TRUE
#' @param log logical. specify whether log2(x+1) transformation should be performed. Default = FALSE
#' @param UQ logical. specify whether upper quantile normalization should be performed. Default = FALSE
#' @param scale logical. specifiy whether standardization should be performed. Default = TRUE
#'
#' @return a processed gene expression matrix
#'
#' @examples
#' train_pcg <- pre.process(train_pcg)
#'
#' @export
pre.process <- function (gex, var.filter=TRUE, log=FALSE, UQ=FALSE, std=TRUE) {
  if (mode(gex)!="numeric" | class(gex)!="matrix") stop ("Error: input data must be a numeric matrix")

  if (var.filter == TRUE) {
    gvar <- apply(gex, 2, var)
    gex <- gex[, which(gvar!=0)]
  }

  if(log == TRUE) gex <- log2(gex + 1)

  if (UQ == TRUE) gex <- apply(gex, 2, function (x) x/quantile(x, 0.75))

  if (std == TRUE) gex <- scale(gex)

  return(gex)
}
