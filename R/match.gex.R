
#' @title Match gene expression matrices by gene ID
#'
#' @description Returns expression matrices with common genes
#'
#' @param train_pcg Gene expression matric with gene IDs as column names
#' @param my_pcg Gene expression matric with gene IDs as column names
#'
#' @return A list containing two matrices with common genes
#'
#' @examples
#' temp <- match.gex(GA.pcg, HS.pcg)
#' Matched.GA.pcg <- data.matrix(temp[[1]])
#' Matched.HS.pcg <- data.matrix(temp[[2]])
#'
#' @export
match.gex <- function (train_pcg, my_pcg) {

  if (mode(train_pcg)!="numeric" | class(train_pcg)!="matrix" |
      mode(my_pcg)!="numeric" | class(my_pcg)!="matrix" ) stop ("Error: input data must be a numeric matrix")

  y <- match (colnames(train_pcg), colnames(my_pcg))
  y1 <- which(!is.na(y))
  y2 <- na.omit(y)

  ncom <- length(y1)
  if (ncom < 0.25*nrow(train_pcg)) stop ("Error: <25% common genes in training and test dataset")
  if (ncom < 200) stop ("Error: <200 common genes in the two datasets")

  new.train_pcg <- train_pcg[, y1]
  new.my_pcg <- my_pcg[, y2]
  cmat <- list(new.train_pcg, new.my_pcg)

  return(cmat)
}
