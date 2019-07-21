#' @title Pre-standardization filter to select miRNA expressed in a specific proportion of samples
#'
#' @description return a subset of the miRNA expression dataset containing miRNAs that are expressed
#' in at least a specified percentage of samples based on a user-defined expression threshold
#'
#' @param train_mir an expression matrix with miRNA in columns, samples in rows
#' @param cutoff percentage of samples in which the miRNA should be expressed
#' @param threshold the numeric threshold defining expression of miRNA. Default threshold = 0.
#' @return filtered miRNA expression matrix
#'
#' @examples
#' New.mir <-  filter.exp(GA.mir, cutoff=95, threshold=0)
#'
#' @export
filter.exp <- function (train_mir, cutoff=75, threshold=0) {

  if (mode(train_mir)!="numeric" | class(train_mir)!="matrix") stop ("Error: input data must be a numeric matrix")

  index <- 0

  cut <- cutoff*nrow(train_mir)/100

  for (i in 1:ncol(train_mir)) {
    if (sum(train_mir[, i] > threshold) > cut)
      index <- append(index, i)
  }

  index <- index[2:length(index)]
  f_train_mir <- train_mir[, index]

  return (f_train_mir)
}
