#' @title Process cross-validation results
#' @description This function returns useful information by organizing the output from cross-validation analysis.
#' Used internally by \link[iMIRAGE]{imirage.cv.loop}
#'
#' @param res The output object from imirage.cv
#'
#' @return a processed matrix containing 3 columns: Spearman's correlation coefficient, P-value and root mean
#' squared error from cross-validation analysis
#' @export CVProc
CVProc <- function (res) {
  df <- matrix(nrow=length(res), ncol=ncol(res[[1]]))
  colnames(df) <- colnames(res[[1]])
  for (i in 1:length(res)) {
    df[i,1] <- mean(res[[i]][,1])
    df[i,2] <- mean(res[[i]][,2])
    df[i,3] <- mean(res[[i]][,3])
  }
  colnames(df) <- c("Coef", "P-value", "RMSE")
  return(df)
}
