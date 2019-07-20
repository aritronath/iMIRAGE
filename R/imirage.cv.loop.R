
#' @title iMIRAGE cross-validation loop function for full miRNA matrix
#' @description Convinient wrapper for \link[iMIRAGE]{imirage.cv} that performs cross-validation analysis for
#' assessing imputation accuracies for all miRNAs using the training datasets
#' @param train_pcg training protein coding dataset. a numeric matrix with with row names indicating
#' samples, and column names indicating protein coding gene IDs.
#' @param train_mir training miRNA expression dataset. a numeric matrix with row names indicating
#' samples, and column names indicating miRNA IDs.
#' @param method method for imputation, either "RF" for random forests, "KNN" for K-nearest neighbor or
#' "SVM" for support vector machines.
#' @param ... optional parameters that can be passed on to the machine-learning method:
#' RF (\link[randomForest]{randomForest}), KNN (\link[FNN]{knn.reg}) or SVM(\link[e1071]{svm})
#'
#' @examples
#'
#' @return a matrix containing Spearman's correlation coefficient, P-value and RMSE from the cross-validation analysis
#' of the complete miRNA training dataset
#' @export imirage.cv.loop
imirage.cv.loop <- function (train_pcg, train_mir, method="KNN", ...) {
  cv.loop <- list()
  for (i in 1:ncol(train_mir)) {
    cv.loop[[i]] <- imirage.cv(train_pcg, train_mir, gene_index=i, method=method, ...)
    print(i)
  }

  cv.results <- CVProc(cv.loop)
  return(cv.results)
}
