
#' @title Cross validation function for iMIRAGE imputation accuracy
#'
#' @description Function to obtain accuracy parameters: correlation coefficient, P-value and RMSE of
#' imputation model
#'
#' @param train_pcg training protein coding dataset. a numeric matrix with row names indicating
#' samples, anSed column names indicating protein coding gene IDs.
#' @param train_mir training miRNA expression dataset. a numeric matrix with row names indicating
#' samples, and column names indicating miRNA IDs
#' @param gene_index either gene name (character) or index (column number) of miRNA to be imputed.
#' @param method method for imputation, either "RF" for random forests, "KNN" for K-nearest neighbor or
#' "SVM" for support vector machines.
#' @param num number of informative protein coding genes to be used in constructing imputation model.
#' Default is 50 genes.
#' @param folds number specifying folds (k) of cross validation to obtain imputation accuracy.
#' Default is k=10.
#' @param target logical, specifying whether protein coding genes should be restricted to predicted
#' targets of the miRNA (from TargetScan) or use all genes as candidates. Default = FALSE.
#' @param ... optional parameters that can be passed on to the machine-learning functions
#' RF (\link[randomForest]{randomForest}), KNN (\link[FNN]{knn.reg}) or SVM(\link[e1071]{svm})
#'
#' @return a matrix with three values corresponding to Spearman's correlation coefficient,
#' P-value of the fit and root mean squared error (RMSE).
#'
#' @examples
#' imirage_cv(train_pcg, train_mir, gene_index="ENSG00000184441", method="RF", num=100)
#' imirage_cv(train_pcg, train_mir, gene_index=25, method="RF", num=100)
#'
#' @import randomForest
#' @import FNN
#' @import e1071
#' @import glmnet
#' @export imirage.cv
imirage.cv <- function (train_pcg, train_mir, gene_index, num=50, method, folds=10, target=FALSE, verbose=TRUE, ...) {

  if (mode(gene_index)!="numeric" & mode(gene_index)!="character") stop ("Error: miRNA not found in training dataset. Please check gene name or rownumber")
  if (mode(gene_index)=="numeric" & gene_index > ncol(train_mir))  stop ("Error: miRNA not found in training dataset. Please check ID or rownumber")
  if (mode(gene_index)=="character" & is.na (match (gene_index, colnames(train_mir)))) stop ("Error: miRNA not found. Please check ID or rownumber")

  cv.res <- matrix (nrow=folds, ncol=3)
  colnames (cv.res) <- c("PCC", "P-Value", "RMSE")

  if (num >= 50 & ncol(train_pcg) >=50 ) x <- corf(train_pcg, train_mir, gene_index, num, target)
  if (ncol(train_pcg) < 50) x <- train_pcg

  if (mode(gene_index)=="numeric") {
    if (num >= 50 & ncol(train_pcg) >=50 ) train_pcg <- corf(train_pcg, train_mir, gene_index, num, target)
    if (ncol(train_pcg) < 50) train_pcg <- train_pcg
  }
  if (mode(gene_index)=="character") {
    if (num >= 50 & ncol(train_pcg) >=50 ) train_pcg <- corf(train_pcg, train_mir, gene_index, num, target)
    if (ncol(train_pcg) < 50) train_pcg <- train_pcg
  }

  cat("\nRunning ", folds,"-folds cross-validation...", sep="")

  kgrp <- split(1:nrow(train_pcg), sample(1:folds, nrow(train_pcg), replace=T))

  for (k in 1:folds) {

    if(verbose==TRUE) {
      cat("\nIteration ", k, sep="")
    } else print("")

    ind <- unlist(kgrp[[k]])
    x <- train_pcg [-ind, ]
    testx <- train_pcg [ind, ]

    if (mode(gene_index)=="numeric") {
      y <- train_mir[-ind, gene_index]
      actual.y <- train_mir [ind, gene_index]
    }

    if (mode(gene_index)=="character") {
      y <- train_mir[-ind, match(gene_index, colnames(train_mir))]
      actual.y <- train_mir[ind, match(gene_index, colnames(train_mir))]
    }

    if(method=="RF") {
      #randomForest
      imp.rf <- randomForest(x, y, ntree=250, ...)
      predict.y <- predict(imp.rf, testx)
      r.rf <- cor.test(predict.y, actual.y, method="spearman")
      rmse.rf <- sqrt(mean((actual.y-predict.y)^2))
      cv.res [k, 1:3] <- c(r.rf$estimate, r.rf$p.value, rmse.rf)
      next()
    }

    if (method=="KNN") {
      #knn regression
      knn.fit <- knn.reg(train = x, test = testx, y = y, k=50, ...)
      r.knn <- cor.test(knn.fit$pred, actual.y, method="spearman")
      rmse.knn <- sqrt(mean((actual.y - knn.fit$pred)^2))
      cv.res [k, 1:3] <- c(r.knn$estimate, r.knn$p.value, rmse.knn)
      next()
    }

    if (method=="SVM") {
      #svm regression
      imp.svm <- svm(x, y, ...)
      predict.y <- predict(imp.svm, testx)
      r.svm <- cor.test(predict.y, actual.y, method="spearman")
      rmse.svm <- sqrt(mean((actual.y-predict.y)^2))
      cv.res [k, 1:3] <- c(r.svm$estimate, r.svm$p.value, rmse.svm)
      next()
    }

  }

  if (verbose==TRUE) cat("\nCross-validation complete\n")
  return (cv.res)
}
