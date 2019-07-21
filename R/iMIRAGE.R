
#' @title iMIRAGE: imputated miRNA Activity from Gene Expression
#'
#' @description Function to impute miRNA expression profile from protein coding expression dataset
#'
#' @param train_pcg training protein coding dataset. a numeric matrix with with row names indicating
#' samples, and column names indicating protein coding gene IDs.
#' @param train_mir training miRNA expression dataset. a numeric matrix with row names indicating
#' samples, and column names indicating miRNA IDs.
#' @param my_pcg test protein coding expression dataset. a numeric matrix with row names indicating
#' samples, and column names indicating protein coding gene IDs.
#' @param gene_index either gene name (character) or index (column number) of miRNA to be imputed.
#' @param method method for imputation, either "RF" for random forests, "KNN" for K-nearest neighbor or
#' "SVM" for support vector machines. Uses KNN by default.
#' @param num number of informative protein coding genes to be used in constructing imputation model.
#' Default is 50 genes.
#' @param target "none" (default), "ts.pairs", or dataframe/matrix/list.
#' this argument accepts character strings to indicate the use of all candidate genes as predictors ("none),
#' or use built-in TargetScan miRNA-gene pairs ("ts.pairs"). also accepts a dataframe , matrix or list object
#' containing a column with names of miRNA and a column with the names of target genes.
#' @param ... optional parameters that can be passed on to the machine-learning method:
#' RF (\link[randomForest]{randomForest}), KNN (\link[FNN]{knn.reg}) or SVM(\link[e1071]{svm})
#'
#' @return a numeric vector with imputed expression levels of the miRNA
#'
#' @examples
#' data(iMIRAGE.datasets)
#' imirage(GA.pcg, GA.mir, HS.pcg, gene_index="hsa-let-7c", method="KNN", num=50)
#' imirage(GA.pcg, GA.mir, HS.pcg, gene_index=25, method="KNN", num=50)
#'
#' @import randomForest
#' @import FNN
#' @import e1071
#'
#' @export imirage
imirage <- function (train_pcg, train_mir, my_pcg, gene_index, method="KNN", num=50, target="none", ...) {

  if (mode(train_pcg)!="numeric" | mode(train_mir)!="numeric" | mode(my_pcg)!="numeric" |
      class(train_pcg)!="matrix" | class(train_mir)!="matrix" | class(my_pcg)!="matrix") stop ("Error: input data must be a numeric matrix")

  if (mode(gene_index)=="numeric" & gene_index > ncol(train_mir))  stop ("Error: miRNA not found in training dataset")
  if (mode(gene_index)=="character" & is.na (match (gene_index, colnames(train_mir)))) stop ("Error: miRNA not found")

  if (mode(gene_index)=="numeric") y <- train_mir[, gene_index]
  if (mode(gene_index)=="character") y <- train_mir[, match(gene_index, colnames(train_mir))]

  if (sd(y)==0) stop ("Error: Standard deviation of miRNA is 0")
  if (nrow(train_pcg)<50) warning("Warning: Sample size is <50")

  if (num >= 50 & ncol(train_pcg) >=50 ) x <- corf(train_pcg, train_mir, gene_index, num, target)
  if (ncol(train_pcg) < 50) x <- train_pcg

  if (method == "RF") {
    rfit <- randomForest(x, y, ntree=250, ...)
    predict.y <-  predict(rfit, my_pcg)
    return(predict.y)
  }

  if (method == "KNN") {
    mX <- match(colnames(x), colnames(my_pcg))
    predict.y <- knn.reg(train = x, test = my_pcg[,mX], y = y, k=50)$pred
    return(predict.y)
  }

  if (method == "SVM") {
    mX <- match(colnames(x), colnames(my_pcg))
    svm.fit <- svm(x = x, y = y)
    predict.y <- predict(svm.fit, my_pcg[,mX])
    return(predict.y)
  }
}


