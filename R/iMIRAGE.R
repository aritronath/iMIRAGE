
#' @title Filter miRNA dataset based on positive expression in a specific percentage of samples
#'
#' @description Function to return miRNA expression subset containing miRNAs that are expressed in a specific percentage of samples
#'
#' @param train_mir a numeric matrix with miRNA expression data
#' @param cutoff percentage of samples in which the miRNA should be expressed
#' @param threshold the numeric threshold that defines whether the miRNA is expressed or not. Defaults to 1 (for RPKM or FPKM)
#' @return filtered miRNA expression matrix
#'
#' @examples posexp(train_mir, cutoff=75, threshold=1)
#' @export
posexp <- function (train_mir, cutoff=75, threshold=1) {

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


#' @title Match matrices by gene ID
#'
#' @description Function to return matched matrices with common genes
#'
#' @param train_pcg a numeric matrix. training protein coding dataset to construct miRNA imputation model. Gene IDs are specified in column names.
#' @param my_pcg a numeric matrix. test protein coding dataset to impute miRNA profile. Gene IDs are specified in column names.
#'
#' @return list containing two matrices with common protein coding genes
#'
#' @examples temp <- match_mat(train_pcg, my_pcg)
#' @examples new_train_pcg <- data.matrix(temp[[1]])
#' @examples new_my_pcg <- data.matrix(temp[[2]])
#'
#' @export
match_mat <- function (train_pcg, my_pcg) {

  if (mode(train_pcg)!="numeric" | class(train_pcg)!="matrix" |
      mode(my_pcg)!="numeric" | class(my_pcg)!="matrix" ) stop ("Error: input data must be a numeric matrix")

  y <- match (colnames(train_pcg), colnames(my_pcg))
  y1 <- which(!is.na(y))
  y2 <- na.omit(y)

  ncom <- length(y1)
  if (ncom < 0.5*nrow(train_pcg)) stop ("Error: <50% common genes in training and test dataset")
  if (ncom < 500) stop ("Error: <500 common genes training and test datasets")

  new.train_pcg <- train_pcg[, y1]
  new.my_pcg <- my_pcg[, y2]
  cmat <- list(new.train_pcg, new.my_pcg)

  return(cmat)
}


#' @title Select informative PCG subset
#'
#' @description Function to filter informative protein coding genes based on correlation with miRNA expression
#'
#' @param train_pcg training protein coding dataset. a numeric matrix with row names indicating samples, and column names indicating protein coding gene IDs.
#' @param train_mir training miRNA expression dataset. a numeric matrix with row names indicating samples, and column names indicating miRNA IDs
#' @param gene_index either gene name (character) or index (column number) of miRNA to be imputed.
#' @param num number of informative protein coding genes to be used in constructing imputation model. Default is 100 genes.
#'
#' @return a numeric matrix. subset of protein coding genes correlated with miRNA of interest.
corf <- function (train_pcg, train_mir, gene_index, num=100) {
  pcor <- abs(cor(train_pcg, train_mir[, gene_index]))
  r_pcor <- rank(-pcor)
  gin <- which (r_pcor < num)
  temp_pcg <- train_pcg[, gin]
  return(temp_pcg)
}


#' @title Cross validation function for imputation accuracy
#'
#' @description Function to obtain accuracy parameters: correlation coefficient, P-value and RMSE of imputation model
#'
#' @param train_pcg training protein coding dataset. a numeric matrix with row names indicating samples, and column names indicating protein coding gene IDs.
#' @param train_mir training miRNA expression dataset. a numeric matrix with row names indicating samples, and column names indicating miRNA IDs
#' @param gene_index either gene name (character) or index (column number) of miRNA to be imputed.
#' @param num number of informative protein coding genes to be used in constructing imputation model. Default is 100 genes.
#' @param folds number specifying folds of cross validation to obtain imputation accuracu. Default is 5.
#'
#' @return a matrix with three values corresponding to Pearson's correlation coefficient, P-value of fit and root mean square error
#'
#' @examples mirage_cv(train_pcg, train_mir, gene_index="ENSG00000184441", num=100)
#' @examples mirage_cv(train_pcg, train_mir, gene_index=25, num=100)
#'
#' @import randomForest
#'
#' @export
mirage_cv <- function (train_pcg, train_mir, gene_index, num=100, folds=5) {

  if (mode(gene_index)!="numeric" & mode(gene_index)!="character") stop ("Error: miRNA not found in training dataset. Please check gene name or rownumber")
  if (mode(gene_index)=="numeric" & gene_index > ncol(train_mir))  stop ("Error: miRNA not found in training dataset. Please check ID or rownumber")
  if (mode(gene_index)=="character" & is.na (match (gene_index, colnames(train_mir)))) stop ("Error: miRNA not found. Please check ID or rownumber")

  cv.res <- matrix (nrow=5,ncol=3)
  colnames (cv.res) <- c("PCC", "P-Value", "RMSE")

  if (mode(gene_index)=="numeric") train_pcg <- scale(corf(train_pcg, train_mir, gene_index, num))
  if (mode(gene_index)=="character") train_pcg <- scale (corf(train_pcg, train_mir, match(gene_index, colnames(train_pcg)), num))

  train_mir <- scale(train_mir)

  cat("\nRunning",folds,"folds cross-validation...")
  kgrp <- split(1:nrwow(train_pcg), sample(1:k, nrow(train_pcg), replace=T))
  for (k in 1:folds) {

    cat("\nIteration",k)
    ind <- unlist(kgrp[[i]])
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

    #randomForest
    imp.rf <- randomForest(x, y, ntree=100)
    predict.y <- predict(imp.rf, testx)
    r.rf <- cor.test(predict.y, actual.y, method="spearman")
    rmse.rf <- sqrt(mean((actual.y-predict.y)^2))
    cv.res [k, 1:3] <- c(r.rf$estimate, r.rf$p.value, rmse.rf)
  }
  cat("\nCross-validation complete\n")
  return (cv.res)
}

#' @title MIRAGE: miRNA Activity from Gene Expression
#'
#' @description Function to impute miRNA expression profile from protein coding expression dataset
#'
#' @param train_pcg training protein coding dataset. a numeric matrix with with row names indicating samples, and column names indicating protein coding gene IDs.
#' @param train_mir training miRNA expression dataset. a numeric matrix with row names indicating samples, and column names indicating miRNA IDs
#' @param my_pcg test protein coding expression dataset. a numberic matrix with row names indicating samples, and column names indicating protein coding gene IDs.
#' @param gene_index either gene name (character) or index (column number) of miRNA to be imputed.
#' @param num number of informative protein coding genes to be used in constructing imputation model. Default is 100 genes.
#'
#' @return a numeric vector containing imputed & standardized expression levels of the miRNA
#'
#' @import randomForest
#'
#' @examples mirage(train_pcg, train_mir, my_pcg, gene_index="ENSG00000228630", num=100)
#' @examples mirage(train_pcg, train_mir, my_pcg, gene_index=25, num=100)
#'
#' @export
mirage <- function (train_pcg, train_mir, my_pcg, gene_index, num=100) {

  if (mode(train_pcg)!="numeric" | mode(train_mir)!="numeric" | mode(my_pcg)!="numeric" |
      class(train_pcg)!="matrix" | class(train_mir)!="matrix" | class(my_pcg)!="matrix") stop ("Error: input data must be numberic matrix")

  if (mode(gene_index)=="numeric" & gene_index > ncol(train_mir))  stop ("Error: miRNA not found in training dataset. Please check ID or rownumber")
  if (mode(gene_index)=="character" & is.na (match (gene_index, colnames(train_mir)))) stop ("Error: miRNA not found. Please check ID or rownumber")
  if (mode(gene_index)=="numeric") y <- scale(train_mir[, gene_index])
  if (mode(gene_index)=="character") y <- scale(train_mir[, match(gene_index, colnames(train_mir))])

  if (sd(y)==0) stop ("Error: Standard deviation of miRNA is 0")
  if (nrow(train_pcg)<200) stop ("Error: Sample size insufficient (<200)")

  x <- scale (corf(train_pcg, train_mir, gene_index, num))

  rfit <- randomForest(x, y, ntree=100)
  predict.y <-  predict(rfit, my_pcg)
  return(predict.y)
}
