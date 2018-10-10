#' @title Pre-processing gene expression matrix 
#' 
#' @description Use this function to perform a number of pre-processing steps, including 
#' log2 transformation, normalization and standardization of gene expression matrix. 
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
#' filter.exp(train_mir, cutoff=75, threshold=0)
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
#' temp <- match.gex(train_pcg, my_pcg)
#' train_pcg <- data.matrix(temp[[1]])
#' my_pcg <- data.matrix(temp[[2]])
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

#' @title Filter training gene expression features  
#'
#' @description Function to filter informative protein coding genes based on correlation with 
#' miRNA expression
#' 
#' @param train_pcg training protein coding expression dataset. a numeric matrix with row names 
#' indicating samples, and column names indicating protein coding gene IDs.
#' @param train_mir training miRNA expression dataset. a numeric matrix with row names indicating 
#' samples, and column names indicating miRNA IDs
#' @param gene_index either gene name (character) or index (column number) of miRNA to be imputed.
#' @param num number of informative protein coding genes to be used in constructing imputation model. 
#' Default is 50 genes.
#'
#' @return a numeric matrix. subset of protein coding genes correlated with miRNA of interest.
corf <- function (train_pcg, train_mir, gene_index, num=50) {
  pcor <- abs(cor(train_pcg, train_mir[, gene_index]))
  r_pcor <- rank(-pcor)
  gin <- which (r_pcor < num)
  temp_pcg <- train_pcg[, gin]
  return(temp_pcg)
}


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
#' @param method method for imputation, either RF or KNN, for random forests and K-nearest neighbor 
#' respectively.
#' @param num number of informative protein coding genes to be used in constructing imputation model. 
#' Default is 50 genes.
#' @param folds number specifying folds (k) of cross validation to obtain imputation accuracy. 
#' Default is k=10.
#'
#' @return a matrix with three values corresponding to Pearson's correlation coefficient,
#' P-value of the fit and root mean squared error (RMSE).
#'
#' @examples 
#' imirage_cv(train_pcg, train_mir, gene_index="ENSG00000184441", method="RF", num=100)
#' imirage_cv(train_pcg, train_mir, gene_index=25, method="RF", num=100)
#' 
#' @import randomForest
#' @import FNN
#' @export
imirage.cv <- function (train_pcg, train_mir, gene_index, num=50, method, folds=10) {

  if (mode(gene_index)!="numeric" & mode(gene_index)!="character") stop ("Error: miRNA not found in training dataset. Please check gene name or rownumber")
  if (mode(gene_index)=="numeric" & gene_index > ncol(train_mir))  stop ("Error: miRNA not found in training dataset. Please check ID or rownumber")
  if (mode(gene_index)=="character" & is.na (match (gene_index, colnames(train_mir)))) stop ("Error: miRNA not found. Please check ID or rownumber")

  cv.res <- matrix (nrow=folds, ncol=3)
  colnames (cv.res) <- c("PCC", "P-Value", "RMSE")

  if (num >= 50 & ncol(train_pcg) >=50 ) x <- scale (corf(train_pcg, train_mir, gene_index, num))
  if (ncol(train_pcg) < 50) x <- scale(train_pcg)
  
  if (mode(gene_index)=="numeric") {
    if (num >= 50 & ncol(train_pcg) >=50 ) train_pcg <- scale(corf(train_pcg, train_mir, gene_index, num))
    if (ncol(train_pcg) < 50) train_pcg <- scale(train_pcg)
  }
  if (mode(gene_index)=="character") {
    if (num >= 50 & ncol(train_pcg) >=50 ) train_pcg <- scale(corf(train_pcg, train_mir, gene_index, num))
    if (ncol(train_pcg) < 50) train_pcg <- scale(train_pcg)
  }
  
  train_mir <- scale(train_mir)

  cat("\nRunning ", folds,"-folds cross-validation...", sep="")

  kgrp <- split(1:nrow(train_pcg), sample(1:folds, nrow(train_pcg), replace=T))

  for (k in 1:folds) {

    cat("\nIteration", k)
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
      imp.rf <- randomForest(x, y, ntree=250)
      predict.y <- predict(imp.rf, testx)
      r.rf <- cor.test(predict.y, actual.y, method="spearman")
      rmse.rf <- sqrt(mean((actual.y-predict.y)^2))
      cv.res [k, 1:3] <- c(r.rf$estimate, r.rf$p.value, rmse.rf)
    }
    
    if (method=="KNN") {
      #knn regression 
      knn.fit <- knn.reg(train = x, test = testx, y = y, k=50)
      r.knn <- cor.test(knn.fit$pred, actual.y, method="spearman")
      rmse.knn <- sqrt(mean((actual.y - knn.fit$pred)^2))
      cv.res [k, 1:3] <- c(r.knn$estimate, r.knn$p.value, rmse.knn)
    }
  }
  cat("\nCross-validation complete\n")
  return (cv.res)
}

#' @title iMIRAGE: miRNA Activity from Gene Expression
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
#' @param method method for imputation, either RF or KNN, for random forests and K-nearest neighbor 
#' respectively.
#' @param num number of informative protein coding genes to be used in constructing imputation model. 
#' Default is 50 genes.
#'
#' @return a numeric vector containing imputed & standardized expression levels of the miRNA.
#'
#' @examples 
#' imirage(train_pcg, train_mir, my_pcg, gene_index="ENSG00000228630", method="KNN", num=50)
#' imirage(train_pcg, train_mir, my_pcg, gene_index=25, method="KNN", num=50)
#' 
#' @import randomForest
#' @import FNN
#' @export
imirage <- function (train_pcg, train_mir, my_pcg, gene_index, method, num=50) {

  if (mode(train_pcg)!="numeric" | mode(train_mir)!="numeric" | mode(my_pcg)!="numeric" |
      class(train_pcg)!="matrix" | class(train_mir)!="matrix" | class(my_pcg)!="matrix") stop ("Error: input data must be a numeric matrix")

  if (mode(gene_index)=="numeric" & gene_index > ncol(train_mir))  stop ("Error: miRNA not found in training dataset")
  if (mode(gene_index)=="character" & is.na (match (gene_index, colnames(train_mir)))) stop ("Error: miRNA not found")
  if (mode(gene_index)=="numeric") y <- scale(train_mir[, gene_index])
  if (mode(gene_index)=="character") y <- scale(train_mir[, match(gene_index, colnames(train_mir))])

  if (sd(y)==0) stop ("Error: Standard deviation of miRNA is 0")
  if (nrow(train_pcg)<50) warning("Warning: Sample size is <50")

  if (num >= 50 & ncol(train_pcg) >=50 ) x <- scale (corf(train_pcg, train_mir, gene_index, num))
  if (ncol(train_pcg) < 50) x <- scale(train_pcg)
  
  if (method=="RF") {
    rfit <- randomForest(x, y, ntree=250)
    predict.y <-  predict(rfit, my_pcg)
    return(predict.y)
  }
  
  if (method=="KNN") {
    knn.fit <- knn.reg(train = x, y = y, k=50)
    return(knn.fit$pred)
  }
}
