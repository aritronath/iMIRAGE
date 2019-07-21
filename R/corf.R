
#' @title Filter training gene expression features
#'
#' @description Internal function used by \link[iMIRAGE]{imirage.cv}, \link[iMIRAGE]{imirage.cv.loop},
#' and \link[iMIRAGE]{imirage.cv} to filter informative protein coding genes based on correlation with
#' miRNA expression.
#'
#' @param train_pcg training protein coding expression dataset. a numeric matrix with row names
#' indicating samples, and column names indicating protein coding gene IDs.
#' @param train_mir training miRNA expression dataset. a numeric matrix with row names indicating
#' samples, and column names indicating miRNA IDs
#' @param gene_index either gene name (character) or index (column number) of the miRNA to be imputed.
#' @param num number of informative protein coding genes to be used in constructing imputation model.
#' Default is 50 genes.
#' @param target "none" (default), "ts.pairs", or dataframe/matrix/list.
#' this argument accepts character strings to indicate the use of all candidate genes as predictors ("none),
#' or use built-in TargetScan miRNA-gene pairs ("ts.pairs"). also accepts a dataframe , matrix or list object
#' containing a column with names of miRNA and a column with the names of target genes.
#'
#' @return a numeric matrix containing subset of protein coding genes correlated with miRNA of interest.
#' @export corf
corf <- function (train_pcg, train_mir, gene_index, num=50, target="none") {
  if (target=="ts.pairs") {
    m1 <- grep(colnames(train_mir)[gene_index], ts.pairs$miRNA) #get all entries in the ts.pair table

    if (sum(!is.na(m1))==0) {
      pcor <- abs(cor(train_pcg, train_mir[, gene_index]))
      r_pcor <- rank(-pcor)
      gin <- which (r_pcor < num)
      temp_pcg <- train_pcg[, gin]
      warning ("miRNA not found in target-pair table. Using all genes")
    }
    if (sum(!is.na(m1)) > 0) {
      m2 <- match(ts.pairs$GeneID[m1], colnames(train_pcg))  #get target IDs
      temp <- train_pcg[, na.omit(m2)]
      pcor <- abs(cor(temp, train_mir[, gene_index]))
      r_pcor <- rank(-pcor)
      gin <- which (r_pcor < num)
      temp_pcg <- train_pcg[, gin]
    }
  } else if (target=="none") {
    pcor <- abs(cor(train_pcg, train_mir[, gene_index]))
    r_pcor <- rank(-pcor)
    gin <- which (r_pcor < num)
    temp_pcg <- train_pcg[, gin]
  } else {
    mirs <- grep("mirna", colnames(target), ignore.case=TRUE)
    gens <- grep("gene", colnames(target), ignore.case=TRUE)
    m1 <- grep(colnames(train_mir)[gene_index], target[,mirs]) #get all entries in the ts.pair table

    if (sum(!is.na(m1))==0) {
      pcor <- abs(cor(train_pcg, train_mir[, gene_index]))
      r_pcor <- rank(-pcor)
      gin <- which (r_pcor < num)
      temp_pcg <- train_pcg[, gin]
      warning ("miRNA not found in target-pair table. Using all genes")
    }

    if (sum(!is.na(m1)) > 0) {
      m2 <- match(target[m1,gens], colnames(train_pcg))  #get target IDs
      temp <- train_pcg[, na.omit(m2)]
      pcor <- abs(cor(temp, train_mir[, gene_index]))
      r_pcor <- rank(-pcor)
      gin <- which (r_pcor < num)
      temp_pcg <- train_pcg[, gin]
    }
  }

  return(temp_pcg)
}
