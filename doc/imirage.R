## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width=7, fig.height=6
)

## ----results="hide"------------------------------------------------------
#Load the package in R
library(iMIRAGE)

## ------------------------------------------------------------------------
# return a pair of matrices with matched columns (protein-coding genes that will be used as training features)

temp <- match.gex(GA.pcg, HS.pcg)
GA.pcg <- temp[[1]]
HS.pcg <- temp[[2]]

## ------------------------------------------------------------------------
#Use **imirage.cv** with the default parameters of using K-nearest neighbors as the prediction algorithm and using 50 protein-coding genes as training features

CV.miRNA <- imirage.cv(train_pcg=GA.pcg, train_mir=GA.mir, gene_index="hsa-let-7c", method="KNN", target="none")

#Print the accuracy metrics from each fold of the cross-validation analysis 
print(CV.miRNA)

#note: this example performs cross-validation analysis for 1 unique miRNA, hsa-let-7c. The name must match with one of the names in the train_mir object's column names

## ----results="hide"------------------------------------------------------
#Perform cross-validation analysis over the entire training dataset 
CV.full <- imirage.cv.loop(train_pcg = GA.pcg, train_mir = GA.mir, method = "KNN", target="none")

#Plot some performance metrics from cross-validation analysis
plot(CV.full[,1], -log10(CV.full[,2]), xlab="Spearman R", ylab="P-value (-log10)", pch=16, main="Cross-validation accuracy")
abline(h=-log10(0.05), lty=3)

## ------------------------------------------------------------------------
#Find out which miRNAs are imputed with good accuracy 
colnames(GA.mir)[which(CV.full[,1] > 0.5)] #arbritarily, Spearman correlation > 0.5 

colnames(GA.mir)[which(CV.full[,2] < 0.05)] #P-value < 0.05

## ------------------------------------------------------------------------
#Use **imirage** with default parameters
Pred.miRNA <- imirage(train_pcg=GA.pcg, train_mir=GA.mir, my_pcg=HS.pcg , gene_index="hsa-let-7c", target="none")

#Display the predicted miRNA expression values
print(Pred.miRNA)

#Compare the predicted miRNA expression values with measured expression 
plot(Pred.miRNA, HS.mir[,"hsa-let-7c"], xlab="Predicted", ylab="Measured", main="hsa-let-7c imputation performance", pch=16)
abline(lm(HS.mir[,"hsa-let-7c"] ~ Pred.miRNA), lty=3)

## ------------------------------------------------------------------------
#Create an empty matrix to store imputed expression 
Pred.full <- matrix(data=NA, ncol=ncol(GA.mir), nrow=nrow(HS.pcg))

#Execute a loop to impute each miRNA using the test protein coding dataset
for (i in 1:ncol(GA.mir)) {
  Pred.full[,i] <- imirage(train_pcg = GA.pcg, train_mir = GA.mir, my_pcg = HS.pcg, gene_index = i, method="KNN", target="none")
}

## ------------------------------------------------------------------------
#Obtain correlation coefficients between imputed and measured miRNA expression 
Pred.Cors <- array(dim=ncol(GA.mir))
for (i in 1:ncol(GA.mir)) {
  Pred.Cors[i] <- cor(HS.mir[,i], Pred.full[,i], method="spearman")
}

#Plot imputation correlations in comparison to cross-validation results
plot(Pred.Cors, CV.full[,1], xlab="Imputation accuracy", ylab="Cross-validation accuracy", main="Imputation performance", pch=16)

## ------------------------------------------------------------------------
load("TCGA_BRCA_Datasets.RData")

## ------------------------------------------------------------------------
#Here, we filter the miRNA datasets to retain all miRNAs that are expressed above a level of 0 in atleast 75% of the samples
ga.mirna <- filter.exp(ga.mirna, cutoff=75, threshold = 0)
hiseq.mirna <- filter.exp(hiseq.mirna, cutoff=75, threshold = 0)


## ------------------------------------------------------------------------
#We also will keep the unprocessed hiseq.gex.1 and ga.gex.1 datasets for subsquent comparisons
temp <- match.gex(hiseq.gex, ga.gex)
hiseq.gex.1 <- temp[[1]]
ga.gex.1 <- temp[[2]]

#For subsquent prediction analyses, we will also select miRNAs are common to both datasets
temp <- match.gex(ga.mirna, hiseq.mirna)
ga.mirna.1 <- temp[[1]]
hiseq.mirna.1 <- temp[[2]]


## ------------------------------------------------------------------------
ga.gex.2 <- pre.process(ga.gex.1, log = TRUE, var.filter = TRUE, UQ = TRUE, std = TRUE)
hiseq.gex.2 <- pre.process(hiseq.gex.1, log = TRUE, var.filter = TRUE, UQ = TRUE, std = TRUE)

#This is how the data looks before and after pre.processing (showing the first 10 samples from each dataset)
par(mfrow=c(2,2), cex=0.75)
boxplot(t(ga.gex.1[1:10,]), main="GA - raw expression")
boxplot(t(ga.gex.2[1:10,]), main="GA - processed expression")
boxplot(t(hiseq.gex.1[1:10,]), main="Hiseq - raw expression")
boxplot(t(hiseq.gex.2[1:10,]), main="Hiseq - processed expression")


## ----results="hide"------------------------------------------------------
#Unprocessed training data
CV.ga.gex1 <- imirage.cv.loop(train_pcg = ga.gex.1, train_mir = ga.mirna.1, method = "KNN", target="none")
#Processed training data
CV.ga.gex2 <- imirage.cv.loop(train_pcg = ga.gex.2, train_mir = ga.mirna.1, method = "KNN", target="none")

#Unprocessed training data
CV.hs.gex1 <- imirage.cv.loop(train_pcg = hiseq.gex.1, train_mir = hiseq.mirna.1, method = "KNN", target="none")
#Processed training data
CV.hs.gex2 <- imirage.cv.loop(train_pcg = hiseq.gex.2, train_mir = hiseq.mirna.1, method = "KNN", target="none")


#Comparison of imputation accuracies between raw and processed data
boxplot(CV.ga.gex1[,1], CV.ga.gex2[,1], names=c("Raw", "Processed"), main="GA CV accuracy")

boxplot(CV.hs.gex1[,1], CV.hs.gex2[,1], names=c("Raw", "Processed"), main="Hiseq CV accuracy")


## ------------------------------------------------------------------------
sessionInfo()

