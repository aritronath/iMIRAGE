# iMIRAGE: imputed microRNA activity from gene expression 

As the name stands, iMIRAGE a set of tools that facilitates imputation of microRNA or miRNA expression using protein-coding genes. 

With a growing repertoire of publicly available transcriptomic data, it is now possible to study miRNA expression in several phenotypes. However, for a large number of transcriptomic studies still do not produce reliable miRNA profiles due to technical limitations of standard library preparation and expression profiling techniques. In case of such studies with reliable protein-coding data, it is possible to *impute* the expression of miRNAs using the iMIRAGE package. Essentially, the iMIRAGE package builds prediction models using protein-coding and miRNA expression data sets using machine learning, and then imputes the miRNA expression in the independent data set of interest using just their protein-coding profiles. 

## Download and Installation in R:
```
#Download using devtools 
devtools::install_github("aritronath/iMIRAGE")

#Load the package in R
library(iMIRAGE)
```
See the [**iMIRAGE vignette**](https://aritronath.github.io/iMIRAGE/articles/imirage.html) for further details, quick start guide and detailed description of the package features. 

See [**Reference**](https://aritronath.github.io/iMIRAGE/reference/index.html) for detailed documentation on various functions included in the package. 

