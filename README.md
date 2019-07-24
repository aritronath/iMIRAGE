# iMIRAGE: imputed microRNA activity from gene expression 

As the name stands, iMIRAGE a set of tools that facilitates imputation of microRNA or miRNA expression using protein-coding genes. 

## Usage

With a growing repertoire of publicly available transcriptomic data, it is now possible to study miRNA expression in several phenotypes. However, for a large number of transcriptomic studies still do not produce reliable miRNA profiles due to technical limitations of standard library preparation and expression profiling techniques. In case of such studies with reliable protein-coding data, it is possible to *impute* the expression of miRNAs using the iMIRAGE package. Essentially, the iMIRAGE package builds prediction models using protein-coding and miRNA expression data sets using machine learning, and then imputes the miRNA expression in the independent data set of interest using just their protein-coding profiles 

See the **iMIRAGE vignette** under **Articles** for further details. 

See **Reference** for detailed documentation on various functions included in the package. 

