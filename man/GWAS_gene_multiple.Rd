\name{GWAS_gene_multiple}
\alias{GWAS_gene_multiple}
\docType{data}
\title{Gene symbol query data for multiple samples 
}
\description{An example of input gene symbol query data for multiple samples tissue-specific enrichment analysis
}
\usage{data("GWAS_gene_multiple")}
\format{
  A data frame with 22134 genes if associated with following 5 neuropsychiatric disorders GWAS traits.
  
  Row is genes symbol and column is sample names.
  
			   ADHD ASD BD MDD SCZ
			   
	A1BG       0   0  0   0   0
	
	A1BG-AS1   0   0  0   0   0
	
	A1CF       0   1  0   0   0
	
	A2M        0   0  0   0   0
	
	A2M-AS1    0   0  0   0   0
	
}
\details{nothing
}
\source{nothing
}
\references{Pei G., Dai Y., Zhao Z., Jia P. (2019) deTS: Tissue-Specific Enrichment Analysis to decode tissue specificity. Bioinformatics, In submission.
}
\examples{
data(GWAS_gene_multiple)
}
\keyword{datasets}
