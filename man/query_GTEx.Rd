\name{query_GTEx}
\alias{query_GTEx}
\docType{data}
\title{GTEx raw query data
}
\description{An example of RNA-Seq query data from GTEx data for tissue-specific enrichment analysis
}
\usage{data("query_GTEx")}
\format{
  A data frame with average expression level of 18067 gene in 47 GTEx tissues.
  
  Row is genes symbol and column is sample names.
  
         Adipose - Subcutaneous Adipose - Visceral (Omentum) Adrenal Gland Artery - Aorta Artery - Coronary
		 
OR4F5                    0.0317                       0.0284        0.0469         0.0133            0.0225

SAMD11                   0.4451                       2.3056        3.8928         3.5822            2.5632

NOC2L                   21.9084                      21.0439       19.4613        19.4929           19.8367

KLHL17                   4.1406                       4.4075        4.4227         5.0840            5.3749

PLEKHN1                  0.4531                       0.3452        1.1795         0.3081            0.3722

}
\details{nothing
}
\source{nothing
}
\references{Pei G., Dai Y., Zhao Z., Jia P. (2019) deTS: Tissue-Specific Enrichment Analysis to decode tissue specificity. Bioinformatics, In submission.
}
\examples{
data(query_GTEx)
}
\keyword{datasets}
