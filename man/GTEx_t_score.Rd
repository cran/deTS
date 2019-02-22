\name{GTEx_t_score}
\alias{GTEx_t_score}
\docType{data}
\title{GTEx t-score to define tissue-specific genes
}
\description{GTEx t-score matrix to define tissue-specific genes
}
\usage{data("GTEx_t_score")}
\format{
  A data frame with t-score of 14725 genes in 47 GTEx tissues.
  
  Row is genes symbol and column is tissue names.
  
         Adipose - Subcutaneous Adipose - Visceral (Omentum) Adrenal Gland Artery - Aorta Artery - Coronary
		 
OR4F5                    -0.524                       -0.597         0.134         -1.109            -0.588

SAMD11                   -9.921                       -1.734         3.633          3.595             0.017

KLHL17                   -6.812                       -4.553        -3.084         -0.744             0.494

PLEKHN1                  -7.785                       -6.882        -3.915         -6.570            -4.892

C1orf170                 -7.113                       -6.257        -4.465         -5.897            -4.004

}
\details{nothing
}
\source{nothing
}
\references{Pei G., Dai Y., Zhao Z., Jia P. (2019) deTS: Tissue-Specific Enrichment Analysis to decode tissue specificity. Bioinformatics, In submission.
}
\examples{
data(GTEx_t_score)
}
\keyword{datasets}
