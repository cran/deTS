\name{ENCODE_z_score}
\alias{ENCODE_z_score}
\docType{data}
\title{ENCODE z-score to define tissue-specific genes
}
\description{ENCODE z-score matrix to define tissue-specific genes
}
\usage{data("ENCODE_z_score")}
\format{
   A data frame with z-score of 14031 genes in 44 ENCODE tissues.
   
   Row is genes symbol and column is tissue names.
  
           Adrenal Gland Body of Pancreas Breast Epithelium Camera-type Eye Cerebellum
		   
C1orf112        -0.674           -0.440            -0.246           3.892      1.333

FGR             -0.078           -0.345             0.159          -0.354     -0.407

CFH             -0.093           -0.365            -0.134          -0.133     -0.160

FUCA2            3.028            1.467             0.040           0.228     -0.601

NFYA            -0.637           -0.872             0.053           2.364      0.619

}
\details{nothing
}
\source{nothing
}
\references{Pei G., Dai Y., Zhao Z., Jia P. (2019) deTS: Tissue-Specific Enrichment Analysis to decode tissue specificity. Bioinformatics, In submission.
}
\examples{
data(ENCODE_z_score)
}
\keyword{datasets}
