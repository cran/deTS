\name{query_ENCODE}
\alias{query_ENCODE}
\docType{data}
\title{ENCODE raw query data
}
\description{An example of RNA-Seq query data from ENCODE data for tissue-specific enrichment analysis
}
\usage{data("query_ENCODE")}
\format{
  A data frame with average expression level of 18661 genes in 44 ENCODE tissues.
  
  Row is genes symbol and column is sample names.
  
             Adrenal Gland Body of Pancreas Breast Epithelium Camera-type Eye Cerebellum
			 
TSPAN6           11.64            5.390             11.04           24.65     13.238

TNMD              0.01            0.147              2.24           12.43      0.090

DPM1             18.82            9.812             14.21           24.02     10.505

SCYL3             3.81            2.593              5.63           10.50      3.783

C1orf112          1.64            2.308              2.86           14.61      7.345

}
\details{nothing
}
\source{nothing
}
\references{Pei G., Dai Y., Zhao Z., Jia P. (2019) deTS: Tissue-Specific Enrichment Analysis to decode tissue specificity. Bioinformatics, In submission.
}
\examples{
data(query_ENCODE)
}
\keyword{datasets}
