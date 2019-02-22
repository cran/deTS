\name{tsea.expression.normalization}
\alias{tsea.expression.normalization}
\title{RNA-Seq expression profiles normalization
}
\description{To avoid the data bias and adapt better data heterogeneity, before tsea.expression.decode() analysis, the raw discrete RPKM value have to normalized to continuous variable meet the normal distribution before t-test.
}
\usage{
tsea.expression.normalization(query_mat, 
correction_factor, normalization = "abundance")
}
\arguments{
  \item{query_mat}{a RNA-seq RPKM object, row name should be gene symbol, and column name should be sample name.
}
  \item{correction_factor}{correction_factor, a gene table object contain genes average expression level and standard variance in GTEx database, can be loaded by data(correction_factor).
}
  \item{normalization}{normalization methods, c("z-score", "abundance")
}
}
\details{As RNA-Seq samples are often heterogeneous, before in-depth analysis, it is necessary to decode tissue heterogeneity to avoid samples with confounding effects. However, the raw discrete RPKM value have to normalized to continuous variable meet the normal distribution before t-test.
}
\value{
  A data frame with normalized RNA-Seq expression profiles.
  
  Rows stand for tissue names and columns stand for sample names.
 
}
\references{Pei G., Dai Y., Zhao Z., Jia P. (2019) deTS: Tissue-Specific Enrichment Analysis to decode tissue specificity. Bioinformatics, In submission.
}
\author{Guangsheng Pei
}
\note{nothing
}

\seealso{https://github.com/bsml320/deTS
}
\examples{
data(query_GTEx)
query_matrix = query_GTEx[,1:2]
data(correction_factor)
query_mat_zscore_nor = tsea.expression.normalization(query_matrix, 
	correction_factor, normalization = "z-score")
}
\keyword{ "Normalization" }
\keyword{ "Profile" }
