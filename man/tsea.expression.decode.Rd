\name{tsea.expression.decode}
\alias{tsea.expression.decode}
\title{Tissue-specific enrichment analysis for RNA-Seq expression profiles
}
\description{Tissue-specific enrichment analysis to decode whether a given RNA-seq sample (RPKM) with potential confounding effects based on expression profiles.
}
\usage{
tsea.expression.decode(query_mat_normalized_score, score, 
ratio = 0.05, p.adjust.method = "BH")
}
\arguments{
  \item{query_mat_normalized_score}{a normalized RNA-seq RPKM object, which produced by "tsea.expression.normalization".
}
  \item{score}{a gene tissue-specific score matrix, c("GTEx_t_score" or "ENCODE_z_score"), can be loaded by data(GTEx) or data(ENCODE), the default value is recommended "GTEx_t_score".
}
  \item{ratio}{the threshold to define tissue-specific genes (with top t-score or z-score), the default value is 0.05.
}
  \item{p.adjust.method}{p.adjust.method, c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
}
}
\details{Tissue-specific enrichment analysis for RNA-Seq expression profiles.
}
\value{
  A data frame with p-value of tissue-specific enrichment result for RNA-Seq expression profiles.
  
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
data(ENCODE_z_score)
query_mat_zscore_nor = tsea.expression.normalization(query_matrix, 
	correction_factor, normalization = "z-score")
tseaed_in_ENCODE = tsea.expression.decode(query_mat_zscore_nor, 
	ENCODE_z_score, 0.05, p.adjust.method = "BH")
}
\keyword{ "RNAseq" }
\keyword{ "Expression" }
