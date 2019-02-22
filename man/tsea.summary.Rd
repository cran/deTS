\name{tsea.summary}
\alias{tsea.summary}
\title{Tissue-specific enrichment analysis result summary
}
\description{Tissue-specific enrichment analysis result summary (list the top 3 most enriched tissues) from the given gene list or RNA-seq expression profiles.
}
\usage{
tsea.summary(tsea_result)
}
\arguments{
  \item{tsea_result}{the result of tissue-specific enrichment analysis, which produced by "tsea.analysis", "tsea.analysis.multiple" or "tsea.expression.decode".
}
}
\details{Tissue-specific enrichment analysis result summary
}
\value{  
  A data frame with summary result of top 3 most enriched tissues.
  
  Rows stand for sample names and columns stand for top 3 most enriched tissues (with p-value).
 
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
query_matrix = query_GTEx
data(correction_factor)
data(ENCODE_z_score)
query_mat_zscore_nor = tsea.expression.normalization(query_matrix, 
	correction_factor, normalization = "z-score")
tseaed_in_ENCODE = tsea.expression.decode(query_mat_zscore_nor, 
	ENCODE_z_score, 0.05, p.adjust.method = "BH")	
tseaed_in_ENCODE_summary = tsea.summary(tseaed_in_ENCODE)
}
\keyword{ "Top tissue" }
\keyword{ "Summary" }
