\name{tsea.plot}
\alias{tsea.plot}
\title{Tissue-specific enrichment analysis result heatmap plot
}
\description{Heat map plot for tissue-specific enrichment analysis result.
}
\usage{
tsea.plot(tsea_result, threshold = 0.05)
}
\arguments{
  \item{tsea_result}{the result of tissue-specific enrichment analysis, which produced by "tsea.analysis", "tsea.analysis.multiple" or "tsea.expression.decode".
}
  \item{threshold}{the p-value threshold to define if the gene list or RNA-seq data enriched in a given tissue, p-value greater than threshold will not be labeled in the plot. The default value is 0.05.
}
}
\details{Heat map plot for tissue-specific enrichment analysis result 
}
\value{
  Heatmap plot

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
data(GWAS_gene_multiple)
data(GTEx_t_score)
query_gene_list = GWAS_gene_multiple
tsea_t_multi = tsea.analysis.multiple(query_gene_list, 
	GTEx_t_score, 0.05, p.adjust.method = "BH")
tsea.plot(tsea_t_multi, 0.05)
}
\keyword{ "Heatmap" }
\keyword{ "Enrichment" }
