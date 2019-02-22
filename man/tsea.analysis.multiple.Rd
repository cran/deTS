\name{tsea.analysis.multiple}
\alias{tsea.analysis.multiple}
\title{Tissue-specific enrichment analysis for multi query gene lists
}
\description{Tissue-specific enrichment analysis by Fisher's Exact Test for multiple gene list.
}
\usage{
tsea.analysis.multiple(query_gene_list, score, ratio = 0.05, 
p.adjust.method = "BH")
}
\arguments{
  \item{query_gene_list}{a 0~1 gene~sample table object, row should be gene symbol, column should be sample name. In the table, gene labeled with 1 indicated it is target gene for a given sample, while 0 indicated it is not target in a given sample.
}
  \item{score}{a gene tissue-specific score matrix, c("GTEx_t_score" or "ENCODE_z_score"), can be loaded by data(GTEx) or data(ENCODE), the default value is recommended "GTEx_t_score".
}
  \item{ratio}{the threshold to define tissue-specific genes (with top t-score or z-score), the default value is 0.05.
}
  \item{p.adjust.method}{p.adjust.method, c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
}
}
\details{Tissue-specific enrichment analysis by Fisher's Exact Test for multiple gene list.
}
\value{
  A data frame with p-value of tissue-specific enrichment result.
  
  Rows stand for tissue names and columns stand for sample names.
  
}
\references{Pei G., Dai Y., Zhao Z. Jia P. (2019) deTS: Tissue-Specific Enrichment Analysis to decode tissue specificity. Bioinformatics, In submission.
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
}
\keyword{ "Multiple" }
\keyword{ "deTS" }
