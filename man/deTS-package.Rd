\name{deTS-package}
\alias{deTS-package}
\alias{deTS}
\docType{package}
\title{Tissue-Specific Enrichment Analysis
\packageTitle{deTS}
}
\description{Tissue-specific enrichment analysis to assess lists of candidate genes and tissue-specific expression decode analysis for RNA-seq data to decode RNA expression matrices tissue heterogeneity.
}
\details{Since disease and physiological condition are often associated with a specific tissue, understanding the tissue-specific genes (TSG) expression patterns will substantially reduce false discoveries in biomedical research. However, due to cell complexity in human system, heterogeneous tissues are frequently collected. Making it difficult to distinguish gene expression variability and mislead result interpretation. Here, we present deTS, an R package that conducts Tissue-Specific Enrichment Analysis (TSEA) using two built-in reference panels: the Genotype-Tissue Expression (GTEx) data and the ENCyclopedia Of DNA Elements (ENCODE) data. We implemented two major functions in TSEA to assess lists of candidate genes or expression matrices.

The DESCRIPTION file:
\packageDESCRIPTION{deTS}
\packageIndices{deTS}
}
\author{
\packageAuthor{deTS}

Maintainer: Guangsheng Pei
}
\references{Pei G., Dai Y., Zhao Z., Jia P. (2019) deTS: Tissue-Specific Enrichment Analysis to decode tissue specificity. Bioinformatics, In submission.
}
\keyword{ deTS; TSEA }
\seealso{https://github.com/bsml320/deTS
}
\examples{
data(GTEx_t_score)
data(ENCODE_z_score)
library(pheatmap)

data(GWAS_gene)
query_gene_list = GWAS_gene
tsea_t = tsea.analysis(query_gene_list, GTEx_t_score, 0.05,
	p.adjust.method = "bonferroni")
tsea_t_summary = tsea.summary(tsea_t)

data(GWAS_gene_multiple)
query_gene_list = GWAS_gene_multiple[,1:2]
tsea_t_multi = tsea.analysis.multiple(query_gene_list, 
	GTEx_t_score, 0.05, p.adjust.method = "BH")

data(query_GTEx)
query_matrix = query_GTEx[,1:2]
data(correction_factor)
query_mat_zscore_nor = tsea.expression.normalization(query_matrix, 
	correction_factor, normalization = "z-score")
tseaed_in_ENCODE = tsea.expression.decode(query_mat_zscore_nor, 
	ENCODE_z_score, 0.05, p.adjust.method = "BH")
tseaed_in_ENCODE_summary = tsea.summary(tseaed_in_ENCODE)
}
