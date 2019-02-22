tsea.analysis <-
function(query_gene_list, score, ratio = 0.05, p.adjust.method = "BH"){
if (mode(query_gene_list)!="character"){
print ("Please input a gene symbol list!")
stop;
}
if (length(query_gene_list) < 20){
print("At least 20 genes are needed for tissue-specific enrichment analysis!")
stop;
}else{
query_list_gene_matched = intersect(query_gene_list, rownames(score))
query.tsea_FET.mat = matrix(1, nrow = ncol(score), ncol = 1);
rownames(query.tsea_FET.mat) = colnames(score)
colnames(query.tsea_FET.mat) = "query";
for(k.tissue in 1:ncol(score)){
which(as.numeric(as.vector(score[,k.tissue])) > quantile(as.numeric(as.vector(score[,k.tissue])), probs = (1 - ratio), na.rm = TRUE )) -> ii
genes.for.test = rownames(score)[ii]
common = length(intersect(query_list_gene_matched, genes.for.test))
unique_query = length(setdiff(query_list_gene_matched, genes.for.test))
unique_tissue = length(setdiff(genes.for.test, query_list_gene_matched))
remain = nrow(score) - common - unique_query - unique_tissue
fisher_test = fisher.test(matrix(c(common, unique_query, unique_tissue, remain), nrow=2), alternative = "greater")$p.value
query.tsea_FET.mat[k.tissue, 1] = fisher_test
cat(".", sep = "")
}
}
query.tsea_FET.mat[,1] = p.adjust(query.tsea_FET.mat[,1], p.adjust.method)
return(query.tsea_FET.mat)
}
