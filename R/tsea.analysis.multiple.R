tsea.analysis.multiple <-
function(query_gene_list, score, ratio = 0.05, p.adjust.method = "BH"){
if (!(mode(query_gene_list) == "numeric" || mode(query_gene_list) == "list")){
print ("Please input a gene symbol list 0~1 matrix!")
stop;
}
if (any(colSums(query_gene_list) < 20)){
print("For sample:")
print(names(which(colSums(query_gene_list) < 20)))
print("At least 20 genes are needed for tissue-specific enrichment analysis, please remove these sample then retry!")
quit;
}else{
query_gene_list_matched = query_gene_list[which(rownames(query_gene_list) %in% intersect(rownames(query_gene_list), rownames(score))),]
query.tsea_FET.mat = matrix(1, nrow = ncol(score), ncol = ncol(query_gene_list_matched));
rownames(query.tsea_FET.mat) = colnames(score)
colnames(query.tsea_FET.mat) = colnames(query_gene_list_matched);
for(k.tissue in 1:ncol(score)){
which(as.numeric(as.vector(score[,k.tissue])) > quantile(as.numeric(as.vector(score[,k.tissue])), probs = (1 - ratio), na.rm =TRUE )) -> ii
genes.for.test = rownames(score)[ii]
for(k.query in 1:ncol(query_gene_list)){
query_gene_list_matched_tissue = rownames(query_gene_list_matched[which(query_gene_list_matched[,k.query] == "1"),])
common = length(intersect(query_gene_list_matched_tissue, genes.for.test))
unique_query = length(setdiff(query_gene_list_matched_tissue, genes.for.test))
unique_tissue = length(setdiff(genes.for.test, query_gene_list_matched_tissue))
remain = nrow(score) - common - unique_query - unique_tissue
fisher_test = fisher.test(matrix(c(common, unique_query, unique_tissue, remain), nrow = 2), alternative = "greater")$p.value
query.tsea_FET.mat[k.tissue, k.query] = fisher_test
}
query.tsea_FET.mat[k.tissue,] = p.adjust(query.tsea_FET.mat[k.tissue,], p.adjust.method)
cat(".", sep = "")
}
}
return(query.tsea_FET.mat)
}
