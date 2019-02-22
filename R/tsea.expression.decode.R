tsea.expression.decode <-
function(query_mat_normalized_score, score, ratio = 0.05, p.adjust.method = "BH"){
query.tsea_t.mat = matrix(1, nrow = ncol(score), ncol = ncol(query_mat_normalized_score));
rownames(query.tsea_t.mat) = colnames(score);
colnames(query.tsea_t.mat) = colnames(query_mat_normalized_score);
for(k.tissue in 1:ncol(score)){
which(as.numeric(as.vector(score[,k.tissue])) > quantile(as.numeric(as.vector(score[,k.tissue])), probs = (1 - ratio), na.rm = TRUE )) -> ii
genes.for.test = rownames(score)[ii] 
match(genes.for.test,rownames(query_mat_normalized_score)) -> idx
idx = idx[!is.na(idx)]
for(k.query in 1:ncol(query_mat_normalized_score)){
p1 = t.test(query_mat_normalized_score[idx, k.query], query_mat_normalized_score[-idx, k.query ], alternative= "greater")$p.value
query.tsea_t.mat[k.tissue, k.query] = p1
}
query.tsea_t.mat[k.tissue,] = p.adjust(query.tsea_t.mat[k.tissue,], p.adjust.method)
cat(".", sep="")
}
return(query.tsea_t.mat)
}
