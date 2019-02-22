tsea.summary <-
function(tsea_result){
query.tsea = matrix("", nrow = ncol(tsea_result), ncol = 6)
rownames(query.tsea) = colnames(tsea_result)
colnames(query.tsea) = c("Top1", "tissue1_p-value", "Top2", "tissue2_p-value", "Top3", "tissue3_p-value")
query.tsea_t.log = -log(tsea_result,10)
for (i in 1:ncol(query.tsea_t.log)){
top1 = order(query.tsea_t.log[, i], decreasing=TRUE)[1]
top2 = order(query.tsea_t.log[, i], decreasing=TRUE)[2]
top3 = order(query.tsea_t.log[, i], decreasing=TRUE)[3]
query.tsea[i,1]=rownames(query.tsea_t.log)[top1]
query.tsea[i,3]=rownames(query.tsea_t.log)[top2]
query.tsea[i,5]=rownames(query.tsea_t.log)[top3]
query.tsea[i,2]=tsea_result[top1, i]
query.tsea[i,4]=tsea_result[top2, i]
query.tsea[i,6]=tsea_result[top3, i]
}
return(query.tsea)
}
