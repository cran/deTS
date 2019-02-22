tsea.plot <-
function(tsea_result, threshold = 0.05){
tsea_result[tsea_result > threshold] <- 1
query.tsea_t.log = -log(tsea_result,10)
query.tsea_t.label = query.tsea_t.log
for (i in 1:ncol(query.tsea_t.log)){
top1 = order(query.tsea_t.log[, i],decreasing=TRUE)[1]
top2 = order(query.tsea_t.log[, i],decreasing=TRUE)[2]
top3 = order(query.tsea_t.log[, i],decreasing=TRUE)[3]
query.tsea_t.label[top1, i]="1"
query.tsea_t.label[top2, i]="2"
query.tsea_t.label[top3, i]="3"
}
query.tsea_t.label[query.tsea_t.label != "1" & query.tsea_t.label != "2" & query.tsea_t.label != "3"] <-""
query.tsea_t.log[query.tsea_t.log == "Inf"] <- max(query.tsea_t.log[query.tsea_t.log != "Inf"])
query.tsea_t.label[tsea_result > threshold] <- ""
pheatmap(query.tsea_t.log, cluster_rows = F, cluster_cols = F, display_numbers = query.tsea_t.label, color = colorRampPalette(c("white","firebrick1"))(50), main = "-log10 p-value")
}
