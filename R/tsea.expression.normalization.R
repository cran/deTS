tsea.expression.normalization <-
function (query_mat, correction_factor, normalization = "abundance") {
query_mat_matched = query_mat[which(rownames(query_mat) %in% rownames(correction_factor)),]
query_mat_matched_sd = correction_factor[which(rownames(correction_factor) %in% rownames(query_mat_matched)),]
query_mat_matched_sd = query_mat_matched_sd[match(rownames(query_mat_matched), rownames(query_mat_matched_sd)),]
query_mat_matched_filter = query_mat_matched[which(as.numeric(as.vector(query_mat_matched_sd[, 2])) > 0),]
query_mat_matched_sd_filter = query_mat_matched_sd[which(as.numeric(as.vector(query_mat_matched_sd[, 2])) > 0),]
if (normalization == "abundance"){
query_mat_normalized = log2((query_mat_matched_filter) + 1)/(log2(as.numeric(as.vector(query_mat_matched_sd_filter[, 1])) + 1) + 1)
}
if (normalization == "z-score"){
query_mat_normalized = (query_mat_matched_filter - as.numeric(as.vector(query_mat_matched_sd_filter[, 1])))/as.numeric(as.vector(query_mat_matched_sd_filter[, 2]))
}
return(query_mat_normalized)
}
