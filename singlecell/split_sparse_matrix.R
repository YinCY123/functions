split_sparse_matrix <- function(mtx, group){
    splice_indices <- split(seq_len(ncol(mtx)), group)
    lapply(splice_indices, function(idx){mtx[, idx, drop = F]})
}