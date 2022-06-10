# Define function to extract marker genes for each annotated cell type in the reference data set.
library(foreach)
library(doParallel)
library(mclust)

rank_markers = function(exp_data, cell_anotations, target_cell_type, test_fun = t.test, sample_size = 500, stratified = T) {
  target_cells = which(cell_anotations == target_cell_type)
  background_cells = cell_anotations[which(cell_anotations != target_cell_type)]
  if (stratified == T) {
    # Use stratified sampling to make sure the background is not dominant by major cell types.
    stratified_background = unlist(tapply(names(background_cells), 
                                          background_cells, sample, size = sample_size, replace = T), use.names = F)
  } else {
    stratified_background = sample(x = names(background_cells), size = sample_size, replace = T)
  }
  ranking.stat = apply(exp_data, 1, FUN = function(x) {
    test_fun(x[target_cells], x[stratified_background],alternative = 'g')$statistic})
  return(ranking.stat)
}

build_marker_db = function(exp_data, cell_anotations, test_fun = t.test, sample_size = 500, n_core = 10, stratified = T) {
  registerDoParallel(n_core)
  cell_type_list = names(table(cell_anotations))
  out_db = data.frame(foreach(x = cell_type_list, .combine = cbind) %dopar% {
    rank_markers(exp_data = exp_data, 
                 cell_anotations = cell_anotations,
                 target_cell_type = x, 
                 test_fun = test_fun,
                 stratified = stratified,
                 sample_size = sample_size)})
  colnames(out_db) = cell_type_list
  return(out_db)
}

contruct_marker_list = function(marker_db, top_n = 30) {
  # Construct a cell marker annotation list (sorted, top N):
  marker_db.list = list()
  for (cell_type in colnames(marker_db)) { 
    tmp_genes = rownames(marker_db)[order(marker_db[,cell_type], decreasing = T)[c(1:top_n)]]
    sapply(names(marker_db.list), function(x){
      if (length(union(tmp_genes,marker_db.list[[x]])) == top_n) {
        warning(paste(cell_type, 'and', x, 'share the same set of markers'), immediate. = T)
      }
    })
    marker_db.list[[cell_type]] = tmp_genes
  }
  return(marker_db.list)
}

