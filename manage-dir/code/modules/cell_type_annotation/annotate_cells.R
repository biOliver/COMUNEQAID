# Functions for cell type annotation.
source('sctype_scoring.R')

select_cutoff = function(x, q = 1, G = 2){
  km <- mclust::Mclust(data = x, G = 2:G, verbose = F)
  cl <- km$classification
  cl_center = km$parameters$mean
  high_cl <- which(cl_center == max(cl_center))
  cutoff <- quantile(x[cl != high_cl], q, na.rm = T)
  return(cutoff)
}

filter_low_score = function(target_cell_type, es.data, q = 1, G = 3, lable_name = 'sctype_label'){
  es.data.target_cells = es.data[,target_cell_type]
  names(es.data.target_cells) = rownames(es.data)
  # Model cell-type score distribution as a mixture of multiple Gaussian (background + true signal).
  # Sometimes it will be more than two, because some cell types are similar to each other, thus produce an extra peak.
  cutoff = select_cutoff(es.data.target_cells, q = q, G = G)
  message(paste0('Threshold for ',target_cell_type, ' is ', round(cutoff, 2), '.'))
  es.data.filtered = es.data
  filtered_cells = which(es.data.target_cells < cutoff & es.data$sctype_label == target_cell_type)
  message(paste0('Filtered ',length(filtered_cells), ' cells.'))
  # For cells with score below the cut-off, set their cell-type to unknown.
  es.data.filtered[names(es.data.target_cells)[filtered_cells],lable_name] = 'Uncertain'
  return(es.data.filtered)
}

cell_level_annotation = function(tmp_object, marker_db.list, exp.var = NA, assay = 'SCT', filter.var = F, label = 'sctype_label'){
  es.max = sctype_score(scRNAseqData = tmp_object[[assay]]@scale.data, exp.var = exp.var, # Using scale normalized data
                        filter.var = filter.var, p = .1, scaling = F, gs = marker_db.list, gene_names_to_uppercase = F) 
  es.max.data = data.frame(t(es.max))
  es.max.data$sctype_label = apply(es.max.data,1,function(x) names(x)[which.max(x)])
  # Setting cells with low annotation score to unknown.
  es.max.data.filtered = es.max.data
  for (x in names(table(es.max.data$sctype_label))) {
    es.max.data.filtered = filter_low_score(x, es.max.data.filtered, q = 1, G = 3)
  }
  tmp_object@meta.data[[label]] = es.max.data$sctype_label
  tmp_object@meta.data[[paste(label,'filtered',sep = '.')]] = es.max.data.filtered$sctype_label
  return(tmp_object)
}

cluster_level_annotation = function(tmp_object, sctype_label = 'sctype_label.filtered') {
  cL_resutls = do.call("rbind", lapply(levels(tmp_object@meta.data$seurat_clusters), function(cl){
    target_cells = which(tmp_object@meta.data$seurat_clusters == cl)
    es.max.cl = sapply(names(table(tmp_object[[sctype_label]])), 
                       function(x) sum(tmp_object[[sctype_label]][target_cells,] == x))
    c(es.max.cl,n_cell = length(target_cells))
  }))
  cL_resutls = data.frame(cL_resutls)
  cL_resutls$cluster = levels(tmp_object@meta.data$seurat_clusters)
  n.cell.type = dim(cL_resutls)[2] - 2
  cL_resutls[,c(1:n.cell.type)] = cL_resutls[,c(1:n.cell.type)]/cL_resutls$n_cell
  cL_resutls$sctype_label = apply(cL_resutls[,c(1:n.cell.type)],1,function(x) names(x)[which.max(x)])
  sctype_label.cluster = paste(sctype_label, 'cluster',sep = '.')
  tmp_object@meta.data[[sctype_label.cluster]] = ""
  for (j in unique(cL_resutls$cluster)) {
    cl_type = cL_resutls[which(cL_resutls$cluster == j),'sctype_label']; 
    tmp_object@meta.data[[sctype_label.cluster]][which(tmp_object@meta.data$seurat_clusters == j)] = as.character(cl_type)
  }
  return(tmp_object)
}

cell_annotation_ref_based = function(seurat_object, marker_db.list, target_cells, assay = 'SCT', 
                                     filter.var = F, filter = T, resolution = 1, n_pc = 15, label = 'sctype_label') {
  tmp_object = subset(seurat_object, cells = target_cells)
  # Define cell clusters ####
  tmp_object <- FindNeighbors(tmp_object, dims = 1:n_pc, verbose = T)
  tmp_object <- FindClusters(tmp_object, resolution = resolution) 
  # For high level annotation, a low resolution + low PC is preferable.
  # For low level annotation, a higher resolution might be better
  if (filter.var) {
    sampled_cells = unlist(tapply(colnames(tmp_object), INDEX = tmp_object$seurat_clusters, sample, size = 100, replace = T))
    exp.var = apply(tmp_object[[assay]]@data[,sampled_cells],1, var)
  } else {exp.var = NA}
  tmp_object = cell_level_annotation(tmp_object, marker_db.list, filter.var = filter.var, exp.var = exp.var, assay = assay, label = label)
  if (filter) {sctype_label = paste(label,'filtered', sep = '.')
  } else {sctype_label = label
  }
  tmp_object = cluster_level_annotation(tmp_object, sctype_label = sctype_label)
  return(tmp_object)
}

