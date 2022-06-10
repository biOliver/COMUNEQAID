# GNU General Public License v3.0 (https://github.com/IanevskiAleksandr/sc-type/blob/master/LICENSE)
# Written by Aleksandr Ianevski <aleksandr.ianevski@helsinki.fi>, June 2021
# See https://www.nature.com/articles/s41467-022-28803-w
# Modified by Bitao Qiu.
# Functions on this page:
# sctype_score: calculate ScType scores and assign cell types
#
# @params: scRNAseqData - input scRNA-seq matrix (rownames - genes, column names - cells), 
# @params: scale - indicates whether the matrix is scaled (TRUE by default)
# @params: gs - list of gene sets positively expressed in the cell type 
# @params: gs2 - list of gene sets that should not be expressed in the cell type (NULL if not applicable)
mus.to.hsa.id <- read.csv('mus_hsa_id.csv',row.names = 1)
Hsa_2_Mus_id <- function(x, available.id){
  tmp.id = mus.to.hsa.id$geneID.mus[which(mus.to.hsa.id$geneID.hsa %in% x)]
  tmp.id.filtered = tmp.id[which(tmp.id %in% available.id)]
  return(tmp.id.filtered)
}

filter_genes = function(x, exp.var, p = .9){
  high.var.genes = exp.var[which(exp.var > quantile(exp.var, p))]
  return(x[which(x %in% names(high.var.genes))])
}

cal_marker_sensitivity <- function(exp.var, filter.var = F, p = .1, gs, gene_names_to_uppercase = T, ...){
  marker_stat = sort(table(unlist(gs)), decreasing = T)    
  if (gene_names_to_uppercase) {
    names(exp.var) = toupper(names(exp.var))
  }
  
  # $BQ: We filter genes with little expression variation among cells, because these genes are not informative 
  if (filter.var) {
    
    # $BQ: Finding genes with high expression variation.
    var.genes = filter_genes(names(exp.var),exp.var = exp.var, p = p) 
    # $BQ: Only retain genes that are informative.
    marker_stat.filtered = marker_stat[which(names(marker_stat) %in% var.genes)] 
    # $BQ: Only retain cell types with at least one marker gene.
    gs.filtered = lapply(gs, function(x) x[which(x %in% var.genes)])  
    gs.filtered = gs.filtered[names(gs.filtered)[unlist(lapply(gs.filtered, function(x) length(x) > 0))]]
  } else {
    gs.filtered = gs
    marker_stat.filtered = marker_stat
  }
  # $BQ: The following examines how often genes are associated with cell types (sensitivity).
  marker_sensitivity = data.frame(score_marker_sensitivity = scales::rescale(as.numeric(marker_stat.filtered), to = c(0,1), from = c(length(gs.filtered),1)),
                                  gene_ = names(marker_stat.filtered), stringsAsFactors = !1) 
  
  # $BQ: The author did not follow their method description. With the above formula, the minimum != 0
  # $BQ: If a gene only expressed in one cell type, sensitivity = 1; If expressed in all cell types (n = length(gs)), sensitivity = 0
  # $BQ: Note that the number of input cell type will influence the sensitivity scores.
  return(list(gs.filtered = gs.filtered, marker_sensitivity = marker_sensitivity))
}

sctype_score <- function(scRNAseqData, exp.var, filter.var = F, p = .1, scaling = F, gs, gs2 = NULL, gene_names_to_uppercase = T, ...){
  # marker sensitivity
  sensitivity.list = cal_marker_sensitivity(exp.var, filter.var = filter.var, p = p, gs, gs2 = gs2, gene_names_to_uppercase = gene_names_to_uppercase)
  gs.filtered = sensitivity.list[[1]]
  marker_sensitivity = sensitivity.list[[2]]
  # convert gene names to Uppercase
  if (gene_names_to_uppercase) {
    rownames(scRNAseqData) = toupper(rownames(scRNAseqData));
  }
  # Select genes only found in data
  names_gs_cp = names(gs.filtered); 
  names_gs_2_cp = names(gs2);
  gs.filtered = lapply(1:length(gs.filtered), function(d_){ 
    GeneIndToKeep = rownames(scRNAseqData) %in% as.character(gs.filtered[[d_]]); 
    rownames(scRNAseqData)[GeneIndToKeep]})
  gs2 = lapply(1:length(gs2), function(d_){ 
    GeneIndToKeep = rownames(scRNAseqData) %in% as.character(gs2[[d_]]); rownames(scRNAseqData)[GeneIndToKeep]})
  names(gs.filtered) = names_gs_cp; 
  gs.filtered = gs.filtered[names(gs.filtered)[unlist(lapply(gs.filtered, function(x) length(x) > 0))]]
  names(gs2) = names_gs_2_cp;
  cell_markers_genes_score = marker_sensitivity[marker_sensitivity$gene_ %in% unique(unlist(gs.filtered)),]
  
  if (scaling) Z <- t(scale(t(scRNAseqData))) else Z <- scRNAseqData
  # subselect only with marker genes
  Z = Z[unique(c(unlist(gs.filtered),unlist(gs2))), ]
  # multiple by marker sensitivity
  # $BQ: Using multiplication for mulitple genes.
  Z[cell_markers_genes_score$gene_, ] = Z[cell_markers_genes_score$gene_, ]*cell_markers_genes_score$score_marker_sensitivity
  # combine scores
  es = do.call('rbind', lapply(names(gs.filtered), function(gss_){ 
    gs_z = Z[gs.filtered[[gss_]],]
    if (length(gs.filtered[[gss_]]) > 1) {
      # $BQ: Using square root will give cell types with more markers a higher score.
      sum_t1 = colSums(gs_z)/sqrt(length(gs.filtered[[gss_]]))
    } else {
      sum_t1 = gs_z
    }
    if (!is.null(gs2[[gss_]])) {
      gs2_z = Z[gs2[[gss_]], ] * -1
      if (length(gs2[[gss_]]) > 1) {
        sum_t2 = colSums(gs2_z)/(length(gs2[[gss_]]))
      } else {
        sum_t2 = sum(gs2_z)
      }} else {
        sum_t2 = 0
      }
    sum_t1 + sum_t2 
  }))
  dimnames(es) = list(names(gs.filtered), colnames(Z))
  es.max <- es[!apply(is.na(es) | es == "", 1, all),] # remove na rows
  return(data.frame(es.max))
}
