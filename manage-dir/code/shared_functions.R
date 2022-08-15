################################################################################
##########                          Functions                         ##########
################################################################################
'%!in%' <- function(x,y)!('%in%'(x,y))

out_put.cex <- 1.5
par(cex = out_put.cex)

load_fry <- function(frydir, which_counts = c('U','S','A'), verbose = FALSE, output_list = F) {
  # read in metadata
  meta_info = rjson::fromJSON(file = file.path(frydir, 'meta_info.json'))
  ng = meta_info[['num_genes']]
  usa_mode = meta_info[['usa_mode']]
  
  if (usa_mode) {
    if (length(which_counts) == 0) {
      stop('Please at least provide one status in \'U\' \'S\' \'A\' ')
    }
    if (verbose) {
      message('processing input in USA mode, will return ', paste(which_counts, collapse = '+'))
    }
  } else if (verbose) {
    message('processing input in standard mode, will return spliced count')
  }
  
  # read in count matrix
  af_raw = readMM(file = file.path(frydir, 'alevin', 'quants_mat.mtx'))
  
  # if usa mode, each gene gets 3 rows, so ng/3
  if (usa_mode) {
    ng = as.integer(ng/3)
  }
  
  # read in gene name file and cell barcode file
  afg = read.csv(file.path(frydir, 'alevin', 'quants_mat_cols.txt'), strip.white = TRUE, header = FALSE, nrows = ng, col.names = c('gene_ids'))
  afc = read.csv(file.path(frydir, 'alevin', 'quants_mat_rows.txt'), strip.white = TRUE, header = FALSE, col.names = c('barcodes'))
  
  # if in usa_mode, sum up counts in different status according to which_counts
  if (output_list) {
    which_counts = c('U','S','A')
    if (usa_mode) {
      rd = list('S' = seq(1, ng), 'U' =  seq(ng + 1, 2*ng), 'A' =  seq(2*ng + 1, 3*ng))
      o = af_raw[, rd[[which_counts[1]]]]
      for (wc in which_counts[-1]) {
        o = o + af_raw[, rd[[wc]]]
      }
    } else {
      o = af_raw
    }
    
    res.rna <- t(o)
    colnames(res.rna) <- afc[['barcodes']]
    rownames(res.rna) <- afg[['gene_ids']]
    
    which_counts = c('S','A')
    if (usa_mode) {
      rd = list('S' = seq(1, ng), 'U' =  seq(ng + 1, 2*ng), 'A' =  seq(2*ng + 1, 3*ng))
      o = af_raw[, rd[[which_counts[1]]]]
      for (wc in which_counts[-1]) {
        o = o + af_raw[, rd[[wc]]]
      }
    } else {
      o = af_raw
    }
    
    res.spl <- t(o)
    colnames(res.spl) <- afc[['barcodes']]
    rownames(res.spl) <- afg[['gene_ids']]
    
    which_counts = c('U')
    if (usa_mode) {
      rd = list('S' = seq(1, ng), 'U' =  seq(ng + 1, 2*ng), 'A' =  seq(2*ng + 1, 3*ng))
      o = af_raw[, rd[[which_counts[1]]]]
      for (wc in which_counts[-1]) {
        o = o + af_raw[, rd[[wc]]]
      }
    } else {
      o = af_raw
    }
    
    res.uns <- t(o)
    colnames(res.uns) <- afc[['barcodes']]
    rownames(res.uns) <- afg[['gene_ids']]
    
    return(list('RNA'         = res.rna,
                'Spliced'     = res.spl,
                'Unspliced'   = res.uns))
  }else {
    if (usa_mode) {
      rd = list('S' = seq(1, ng), 'U' =  seq(ng + 1, 2*ng), 'A' =  seq(2*ng + 1, 3*ng))
      o = af_raw[, rd[[which_counts[1]]]]
      for (wc in which_counts[-1]) {
        o = o + af_raw[, rd[[wc]]]
      }
    } else {
      o = af_raw
    }
    res.counts <- t(o)
    colnames(res.counts) <- afc[['barcodes']]
    rownames(res.counts) <- afg[['gene_ids']]
    
    return(res.counts)
  }
}

annotateMat <- function(mat, orga){
  
  if (orga == 'Mouse') {
    require(org.Mm.eg.db, quietly = T)
    orgDb <- org.Mm.eg.db
  }
  if (orga == 'Human') {
    require(org.Hs.eg.db, quietly = T)
    orgDb <- org.Hs.eg.db
  }
  if (orga == 'Rat') {
    require(org.Rn.eg.db, quietly = T)
    orgDb <- org.Rn.eg.db
  }
  if (orga == 'Rhesus') {
    require(org.Mmu.eg.db, quietly = T)
    orgDb <- org.Mmu.eg.db
  }
  
  # Annotate genes
  ann <- suppressMessages(data.frame(mapIds(orgDb,
                                            sapply(strsplit(rownames(mat), '[.]'), '[', 1),
                                            column = 'SYMBOL', keytype = 'ENSEMBL')))
  ann[['ens']] <- as.character(rownames(ann))
  colnames(ann)[1] <- 'gene'
  ann[['gene']] <- as.character(ann[['gene']])
  converted <- ifelse(is.na(ann[['gene']]), ann[['ens']], ann[['gene']])
  rownames(mat) <- converted
  
  return(mat)
}

emptyDropsmodal <- function(q, verbose = T, plot = T, format = 'save', skipModCheck = F){
  if (verbose) {
    cat('#\tComputing knee & inflection from barcode ranks..\n',
        sep = '')
  }
  sum.col <- Matrix::colSums(q)
  bc.calls <- barcodeRanks(q)
  bcs.knee <- names(sum.col)[sum.col > metadata(bc.calls)[['knee']]]
  bcs.infl <- names(sum.col)[sum.col > metadata(bc.calls)[['inflection']]]
  if (verbose) {
    cat('#\t-\tcell barcodes above knee\t=\t',length(bcs.knee),'\n',
        '#\t-\tcell barcodes above inflection\t=\t',length(bcs.infl),'\n',
        sep = '')
  }
  
  if (skipModCheck == F) {
    cat('#\t-\tPerforming Hartigan\'s dip test on cell barcodes above knee..\n',
        '#\t-\thypothesis: data is unimodal\n',
        sep = '')
    
    dip <- diptest::dip.test(sum.col[bcs.knee])
  }else {
    dip <- diptest::dip.test(sum.col[bcs.knee])
    dip[['p.value']] <- 1
  }
  
  if (verbose) {
    cat('#\t-\t\tP value = ',dip[['p.value']],'\n',
        sep = '')
  }
  
  if (dip[['p.value']] < 0.05) {
    if (verbose) {
      cat('#\t-\tData is multimodal..\n',
          '#\t-\t-\tfitting mixture model on cell barcodes above inflection..\n',
          '#\t-\t-\tfitting to: bimodal\n',
          sep = '')
    }
    
    tmp.fit <- mclust::Mclust(data = sum.col[bcs.infl], G = 1:2, verbose = F)
    bcs.infl <- bcs.infl[tmp.fit[['classification']] == 2]
    
    if (verbose) {
      cat('#\t-\t-\t\tmean\n',
          '#\t-\t-\tmod 1\t',round(tmp.fit[['parameters']][['mean']][[1]]),'\n',
          '#\t-\t-\tmod 2\t',round(tmp.fit[['parameters']][['mean']][[2]]),'\n',
          '#\t-\t-\tretaining ',length(bcs.infl),' cell barcodes assigned to mod 2','\n',
          sep = '')
      
      cat('#\t-\t-\tcomputing knee & inflection from barcode ranks..\n',
          sep = '')
    }
    bc.calls.new <- barcodeRanks(q[,bcs.infl])
    bcs.knee.new <- names(sum.col)[sum.col > metadata(bc.calls.new)[['knee']]]
    bcs.infl.new <- names(sum.col)[sum.col > metadata(bc.calls.new)[['inflection']]]
    if (verbose) {
      cat('#\t-\t-\tcell barcodes above knee\t=\t',length(bcs.knee.new),'\n',
          '#\t-\t-\tcell barcodes above inflection\t=\t',length(bcs.infl.new),'\n',
          '#\t-\t-\tretaining cells above inflection..\n',
          sep = '')
    }
    if (plot) {
      if (format == 'save') {
        png(paste0(dir.outs.qc.plots,'/',tmp.pool[['Index (10x)']],'_barcode-ranks-plot.png'), width = 1500, height = 500)
      }
      par(mfrow = c(1,2), cex = out_put.cex)
      r <- rank(-bc.calls[['total']])
      suppressWarnings(plot(r, bc.calls[['total']], log = 'xy', xlab = 'Rank', ylab = 'Total UMI count', main = ''))
      abline(h = metadata(bc.calls.new)[['knee']], col = '#ff6c00', lty = 2, lwd = 3)
      abline(h = metadata(bc.calls.new)[['inflection']], col = '#04c2c4', lty = 2, lwd = 3)
      abline(h = metadata(bc.calls)[['knee']], col = '#d6aa89', lty = 3, lwd = 3)
      abline(h = metadata(bc.calls)[['inflection']], col = '#9fc9c9', lty = 3, lwd = 3)
      abline(v = 40000, col = '#42e373', lty = 2, lwd = 3)
      abline(v = 80000, col = '#e632b9', lty = 2, lwd = 3)
      legend('bottomleft', bty = "n", lty = c(2,2,3,3), lwd = 3, col = c('#ff6c00', '#04c2c4','#42e373','#e632b9','#d6aa89', '#9fc9c9'), 
             legend = c('knee (updated)', 'inflection (updated)','rank 40k','rank 80k','knee (original)', 'inflection (original)'))
      hist(log10(bc.calls[['total']]), xlab = 'Log[10] UMI count', main = '')
      abline(v = log10(metadata(bc.calls.new)[['knee']]), col = '#ff6c00', lty = 2, lwd = 3)
      abline(v = log10(metadata(bc.calls.new)[['inflection']]), col = '#04c2c4', lty = 2, lwd = 3)
      abline(v = log10(metadata(bc.calls)[['knee']]), col = '#d6aa89', lty = 3, lwd = 3)
      abline(v = log10(metadata(bc.calls)[['inflection']]), col = '#9fc9c9', lty = 3, lwd = 3)
      if (format == 'save') {
        dev.off()
      }
    }
  }else{
    if (verbose) {
      cat('#\t-\tretaining cells above inflection..\n',
          sep = '')
    }
    if (plot) {
      if (format == 'save') {
        png(paste0(dir.outs.qc.plots,'/',tmp.pool[['Index (10x)']],'_barcode-ranks-plot.png'), width = 1500, height = 500)
      }
      par(mfrow = c(1,2), cex = out_put.cex)
      r <- rank(-bc.calls[['total']])
      suppressWarnings(plot(r, bc.calls[['total']], log = 'xy', xlab = 'Rank', ylab = 'Total UMI count', main = ''))
      abline(h = metadata(bc.calls)[['knee']], col = '#ff6c00', lty = 2, lwd = 3)
      abline(h = metadata(bc.calls)[['inflection']], col = '#04c2c4', lty = 2, lwd = 3)
      abline(v = 40000, col = '#42e373', lty = 2, lwd = 3)
      abline(v = 80000, col = '#e632b9', lty = 2, lwd = 3)
      legend('bottomleft', bty = "n", lty = c(2,2), lwd = 3, col = c('#ff6c00', '#04c2c4','#42e373','#e632b9'), 
             legend = c('knee','inflection','rank 40k','rank 80k'))
      hist(log10(bc.calls[['total']]), xlab = 'Log[10] UMI count', main = '')
      abline(v = log10(metadata(bc.calls)[['knee']]), col = '#ff6c00', lty = 2, lwd = 3)
      abline(v = log10(metadata(bc.calls)[['inflection']]), col = '#04c2c4', lty = 2, lwd = 3)
      abline(h = 40000, col = '#42e373', lty = 2, lwd = 3)
      abline(h = 80000, col = '#e632b9', lty = 2, lwd = 3)
      if (format == 'save') {
        dev.off()
      }
    }
  }
  if (verbose) {
    cat('#\n',
        '#\n',
        sep = '')
  }
  return(bcs.infl)
}

emptyDropsManu <- function(q, manualInfl){
  
  cat('#\tComputing knee & inflection from barcode ranks..\n',
      sep = '')
  sum.col <- Matrix::colSums(q)
  bc.calls <- barcodeRanks(q)
  bcs.knee <- names(sum.col)[sum.col > metadata(bc.calls)[['knee']]]
  bcs.infl <- names(sum.col)[sum.col > metadata(bc.calls)[['inflection']]]
  bcs.manu <- names(sum.col)[sum.col > manualInfl]
  cat('#\t-\tcell barcodes above knee\t=\t',length(bcs.knee),'\n',
      '#\t-\tcell barcodes above inflection\t=\t',length(bcs.infl),'\n',
      '#\t-\tcell barcodes above manual TH\t=\t',length(bcs.manu),'\n',
      sep = '')
  
  
  
  cat('#\t-\tretaining cells above inflection..\n',
      sep = '')
  
  png(paste0(dir.outs.qc.plots,'/',tmp.pool[['Index (10x)']],'_barcode-ranks-plot.png'), width = 1500, height = 500)
  par(mfrow = c(1,2), cex = out_put.cex)
  r <- rank(-bc.calls[['total']])
  suppressWarnings(plot(r, bc.calls[['total']], log = 'xy', xlab = 'Rank', ylab = 'Total gene count', main = ''))
  abline(h = metadata(bc.calls)[['knee']], col = '#ff6c00', lty = 2, lwd = 3)
  abline(h = metadata(bc.calls)[['inflection']], col = '#04c2c4', lty = 2, lwd = 3)
  abline(h = manualInfl, col = '#fa55fa', lty = 2, lwd = 3)
  abline(v = 40000, col = '#42e373', lty = 2, lwd = 3)
  abline(v = 80000, col = '#e632b9', lty = 2, lwd = 3)
  legend('bottomleft', bty = "n",lty = c(2,2,2), lwd = 3, col = c('#ff6c00', '#04c2c4','#fa55fa','#42e373','#e632b9'), 
         legend = c('knee','inflection','Manual threshold','rank 40k','rank 80k'))
  hist(log10(bc.calls[['total']]), xlab = 'Log[10] gene count', main = '')
  abline(v = log10(metadata(bc.calls)[['knee']]), col = '#ff6c00', lty = 2, lwd = 3)
  abline(v = log10(metadata(bc.calls)[['inflection']]), col = '#04c2c4', lty = 2, lwd = 3)
  abline(v = log10(manualInfl), col = '#fa55fa', lty = 2, lwd = 3)
  
  dev.off()
  
  return(bcs.manu)
}

is_multimodal <- function(x, p.cutoff = 1e-2) {
  # Test if the expression distribution violates unimodal distribution.
  p = diptest::dip.test(x)$p.value
  return(p < p.cutoff)
}

select_hash_cutoff_mcl <- function(x, q_l = 1, q_h = 0.01, seed = 42) {
  # Model HTO data as a mixture of two Gaussian distributions (for normalized [across cells] data)
  # And select HTO cutoff based on mclust (Model based clustering).
  assertthat::assert_that(class(x) == "numeric")
  assertthat::is.number(seed)
  assertthat::assert_that(length(seed) == 1)
  set.seed(seed)
  km <- mclust::Mclust(data = x, G = 2, verbose = F)
  cl <- km$classification
  cl_center = km$parameters$mean
  high_cl <- which(cl_center == max(cl_center))
  low_cl <- which(cl_center != max(cl_center))
  # q_l and q_h are the quantiles for negative and postive cluster, respectively. 
  cutoff <- max(quantile(x[cl == low_cl], q_l), quantile(x[cl == high_cl], q_h))
  # The higher the cut off, the less false positive (the more false negative).
  return(cutoff)
}

hash_mcl_p <- function(x, seed = 3030, q_l = 1, q_h = 0.001) {
  assertthat::assert_that(class(x) == "numeric")
  assertthat::is.number(seed)
  assertthat::assert_that(length(seed) == 1)
  set.seed(seed)
  km <- mclust::Mclust(data = x, G = 2, verbose = F)
  cl <- km$classification
  cl_center = km$parameters$mean
  high_cl <- which(cl_center == max(cl_center))
  low_cl <- which(cl_center != max(cl_center))
  p.high_cl <- km$z[,high_cl]
  # Correct assignment error from Mclust
  p.high_cl[which(x < max(quantile(cl_center[low_cl], q_l), quantile(cl_center[low_cl], q_h)))] = 0
  names(p.high_cl) = names(x)
  return(p.high_cl)
}

HTO_classifcation = function(discrete, hto_mcl.p, assay){
  # Based on HTODemux (Seurat)
  npositive <- colSums(x = discrete)
  classification.global <- npositive
  classification.global[npositive == 0] <- "Negative"
  classification.global[npositive == 1] <- "Singlet"
  classification.global[npositive > 1] <- "Doublet"
  donor.id = rownames(x = discrete)
  hash.max <- apply(X = hto_mcl.p, MARGIN = 2, FUN = max) # This returns the probability of the most likely HashID (based on the Hashtag distribution among cells)
  hash.maxID <- as.character(donor.id[apply(X = hto_mcl.p, MARGIN = 2, FUN = which.max)])
  hash.second <- apply(X = hto_mcl.p, MARGIN = 2, FUN = function(x) sort(x,decreasing = T)[2])
  hash.secondID <- as.character(donor.id[apply(X = hto_mcl.p, MARGIN = 2, FUN = function(x) order(x,decreasing = T)[2])])
  hash.margin <- hash.max - hash.second
  doublet_id <- sapply(X = 1:length(x = hash.maxID), FUN = function(x) {
    return(paste(sort(x = c(hash.maxID[x], hash.secondID[x])), 
                 collapse = "_"))
  })
  classification <- classification.global
  classification[classification.global == "Negative"] <- "Negative"
  classification[classification.global == "Singlet"] <- hash.maxID[which(x = classification.global == "Singlet")]
  classification[classification.global == "Doublet"] <- doublet_id[which(x = classification.global == "Doublet")]
  classification.metadata <- data.frame(hash.maxID, hash.secondID, hash.margin, classification, classification.global)
  colnames(x = classification.metadata) <- paste(assay, 'mcl', c("maxID", "secondID", "margin", "classification", "classification.global"), sep = "_")
  return(classification.metadata)
}

HTODemux.mcl <- function(object, assay = "HTO", q_l = 1, q_h = 0.005, seed = 42){
  # A function to find the threshold for each hastag, the P, and singlet, doublets and negative.
  # The input is the HTO data matrix (normalized across cells).
  
  assay <- assay %||% DefaultAssay(object = object)
  data <- GetAssayData(object = object, assay = assay, slot = 'data')
  hto_mcl.cutoff = data.frame(cut_off = future.apply::future_apply(data,1,function(x) select_hash_cutoff_mcl(x, q_l = q_l, q_h = q_h), future.seed = T))
  hto_mcl.cutoff$Multi_modal = apply(data,1,function(x) is_multimodal(x))
  print(hto_mcl.cutoff)
  hto_mcl.p = t(apply(data,1,function(x) hash_mcl_p(x, seed = seed, q_l = q_l, q_h = q_h)))
  discrete <- data
  discrete[discrete > 0] <- 0
  for (iter in rownames(x = data)) {
    values <- data[iter, ]
    cutoff <- hto_mcl.cutoff[iter,'cut_off']
    discrete[iter, names(x = which(x = values > cutoff))] <- 1
  }
  classification.metadata <- HTO_classifcation(discrete, hto_mcl.p, assay)
  object <- AddMetaData(object = object, metadata = classification.metadata)
  Idents(object) <- paste(assay, "mcl", "classification", sep = '_')
  doublets <- rownames(x = object[[]])[which(object[[paste(assay, "mcl","classification.global", sep = "_")]] == "Doublet")]
  Idents(object = object, cells = doublets) <- "Doublet"
  Idents(object) = factor(Idents(object), levels = c('Doublet', 'Negative', rownames(object@assays$HTO)))
  object$hash.mcl.ID <- Idents(object = object)
  return(object)
}

HTODemux_mcl.visualization <- function(object, assay = "HTO", q_l = 1, q_h = 0.005, seed = 42){
  assay <- assay %||% DefaultAssay(object = object)
  data <- GetAssayData(object = object, assay = assay, slot = 'data')
  hto.data.wide = data.frame(t(data))
  hto.data.long = data.table::melt(data.table::setDT(hto.data.wide,keep.rownames = T), id.vars = 'rn', value.name = 'Expression',variable.name = 'hto')
  hto_mcl.cutoff = data.frame(cut_off = future.apply::future_apply(data,1,function(x) select_hash_cutoff_mcl(x, q_l = q_l, q_h = q_h), future.seed = T), hto = colnames(hto.data.wide)[-1])
  
  p <- ggplot(hto.data.long, aes(x = Expression)) +
    geom_histogram(bins = 100) +
    geom_vline(data = hto_mcl.cutoff, aes(xintercept = cut_off), col = 'red') +
    facet_wrap(~hto,scales = 'free',ncol = 2) +
    xlab('Expression') +
    ylab('Counts') +
    # scale_y_sqrt() +
    ggtitle('Individual HTO distributions', subtitle = paste0('Q_negative = ', q_l, '; Q_positive = ', q_h)) +
    theme_minimal(base_size = 20)
  
  ggsave(filename = paste0(tmp.pool[['Index (10x)']],'_HTO-hist.png'),
         plot = p,
         path = dir.outs.qc.plots,
         width = 8,
         height = 8)
}

dub_cutoff <- function(x) {
  
  # calculate number of neighbors at each proportion that are doublets
  data.frame('prop' = x[['proportion_dub_neighbors']]) %>%
    group_by(prop) %>% 
    summarize(n = n()) %>% 
    mutate(pct = n/sum(n)) -> data 
  # find point at which we gain very few doublets as proportion increases 
  cut <- data[['prop']][PCAtools::findElbowPoint(variance = sort(data[['n']], decreasing = T)) + 1]
  vec <- if_else(x[['proportion_dub_neighbors']] <= cut, F, T)
  return(vec)
}

makeUmiPlot <- function(feat.RNA, feat.HTO, CBs.called, title = 'Please give this plot a title...'){
  
  common.CBs <- intersect(feat.RNA$CB,feat.HTO$CB)
  
  feat.RNA %>% 
    filter(CB %in% common.CBs) %>% 
    arrange(CB) -> feat.RNA
  
  feat.HTO %>% 
    filter(CB %in% common.CBs) %>% 
    arrange(CB) -> feat.HTO
  
  p <- ggplot(mapping = aes(x = feat.RNA$DeduplicatedReads,
                            y = feat.HTO$DeduplicatedReads,
                            col = ifelse(sort(common.CBs) %in% CBs.called, 'Called', 'Uncalled'))) +
    geom_point() + scale_x_log10() + scale_y_log10() +
    labs(x = 'UMIs (RNA)',
         y = 'UMIs (HTO)',
         col = '') +
    theme_minimal(base_size = 25) + 
    ggtitle(title)
  
  return(p)
}

