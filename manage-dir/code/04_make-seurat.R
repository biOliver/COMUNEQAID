################################################################################
##########                          Packages                          ##########
################################################################################
suppressPackageStartupMessages({
  library(foreach)
  library(parallel)
  library(doParallel)
  library(Seurat)
  library(tximport)       
  library(SingleCellExperiment)
  library(future)
  library(tidyverse)
  library(scater)
  library(scran)
  library(sctransform)
  library(DropletUtils)
  library(Matrix)
  library(mclust)
  library(AnnotationDbi)
})

################################################################################
##########                          Functions                         ##########
################################################################################
'%!in%' <- function(x,y)!('%in%'(x,y))

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
      par(mfrow = c(1,2))
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
      par(mfrow = c(1,2))
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
  par(mfrow = c(1,2))
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

select_hash_cutoff_mcl <- function(x, q = 1, seed = 42) {
  # Model HTO data as a mixture of two gaussian distributions (for normalized [across cells] data)
  # And select HTO cutoff based on mclust (Model based clustering).
  assertthat::assert_that(class(x) == "numeric")
  assertthat::is.number(seed)
  assertthat::assert_that(length(seed) == 1)
  set.seed(seed)
  # For each HTO, the data can be modeled as 
  km <- mclust::Mclust(data = x, G = 2, verbose = F)
  cl <- km$classification
  cl_center = km$parameters$mean
  high_cl <- which(cl_center == max(cl_center))
  low_cl <- which(cl_center != max(cl_center))
  #high_center <- cl_center[high_cl]
  #low_center <- cl_center[low_cl]
  #fc <- high_center / low_center
  #print(fc)
  cutoff <- quantile(x[cl == low_cl], q)
  return(cutoff)
}

hash_mcl_p <- function(x, seed = 3030) {
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
  p.high_cl[which(x < median(cl_center[low_cl]))] = 0
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

HTODemux.mcl <- function(object, assay = "HTO", q = 1, seed = 42){
  # A function to find the threshold for each hastag, the P, and singlet, doublets and negative.
  # The input is the HTO data matrix (normalized across cells).
  
  assay <- assay %||% DefaultAssay(object = object)
  data <- GetAssayData(object = object, assay = assay, slot = 'data')
  hto_mcl.cutoff = data.frame(cut_off = future.apply::future_apply(data,1,function(x) select_hash_cutoff_mcl(x,q = q), future.seed = T))
  hto_mcl.cutoff$Multi_modal = apply(data,1,function(x) is_multimodal(x))
  print(hto_mcl.cutoff)
  hto_mcl.p = t(apply(data,1,function(x) hash_mcl_p(x, seed = seed)))
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
  object$hash.mcl.ID <- Idents(object = object)
  return(object)
}

HTODemux.mcl.visualization <- function(object, assay = "HTO", q = 1, seed = 42){
  assay <- assay %||% DefaultAssay(object = object)
  data <- GetAssayData(object = object, assay = assay, slot = 'data')
  hto.data.wide = data.frame(t(data))
  hto.data.long = data.table::melt(data.table::setDT(hto.data.wide,keep.rownames = T), id.vars = 'rn', value.name = 'Expression',variable.name = 'hto')
  hto_mcl.cutoff = data.frame(cut_off = future.apply::future_apply(data,1,function(x) select_hash_cutoff_mcl(x,q = q), future.seed = T), hto = colnames(hto.data.wide)[-1])
  
  p <- ggplot(hto.data.long, aes(x = Expression))+
    geom_histogram(bins = 100)+
    geom_vline(data = hto_mcl.cutoff, aes(xintercept = cut_off), col = 'red')+
    facet_wrap(~hto,scales = 'free',ncol = 2) +
    xlab('Expression') +
    ylab('Counts') +
    ggtitle('Individual HTO distributions') +
    theme_minimal()
  
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
    theme_minimal() +
    ggtitle(title)
  
  return(p)
}


################################################################################
##########                            Init                            ##########
################################################################################
scop.ID     <- snakemake@params[['scopID']]
com.ID      <- snakemake@params[['comID']]
com.ID.list <- str_split(com.ID, ',', simplify = T)[1,]
user.ID     <- snakemake@params[['userID']]

dir.proj <- paste0('/', file.path('projects',user.ID,'COMUNEQAID','outs',scop.ID))
dir.outs.qc <- file.path(dir.proj, 'scRNAseq', '00_QC')
dir.bcls <- file.path(dir.proj, 'scRNAseq', '01_BCL')

for (com.ID in com.ID.list) {
  unique.ID <- paste(scop.ID,com.ID,sep = '_')
  
  dir.data <- paste0('/', file.path('projects',user.ID,'COMUNEQAID','manage-dir','tmp-data',unique.ID))
  dir.outs.log <- file.path(dir.proj, 'scRNAseq', '04_Log', com.ID)
  dir.outs.qc.plots <- file.path(dir.outs.qc, com.ID, 'plots')
  
  dir.outs.indi <- file.path(dir.proj, 'Output', 'data', 'reaction-unfiltered', com.ID)
  
  var.orga <- scan(file.path(dir.data,'tmp_organism.txt'), what = '', sep = '\n', quiet = T)
  var.wofl <- scan(file.path(dir.data,'tmp_workflow.txt'), what = '', sep = '\n', quiet = T)
  
  pool.table <- read.table(file.path(dir.data,'poolTable-updated.csv'))
  
  if (var.wofl == '10x') {
    colnames(pool.table) <- c('Index (10x)','Lane','Loaded Cells','BCL PIN (10x)',"SEQ NAMES (10x)","READS (10x)")
  }
  if (var.wofl == '10x + HTO') {
    colnames(pool.table) <- c('Index (10x)','Index (HTO)','Lane','Loaded Cells','BCL PIN (10x)','BCL PIN (HTO)',"SEQ NAMES (HTO)","READS (HTO)","SEQ NAMES (10x)","READS (10x)")
  }
  
  
  nPools <- dim(pool.table)[1]
  
  dir.create(dir.outs.indi, recursive = T, showWarnings = F)
  dir.create(dir.outs.qc.plots, recursive = T, showWarnings = F)
  
  ################################################################################
  ##########       Processing single nucleus RNA sequencing data        ##########
  ################################################################################
  
  ####################  Iterate over pools  ####################
  cl <- makeCluster(nPools) 
  registerDoParallel(cl)
  foreach(i = seq(nPools),
          .combine = 'c',
          .packages = c('Seurat','tximport','SingleCellExperiment','future',
                        'tidyverse','scater','scran','sctransform',
                        'DropletUtils','Matrix','mclust')) %dopar% {
    
    ####################  Prep processing of sample i  ####################
    tmp.pool        <- pool.table[i,]
    
    pool.10x      <- tmp.pool[['Index (10x)']]
    pins.10x      <- str_split(tmp.pool[['BCL PIN (10x)']], ',')
    
    if (var.wofl == '10x + HTO') {
      pool.hto      <- tmp.pool[['Index (HTO)']]
      pins.hto      <- str_split(tmp.pool[['BCL PIN (HTO)']], ',')
      
      hto.string <- paste0('#\tand pairing with HTO index:\t',pool.hto,'\t\tfrom sequencing run(s):\t',paste(tmp.pool[['SEQ NAMES (HTO)']], collapse = ','),'\n')
      }else{
        hto.string <- NULL
      }
    
    sink.file <- file(paste0('logs/04_make-seurat_',pool.10x,'.log'), 'a')
    sink(file = sink.file, append = T,
         type = 'output', split = T)
    cat(as.character(lubridate::now()),'\n', sep = '')
    
    cat(rep('#',80),'\n',
        '#####                                                                      #####\n',
        '#####                     Summarizing sequencing lanes                     #####\n',
        '#####                                                    **** / ******     #####\n',
        rep('#',80),'\n',
        '#####\n',
        '##\n',
        '##\n',
        '##\tPROJECT ID:     ',scop.ID,'\n',
        '##\n',
        '##\tCOM ID:         ',com.ID,'\n',
        '##\n',
        '##\tORGANISM:       ',var.orga,'\n',
        '##\n',
        '##\tWORKFLOW:       ',var.wofl,'\n',
        '##\n',
        '##\n',
        '###\n',
        rep('#',80),'\n',
        sep = '')
    
    cat('##\n',
        '#\n',
        '#\n',
        '#\tProccessing 10x index:\t\t',pool.10x,'\tfrom sequencing run(s):\t',paste(tmp.pool[['SEQ NAMES (10x)']], collapse = ','),'\n',
        hto.string,
        '#\n',
        '#\t(pool ',i,' out of ',nPools,')\n',
        '#\n',
        '#\n',
        '###\n',
        rep('#',80),'\n',
        '##\n',
        '#\n',
        '#\n',
        sep = '')
    
    
    
    mat.files.10x <- file.path(dir.proj,'scRNAseq','03_PipelineOut',com.ID,'10x',pool.10x,'res')
    if (!file.exists(mat.files.10x)) { stop('Not able to locate quants_mat.gz (gene expression)')}
    
    if (var.wofl == '10x + HTO') {
      mat.files.hto <- file.path(dir.proj,'scRNAseq','03_PipelineOut',com.ID,'hto',pool.hto,'res')
      if (var.wofl == '10x + HTO') {if (!file.exists(mat.files.hto)) { stop('Not able to locate quants_mat.gz (HTO)')}}
    }
    
    ####################  Read RNA matrices ####################
    cat('#\tReading count matrices..\n',
        sep = '')
    counts.all <- load_fry(frydir = mat.files.10x, output_list = T)
    counts.rna <- counts.all[['RNA']]
    counts.spl <- counts.all[['Spliced']]
    counts.uns <- counts.all[['Unspliced']]
    
    if (var.wofl == '10x + HTO') {
      cat('#\t-\tHTO\n',
          sep = '')
      counts.hto <- load_fry(frydir = mat.files.hto)
    }
    cat('#\n',
        '#\n',
        sep = '')
    
    ####################  Filter low rank barcodes  ####################
    cells.keep <- emptyDropsmodal(counts.uns)
    
    ############ MAKE UMI PLOT  (MAYBE NEW VERSION) ############
    if (var.wofl == '10x + HTO') {
      featDump.10x <- suppressMessages(read_delim(file.path(mat.files.10x,'featureDump.txt'), delim = '\t'))
      featDump.hto <- suppressMessages(read_delim(file.path(mat.files.hto,'featureDump.txt'), delim = '\t'))
      
      p <- makeUmiPlot(feat.RNA = featDump.10x,
                       feat.HTO = featDump.hto,
                       CBs.called = cells.keep,
                       title = paste0('UMI-UMI plot - ', tmp.pool[['Index (10x)']]))
      
      ggsave(filename = paste0(tmp.pool[['Index (10x)']],'_UMI-UMI-plot.png'),
             plot = p,
             path = dir.outs.qc.plots,
             width = 8,
             height = 8)
    }
    ############################################################
    
    counts.rna <- counts.rna[,cells.keep]
    counts.spl <- counts.spl[,cells.keep]
    counts.uns <- counts.uns[,cells.keep]
    
    ####################  Annotation  ####################
    cat('#\tAnnotating..\n', 
        sep = '')
    counts.rna <- annotateMat(counts.rna, var.orga)
    cat('#\t-\tRNA assay\n',
        sep = '')
    counts.spl <- annotateMat(counts.spl, var.orga)
    cat('#\t-\tspliced assay\n',
        sep = '')
    counts.uns <- annotateMat(counts.uns, var.orga)
    cat('#\t-\tunspliced assay\n',
        '#\n',
        '#\n',
        sep = '')
    
    if (var.wofl == '10x + HTO') {
      cat('#\tProcessing HTOs..\n',
          sep = '')
      
      cells.keep <- intersect(cells.keep,colnames(counts.hto))
      cat('#\t-\tretaining ',length(cells.keep),' cells overlapping between RNA & HTO\n',
          '#\n',
          '#\n',
          sep = '')
      
      cat('#\tConverting matrices to Seurat object(s)..\n',
          sep = '')
      cat('#\t-\tRNA\n')
      seur.full <- CreateSeuratObject(counts.rna[,cells.keep], project = pool.10x)
      seur.full[['Index.10x']] <- pool.10x
      cat('#\t-\tHTO\n')
      seur.full[['HTO']] <- CreateAssayObject(counts.hto[,cells.keep])
      seur.full[['Index.HTO']] <- pool.hto
      cat('#\t-\tSpliced\n')
      seur.full[['spliced']] <- CreateAssayObject(counts.spl[,cells.keep])
      cat('#\t-\tUnspliced\n')
      seur.full[['unspliced']] <- CreateAssayObject(counts.uns[,cells.keep])
      
      cat('#\n',
          '#\n',
          sep = '')
      
      cat('#\tHTO demultiplexing..\n')
      cat('#\t-\tnormalizing..\n')
      #seur.full <- NormalizeData(seur.full, assay = 'HTO', normalization.method = 'CLR', margin = 1)
      seur.full <- NormalizeData(seur.full, assay = 'HTO', normalization.method = 'CLR', margin = 2, verbose = F)
      
      #VariableFeatures(seur.full, assay = 'HTO') <- rownames(seur.full[['HTO']]@counts)
      
      #cat('#\t-\tscaling..\n',
      #    sep = '')
      #seur.full <- ScaleData(seur.full, assay = 'HTO', verbose = F)
      
      #cat('#\t-\toptimizing HTO-demultiplexing parameters for singlet abundance..\n',
      #    sep = '')
      
      #x <- seur.full
      #thresh <- .95
      #sing.max <- 0
      
      #for (q in seq(.99, .5, -.01)) {
      #  x <- HTODemux(x,
      #                assay = 'HTO',
      #                kfunc = 'kmeans',
      #                positive.quantile = q,
      #                verbose = F)
      #  
      #  doub.tmp <- table(x[['HTO_classification.global']])['Doublet']
      #  nega.tmp <- table(x[['HTO_classification.global']])['Negative']
      #  sing.tmp <- table(x[['HTO_classification.global']])['Singlet']
      #  
      #  if (sing.tmp > sing.max) {
      #    doub.max <- doub.tmp
      #    nega.max <- nega.tmp
      #    sing.max <- sing.tmp
      #    thresh <- q
      #  }
      #}
      #cat('#\t-\t\tOptimal quantile:\t',thresh,'\n',
      #    '#\t-\t\tDoublets:\t\t',doub.max,'\n',
      #    '#\t-\t\tNegatives:\t\t',nega.max,'\n',
      #    '#\t-\t\tSinglets:\t\t',sing.max,'\n',
      #    sep = '')
      #
      #remove(list = c('x'))
      #
      cat('#\t-\tdemultiplexing HTOs..\n',
          sep = '')
      #seur.full <- HTODemux(seur.full, assay = 'HTO', kfunc = 'kmeans', positive.quantile = thresh, verbose = F)
      seur.full <- HTODemux.mcl(seur.full, q = 1)
      seur.full[,seur.full$HTO_mcl_classification.global == 'Doublet']
      
      table(seur.full$HTO_mcl_classification.global)[names(table(seur.full$HTO_mcl_classification.global)) == 'Doublet']
      
      cat('#\t-\t\tDoublets:\t\t',table(seur.full$HTO_mcl_classification.global)[names(table(seur.full$HTO_mcl_classification.global)) == 'Doublet'],'\n',
          '#\t-\t\tNegatives:\t\t',table(seur.full$HTO_mcl_classification.global)[names(table(seur.full$HTO_mcl_classification.global)) == 'Negative'],'\n',
          '#\t-\t\tSinglets:\t\t',table(seur.full$HTO_mcl_classification.global)[names(table(seur.full$HTO_mcl_classification.global)) == 'Singlet'],'\n',
          sep = '')
      
      HTODemux.mcl.visualization(seur.full)
      
      cat('#\t-\tperforming tSNE\n',
          sep = '')
      seur.full <- RunTSNE(seur.full,
                           distance.matrix = as.matrix(dist(t(GetAssayData(object = seur.full, assay = 'HTO')))),
                           reduction.name = 'hto.tsne')
      
      cat('#\t-\tcreating plots\n')
      p <- VlnPlot(seur.full, features = c('nCount_HTO'), log = T)
      ggsave(filename = paste0(tmp.pool[['Index (10x)']],'_HTO-violin-plot.png'),
             plot = p,
             path = dir.outs.qc.plots,
             width = 8,
             height = 8)

      p <- DimPlot(seur.full, group.by = 'hash.mcl.ID') + ggtitle(paste0('HTO demultiplexing'))
      ggsave(filename = paste0(tmp.pool[['Index (10x)']],'_HTO-tSNE-plot.png'),
             plot = p,
             path = dir.outs.qc.plots,
             width = 8,
             height = 8)
      
      #Idents(seur.full) <- 'HTO_maxID'
      #p <- RidgePlot(seur.full, features = rownames(seur.full[['HTO']]@counts), assay = 'HTO', ncol = 1)
      #suppressMessages(
      #  ggsave(filename = paste0(tmp.pool[['Index (10x)']],'_HTO-ridge-plot.png'),
      #         plot = p,
      #         path = dir.outs.qc.plots,
      #         width = 8,
      #         height = (2*length(table(seur.full$HTO_maxID))))
      #)
      
      cat('#\n',
          '#\n',
          sep = '')
      
      cat('#\tConverting matrices to SinglecellExperiment object(s)..\n',
          sep = '')
      sce.full <- as.SingleCellExperiment(seur.full)
      
      
      # identify intra-hash doublets
      cat('#\t-\tidentifying intra-hash doublets','\n', sep = '')
      sce.full <- logNormCounts(sce.full)
      dec.hash <- modelGeneVar(sce.full)
      top.hash <- getTopHVGs(dec.hash, n = 1000)
      set.seed(1011110)
      sce.full <- runPCA(sce.full, subset_row = top.hash, ncomponents = 20)
      sce.full[['doublet']] <- if_else(sce.full[['hash.mcl.ID']] == 'Doublet', true = T, false = F)
      
      # Recovering the intra-sample doublets:
      cat('#\t-\trecovering intra-sample doublets\n',
          sep = '')
      hashed.doublets <- scDblFinder::recoverDoublets(sce.full,
                                                      use.dimred = 'PCA',
                                                      doublets = sce.full[['doublet']],
                                                      samples = table(sce.full[['hash.mcl.ID']]))
      
      sce.full[['proportion_dub_neighbors']] <- hashed.doublets[['proportion']]
      sce.full[['predicted_dub_std']] <- hashed.doublets[['predicted']]
      sce.full[['predicted_dub_cut']] <- dub_cutoff(sce.full)
      
      seur.full[['doublet']] <- sce.full[['doublet']]
      seur.full[['predicted_dub_std']] <- sce.full[['predicted_dub_std']]
      seur.full[['predicted_dub_cut']] <- sce.full[['predicted_dub_cut']]
    }
    
    if (var.wofl == '10x') {
      cat('#\tConverting matrices to Seurat object(s)..\n',
          sep = '')
      cat('#\t-\tRNA\n')
      seur.full <- CreateSeuratObject(counts.rna[,cells.keep], project = pool.10x)
      seur.full[['Index (10x)']] <- pool.10x
      cat('#\t-\tSpliced\n')
      seur.full[['spliced']] <- CreateAssayObject(counts.spl[,cells.keep])
      cat('#\t-\tUnspliced\n')
      seur.full[['unspliced']] <- CreateAssayObject(counts.uns[,cells.keep])
      
      sce.full <- as.SingleCellExperiment(seur.full)
      
      cat('#\t-\tidentifying intra-hash doublets','\n', sep = '')
      sce.full <- logNormCounts(sce.full)
      top.rna <- getTopHVGs(modelGeneVar(sce.full))
      set.seed(1011110)
      sce.full <- runPCA(sce.full, subset_row = top.rna, n = 1000, ncomponents = 20)
      
      dbl.dens <- scDblFinder::computeDoubletDensity(sce.full, subset.row = top.rna,
                                                     d = ncol(reducedDim(sce.full)))
      
      seur.full[['DoubletScore']] <- dbl.dens
      
      
    }
    
    p <- VlnPlot(seur.full, c('nCount_RNA','nFeature_RNA'), log = T)
    ggsave(filename = paste0(tmp.pool[['Index (10x)']],'_RNA-violin-plot.png'),
           plot = p,
           path = dir.outs.qc.plots,
           width = 8,
           height = 8)
    
    p <- FeatureScatter(seur.full, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA')
    ggsave(filename = paste0(tmp.pool[['Index (10x)']],'_RNA-feature-scatter-plot.png'),
           plot = p,
           path = dir.outs.qc.plots,
           width = 8,
           height = 8)
    
    cat('#\n',
        '#\n',
        '#\tSaving object..\n',
        sep = '')
    
    out.name.seur <- paste0(com.ID,'_',tmp.pool[['Index (10x)']],'_seurat.rds')
    out.path.seur <- file.path(dir.outs.indi,out.name.seur)
    
    cat('#\t-\tsaving Seurat object..\n',
        sep = '')
    saveRDS(seur.full, out.path.seur)
    
    # Memory cleanup
    remove(
      list = c('sce.full','seur.full','counts.all','counts.rna','counts.spl','counts.uns')
    )
    
    cat('#\n',
        '#\n',
        '###\n',
        rep('#',80),'\n',
        sep = '')
    
    sink()
    
    NULL
          }
}

file.create(snakemake@output[[1]])