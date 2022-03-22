################################################################################
##########                          Packages                          ##########
################################################################################
suppressPackageStartupMessages({
  library(Seurat)
  library(tximport)       
  library(AnnotationDbi)
  library(SingleCellExperiment)
  library(future)
  library(tidyverse)
  library(scater)
  library(scran)
  library(sctransform)
  library(DropletUtils)
  library(rjson)
  library(jsonlite)
  library(Matrix)
  library(mclust)
  library(kableExtra)
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
        png(paste0(dir.outs.qc.plots,'/',tmp.pool[['Index (10x)']],'_KneePlot.png'), width = 1500, height = 500)
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
        png(paste0(dir.outs.qc.plots,'/',tmp.pool[['Index (10x)']],'_KneePlot.png'), width = 1500, height = 500)
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
  
  png(paste0(dir.outs.qc.plots,'/',tmp.pool[['Index (10x)']],'_KneePlot.png'), width = 1500, height = 500)
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

add_lane_to_df <- function(i, l){
  df <- l[[i]]
  df[['Lane']] <- i
  return(df)
}

collapse_demux_results <- function(l){
  res <- lapply(seq_along(l), add_lane_to_df, l = l)
  df <- do.call('rbind', res)
  return(df)
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
scop.ID <- snakemake@params[['scopID']]
com.ID  <- snakemake@params[['comID']]
user.ID  <- snakemake@params[['userID']]

unique.ID <- paste(scop.ID,com.ID,sep = '_')

dir.proj <- paste0('/', file.path('projects',user.ID,'COMUNEQAID','outs',scop.ID)) # <- change to your own WD
dir.data <- paste0('/', file.path('projects','SCOP','pipelines','COMUNEQAID','COMUNEQAID-app','app-data','tmp-data',unique.ID))

dir.outs.qc <- file.path(dir.proj, 'scRNAseq', '00_QC')
dir.bcls <- file.path(dir.proj, 'scRNAseq', '01_BCL')
dir.outs.log <- file.path(dir.proj, 'scRNAseq', '04_Log', com.ID)
dir.outs.qc.plots <- file.path(dir.outs.qc, 'plots',com.ID)
dir.outs.qc.summa <- file.path(dir.outs.qc, 'summary',com.ID)

dir.outs.indi <- file.path(dir.proj, 'Output', 'data', 'indi', com.ID)
dir.outs.comb <- file.path(dir.proj, 'Output', 'data', 'comb', com.ID)

var.type <- scan(file.path(dir.data,'tmp_seqType.txt'), what = '', sep = '\n', quiet = T)
var.orga <- scan(file.path(dir.data,'tmp_organism.txt'), what = '', sep = '\n', quiet = T)
var.wofl <- scan(file.path(dir.data,'tmp_workflow.txt'), what = '', sep = '\n', quiet = T)



################################################################################
##########                    Update Pool Table                       ##########
################################################################################

tmp.seqs.10x  <- list()
tmp.seqs.hto  <- list()

read.df.all   <- list()
  
if (var.wofl == '10x') {
  pool.table <- read.csv(file.path(dir.data,'poolTable.csv'),
                         colClasses = c('BCL.PIN..10x.' = 'character'))
  colnames(pool.table) <- c('Index (10x)','Lane','Loaded Cells','BCL PIN (10x)')
}
if (var.wofl == '10x + HTO') {
  pool.table <- read.csv(file.path(dir.data,'poolTable.csv'),
                         colClasses = c('BCL.PIN..10x.' = 'character',
                                        'BCL.PIN..HTO.' = 'character'))
  colnames(pool.table) <- c('Index (10x)','Index (HTO)','Lane','Loaded Cells','BCL PIN (10x)','BCL PIN (HTO)')
  
  pool.table[['SEQ NAMES (HTO)']] <- ''
  pool.table[['READS (HTO)']] <- 0
}

pool.table[['SEQ NAMES (10x)']] <- ''
pool.table[['READS (10x)']] <- 0
  
n.pools <- dim(pool.table)[1]

################################################################################
##########                    Iterate over lane                       ##########
################################################################################

for (pool.i in seq(n.pools)) {
  tmp.pool      <- pool.table[pool.i,]
  
  pins.10x      <- str_split(tmp.pool[['BCL PIN (10x)']], ',', simplify = T)[1,]
  pool.10x      <- str_split(tmp.pool[['Index (10x)']], ',', simplify = T)[1,]
  seqs.10x      <- listLen(pins.10x)
  
  for (pin.i in seq(pins.10x)) {
    tmp.seq <- grep(list.files(dir.bcls), pattern = pins.10x[pin.i], value = T)
    seqs.10x[pin.i] <- tmp.seq
    tmp.seqs.10x <- c(tmp.seqs.10x, tmp.seq)
  }
  
  year <- substr(tmp.seq, start = 1, stop = 2)
  tmp.seqs.10x.unique <- unique(tmp.seqs.10x)
  tmp.pool[['SEQ NAMES (10x)']] <- paste(tmp.seqs.10x.unique, collapse = ',')
  tmp.seqs.all <- tmp.seqs.10x.unique
  
  if (var.wofl == '10x + HTO') {
    pins.hto      <- str_split(tmp.pool[['BCL PIN (HTO)']], ',', simplify = T)[1,]
    pool.hto      <- str_split(tmp.pool[['Index (HTO)']], ',', simplify = T)[1,]
    seqs.hto      <- listLen(pins.hto)
    
    for (pin.i in seq(pins.hto)) {
      tmp.seq <- grep(list.files(dir.bcls), pattern = pins.hto[pin.i], value = T)
      seqs.hto[pin.i] <- tmp.seq
      tmp.seqs.hto <- c(tmp.seqs.hto,tmp.seq)
    }
    
    tmp.seqs.hto.unique <- unique(tmp.seqs.hto)
    tmp.pool[['SEQ NAMES (HTO)']] <- paste(tmp.seqs.hto.unique, collapse = ',')
    tmp.seqs.all <- unique(tmp.seqs.10x.unique,tmp.seqs.hto.unique)
  }
  
  for (bcl in rev(tmp.seqs.all)) {
    stats.json.path <- file.path(dir.proj,'scRNAseq','02_FASTQ',bcl,'fastq-path/Stats/Stats.json')
    js <- jsonlite::read_json(stats.json.path, simplifyVector = TRUE)
    lane_data <- collapse_demux_results(js[['ConversionResults']][['DemuxResults']])
    lane_data[['Lane']] <- factor(lane_data[['Lane']])
    known_barcodes <- lane_data[,c('Lane', 'SampleId', 'NumberReads')]
    
    if (year <= 20) {
      known_barcodes[['SampleId']] <- substr(known_barcodes[['SampleId']],1,nchar(known_barcodes[['SampleId']]) - 2)
    }
    
    known_barcodes %>%
      group_by(SampleId) %>%
      summarise(Counts = sum(NumberReads)) -> read.df
    
    for (pool.10x.i in seq(pool.10x)) {
        
        if (year <= 20) {
          tmp.index <- paste0('SI-GA-', pool.10x[pool.10x.i])
          tmp.pool[['Index (10x)']] <- paste(paste0('SI-GA-', pool.10x), collapse = ',')
        }
        if (year >= 21) {
          tmp.index <- paste0('SI-TT-', pool.10x[pool.10x.i])
          tmp.pool[['Index (10x)']] <- paste(paste0('SI-TT-', pool.10x), collapse = ',')
        }
        if (tmp.index %in% read.df[['SampleId']]) {
          tmp.pool[['READS (10x)']] <- tmp.pool[['READS (10x)']] + read.df[read.df[['SampleId']] == tmp.index,] [['Counts']]
        }
    }
    if (var.wofl == '10x + HTO') {
      for (pool.hto.i in seq(pool.hto)) {
        tmp.index <- pool.hto[pool.hto.i]
        
        if (tmp.index %in% read.df[['SampleId']]) {
          tmp.pool[['READS (HTO)']] <- tmp.pool[['READS (HTO)']] + read.df[read.df[['SampleId']] == tmp.index,][['Counts']]
        }
      }
    }
  }
  pool.table[pool.i,] <- tmp.pool
}

nPools        <- dim(pool.table)[1]

dir.create(dir.outs.indi, recursive = T, showWarnings = F)
dir.create(dir.outs.comb, recursive = T, showWarnings = F)
dir.create(dir.outs.qc.plots, recursive = T, showWarnings = F)
dir.create(dir.outs.qc.summa, recursive = T, showWarnings = F)
dir.create(dir.outs.log, recursive = T, showWarnings = F)


################################################################################
##########                          QC plots                          ##########
################################################################################

# Top n BCs
bc.obs.tib <- tibble(
  'Barcode' = character(),
  'n_obs_comb' = numeric(),
  'Seq' = character()
)

if (var.wofl == '10x') {
  tmp.seqs.unique <- tmp.seqs.10x.unique
}
if (var.wofl == '10x + HTO') {
  tmp.seqs.unique <- unique(tmp.seqs.10x.unique,tmp.seqs.hto.unique)
}

for (bcl in tmp.seqs.unique) {
  tmp.tib <- tibble(
    'Barcode' = character(),
    'n_obs' = numeric()
  )
  
  stats.json.path <- file.path(dir.proj,'scRNAseq','02_FASTQ',bcl,'fastq-path/Stats/Stats.json')
  
  js.nim <- jsonlite::read_json(stats.json.path, simplifyVector = F)
  
  for (i in seq(js.nim[['UnknownBarcodes']])) {
    tmp.bc.list <- js.nim[['UnknownBarcodes']][[i]][['Barcodes']]
    
    for (j in seq(tmp.bc.list)) {
      tmp.tib <- add_row(tmp.tib,
                         'Barcode' = names(tmp.bc.list[j]),
                         'n_obs' = tmp.bc.list[[j]]
      )
    }
  }
  group_by(tmp.tib, Barcode) %>% 
    summarise(n_obs_comb = sum(n_obs)) %>% 
    arrange(desc(n_obs_comb)) -> tmp.tib.sum
  
  tmp.tib.sum.sub <- tmp.tib.sum[1:10,]
  
  p <- ggplot(tmp.tib.sum.sub) +
    geom_bar(mapping = aes(x = reorder(Barcode, n_obs_comb), y = n_obs_comb),stat = 'identity') +
    coord_flip() +
    xlab('Bardcode') +
    ylab('Count') +
    ggtitle(paste0('Top 10 unknown barcodes - ',bcl)) +
    theme_minimal()
  
  ggsave(filename = paste0('top10UndeterminedBCs_',bcl,'.png'),
         plot = p,
         path = dir.outs.qc.plots,
         width = 12,
         height = 8)
  
}

read.df.all <- list()

for (bcl in tmp.seqs.unique) {
  stats.json.path <- file.path(dir.proj,'scRNAseq','02_FASTQ',bcl,'fastq-path/Stats/Stats.json')
  js.sim <- jsonlite::read_json(stats.json.path, simplifyVector = T)
  lane_data <- collapse_demux_results(js.sim[['ConversionResults']][['DemuxResults']])
  unde_data <- js.sim[['ConversionResults']][['Undetermined']]
  lane_data[['Lane']] <- factor(lane_data[['Lane']])
  known_barcodes <- lane_data[,c('Lane', 'SampleId', 'NumberReads')]
  
  if (year <= 20) {
    known_barcodes[['SampleId']] <- substr(known_barcodes[['SampleId']],1,nchar(known_barcodes[['SampleId']]) - 2)
  }
  
  known_barcodes %>%
    group_by(SampleId) %>%
    summarise(Counts = sum(NumberReads)) -> read.df
  
  und.df <- list(as.character('Undetermined'),
                 as.integer(sum(unde_data[['NumberReads']])))
  read.df <- rbind(read.df,und.df)
  read.df[['Seq']] <- bcl
  read.df.all <- rbind(read.df.all,read.df)
  
}

p <- ggplot(read.df.all, aes(x = SampleId, y = Counts)) +
  geom_bar(stat = 'identity') +
  coord_flip() +
  facet_wrap(~Seq)

ggsave(filename = paste0('readDistribution_bcl2fastq.png'),
       plot = p,
       path = dir.outs.qc.plots,
       width = 12,
       height = 12)

################################################################################
##########       Processing single nucleus RNA sequencing data        ##########
################################################################################

####################  Start sink  ####################
sink.file <- file(file.path(dir.outs.log,'03_summarize.log'), 'a')
sink(file = sink.file, append = T,
     type = 'output', split = T)
cat(as.character(lubridate::now()),'\n', sep = '')

cat(rep('#',80),'\n',
    '#####                                                                      #####\n',
    '#####                     Summarizing sequencing lanes                     #####\n',
    '#####                                                     **** / *****     #####\n',
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

####################  Iterate over pools  ####################
for (pool.i in seq(n.pools)) {
  ####################  Prep processing of sample i  ####################
    tmp.pool        <- pool.table[pool.i,]
    
    pool.10x      <- tmp.pool[['Index (10x)']]
    pins.10x      <- str_split(tmp.pool[['BCL PIN (10x)']], ',')
    
    if (var.wofl == '10x + HTO') {
      pool.hto      <- tmp.pool[['Index (HTO)']]
      pins.hto      <- str_split(tmp.pool[['BCL PIN (HTO)']], ',')

      hto.string <- paste0('#\tand pairing with HTO index:\t',pool.hto,'\t\tfrom sequencing run(s):\t',paste(tmp.pool[['SEQ NAMES (HTO)']], collapse = ','),'\n')
    }else{
      hto.string <- NULL
    }
    
    cat('##\n',
        '#\n',
        '#\n',
        '#\tProccessing 10x index:\t\t',pool.10x,'\tfrom sequencing run(s):\t',paste(tmp.pool[['SEQ NAMES (10x)']], collapse = ','),'\n',
        hto.string,
        '#\n',
        '#\t(pool ',pool.i,' out of ',nPools,')\n',
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
      cat('#\t-\tHTO\n')
      seur.full[['HTO']] <- CreateAssayObject(counts.hto[,cells.keep])
      cat('#\t-\tSpliced\n')
      seur.full[['spliced']] <- CreateAssayObject(counts.spl[,cells.keep])
      cat('#\t-\tUnspliced\n')
      seur.full[['unspliced']] <- CreateAssayObject(counts.uns[,cells.keep])
      
      cat('#\n',
          '#\n',
          sep = '')
      
      cat('#\tHTO demultiplexing..\n')
      cat('#\t-\tnormalizing..\n')
      seur.full <- NormalizeData(seur.full, assay = 'HTO', normalization.method = 'CLR', margin = 1, verbose = F)
      
      VariableFeatures(seur.full, assay = 'HTO') <- rownames(seur.full[['HTO']]@counts)
      
      cat('#\t-\tscaling..\n',
          sep = '')
      seur.full <- ScaleData(seur.full, assay = 'HTO', verbose = F)
      
      cat('#\t-\toptimizing HTO-demultiplexing parameters for singlet abundance..\n',
          sep = '')
      
      x <- seur.full
      thresh <- .95
      sing.max <- 0
      
      for (q in seq(.99, .5, -.01)) {
        x <- HTODemux(x,
                      assay = 'HTO',
                      kfunc = 'kmeans',
                      positive.quantile = q,
                      verbose = F)
        
        doub.tmp <- table(x[['HTO_classification.global']])['Doublet']
        nega.tmp <- table(x[['HTO_classification.global']])['Negative']
        sing.tmp <- table(x[['HTO_classification.global']])['Singlet']
        
        if (sing.tmp > sing.max) {
          doub.max <- doub.tmp
          nega.max <- nega.tmp
          sing.max <- sing.tmp
          thresh <- q
          }
      }
      cat('#\t-\t\tOptimal quantile:\t',thresh,'\n',
          '#\t-\t\tDoublets:\t\t',doub.max,'\n',
          '#\t-\t\tNegatives:\t\t',nega.max,'\n',
          '#\t-\t\tSinglets:\t\t',sing.max,'\n',
          sep = '')
      
      remove(list = c('x'))
      
      cat('#\t-\tdemultiplexing HTOs..\n',
          sep = '')
      seur.full <- HTODemux(seur.full, assay = 'HTO', kfunc = 'kmeans', positive.quantile = thresh, verbose = F)

      cat('#\t-\tperforming tSNE\n',
          sep = '')
      seur.full <- RunTSNE(seur.full,
                           distance.matrix = as.matrix(dist(t(GetAssayData(object = seur.full, assay = 'HTO')))),
                           reduction.name = 'hto.tsne')

      cat('#\t-\tcreating plots\n')
      p <- VlnPlot(seur.full, features = c('nCount_HTO'), log = T)
      ggsave(filename = paste0(tmp.pool[['Index (10x)']],'_HTO_VlnPlot.png'),
             plot = p,
             path = dir.outs.qc.plots,
             width = 8,
             height = 8)
      
      p <- DimPlot(seur.full, group.by = 'hash.ID') + ggtitle(paste0('HTO demultiplexing (quantile: ',thresh,')'))
      ggsave(filename = paste0(tmp.pool[['Index (10x)']],'_HTO_tSNEPlot.png'),
             plot = p,
             path = dir.outs.qc.plots,
             width = 8,
             height = 8)
      
      Idents(seur.full) <- 'HTO_maxID'
      p <- RidgePlot(seur.full, features = rownames(seur.full[['HTO']]@counts), assay = 'HTO', ncol = 1)
      suppressMessages(
        ggsave(filename = paste0(tmp.pool[['Index (10x)']],'_HTO_ridgePlot.png'),
               plot = p,
               path = dir.outs.qc.plots,
               width = 8,
               height = (2*length(table(seur.full$HTO_maxID))))
      )
      
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
      sce.full[['doublet']] <- if_else(sce.full[['hash.ID']] == 'Doublet', true = T, false = F)
      
      # Recovering the intra-sample doublets:
      cat('#\t-\trecovering intra-sample doublets\n',
          sep = '')
      hashed.doublets <- scDblFinder::recoverDoublets(sce.full,
                                                      use.dimred = 'PCA',
                                                      doublets = sce.full[['doublet']],
                                                      samples = table(sce.full[['hash.ID']]))
      
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
    ggsave(filename = paste0(tmp.pool[['Index (10x)']],'_RNA_VlnPlot.png'),
           plot = p,
           path = dir.outs.qc.plots,
           width = 8,
           height = 8)
    
    p <- FeatureScatter(seur.full, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA')
    ggsave(filename = paste0(tmp.pool[['Index (10x)']],'_RNA_FeaScaPlot.png'),
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

}

################################################################################
##########                 Preparing aggregated output                ##########
################################################################################

cat('#####                                                                      #####\n',
    '#####                     Summarizing sequencing lanes                     #####\n',
    '#####                                                    ***** / *****     #####\n',
    rep('#',80),'\n',
    '#####\n',
    '##\n',
    '##\n',
    sep = '')
  
cat('#\tReading the following files: \n',
    sep = '')
indi.files <- grep(list.files(dir.outs.indi), pattern = '_seurat.rds', value = T)
for (fname in indi.files) {
  cat('#\t-\t',fname,'\n',
      sep = '')
}
  
seur.list <- lapply(file.path(dir.outs.indi,indi.files), readRDS)
cat('#\n',
    '#\n',
    sep = '')
  
cat('#\tFinding overlapping genes..\n',
    sep = '')
  
# Extract an overlapping set of rownames for each object
common.genes <- purrr::map(seur.list, rownames) %>% 
  Reduce(intersect, .)
cat('#\t-\tfound:\t',length(common.genes),'\n',
    '#\n',
    '#\n',
    sep = '')
  
cat('#\tFiltering objects to contain only common genes and binding to final object..\n',
    '#\n',
    '#\n',
    sep = '')


# Merge objects
if (length(indi.files) > 1) {
  seur.comb <- merge(seur.list[[1]],
                     seur.list[2:length(seur.list)])
}else {
  seur.comb <- seur.list[[1]]
}
cat('#\tFiltering genes with no detection..\n',
    sep = '')

# Filter genes with no detection
rowsums <- Matrix::rowSums(seur.comb[['RNA']]@counts)
num.genes <- length(rownames(seur.comb))
seur.comb <- seur.comb[!rowsums == 0,]
num.genes.red <- length(rownames(seur.comb))

cat('#\t-\tthrowing out ',num.genes - num.genes.red,' genes\n',
    '#\t-\tretaining ',num.genes.red,' genes\n',
    '#\n',
    '#\n',
    sep = '')

if (var.wofl == '10x + HTO') {
  Idents(seur.comb) <- 'predicted_dub_std'
  cells.remove.predDubStd <- WhichCells(seur.comb, idents = TRUE)
  Idents(seur.comb) <- 'predicted_dub_cut'
  cells.remove.predDubCut <- WhichCells(seur.comb, idents = TRUE)
  
  cells.remove.predDub <- union(cells.remove.predDubStd,cells.remove.predDubCut)
  seur.comb[['predicted_dub_all']] <- colnames(seur.comb) %in% cells.remove.predDub
  
  p <- VlnPlot(seur.comb, 'nCount_RNA', group.by = 'predicted_dub_all', log = T) + ggtitle('Inferential doublets')
  ggsave(filename = paste0(com.ID,'_cells_remove_infDoubl.png'),
         plot = p,
         path = dir.outs.qc.plots,
         width = 8,
         height = 8)
}

cat('#\tRemoving blacklisted cells..\n',
    '#\t(starting with ',length(colnames(seur.comb)),' cells)\n',
    sep = '')
if (var.wofl == '10x + HTO') {
  cat('#\t-\tremoving negatives and doublets\n',
      sep = '')
  # retain only cells classified as singlets
  Idents(seur.comb) <- 'HTO_classification.global'
  seur.comb <- subset(seur.comb, idents = 'Singlet')
  cat('#\t\t(',length(colnames(seur.comb)),' cells remaining)\n',
      sep = '')
  cat('#\t-\tremoving inferential doublets\n',
      sep = '')
  # remove manually called doublets
  Idents(seur.comb) <- 'predicted_dub_all'
  seur.comb <- subset(seur.comb, idents = FALSE)
    
  cat('#\t\t(',length(colnames(seur.comb)),' cells remaining)\n',
      sep = '')
}

cat('#\n',
    '#\n',
    sep = '')

Idents(seur.comb) <- 'hash.ID'
cat('#\tAnalyzing object..\n',
    sep = '')

# Analysis
cat('#\t-\tperforming SCTransform','\n',
    sep = '')
seur.comb <- SCTransform(seur.comb, method = 'qpoisson', verbose = FALSE)
cat('#\t-\tperforming PCA\n',
    sep = '')
seur.comb <- RunPCA(seur.comb, verbose = FALSE)
cat('#\t-\testimating dimensionality','\n',
    sep = '')
dims <- round(
  as.numeric(
    intrinsicDimension::maxLikGlobalDimEst(
      data = seur.comb@reductions[['pca']][, 1:50],
      k = 20)))
cat('#\t-\tperforming UMAP\n',
    sep = '')
seur.comb <- suppressWarnings(RunUMAP(seur.comb, dims = seq(dims), verbose = FALSE))
cat('#\t-\tfinding neighbors\n',
    sep = '')
seur.comb <- FindNeighbors(seur.comb, dims = seq(dims), verbose = FALSE)
cat('#\t-\tfinding clusters\n',
    sep = '')
seur.comb <- FindClusters(seur.comb, resolution = .8, verbose = FALSE)

cat('#\t-\tsaving dimentional reduction plot\n',
    sep = '')
if (var.wofl == '10x + HTO') {
  p <- DimPlot(seur.comb, group.by = 'hash.ID') + ggtitle(paste0('Samples (',dims,' PCs)'))
}else {
  p <- DimPlot(seur.comb, group.by = 'orig.ident') + ggtitle(paste0('Samples (',dims,' PCs)'))
}
  
ggsave(filename = paste0(com.ID,'_RNA_UMAP.png'),
       plot = p,
       path = dir.outs.qc.plots,
       width = 8,
       height = 8)


cat('#\n',
    '#\n',
    '#\tSaving objects..\n',
    sep = '')
  
# Output objects
out.name.rna <- paste0(com.ID,'_rna_seurat.rds')
out.name.full <- paste0(com.ID,'_full_seurat.rds')

out.path.rna <- file.path(dir.outs.comb,out.name.rna)
out.path.full <- file.path(dir.outs.comb,out.name.full)  

cat('#\t-\tsaving full object\n',
    sep = '')
saveRDS(seur.comb, out.path.full)

cat('#\t-\tsaving RNA object\n',
    sep = '')
seur.comb[['spliced']] <- NULL
seur.comb[['unspliced']] <- NULL

saveRDS(seur.comb, out.path.rna)

cat('#\n',
    '#\n',
    sep = '')
  
################################################################################
##########                          QC plots                          ##########
################################################################################

df <- data.frame(Type = character(),
                 CountSum = double(),
                 CountPer = double(),
                 CountTot = double(),
                 Pool = character())

for (i in seq(seur.list)) {
  
  tmp.pool <- pool.table[i,]
  
  pool.10x <- tmp.pool[['Index (10x)']]
  pool.hto <- tmp.pool[['Index (HTO)']]
  
  mat.files.10x <- file.path(dir.proj,'scRNAseq','03_PipelineOut',com.ID,'10x',pool.10x,'res')
  mat.files.hto <- file.path(dir.proj,'scRNAseq','03_PipelineOut',com.ID,'hto',pool.hto,'res')

  counts.all <- load_fry(frydir = mat.files.10x, output_list = T)
  counts.spl <- counts.all[['Spliced']]
  counts.uns <- counts.all[['Unspliced']]
  #counts.hto <- load_fry(frydir = mat.files.hto)

  tmp.counts.spl.called <- sum(Matrix::colSums(counts.spl[,colnames(seur.list[[i]])]))
  tmp.counts.spl.unlled <- sum(Matrix::colSums(counts.spl[,colnames(counts.spl) %!in% colnames(seur.list[[i]])]))
  tmp.counts.uns.called <- sum(Matrix::colSums(counts.uns[,colnames(seur.list[[i]])]))
  tmp.counts.uns.unlled <- sum(Matrix::colSums(counts.uns[,colnames(counts.uns) %!in% colnames(seur.list[[i]])]))
  tmp.counts.tot.called <- tmp.counts.spl.called + tmp.counts.uns.called
  tmp.counts.tot.unlled <- tmp.counts.spl.unlled + tmp.counts.uns.unlled
  
  tmp.percen.spl.called <- round(tmp.counts.spl.called/tmp.counts.tot.called*100)
  tmp.percen.spl.unlled <- round(tmp.counts.spl.unlled/tmp.counts.tot.unlled*100)
  tmp.percen.uns.called <- round(tmp.counts.uns.called/tmp.counts.tot.called*100)
  tmp.percen.uns.unlled <- round(tmp.counts.uns.unlled/tmp.counts.tot.unlled*100)
  
  tmp.df <- data.frame(ReadType = c('Spliced','Unspliced','Spliced','Unspliced'),
                       CellType = c('Called','Called','Uncalled','Uncalled'),
                       CountSum = c(tmp.counts.spl.called,tmp.counts.uns.called,tmp.counts.spl.unlled,tmp.counts.uns.unlled),
                       CountPer = c(tmp.percen.spl.called,tmp.percen.uns.called,tmp.percen.spl.unlled,tmp.percen.uns.unlled),
                       CountTot = c(tmp.counts.tot.called,tmp.counts.tot.called,tmp.counts.tot.unlled,tmp.counts.tot.unlled),
                       Pool = c(pool.10x,pool.10x,pool.10x,pool.10x))
  
  df <- rbind(df,tmp.df)
  
}

df[['Labelpos']] <- ifelse(df[['ReadType']] == 'Unspliced',
                           df[['CountSum']]/2, df[['CountTot']] - df[['CountSum']]/2)

p <- ggplot(data = df, aes(x = Pool, y = CountSum, fill = ReadType)) +
  geom_bar(stat = 'identity') + 
  geom_text(aes(label = paste0(CountPer,'%'),y = Labelpos),size = 3) +
  coord_flip() + ggtitle('Spliced/Unspliced - Counts') +
  xlab('Index') + ylab('Summed counts') + scale_fill_discrete(name = 'Read type') +
  facet_wrap(~CellType, nrow = 2) +
  theme_minimal()

ggsave(filename = paste0('readDistribution_alevin.png'),
       plot = p,
       path = dir.outs.qc.plots,
       width = 12,
       height = 8)


################################################################################
##########                        Create HTML                         ##########
################################################################################

cat('#\tCreating summary..\n',
    sep = '')
  
# Merge objects - Again
if (length(indi.files) > 1) {
  seur.comb <- merge(seur.list[[1]],
                     seur.list[2:length(seur.list)])
}else{
  seur.comb <- seur.list[[1]]
}

remove(list = c('seur.list'))

# Filter genes with no detection
rowsums <- Matrix::rowSums(seur.comb[['RNA']]@counts)
num.genes <- length(rownames(seur.comb))
seur.comb <- seur.comb[!rowsums == 0,]
num.genes.red <- length(rownames(seur.comb))
  
if (var.wofl == '10x + HTO') {
  stat.tib.10x <- tibble(
    # Raw - reads from fastqs
    'Cells (Loaded)' = double(),
    'Reads (Raw)' = double(),
      
    # Barcode ranks plot
    'Plot (Barcode Ranks)' = ' ',
      
    # Called cells - (Barcode ranks)
    'Cells (Called)' = double(),  'Cells % Loaded (Called)' = character(),
    'Reads (Called)' = double(),  'Reads % Raw (Called)' = character(), 'Reads/Cell (Called)' = double(),
    'UMIs (Called)' = double(), 'UMIs/Cell (Called)' = double(),
    
    # Called singlets (inter-HTO)
    'Cells (inter-HTO)' = double(),  'Cells % Called (inter-HTO)' = character(),
    'Reads (inter-HTO)' = double(),  'Reads % Called (inter-HTO)' = character(), 'Reads/Cell (inter-HTO)' = double(),
    'UMIs (inter-HTO)' = double(), 'UMIs/Cell (inter-HTO)' = double(),
    
    # Called singlets (intra-HTO)
    'Cells (intra-HTO)' = double(),  'Cells % inter-HTO (intra-HTO)' = character(),
    'Reads (intra-HTO)' = double(),  'Reads % inter-HTO (intra-HTO)' = character(), 'Reads/Cell (intra-HTO)' = double(),
    'UMIs (intra-HTO)' = double(), 'UMIs/Cell (intra-HTO)' = double()
    )
  
  
  knee.df.10x.x <- list()
  knee.df.10x.y <- list()
  
  
  stat.tib.hto <- tibble(
    # Raw - reads from fastqs
    'Cells (Loaded)' = double(),
    'Reads (Raw)' = double(),
    
    # UMI vs UMI plot
    'Plot (UMIvsUMI) - cal' = ' ',
    'Plot (UMIvsUMI) - unc' = ' ',
      
    # Called cells - (Barcode ranks)
    'Cells (Called)' = double(),  'Cells % Loaded (Called)' = character(),
    'Reads (Called)' = double(),  'Reads % Raw (Called)' = character(), 'Reads/Cell (Called)' = double(),
    'UMIs (Called)' = double(), 'UMIs/Cell (Called)' = double(),
    
    # Called singlets (inter-HTO)
    'Cells (inter-HTO)' = double(),  'Cells % Called (inter-HTO)' = character(), 
    'Reads (inter-HTO)' = double(),  'Reads % Called (inter-HTO)' = character(), 'Reads/Cell (inter-HTO)' = double(),
    'UMIs (inter-HTO)' = double(), 'UMIs/Cell (inter-HTO)' = double(),
    'Cells (Negative %)' = character(),    
    'Cells (Doublet %)' = character(),    
    'Cells (Singlet %)' = character(),
    'Reads - HTO (Negative %)' = character(),    
    'Reads - HTO (Doublet %)' = character(),    
    'Reads - HTO (Singlet %)' = character(),
    'Reads - RNA (Negative %)' = character(),    
    'Reads - RNA (Doublet %)' = character(),    
    'Reads - RNA (Singlet %)' = character()
    )
    
  umi.df.hto.cal.x <- list()
  umi.df.hto.cal.y <- list()
  umi.df.hto.unc.x <- list()
  umi.df.hto.unc.y <- list()
   
  for (pool.i in seq(n.pools)) {
    tmp.pool      <- pool.table[pool.i,]

      #pool.10x <- pool.table[i,][['Pool..10x.']]
      #pool.hto <- pool.table[i,][['Pool..HTO.']]
      #load.cells <- pool.table[i,][['Loaded.Cells']]
      
      
      #pins.10x      <- str_split(tmp.pool[['BCL PIN (10x)']], ',', simplify = T)[1,]
      
    pool.10x      <- str_split(tmp.pool[['Index (10x)']], ',', simplify = T)[1,]
    pool.hto      <- str_split(tmp.pool[['Index (HTO)']], ',', simplify = T)[1,]
      
      #seqs.10x      <- listLen(pins.10x)
      
    
      # Init
    stats.10x <- list()
    stats.hto <- list()
    
    mat.files.10x <- file.path(dir.proj,'scRNAseq','03_PipelineOut',com.ID,'10x',pool.10x,'res')
    mat.files.hto <- file.path(dir.proj,'scRNAseq','03_PipelineOut',com.ID,'hto',pool.hto,'res')
      
    # Collect barcode ranks data
    counts.10x <- load_fry(frydir = mat.files.10x, which_counts = c('U'))
    bc.calls.10x <- barcodeRanks(counts.10x)
    uniq.10x <- !duplicated(bc.calls.10x[['rank']])
      
    knee.df.10x.x[[pool.10x]] <- log(bc.calls.10x[['rank']][uniq.10x])
    knee.df.10x.y[[pool.10x]] <- log(bc.calls.10x[['total']][uniq.10x])
      
    knee.df.10x.x[[pool.10x]][is.infinite(knee.df.10x.x[[pool.10x]])] <- NA
    knee.df.10x.y[[pool.10x]][is.infinite(knee.df.10x.y[[pool.10x]])] <- NA
      
    stats.10x[['Plot (Barcode Ranks)']] <- ' '
      
      ###
      
    cells.ranks <- emptyDropsmodal(q = counts.10x,
                                   verbose = F,
                                   plot = F,
                                   format = 'noSave', skipModCheck = T)
      
      ###
      
    # Read featureDump.txt
    featDump.10x <- suppressMessages(read_delim(file.path(mat.files.10x,'featureDump.txt'), delim = '\t'))
    featDump.hto <- suppressMessages(read_delim(file.path(mat.files.hto,'featureDump.txt'), delim = '\t'))
      
    # Collect UMI vs UMI data
    common.CBs <- intersect(featDump.10x[['CB']],featDump.hto[['CB']])
    
    common.CBs.cal <- intersect(common.CBs,cells.ranks)
    common.CBs.unc <- setdiff(common.CBs,cells.ranks)
    
    
    featDump.10x %>% 
      filter(CB %in% common.CBs) %>% 
      arrange(CB) -> featDump.10x.sub
    
    featDump.10x.sub.cal <- filter(featDump.10x.sub, CB %in% common.CBs.cal)
    featDump.10x.sub.unc <- filter(featDump.10x.sub, CB %in% common.CBs.unc)
    
    uniq.10x.cal <- !duplicated(featDump.10x.sub.cal[['DeduplicatedReads']])
    uniq.10x.unc <- !duplicated(featDump.10x.sub.unc[['DeduplicatedReads']])
    
    featDump.hto %>% 
      filter(CB %in% common.CBs) %>% 
      arrange(CB) -> featDump.hto.sub
    
    featDump.hto.sub.cal <- filter(featDump.hto.sub, CB %in% common.CBs.cal)
    featDump.hto.sub.unc <- filter(featDump.hto.sub, CB %in% common.CBs.unc)
    
    umi.df.hto.cal.x[[pool.10x]] <- log(featDump.10x.sub.cal[['DeduplicatedReads']][uniq.10x.cal])
    umi.df.hto.cal.y[[pool.10x]] <- log(featDump.hto.sub.cal[['DeduplicatedReads']][uniq.10x.cal])
    umi.df.hto.unc.x[[pool.10x]] <- log(featDump.10x.sub.unc[['DeduplicatedReads']][uniq.10x.unc])
    umi.df.hto.unc.y[[pool.10x]] <- log(featDump.hto.sub.unc[['DeduplicatedReads']][uniq.10x.unc])
    
    umi.df.hto.cal.x[[pool.10x]][is.infinite(umi.df.hto.cal.x[[pool.10x]])] <- NA
    umi.df.hto.cal.y[[pool.10x]][is.infinite(umi.df.hto.cal.y[[pool.10x]])] <- NA
    umi.df.hto.unc.x[[pool.10x]][is.infinite(umi.df.hto.unc.x[[pool.10x]])] <- NA
    umi.df.hto.unc.y[[pool.10x]][is.infinite(umi.df.hto.unc.y[[pool.10x]])] <- NA
    
    stats.hto[['Plot (UMIvsUMI) - cal']] <- ' '
    stats.hto[['Plot (UMIvsUMI) - unc']] <- ' '
    
    #stats.10x[['Reads (Raw)']] <- read.df[read.df[['SampleId']] == pool.10x,][[2]]
    #stats.hto[['Reads (Raw)']] <- read.df.hto[read.df.hto[['SampleId']] == pool.hto,][[2]]
    stats.10x[['Reads (Raw)']] <- tmp.pool[['READS (10x)']]
    stats.hto[['Reads (Raw)']] <- tmp.pool[['READS (HTO)']]

    #stats.10x[['Cells (Loaded)']] <- load.cells
    #stats.hto[['Cells (Loaded)']] <- load.cells
    stats.10x[['Cells (Loaded)']] <- tmp.pool[['Loaded Cells']]
    stats.hto[['Cells (Loaded)']] <- tmp.pool[['Loaded Cells']]
    
    ## Read featureDump.txt
    #featDump.10x <- suppressMessages(read_delim(file.path(mat.files.10x,'featureDump.txt'), delim = '\t'))
    #featDump.hto <- suppressMessages(read_delim(file.path(mat.files.hto,'featureDump.txt'), delim = '\t'))
    
    # CALLED CELLS - BARCODE RANKS
    Idents(seur.comb) <- 'orig.ident'
    seur.comb.sub <- subset(seur.comb, idents = pool.10x)
    #cells.ranks <- sapply(strsplit(colnames(seur.comb.sub),'[_]'), '[', 1)
    
    # 10x
    stats.10x[['Cells (Called)']] <- length(unique(featDump.10x[featDump.10x[['CB']] %in% cells.ranks,][['CB']]))
    stats.10x[['Cells % Loaded (Called)']] <- paste(as.character(round(stats.10x[['Cells (Called)']] / stats.10x[['Cells (Loaded)']] * 100, 1)),'%')
    
    stats.10x[['Reads (Called)']] <- sum(featDump.10x[featDump.10x[['CB']] %in% cells.ranks,][['MappedReads']])
    stats.10x[['Reads % Raw (Called)']] <- paste(as.character(round(stats.10x[['Reads (Called)']] / stats.10x[['Reads (Raw)']] * 100, 1)),'%')
    stats.10x[['Reads/Cell (Called)']] <- round(median(featDump.10x[featDump.10x[['CB']] %in% cells.ranks,][['MappedReads']]))
    
    stats.10x[['UMIs (Called)']] <- round(sum(featDump.10x[featDump.10x[['CB']] %in% cells.ranks,][['DeduplicatedReads']]))
    stats.10x[['UMIs/Cell (Called)']] <- round(median(featDump.10x[featDump.10x[['CB']] %in% cells.ranks,][['DeduplicatedReads']]))
    
    # HTO
    stats.hto[['Cells (Called)']] <- length(unique(featDump.hto[featDump.hto[['CB']] %in% cells.ranks,][['CB']]))
    stats.hto[['Cells % Loaded (Called)']] <- paste(as.character(round(stats.hto[['Cells (Called)']] / stats.hto[['Cells (Loaded)']] * 100, 1)),'%')
    
    stats.hto[['Reads (Called)']] <- sum(featDump.hto[featDump.hto[['CB']] %in% cells.ranks,][['MappedReads']])
    stats.hto[['Reads % Raw (Called)']] <- paste(as.character(round(stats.hto[['Reads (Called)']] / stats.hto[['Reads (Raw)']] * 100, 1)),'%')
    stats.hto[['Reads/Cell (Called)']] <- round(median(featDump.hto[featDump.hto[['CB']] %in% cells.ranks,][['MappedReads']]))
    
    stats.hto[['UMIs (Called)']] <- round(sum(featDump.hto[featDump.hto[['CB']] %in% cells.ranks,][['DeduplicatedReads']]))
    stats.hto[['UMIs/Cell (Called)']] <- round(median(featDump.hto[featDump.hto[['CB']] %in% cells.ranks,][['DeduplicatedReads']]))
    
    
    # CALLED SINGLETS (INTER-HTO)
    Idents(seur.comb.sub) <- 'HTO_classification.global'
    seur.comb.neg <- subset(seur.comb.sub, idents = 'Negative')
    seur.comb.dou <- subset(seur.comb.sub, idents = 'Doublet')
    seur.comb.sub <- subset(seur.comb.sub, idents = 'Singlet')
    
    cells.negat.inter <- sapply(strsplit(colnames(seur.comb.neg),'[_]'), '[', 1)
    cells.doubl.inter <- sapply(strsplit(colnames(seur.comb.dou),'[_]'), '[', 1)
    cells.singl.inter <- sapply(strsplit(colnames(seur.comb.sub),'[_]'), '[', 1)
    
    stats.hto[['Cells (Negative %)']] <- paste(as.character(round(length(cells.negat.inter) / length(cells.ranks) * 100, 1)),'%')
    stats.hto[['Cells (Doublet %)']] <- paste(as.character(round(length(cells.doubl.inter) / length(cells.ranks) * 100, 1)),'%')
    stats.hto[['Cells (Singlet %)']] <- paste(as.character(round(length(cells.singl.inter) / length(cells.ranks) * 100, 1)),'%')
    
    stats.hto[['Reads - HTO (Negative %)']] <- paste(as.character(round(sum(featDump.hto[featDump.hto[['CB']] %in% cells.negat.inter,][['MappedReads']]) / stats.hto[['Reads (Called)']] * 100, 1)),'%')
    stats.hto[['Reads - HTO (Doublet %)']] <- paste(as.character(round(sum(featDump.hto[featDump.hto[['CB']] %in% cells.doubl.inter,][['MappedReads']]) / stats.hto[['Reads (Called)']] * 100, 1)),'%')
    stats.hto[['Reads - HTO (Singlet %)']] <- paste(as.character(round(sum(featDump.hto[featDump.hto[['CB']] %in% cells.singl.inter,][['MappedReads']]) / stats.hto[['Reads (Called)']] * 100, 1)),'%')
    
    stats.hto[['Reads - RNA (Negative %)']] <- paste(as.character(round(sum(featDump.10x[featDump.10x[['CB']] %in% cells.negat.inter,][['MappedReads']]) / stats.10x[['Reads (Called)']] * 100, 1)),'%')
    stats.hto[['Reads - RNA (Doublet %)']] <- paste(as.character(round(sum(featDump.10x[featDump.10x[['CB']] %in% cells.doubl.inter,][['MappedReads']]) / stats.10x[['Reads (Called)']] * 100, 1)),'%')
    stats.hto[['Reads - RNA (Singlet %)']] <- paste(as.character(round(sum(featDump.10x[featDump.10x[['CB']] %in% cells.singl.inter,][['MappedReads']]) / stats.10x[['Reads (Called)']] * 100, 1)),'%')
    
    # 10x
    stats.10x[['Cells (inter-HTO)']] <- length(unique(featDump.10x[featDump.10x[['CB']] %in% cells.singl.inter,][['CB']]))
    stats.10x[['Cells % Called (inter-HTO)']] <- paste(as.character(round(stats.10x[['Cells (inter-HTO)']] / stats.10x[['Cells (Called)']] * 100, 1)),'%')
    
    stats.10x[['Reads (inter-HTO)']] <- sum(featDump.10x[featDump.10x[['CB']] %in% cells.singl.inter,][['MappedReads']])
    stats.10x[['Reads % Called (inter-HTO)']] <- paste(as.character(round(stats.10x[['Reads (inter-HTO)']] / stats.10x[['Reads (Called)']] * 100, 1)),'%')
    stats.10x[['Reads/Cell (inter-HTO)']] <- round(median(featDump.10x[featDump.10x[['CB']] %in% cells.singl.inter,][['MappedReads']]))
    
    stats.10x[['UMIs (inter-HTO)']] <- round(sum(featDump.10x[featDump.10x[['CB']] %in% cells.singl.inter,][['DeduplicatedReads']]))
    stats.10x[['UMIs/Cell (inter-HTO)']] <- round(median(featDump.10x[featDump.10x[['CB']] %in% cells.singl.inter,][['DeduplicatedReads']]))
    
    # HTO
    stats.hto[['Cells (inter-HTO)']] <- length(unique(featDump.hto[featDump.hto[['CB']] %in% cells.singl.inter,][['CB']]))
    stats.hto[['Cells % Called (inter-HTO)']] <- paste(as.character(round(stats.hto[['Cells (inter-HTO)']] / stats.hto[['Cells (Called)']] * 100, 1)),'%')
    
    stats.hto[['Reads (inter-HTO)']] <- sum(featDump.hto[featDump.hto[['CB']] %in% cells.singl.inter,][['MappedReads']])
    stats.hto[['Reads % Called (inter-HTO)']] <- paste(as.character(round(stats.hto[['Reads (inter-HTO)']] / stats.hto[['Reads (Called)']] * 100, 1)),'%')
    stats.hto[['Reads/Cell (inter-HTO)']] <- round(median(featDump.hto[featDump.hto[['CB']] %in% cells.singl.inter,][['MappedReads']]))
    
    stats.hto[['UMIs (inter-HTO)']] <- round(sum(featDump.hto[featDump.hto[['CB']] %in% cells.singl.inter,][['DeduplicatedReads']]))
    stats.hto[['UMIs/Cell (inter-HTO)']] <- round(median(featDump.hto[featDump.hto[['CB']] %in% cells.singl.inter,][['DeduplicatedReads']]))
    
    
    # CALLED SINGLETS (INTRA-HTO)
    
    # Collect all inferential doublets
    Idents(seur.comb.sub) <- 'predicted_dub_std'
    cells.remove.predDubStd <- WhichCells(seur.comb.sub, idents = TRUE)
    Idents(seur.comb.sub) <- 'predicted_dub_cut'
    cells.remove.predDubCut <- WhichCells(seur.comb.sub, idents = TRUE)
    cells.remove.predDub <- union(cells.remove.predDubStd,cells.remove.predDubCut)
    seur.comb.sub[['predicted_dub_all']] <- colnames(seur.comb.sub) %in% cells.remove.predDub
    
    Idents(seur.comb.sub) <- 'predicted_dub_all'
    seur.comb.sub <- subset(seur.comb.sub, idents = FALSE)
    cells.singl.intra <- sapply(strsplit(colnames(seur.comb.sub),'[_]'), '[', 1)
    
    # 10x
    stats.10x[['Cells (intra-HTO)']] <- length(unique(featDump.10x[featDump.10x[['CB']] %in% cells.singl.intra,][['CB']]))
    stats.10x[['Cells % inter-HTO (intra-HTO)']] <- paste(as.character(round(stats.10x[['Cells (intra-HTO)']] / stats.10x[['Cells (inter-HTO)']] * 100, 1)),'%')
    
    stats.10x[['Reads (intra-HTO)']] <- sum(featDump.10x[featDump.10x[['CB']] %in% cells.singl.intra,][['MappedReads']])
    stats.10x[['Reads % inter-HTO (intra-HTO)']] <- paste(as.character(round(stats.10x[['Reads (intra-HTO)']] / stats.10x[['Reads (inter-HTO)']] * 100, 1)),'%')
    stats.10x[['Reads/Cell (intra-HTO)']] <- round(median(featDump.10x[featDump.10x[['CB']] %in% cells.singl.intra,][['MappedReads']]))
    
    stats.10x[['UMIs (intra-HTO)']] <- round(sum(featDump.10x[featDump.10x[['CB']] %in% cells.singl.intra,][['DeduplicatedReads']]))
    stats.10x[['UMIs/Cell (intra-HTO)']] <- round(median(featDump.10x[featDump.10x[['CB']] %in% cells.singl.intra,][['DeduplicatedReads']]))

    # Bind to collective df
    stat.tib.10x <- bind_rows(stat.tib.10x,stats.10x)
    stat.tib.hto <- bind_rows(stat.tib.hto,stats.hto)
  }
    

  stat.df.10x <- as.data.frame(stat.tib.10x)
  rownames(stat.df.10x) <- pool.table[['Index (10x)']]
    
  tmp.df.10x <- stat.df.10x
    
    
  tmp.df.10x %>%
    kbl(escape = F, col.names = c('','','',
                                  'remaining','% loaded','remaining','% sequenced','per cell','remaining','per cell',
                                  'remaining','% called','remaining','% called','per cell','remaining','per cell',
                                  'remaining','% singlets','remaining','% singlets','per cell','remaining','per cell'),
        format.args = list(big.mark = ','),align = rep('c',25),
        booktabs = T) %>% 
    row_spec(0, color = '#dedede', angle = 25, align = 'center', font_size = '12') %>%
    kable_material_dark(full_width = F, html_font = 'helvetica', lightable_options = 'basic') %>%
    column_spec(1, bold = T, color = '#dedede') %>% 
    column_spec(4, image = spec_plot(x = knee.df.10x.x,
                                     y = knee.df.10x.y,
                                     width = 250,
                                     height = 250,
                                     type = "p")) %>% 
    column_spec(column = c(5,6,7,8,9,10,11,19,20,21,22,23,24,25), color = '#ffbd69') %>% 
    column_spec(column = c(12,13,14,15,16,17,18), color = '#adffff') %>% 
    add_header_above(c(
      ' ' = 1,
      'Cells' = 1,
      'Reads' = 1,
      'Barcode Rank Plot' = 1,
      'Cells' = 2, 'Reads' = 3, 'UMIs' = 2,
      'Cells' = 2, 'Reads' = 3, 'UMIs' = 2,
      'Cells' = 2, 'Reads' = 3, 'UMIs' = 2),
      color = '#dedede') %>% 
    add_header_above(c(' ' = 1,
                       'Loaded' = 1,
                       'Sequenced' = 1,
                       ' ' = 1,
                       'Called Cells' = 7,
                       'Between-HTO doublets and negatives removed' = 7,
                       'Within-HTO doublets removed' = 7),
                     align = 'center',
                     bold = T, color = c('#dedede','#dedede','#dedede','#dedede',
                                         '#ffbd69','#adffff','#ffbd69','#adffff','#ffbd69','#adffff')) %>% 
    add_header_above(c('Summary stats (RNA expression)' = 4,
                       ' ' = 21),
                     align = 'center',
                     color = '#dedede',
                     font_size = 20,
                     bold = T) %>% 
    save_kable(file = file.path(dir.outs.qc.summa,'rna_summary.html'), self_contained = T)
  
  
  # HTO time
  stat.df.hto <- as.data.frame(stat.tib.hto)
  rownames(stat.df.hto) <- paste0(pool.table[['Index (HTO)']],' (',pool.table[['Index (10x)']],')')

  tmp.df.hto <- stat.df.hto

  
  tmp.df.hto %>%
    kbl(escape = F, col.names = c('','','filtered CBs','called CBs',
                                  'remaining','% loaded','remaining','% sequenced','per cell','remaining','per cell',
                                  'remaining','% called','remaining','% called','per cell','remaining','per cell',
                                  'Negative %','Doublet %','Singlet %',
                                  'Negative %','Doublet %','Singlet %',
                                  'Negative %','Doublet %','Singlet %'),
        format.args = list(big.mark = ','),
        booktabs = T) %>% 
    row_spec(0, color = '#dedede', angle = 25,align = 'center', font_size = '12') %>%
    kable_material_dark(full_width = F, html_font = 'helvetica', lightable_options = 'basic') %>%
    column_spec(1, bold = T, color = '#dedede') %>%
    column_spec(4, image = spec_plot(x = umi.df.hto.unc.x,
                                     y = umi.df.hto.unc.y,
                                     width = 125,
                                     height = 250,
                                     type = "p")) %>%
    column_spec(5, image = spec_plot(x = umi.df.hto.cal.x,
                                     y = umi.df.hto.cal.y,
                                     width = 125,
                                     height = 250,
                                     type = "p")) %>% 
    column_spec(column = c(6,7,8,9,10,11,12,20,21,22,23,24,25,26,27,28), color = '#ffbd69') %>% 
    column_spec(column = c(13,14,15,16,17,18,19), color = '#adffff') %>%
    add_header_above(c(' ' = 1,
                       'Cells' = 1,
                       'Reads' = 1,
                       'UMI Plot (RNA/HTO)' = 2,
                       'Cells' = 2, 'Reads' = 3, 'UMIs' = 2,
                       'Cells' = 2, 'Reads' = 3, 'UMIs' = 2,
                       'Cells' = 3, 'Reads (HTO)' = 3, 'Reads (RNA)' = 3),
                     color = '#dedede') %>% 
    add_header_above(c(' ' = 1,
                       'Loaded' = 1,
                       'Sequenced' = 1,
                       ' ' = 2,
                       'Called Cells' = 7,
                       'Between-HTO doublets and negatives removed' = 7,
                       'Classification proportions' = 9),
                     align = 'center',
                     bold = T, color = c('#dedede','#dedede','#dedede','#dedede',
                                         '#ffbd69','#adffff','#ffbd69')) %>% 
    add_header_above(c('Summary stats (HTO expression)' = 4,
                       ' ' = 24),
                     align = 'center',
                     color = '#dedede',
                     font_size = 20,
                     bold = T) %>% 
    save_kable(file = file.path(dir.outs.qc.summa,'hto_summary.html'), self_contained = T)
}
  
if (var.wofl == '10x') {
  stat.tib.10x <- tibble(
    # Raw - reads from fastqs
    'Cells (Loaded)' = double(),
    'Reads (Raw)' = double(),
    
    # Barcode ranks plot
    'Plot (Barcode Ranks)' = ' ',
    
    # Called cells - (Barcode ranks)
    'Cells (Called)' = double(),  'Cells % Loaded (Called)' = character(),
    'Reads (Called)' = double(),  'Reads % Raw (Called)' = character(), 'Reads/Cell (Called)' = double(),
    'UMIs (Called)' = double(), 'UMIs/Cell (Called)' = double(),
    )
  
  knee.df.10x.x <- list()
  knee.df.10x.y <- list()
  
  for (pool.i in seq(n.pools)) {
    tmp.pool      <- pool.table[pool.i,]
    pool.10x      <- str_split(tmp.pool[['Index (10x)']], ',', simplify = T)[1,]
    
    # Init
    stats.10x <- list()
    
    mat.files.10x <- file.path(dir.proj,'scRNAseq','03_PipelineOut',com.ID,'10x',pool.10x,'res')
    
    # Collect barcode ranks data
    counts.10x <- load_fry(frydir = mat.files.10x, which_counts = c('U'))
    bc.calls.10x <- barcodeRanks(counts.10x)
    uniq.10x <- !duplicated(bc.calls.10x[['rank']])
    
    knee.df.10x.x[[pool.10x]] <- log(bc.calls.10x[['rank']][uniq.10x])
    knee.df.10x.y[[pool.10x]] <- log(bc.calls.10x[['total']][uniq.10x])
    
    knee.df.10x.x[[pool.10x]][is.infinite(knee.df.10x.x[[pool.10x]])] <- NA
    knee.df.10x.y[[pool.10x]][is.infinite(knee.df.10x.y[[pool.10x]])] <- NA
    
    stats.10x[['Plot (Barcode Ranks)']] <- ' '
    
    cells.ranks <- emptyDropsmodal(q = counts.10x,
                                   verbose = F,
                                   plot = F,
                                   format = 'noSave')
    
    # Read featureDump.txt
    featDump.10x <- suppressMessages(read_delim(file.path(mat.files.10x,'featureDump.txt'), delim = '\t'))
    
    stats.10x[['Reads (Raw)']] <- tmp.pool[['READS (10x)']]
    stats.10x[['Cells (Loaded)']] <- tmp.pool[['Loaded Cells']]
    
    # CALLED CELLS - BARCODE RANKS
    Idents(seur.comb) <- 'orig.ident'
    seur.comb.sub <- subset(seur.comb, idents = pool.10x)
    
    # 10x
    stats.10x[['Cells (Called)']] <- length(unique(featDump.10x[featDump.10x[['CB']] %in% cells.ranks,][['CB']]))
    stats.10x[['Cells % Loaded (Called)']] <- paste(as.character(round(stats.10x[['Cells (Called)']] / stats.10x[['Cells (Loaded)']] * 100, 1)),'%')
    
    stats.10x[['Reads (Called)']] <- sum(featDump.10x[featDump.10x[['CB']] %in% cells.ranks,][['MappedReads']])
    stats.10x[['Reads % Raw (Called)']] <- paste(as.character(round(stats.10x[['Reads (Called)']] / stats.10x[['Reads (Raw)']] * 100, 1)),'%')
    stats.10x[['Reads/Cell (Called)']] <- round(median(featDump.10x[featDump.10x[['CB']] %in% cells.ranks,][['MappedReads']]))
    
    stats.10x[['UMIs (Called)']] <- round(sum(featDump.10x[featDump.10x[['CB']] %in% cells.ranks,][['DeduplicatedReads']]))
    stats.10x[['UMIs/Cell (Called)']] <- round(median(featDump.10x[featDump.10x[['CB']] %in% cells.ranks,][['DeduplicatedReads']]))

    # Bind to collective df
    stat.tib.10x <- bind_rows(stat.tib.10x,stats.10x)
  }
  
  stat.df.10x <- as.data.frame(stat.tib.10x)
  rownames(stat.df.10x) <- pool.table[['Index (10x)']]
  
  tmp.df.10x <- stat.df.10x
  
  
  tmp.df.10x %>%
    kbl(escape = F, col.names = c('','','',
                                  'remaining','% loaded','remaining','% sequenced','per cell','remaining','per cell'),
        format.args = list(big.mark = ','),align = rep('c',11),
        booktabs = T) %>% 
    row_spec(0, color = '#dedede', angle = 25, align = 'center', font_size = '12') %>%
    kable_material_dark(full_width = F, html_font = 'helvetica', lightable_options = 'basic') %>%
    column_spec(1, bold = T, color = '#dedede') %>% 
    column_spec(4, image = spec_plot(x = knee.df.10x.x,
                                     y = knee.df.10x.y,
                                     width = 250,
                                     height = 250,
                                     type = "p")) %>% 
    column_spec(column = c(5,6,7,8,9,10,11), color = '#ffbd69') %>% 
    #column_spec(column = c(12,13,14,15,16,17,18), color = '#adffff') %>% 
    add_header_above(c(
      ' ' = 1,
      'Cells' = 1,
      'Reads' = 1,
      'Barcode Rank Plot' = 1,
      'Cells' = 2, 'Reads' = 3, 'UMIs' = 2),
      color = '#dedede') %>% 
    add_header_above(c(' ' = 1,
                       'Loaded' = 1,
                       'Sequenced' = 1,
                       ' ' = 1,
                       'Called Cells' = 7),
                     align = 'center',
                     bold = T, color = c('#dedede','#dedede','#dedede','#dedede',
                                         '#ffbd69','#adffff','#ffbd69','#adffff','#ffbd69','#adffff')) %>% 
    add_header_above(c('Summary stats (RNA expression)' = 4,
                       ' ' = 7),
                     align = 'center',
                     color = '#dedede',
                     font_size = 20,
                     bold = T) %>% 
    save_kable(file = file.path(dir.outs.qc.summa,'rna_summary.html'), self_contained = T)
    
  }
  
  cat('#\n',
      '#\n',
      '###\n',
      rep('#',80),'\n',
      sep = '')
  
  sink()




################################################################################
################################################################################
################################################################################



file.create(snakemake@output[[1]])
