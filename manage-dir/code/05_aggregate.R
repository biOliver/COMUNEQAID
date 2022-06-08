################################################################################
##########                          Packages                          ##########
################################################################################
suppressPackageStartupMessages({
  library(foreach)
  library(parallel)
  library(doParallel)
  library(Seurat)
  library(tidyverse)
  library(sctransform)
  library(DropletUtils)
  library(Matrix)
  library(mclust)
  library(kableExtra)
})

################################################################################
##########                            Init                            ##########
################################################################################
scop.ID     <- snakemake@params[['scopID']]
com.ID      <- snakemake@params[['comID']]
com.ID.list <- str_split(com.ID, ',', simplify = T)[1,]
user.ID     <- snakemake@params[['userID']]

source(file.path('/projects', user.ID, 'COMUNEQAID/manage-dir/code/shared_functions.R'))

dir.proj <- paste0('/', file.path('projects',user.ID,'COMUNEQAID','outs',scop.ID))
dir.outs.qc <- file.path(dir.proj, 'scRNAseq', '00_QC')

cl <- makeCluster(length(com.ID.list)) 
registerDoParallel(cl)
foreach(i = seq(com.ID.list),
        .combine = 'c',
        .packages = c('Seurat','tidyverse','sctransform','DropletUtils',
                      'Matrix','mclust','kableExtra')) %dopar% {

  com.ID <- com.ID.list[i]
  unique.ID <- paste(scop.ID,com.ID,sep = '_')
  
  dir.data <- paste0('/', file.path('projects',user.ID,'COMUNEQAID','manage-dir','tmp-data',unique.ID))
  dir.outs.log <- file.path(dir.proj, 'scRNAseq', '04_Log', com.ID)
  dir.outs.qc.plots <- file.path(dir.outs.qc, com.ID, 'plots')
  dir.outs.qc.summa <- file.path(dir.outs.qc, com.ID, 'summary')
  
  dir.outs.indi <- file.path(dir.proj, 'Output', 'data', 'reaction-unfiltered', com.ID)
  dir.outs.comb <- file.path(dir.proj, 'Output', 'data', 'aggregated-filtered', com.ID)
  
  var.wofl <- scan(file.path(dir.data,'tmp_workflow.txt'), what = '', sep = '\n', quiet = T)

  pool.table <- read.table(file.path(dir.data,'poolTable-updated.csv'))
  
  if (var.wofl == '10x') {
    colnames(pool.table) <- c('Index (10x)','Lane','Loaded Cells','BCL PIN (10x)',"SEQ NAMES (10x)","READS (10x)")
  }
  if (var.wofl == '10x + HTO') {
    colnames(pool.table) <- c('Index (10x)','Index (HTO)','Lane','Loaded Cells','BCL PIN (10x)','BCL PIN (HTO)',"SEQ NAMES (HTO)","READS (HTO)","SEQ NAMES (10x)","READS (10x)")
  }
  n.pools <- dim(pool.table)[1]
  
  dir.create(dir.outs.comb, recursive = T, showWarnings = F)
  dir.create(dir.outs.qc.summa, recursive = T, showWarnings = F)
  
  ################################################################################
  ##########                 Preparing aggregated output                ##########
  ################################################################################
  
  sink.file <- file(paste0('logs/05_aggregate_',com.ID,'.log'), 'a')
  sink(file = sink.file, append = T,
       type = 'output', split = T)
  cat(as.character(lubridate::now()),'\n', sep = '')
  cat(rep('#',80),'\n',
      '#####                                                                      #####\n',
      '#####                     Summarizing sequencing lanes                     #####\n',
      '#####                                                   ***** / ******     #####\n',
      rep('#',80),'\n',
      '#####\n',
      '##\n',
      '##\tCOM ID:         ',com.ID,'\n',
      '##\n',
      '###\n',
      rep('#',80),'\n',
      '####\n',
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
  seur.list <- lapply(file.path(dir.outs.indi, indi.files), readRDS)
  
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
    ggsave(filename = paste0(com.ID,'_cells_remove_inf-doubl.png'),
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
    Idents(seur.comb) <- 'HTO_mcl_classification.global'
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
  
  Idents(seur.comb) <- 'hash.mcl.ID'
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
  #dims <- round(
   # as.numeric(
   #   intrinsicDimension::maxLikGlobalDimEst(
    #    data = seur.comb@reductions[['pca']][, 1:50],
     #   k = 20)))
  dims <- 50
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
    p <- DimPlot(seur.comb, group.by = 'hash.mcl.ID') + ggtitle(paste0('Samples (',dims,' PCs)'))
  }else {
    p <- DimPlot(seur.comb, group.by = 'orig.ident') + ggtitle(paste0('Samples (',dims,' PCs)'))
  }
  
  ggsave(filename = paste0(com.ID,'_RNA-UMAP.png'),
         plot = p,
         path = dir.outs.qc.plots,
         width = 8,
         height = 8)
  
  
  cat('#\n',
      '#\n',
      '#\tSaving objects..\n',
      sep = '')
  
  # Output objects
  out.name.rna <- paste0(com.ID,'_rna-seurat.rds')
  out.name.full <- paste0(com.ID,'_full-seurat.rds')
  
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
      
      pool.10x      <- str_split(tmp.pool[['Index (10x)']], ',', simplify = T)[1,]
      pool.hto      <- str_split(tmp.pool[['Index (HTO)']], ',', simplify = T)[1,]
      
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
      
      cells.ranks <- emptyDropsmodal(q = counts.10x,
                                     verbose = F,
                                     plot = F,
                                     format = 'noSave')
      
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
      
      stats.10x[['Reads (Raw)']] <- tmp.pool[['READS (10x)']]
      stats.hto[['Reads (Raw)']] <- tmp.pool[['READS (HTO)']]
      
      stats.10x[['Cells (Loaded)']] <- tmp.pool[['Loaded Cells']]
      stats.hto[['Cells (Loaded)']] <- tmp.pool[['Loaded Cells']]
      
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
      
      # HTO
      stats.hto[['Cells (Called)']] <- length(unique(featDump.hto[featDump.hto[['CB']] %in% cells.ranks,][['CB']]))
      stats.hto[['Cells % Loaded (Called)']] <- paste(as.character(round(stats.hto[['Cells (Called)']] / stats.hto[['Cells (Loaded)']] * 100, 1)),'%')
      
      stats.hto[['Reads (Called)']] <- sum(featDump.hto[featDump.hto[['CB']] %in% cells.ranks,][['MappedReads']])
      stats.hto[['Reads % Raw (Called)']] <- paste(as.character(round(stats.hto[['Reads (Called)']] / stats.hto[['Reads (Raw)']] * 100, 1)),'%')
      stats.hto[['Reads/Cell (Called)']] <- round(median(featDump.hto[featDump.hto[['CB']] %in% cells.ranks,][['MappedReads']]))
      
      stats.hto[['UMIs (Called)']] <- round(sum(featDump.hto[featDump.hto[['CB']] %in% cells.ranks,][['DeduplicatedReads']]))
      stats.hto[['UMIs/Cell (Called)']] <- round(median(featDump.hto[featDump.hto[['CB']] %in% cells.ranks,][['DeduplicatedReads']]))
      
      
      # CALLED SINGLETS (INTER-HTO)
      Idents(seur.comb.sub) <- 'HTO_mcl_classification.global'
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
  
  ################################################################################
  ##########                          QC plots                          ##########
  ################################################################################
  #remove(list = c('seur.list'))
  df <- data.frame(Type = character(),
                   CountSum = double(),
                   CountPer = double(),
                   CountTot = double(),
                   Pool = character())
  
  df.hto <- data.frame(Library = character(),
                       Class = character(),
                       CountPer = double(),
                       Labelpos = double(),
                       Pool = character())
  
  for (i in seq(seur.list)) {
    
    tmp.pool <- pool.table[i,]
    
    pool.10x <- tmp.pool[['Index (10x)']]
    
    mat.files.10x <- file.path(dir.proj,'scRNAseq','03_PipelineOut',com.ID,'10x',pool.10x,'res')
    
    counts.all <- load_fry(frydir = mat.files.10x, output_list = T)
    counts.spl <- counts.all[['Spliced']]
    counts.uns <- counts.all[['Unspliced']]
    
    tmp.counts.spl.called <- sum(Matrix::colSums(counts.spl[,colnames(seur.list[[i]])]))
    tmp.counts.spl.uncalled <- sum(Matrix::colSums(counts.spl[,colnames(counts.spl) %!in% colnames(seur.list[[i]])]))
    tmp.counts.uns.called <- sum(Matrix::colSums(counts.uns[,colnames(seur.list[[i]])]))
    tmp.counts.uns.uncalled <- sum(Matrix::colSums(counts.uns[,colnames(counts.uns) %!in% colnames(seur.list[[i]])]))
    tmp.counts.tot.called <- tmp.counts.spl.called + tmp.counts.uns.called
    tmp.counts.tot.uncalled <- tmp.counts.spl.uncalled + tmp.counts.uns.uncalled
    
    tmp.percen.spl.called <- round(tmp.counts.spl.called/tmp.counts.tot.called*100)
    tmp.percen.spl.uncalled <- round(tmp.counts.spl.uncalled/tmp.counts.tot.uncalled*100)
    tmp.percen.uns.called <- round(tmp.counts.uns.called/tmp.counts.tot.called*100)
    tmp.percen.uns.uncalled <- round(tmp.counts.uns.uncalled/tmp.counts.tot.uncalled*100)
    
    tmp.df <- data.frame(ReadType = c('Spliced','Unspliced','Spliced','Unspliced'),
                         CellType = c('Called','Called','Uncalled','Uncalled'),
                         CountSum = c(tmp.counts.spl.called,tmp.counts.uns.called,tmp.counts.spl.uncalled,tmp.counts.uns.uncalled),
                         CountPer = c(tmp.percen.spl.called,tmp.percen.uns.called,tmp.percen.spl.uncalled,tmp.percen.uns.uncalled),
                         CountTot = c(tmp.counts.tot.called,tmp.counts.tot.called,tmp.counts.tot.uncalled,tmp.counts.tot.uncalled),
                         Pool = c(pool.10x,pool.10x,pool.10x,pool.10x))
    
    df <- rbind(df,tmp.df)
    
    if (var.wofl == '10x + HTO') {
      pool.hto <- tmp.pool[['Index (HTO)']]
      mat.files.hto <- file.path(dir.proj,'scRNAseq','03_PipelineOut',com.ID,'hto',pool.hto,'res')
      counts.hto <- load_fry(frydir = mat.files.hto)
      counts.rna <- counts.all[['RNA']]
      
      Idents(seur.list[[i]]) <- 'HTO_mcl_classification.global'
      seur.comb.neg <- subset(seur.list[[i]], idents = 'Negative')
      seur.comb.sub <- subset(seur.list[[i]], idents = 'Singlet')
      seur.comb.dou <- subset(seur.list[[i]], idents = 'Doublet')
      
      bcs.negative <- sapply(strsplit(colnames(seur.comb.neg),'[_]'), '[', 1)
      bcs.singlet <- sapply(strsplit(colnames(seur.comb.sub),'[_]'), '[', 1)
      bcs.doublet <- sapply(strsplit(colnames(seur.comb.dou),'[_]'), '[', 1)
      
      tmp.counts.rna.called <- sum(Matrix::colSums(counts.rna[,colnames(seur.list[[i]])]))
      tmp.counts.rna.uncalled <- sum(Matrix::colSums(counts.rna[,colnames(counts.rna) %!in% colnames(seur.list[[i]])]))
      tmp.counts.rna.negative <- sum(Matrix::colSums(counts.rna[,bcs.negative]))
      tmp.counts.rna.singlet <- sum(Matrix::colSums(counts.rna[,bcs.singlet]))
      tmp.counts.rna.doublet <- sum(Matrix::colSums(counts.rna[,bcs.doublet]))
      tmp.counts.rna.total <- tmp.counts.rna.called + tmp.counts.rna.uncalled
      
      tmp.counts.hto.called <- sum(Matrix::colSums(counts.hto[,colnames(seur.list[[i]])]))
      tmp.counts.hto.uncalled <- sum(Matrix::colSums(counts.hto[,colnames(counts.hto) %!in% colnames(seur.list[[i]])]))
      tmp.counts.hto.negative <- sum(Matrix::colSums(counts.hto[,bcs.negative]))
      tmp.counts.hto.singlet <- sum(Matrix::colSums(counts.hto[,bcs.singlet]))
      tmp.counts.hto.doublet <- sum(Matrix::colSums(counts.hto[,bcs.doublet]))
      tmp.counts.hto.total <- tmp.counts.hto.called + tmp.counts.hto.uncalled
      
      tmp.percen.rna.uncalled <- round(tmp.counts.rna.uncalled/tmp.counts.rna.total*100)
      tmp.percen.rna.negative <- round(tmp.counts.rna.negative/tmp.counts.rna.total*100)
      tmp.percen.rna.singlet <- round(tmp.counts.rna.singlet/tmp.counts.rna.total*100)
      tmp.percen.rna.doublet <- round(tmp.counts.rna.doublet/tmp.counts.rna.total*100)
      
      tmp.percen.hto.uncalled <- round(tmp.counts.hto.uncalled/tmp.counts.hto.total*100)
      tmp.percen.hto.negative <- round(tmp.counts.hto.negative/tmp.counts.hto.total*100)
      tmp.percen.hto.singlet <- round(tmp.counts.hto.singlet/tmp.counts.hto.total*100)
      tmp.percen.hto.doublet <- round(tmp.counts.hto.doublet/tmp.counts.hto.total*100)
      
      tmp.labpos.rna.uncalled <- tmp.counts.rna.uncalled/2
      tmp.labpos.rna.negative <- tmp.counts.rna.uncalled + tmp.counts.rna.negative/2
      tmp.labpos.rna.singlet <- tmp.counts.rna.uncalled + tmp.counts.rna.negative + tmp.counts.rna.singlet/2
      tmp.labpos.rna.doublet <- tmp.counts.rna.uncalled + tmp.counts.rna.negative + tmp.counts.rna.singlet + tmp.counts.rna.doublet/2

      tmp.labpos.hto.uncalled <- tmp.counts.hto.uncalled/2
      tmp.labpos.hto.negative <- tmp.counts.hto.uncalled + tmp.counts.hto.negative/2
      tmp.labpos.hto.singlet <- tmp.counts.hto.uncalled + tmp.counts.hto.negative + tmp.counts.hto.singlet/2
      tmp.labpos.hto.doublet <- tmp.counts.hto.uncalled + tmp.counts.hto.negative + tmp.counts.hto.singlet + tmp.counts.hto.doublet/2
      
      tmp.df.hto <- data.frame(Library = c('RNA','RNA','RNA','RNA',
                                           'HTO','HTO','HTO','HTO'),
                               Class = c('Uncalled','Negative','Singlet','Doublet',
                                         'Uncalled','Negative','Singlet','Doublet'),
                               CountSum = c(tmp.counts.rna.uncalled,tmp.counts.rna.negative,tmp.counts.rna.singlet,tmp.counts.rna.doublet,
                                            tmp.counts.hto.uncalled,tmp.counts.hto.negative,tmp.counts.hto.singlet,tmp.counts.hto.doublet),
                               CountPer = c(tmp.percen.rna.uncalled,tmp.percen.rna.negative,tmp.percen.rna.singlet,tmp.percen.rna.doublet,
                                            tmp.percen.hto.uncalled,tmp.percen.hto.negative,tmp.percen.hto.singlet,tmp.percen.hto.doublet),
                               Labelpos = c(tmp.labpos.rna.uncalled,tmp.labpos.rna.negative,tmp.labpos.rna.singlet,tmp.labpos.rna.doublet,
                                            tmp.labpos.hto.uncalled,tmp.labpos.hto.negative,tmp.labpos.hto.singlet,tmp.labpos.hto.doublet),
                               Pool = c(pool.10x,pool.10x,pool.10x,pool.10x,
                                        pool.hto,pool.hto,pool.hto,pool.hto))
      
      df.hto <- rbind(df.hto,tmp.df.hto)
    }
  }
  
  df[['Labelpos']] <- ifelse(df[['ReadType']] == 'Unspliced',
                             df[['CountSum']]/2, df[['CountTot']] - df[['CountSum']]/2)
  
  
  p <- ggplot(data = df, aes(x = Pool, y = CountSum, fill = ReadType)) +
    geom_bar(stat = 'identity') + 
    geom_text(aes(label = paste0(CountPer,'%'),y = Labelpos),size = 3) +
    coord_flip() + ggtitle('Spliced/Unspliced - Counts') +
    xlab('Index') + ylab('Summed counts') + scale_fill_discrete(name = 'Read type') +
    facet_wrap(~CellType, nrow = 2) +
    theme_minimal(base_size = 20)
  
  ggsave(filename = paste0('read-distribution_alevin.png'),
         plot = p,
         path = dir.outs.qc.plots,
         width = 12,
         height = 8)
  
  if (var.wofl == '10x + HTO') {
    
    p.rna <- ggplot(data = df.hto[df.hto[['Library']] == 'RNA',],
                    aes(x = Pool, y = CountSum, fill=factor(Class, levels = c('Doublet','Singlet','Negative','Uncalled')))) +
      geom_bar(stat = 'identity') + 
      geom_text(aes(label = paste0(CountPer, '%'),
                    y = Labelpos),
                size = 3) +
      coord_flip() + ggtitle('Read distribution by HTO class', subtitle = 'RNA reads') +
      xlab('Index') + ylab('Summed counts') + scale_fill_discrete(name = 'Classification') +
      scale_fill_manual(values = c('#D34F73','#ff66b3','#ffb366','grey')) +
      labs(fill='Classification') + 
      theme_minimal(base_size = 20)
    
    
    p.hto <- ggplot(data = df.hto[df.hto[['Library']] == 'HTO',],
                    aes(x = Pool, y = CountSum, fill=factor(Class, levels = c('Doublet','Singlet','Negative','Uncalled')))) +
      geom_bar(stat = 'identity') + 
      geom_text(aes(label = paste0(CountPer, '%'),
                    y = Labelpos),
                size = 3) +
      coord_flip() + ggtitle('', subtitle = 'HTO reads') +
      xlab('') + ylab('Summed counts') + scale_fill_discrete(name = 'Classification') +
      scale_fill_manual(values = c('#D34F73','#ff66b3','#ffb366','grey')) +
      labs(fill='Classification') + 
      theme_minimal(base_size = 20)
    
    
    plot.comb <- cowplot::plot_grid(p.rna+NoLegend(),p.hto,ncol = 2,rel_widths = c(2,1.5))
    ggsave(filename = paste0('read-distribution_HTO-class.png'),
           plot = plot.comb,
           path = dir.outs.qc.plots,
           width = 16,
           height = length(df.hto$Pool)/2)
  }
  
  cat('#\n',
      '#\n',
      '###\n',
      rep('#',80),'\n',
      sep = '')
}

file.create(snakemake@output[[1]])
