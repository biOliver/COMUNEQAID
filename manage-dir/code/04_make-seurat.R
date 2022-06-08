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
##########                            Init                            ##########
################################################################################
scop.ID     <- snakemake@params[['scopID']]
com.ID      <- snakemake@params[['comID']]
com.ID.list <- str_split(com.ID, ',', simplify = T)[1,]
user.ID     <- snakemake@params[['userID']]

source(file.path('/projects', user.ID, 'COMUNEQAID/manage-dir/code/shared_functions.R'))

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
    
    ############# MAKE UMI PLOT  (MAYBE NEW VERSION) ############
    #if (var.wofl == '10x + HTO') {
    #  featDump.10x <- suppressMessages(read_delim(file.path(mat.files.10x,'featureDump.txt'), delim = '\t'))
    #  featDump.hto <- suppressMessages(read_delim(file.path(mat.files.hto,'featureDump.txt'), delim = '\t'))
    #  
    #  p <- makeUmiPlot(feat.RNA = featDump.10x,
    #                   feat.HTO = featDump.hto,
    #                   CBs.called = cells.keep,
    #                   title = paste0('UMI-UMI plot - ', tmp.pool[['Index (10x)']]))
    #  
    #  ggsave(filename = paste0(tmp.pool[['Index (10x)']],'_UMI-UMI-plot.png'),
    #         plot = p,
    #         path = dir.outs.qc.plots,
    #         width = 8,
    #         height = 8)
    #}
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
      seur.full <- NormalizeData(seur.full, assay = 'HTO', normalization.method = 'CLR', margin = 2, verbose = F)
      
      cat('#\t-\tdemultiplexing HTOs..\n',
          sep = '')
      # Setting the threshold based on the quantiles of the negative and the positive clusters.
      q_l = 1
      q_h = 0.001
      seur.full <- HTODemux.mcl(seur.full, q_l =  q_l, q_h = q_h)
      seur.full[,seur.full$HTO_mcl_classification.global == 'Doublet']
      
      table(seur.full$HTO_mcl_classification.global)[names(table(seur.full$HTO_mcl_classification.global)) == 'Doublet']
      
      cat('#\t-\t\tDoublets:\t\t',table(seur.full$HTO_mcl_classification.global)[names(table(seur.full$HTO_mcl_classification.global)) == 'Doublet'],'\n',
          '#\t-\t\tNegatives:\t\t',table(seur.full$HTO_mcl_classification.global)[names(table(seur.full$HTO_mcl_classification.global)) == 'Negative'],'\n',
          '#\t-\t\tSinglets:\t\t',table(seur.full$HTO_mcl_classification.global)[names(table(seur.full$HTO_mcl_classification.global)) == 'Singlet'],'\n',
          sep = '')
      
      HTODemux_mcl.visualization(seur.full, q_l =  q_l, q_h = q_h)
      
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

      p <- DimPlot(seur.full, group.by = 'hash.mcl.ID') + 
        ggtitle(paste0('HTO demultiplexing'))
      
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
      
      ############ MAKE UMI PLOT  (MAYBE NEW VERSION) ############
      featDump.10x <- suppressMessages(read_delim(file.path(mat.files.10x,'featureDump.txt'), delim = '\t'))
      featDump.hto <- suppressMessages(read_delim(file.path(mat.files.hto,'featureDump.txt'), delim = '\t'))
      
      Idents(seur.full) <- 'HTO_mcl_classification.global'
      seur.full.neg <- subset(seur.full, idents = 'Negative')
      seur.full.dou <- subset(seur.full, idents = 'Doublet')
      seur.full.sub <- subset(seur.full, idents = 'Singlet')
      
      common.CBs <- intersect(featDump.10x$CB,featDump.hto$CB)
      
      featDump.10x %>% 
        filter(CB %in% common.CBs) %>% 
        arrange(CB) -> featDump.10x
      
      featDump.hto %>% 
        filter(CB %in% common.CBs) %>% 
        arrange(CB) -> featDump.hto
      
      common.CBs.sort <- sort(common.CBs)
      
      tmpdf <- ifelse(common.CBs.sort %in% colnames(seur.full.neg), 'Negative',
                      ifelse(common.CBs.sort %in% colnames(seur.full.dou), 'Doublet',
                             ifelse(common.CBs.sort %in% colnames(seur.full.sub), 'Singlet', 'Uncalled')))

      p <- ggplot(mapping = aes(x = featDump.10x$DeduplicatedReads,
                                y = featDump.hto$DeduplicatedReads,
                                col = tmpdf)) +
        geom_point() + scale_x_log10() + scale_y_log10() +
        labs(x = 'UMIs (RNA)',
             y = 'UMIs (HTO)',
             col = '') +
        theme_minimal() +
        ggtitle(paste0('UMI-UMI plot - ', tmp.pool[['Index (10x)']]))
      
      ggsave(filename = paste0(tmp.pool[['Index (10x)']],'_UMI-UMI-plot.png'),
             plot = p,
             path = dir.outs.qc.plots,
             width = 8,
             height = 8)
      
      ############################################################
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