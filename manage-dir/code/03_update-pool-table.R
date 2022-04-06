####################  Start sink  ####################
sink.file <- file(file.path('logs/03_update-pool-table.log'), 'w')
sink(file = sink.file, append = T,
     type = 'output', split = T)
cat(as.character(lubridate::now()),'\n', sep = '')

################################################################################
##########                          Packages                          ##########
################################################################################
cat(rep('#',80),'\n',
    '#####                                                                      #####\n',
    '#####                     Updating pool table & fastqQC                    #####\n',
    '#####                                                    **** / ******     #####\n',
    rep('#',80),'\n',
    '#####\n',
    '##\n',
    '##\n',
    '##\tPROJECT ID:     ',snakemake@params[['scopID']],'\n',
    '##\n',
    '##\n',
    '###\n',
    rep('#',80),'\n',
    sep = '')

cat('#\tLoading packages..\n',
    '#\n',
    '#\n',
    sep = '')

suppressPackageStartupMessages({
  library(future)
  library(tidyverse)
  library(Biobase)
})
  
################################################################################
##########                          Functions                         ##########
################################################################################
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

################################################################################
##########                            Init                            ##########
################################################################################
cat('#\tLoading info..\n',
    '#\n',
    '#\n',
    sep = '')

scop.ID     <- snakemake@params[['scopID']]
com.ID      <- snakemake@params[['comID']]
com.ID.list <- str_split(com.ID, ',', simplify = T)[1,]
user.ID     <- snakemake@params[['userID']]

dir.proj <- paste0('/', file.path('projects',user.ID,'COMUNEQAID','outs',scop.ID))
dir.outs.qc <- file.path(dir.proj, 'scRNAseq', '00_QC')
dir.bcls <- file.path(dir.proj, 'scRNAseq', '01_BCL')

indices.single <- read_csv('/projects/SCOP/resources/non-species-specific/indices/chromium-single-indices.csv',
                           col_names = c('Name','X1','X2','X3','X4')) %>% 
  gather(key = Var,
         value = Barcode,
         X1,X2,X3,X4)
indices.dual <- read_csv('/projects/SCOP/resources/non-species-specific/indices/chromium-dual-indices.csv',
                         col_names = c('Name','X1','X2')) %>% 
  mutate(Barcode = paste0(X1,'+',X2), .keep = 'unused')

cat('#\tLooping over runs..\n',
    '#\n',
    '#\n',
    sep = '')

for (com.ID in com.ID.list) {
  cat(rep('#',80),'\n',
      '#####\n',
      '##\n',
      '##\n',
      '##\tCOMUNEQA ID:     ',com.ID,'\n',
      '##\n',
      '##\n',
      '###\n',
      rep('#',80),'\n',
      '###\n',
      '#\n',
      '#\n',
      sep = '')
  
  unique.ID <- paste(scop.ID, com.ID,sep = '_')
  
  dir.data <- paste0('/', file.path('projects',user.ID,'COMUNEQAID','manage-dir','tmp-data',unique.ID))
  
  dir.outs.log <- file.path(dir.proj, 'scRNAseq', '04_Log', com.ID)
  dir.outs.qc.plots <- file.path(dir.outs.qc, com.ID, 'plots')
  dir.create(dir.outs.log, recursive = T, showWarnings = F)
  dir.create(dir.outs.qc.plots, recursive = T, showWarnings = F)
  
  var.wofl <- scan(file.path(dir.data,'tmp_workflow.txt'), what = '', sep = '\n', quiet = T)
  var.vers <- scan(file.path(dir.data,'tmp_versionRNA.txt'), what = '', sep = '\n', quiet = T)
  
  ################################################################################
  ##########                    Update Pool Table                       ##########
  ################################################################################
  
  cat('#\tReading pool table..\n',
      '#\n',
      '#\n',
      sep = '')
  
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
  cat('#\tUpdating pool table..\n',
      '#\n',
      '#\n',
      sep = '')
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
      
      if (var.vers == 'RNA v3.0') {
        known_barcodes.10x <- filter(known_barcodes, str_detect(SampleId, 'SI-GA-'))
        known_barcodes.10x[['SampleId']] <- substr(known_barcodes.10x[['SampleId']],1,nchar(known_barcodes.10x[['SampleId']]) - 2)

        known_barcodes.hto <- filter(known_barcodes, str_detect(SampleId, 'D7'))
        
        known_barcodes <- bind_rows(known_barcodes.10x,known_barcodes.hto)
      }
      
      known_barcodes %>%
        group_by(SampleId) %>%
        summarise(Counts = sum(NumberReads)) -> read.df
      
      for (pool.10x.i in seq(pool.10x)) {
        
        if (var.vers == 'RNA v3.0') {
          tmp.index <- paste0('SI-GA-', pool.10x[pool.10x.i])
          tmp.pool[['Index (10x)']] <- paste(paste0('SI-GA-', pool.10x), collapse = ',')
        }
        if (var.vers == 'RNA v3.1') {
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
  
  write.table(x = pool.table,
              file = file.path(dir.data,'poolTable-updated.csv'))
  
  
  nPools        <- dim(pool.table)[1]
  

  ################################################################################
  ##########                          QC plots                          ##########
  ################################################################################
  cat('#\tMaking QC plots..\n',
      '#\n',
      '#\n',
      sep = '')
  
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
    
    ggsave(filename = paste0('undetermined-BCs-top10_',bcl,'.png'),
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
    
    if (var.vers == 'RNA v3.0') {
      known_barcodes.10x <- filter(known_barcodes, str_detect(SampleId, 'SI-GA-'))
      known_barcodes.10x[['SampleId']] <- substr(known_barcodes.10x[['SampleId']],1,nchar(known_barcodes.10x[['SampleId']]) - 2)
      
      known_barcodes.hto <- filter(known_barcodes, str_detect(SampleId, 'D7'))
      
      known_barcodes <- bind_rows(known_barcodes.10x,known_barcodes.hto)
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
  
  ggsave(filename = paste0('read-distribution_bcl2fastq.png'),
         plot = p,
         path = dir.outs.qc.plots,
         width = 12,
         height = 12)
}
  
sink()

file.create(snakemake@output[[1]])
