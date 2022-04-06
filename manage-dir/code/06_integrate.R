################################################################################
##########                         Integrate                          ##########
################################################################################

seur.list <- list()

cat('#\tReading the following files: \n',
    sep = '')
for (com.ID_i in seq(com.ID.list)) {
  
  com.ID <- com.ID.list[[com.ID_i]]
  
  dir.outs.comb <- file.path(dir.proj, 'Output', 'data', 'comb', com.ID)
  
  comb.files <- grep(list.files(file.path(dir.outs.comb)), pattern = 'rna_seurat.rds', value = T)
  
  for (fname in comb.files) {
    cat('#\t-\t',fname,'\n',
        sep = '')
  }
  seur.list[com.ID_i] <- lapply(file.path(dir.outs.comb, comb.files), readRDS)
}

# Merge objects
if (length(seur.list) > 1) {
  seur.comb <- merge(seur.list[[1]],
                     seur.list[2:length(seur.list)])
}else {
  seur.comb <- seur.list[[1]]
}

seur.comb.list <- SplitObject(seur.comb, split.by = 'orig.ident')
seur.comb.list <- lapply(X = seur.comb.list, FUN = SCTransform)
features <- SelectIntegrationFeatures(object.list = seur.comb.list, nfeatures = 3000)
seur.comb.list <- PrepSCTIntegration(object.list = seur.comb.list, anchor.features = features)

seur.comb.anchors <- FindIntegrationAnchors(object.list = seur.comb.list, normalization.method = "SCT",
                                            anchor.features = features)
seur.comb.sct <- IntegrateData(anchorset = seur.comb.anchors, normalization.method = "SCT")
seur.comb.sct <- RunPCA(seur.comb.sct, verbose = FALSE)
dims <- round(
  as.numeric(
    intrinsicDimension::maxLikGlobalDimEst(
      data = seur.comb.sct@reductions[['pca']][, 1:50],
      k = 20)))
seur.comb.sct <- RunUMAP(seur.comb.sct, reduction = "pca", dims = 1:dims)


atlas.ID <- paste0(com.ID.list, collapse = '_')

dir.outs.qc.plots <- file.path(dir.outs.qc, atlas.ID, 'plots')
dir.outs.qc.summa <- file.path(dir.outs.qc, atlas.ID, 'summary')
dir.outs.atlas <- file.path(dir.proj, 'Output', 'data', 'integrate', atlas.ID)

dir.create(dir.outs.qc.plots, recursive = T, showWarnings = F)
dir.create(dir.outs.qc.summa, recursive = T, showWarnings = F)
dir.create(dir.outs.atlas, recursive = T, showWarnings = F)

# Output objects
out.name.atl <- paste0(atlas.ID,'_atl_seurat.rds')

out.path.atl <- file.path(dir.outs.atlas,out.name.atl)

cat('#\t-\tsaving atlas\n',
    sep = '')
saveRDS(seur.comb.sct, out.path.atl)

cat('#\n',
    '#\n',
    sep = '')
################################################################################
################################################################################
################################################################################



