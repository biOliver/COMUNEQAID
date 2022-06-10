# A script to extract marker genes in a hierarchy manner.
# We will first extract marker genes at the higher level (e.g., within hypothalamus, markers for neurons).
# Then within each category, extract marker genes at the lower level (e.g., neuron cell types).
# Extracting markers in a hierarchy manner should provide better markers. [relative markers]
library(Matrix)
library(Seurat)
library(tidyverse)
source('sctype_scoring.R')
source('find_markers.R')

# Load ref data set, should have been normalized by sequencing depth.
hypo.ref <- readRDS('ref.cam.hypo.rds') 
table(hypo.ref$taxonomy_lvl2, hypo.ref$cell_type_all_lvl2)

# Step 1. Extract markers at the higher level.
n_cell_types = length(table(hypo.ref$taxonomy_lvl2))
marker_db_lvl1 = build_marker_db(hypo.ref@assays$SCT@data,
                                        hypo.ref$taxonomy_lvl2, 
                                        sample_size = 500, 
                                        # Sampling N cells from each cell type, making sure that the differences are not dictated by the dominant cell type(s).
                                        test_fun = wilcox.test, 
                                        # Using Wilcox.text provides better performance for extracting marker genes, because gene expression may not follow normal distribution.
                                        n_core = n_cell_types)

# Step 2. Extract markers at lower level.
marker_db_lvl2 = list()

for (i in names(table(hypo.ref$taxonomy_lvl2))) {
  target_cells = which(hypo.ref$taxonomy_lvl2 == i)
  n_cell_types = length(table(hypo.ref$cell_type_all_lvl2[target_cells]))
  if (number_cell_types > 1) {marker_db_lvl2[[i]] = build_marker_db(hypo.ref@assays$SCT@data[,target_cells],
                                                                    hypo.ref$cell_type_all_lvl2[target_cells], 
                                                                    sample_size = 100, 
                                                                    test_fun = wilcox.test,
                                                                    n_core = n_cell_types)}
}

# Combine lvl 1 and lvl2 marker db into a single list.
marker_db_cam_hypo = list(lv1 = marker_db_lvl1, lv2 = marker_db_lvl2)
save(marker_db_cam_hypo, file = 'marker_db_cam_hypo.Rdata')
