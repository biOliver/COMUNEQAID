# Annotation with hierarchy marker database.
source('find_markers.R')
source('annotate_cells.R')
library(future)
library(styler)
library(Matrix)
library(Seurat)
library(tidyverse)

# Load query data.
load('test.data.B10_ref.based.Rdata') 
# Load reference RNAseq data
hypo.ref <- readRDS('ref.cam.hypo.rds') 
# Load markers, extracted by extracter_marker_hierarcy.R
load('marker_db_cam_hypo.Rdata')

# Level 1 annotation ####

# Construct marker gene list:
# Users can also provide their own marker gene list: list(cell_type: markers)
# Select the top n DEGs (target cells vs. background) from each cell type as marker genes. 
# For high level annotation, n within 100 and 200, for low level annotation, n = 30 - 50 should be fine.
top_n = 200
marker_gene_lv1 = contruct_marker_list(marker_db_cam_hypo$lv1, top_n = top_n)

# Annotation:
test.data.B10.lv1 =  cell_annotation_ref_based(test.data.B10, marker_gene_lv1, target_cells = colnames(test.data.B10), label = 'sctype_label.lv1', filter = T)

# Visualization:
DimPlot(test.data.B10.lv1, group.by = c('sctype_label.lv1.filtered','sctype_label.lv1.filtered.cluster'),label = T) & NoLegend()

# Level 2 annotation ####
# Neuron cell types as an example.
# Construct marker gene list:
top_n = 50
marker_gene_lv2.neuron = contruct_marker_list(marker_db_cam_hypo$lv2$Neuron, top_n = top_n)

# Annotation:
n_pc = 50
neurone_cells = colnames(test.data.B10.lv1)[which(test.data.B10.lv1$sctype_label.lv1.filtered.cluster %in% c('Neuron'))]
test.data.B10.neuron =  cell_annotation_ref_based(test.data.B10.lv1, marker_gene_lv2.neuron, n_pc = n_pc, 
                                        target_cells = neurone_cells,
                                        label = 'sctype_label.lv2')
# Visualization: 
test.data.B10.neuron = RunPCA(test.data.B10.neuron)
test.data.B10.neuron = RunUMAP(test.data.B10.neuron, dims = c(1:n_pc))
DimPlot(test.data.B10.neuron, group.by = c('sctype_label.lv2.filtered.cluster', 'seurat_clusters'),label = T) & NoLegend()

# Comparison with Seurat annotation (anchor transfer) ####
hypo.ref.neuron = subset(hypo.ref, subset = `taxonomy_lvl2` == 'Neuron')
hypo.ref.neuron = RunPCA(hypo.ref.neuron)
hypo.ref.neuron = RunUMAP(hypo.ref.neuron, dims = c(1:n_pc))
# Seurat cell level annotation:
hypo.ref.neuron.anchors <- FindTransferAnchors(reference = hypo.ref.neuron, query = test.data.B10.neuron, recompute.residuals = F,
                                               dims = 1:n_pc, reference.reduction = "pca")
predictions.neuron <- TransferData(anchorset = hypo.ref.neuron.anchors, refdata = hypo.ref.neuron$cell_type_all_lvl2, dims = 1:n_pc)
test.data.B10.neuron <- AddMetaData(test.data.B10.neuron, metadata = predictions.neuron)
# Seurat cluster level annotation.
test.data.B10.neuron = cluster_level_annotation(test.data.B10.neuron, sctype_label = 'predicted.id')

# Visualization:
DimPlot(test.data.B10.neuron, group.by = c('sctype_label.lv2.filtered.cluster','sctype_label.lv2.filtered','seurat_clusters','predicted.id.cluster'),label = T) & NoLegend()
# Plotting candidate genes:
FeaturePlot(test.data.B10.neuron, features = c('Pomc','Anxa2','Ttr','Qrfp'))
