library(Seurat)
library(tidyverse)
source('detect_DEGs.R')
# Load Seurat object; Expression data should have been probably normalized.
load('test_object.Rdata')
test_data <- FindNeighbors(test_data, dims = 1:50)
test_data <- FindClusters(test_data, resolution = 0.8)
DimPlot(test_data, group.by = c('labels','treatment','time'), split.by = 'geno') & NoLegend()
DimPlot(test_data, group.by = c('seurat_clusters'), split.by = 'geno') & NoLegend()

# Choose which cell population (clusters) to perform differential expression analysis.
selected_cells = which(test_data$labels %in% 'Agrp' & test_data$time == '3')
test_data.sub = subset(test_data, cells = selected_cells)
test_data.sub_exp = data.frame(t(as.matrix(test_data.sub@assays$SCT@data)), test_data.sub@meta.data)

# Example for running glm for a single gene:
run_glm(test_data.sub_exp$Stat3, test_data.sub@meta.data, model = '(exp) ~ ((treatment*geno))', mixed_model = F)
run_glm(test_data.sub_exp$Stat3, test_data.sub@meta.data, model = '(exp) ~ treatment*geno + (1|hash_id)', mixed_model = T)

agrp_deg.glmr = snRNA_DEGs(test_data.sub, fixed_effects = 'treatment*geno',random_effect = '(1|hash_id)', deg_type = 'abs', n_cores = 50)
# agrp_deg.glm = snRNA_DEGs(test_data.sub, fixed_effects = 'treatment*geno', mixed_model = F, deg_type = 'abs', n_cores = 50) # Typical GLM. Fast!

# Visualization by plotting the top candidate gene.
head(agrp_deg.glmr[order(agrp_deg.glmr$treatmentSal.genowt, decreasing = T),], n = 20)
# head(agrp_deg.glm[order(agrp_deg.glm$treatmentSal.genowt, decreasing = T),], n = 20)
ggplot(test_data.sub_exp,aes(y = Epb41l4a, x = treatment, fill = hash_id)) +
  geom_boxplot()+
  facet_wrap(~geno)+
  geom_jitter(position=position_dodge(0.75))+
  theme_bw()+
  theme(legend.position = 'none', 
        axis.text.x = element_text(size = 20), 
        strip.text.x = element_text(size = 20),
        axis.title = element_text(size = 20))
FeaturePlot(test_data, features = 'Epb41l4a', split.by = 'geno') & NoLegend()
