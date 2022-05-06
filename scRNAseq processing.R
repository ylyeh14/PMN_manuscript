#Instant the packages from Bioconductor first
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("SingleR")
#BiocManager::install("SingleCellExperiment")
#BiocManager::install("celldex")
#if (!requireNamespace("remotes", quietly = TRUE)) {
#  install.packages("remotes")
#}
#remotes::install_github("mojaveazure/seurat-disk")
#if (!requireNamespace("remotes", quietly = TRUE)) {
#  install.packages("remotes")
#}
#remotes::install_github("satijalab/seurat-wrappers")



#rm(list = ls())
gc()


#Load all the required libraries for analysis
library(Seurat)
library(ggplot2)
library(SingleR)
library(dplyr)
library(celldex)
library(RColorBrewer)
library(SingleCellExperiment)
library(patchwork)
library(harmony)
library(dplyr)
library(SeuratDisk)
library(SeuratWrappers)
library(cowplot)

# Load all datasets
Control_1 <- Read10X(data.dir = "~/Control_1")
Control_2 <- Read10X(data.dir = "~/Control_2")
Cancer_1 <- Read10X(data.dir = "~/Cancer_1")
Cancer_2 <- Read10X(data.dir = "~/Cancer_2")

# Initialize all Seurat objects with the raw (non-normalized data).
s_control_1 <- CreateSeuratObject(counts =Control_1, project = "control_1", min.cells = 3, min.features = 200)
s_control_2 <- CreateSeuratObject(counts =Control_2, project = "control_2", min.cells = 3, min.features = 200)
s_cancer_1 <- CreateSeuratObject(counts =Cancer_1, project = "cancer_1", min.cells = 3, min.features = 200)
s_cancer_2 <- CreateSeuratObject(counts =Cancer_2, project = "cancer_2", min.cells = 3, min.features = 200)

# Set all dgCMatrix to be NULL to save memory
Control_1 <- NULL
Control_2 <- NULL
Cancer_1 <- NULL
Cancer_2 <- NULL

#Assign the batch number to each Seurat object
s_control_1$batch <- 1
s_control_2$batch <- 2
s_cancer_1$batch <- 3
s_cancer_2$batch <- 3


#Calcultate mitochondria gene percentage
s_control_1[["percent.mt"]] <- PercentageFeatureSet(s_control_1, pattern = "mt-")
s_control_2[["percent.mt"]] <- PercentageFeatureSet(s_control_2, pattern = "mt-")
s_cancer_1[["percent.mt"]] <- PercentageFeatureSet(s_cancer_1, pattern = "mt-")
s_cancer_2[["percent.mt"]] <- PercentageFeatureSet(s_cancer_2, pattern = "mt-")


#Violin plots for individual sample groups 
VlnPlot(s_control_1, features = c("nFeature_RNA","nCount_RNA","percent.mt"),ncol = 3,pt.size = 0.1) &
  theme(plot.title = element_text(size=10))
VlnPlot(s_control_2, features = c("nFeature_RNA","nCount_RNA","percent.mt"),ncol = 3,pt.size = 0.1) &
  theme(plot.title = element_text(size=10))
VlnPlot(s_cancer_1, features = c("nFeature_RNA","nCount_RNA","percent.mt"),ncol = 3,pt.size = 0.1) &
  theme(plot.title = element_text(size=10))
VlnPlot(s_cancer_2, features = c("nFeature_RNA","nCount_RNA","percent.mt"),ncol = 3,pt.size = 0.1) &
  theme(plot.title = element_text(size=10))

#remove cells with too high or too low UMIs or high mitochondrial gene count
s_control_1 <- subset(s_control_1, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 15)
s_control_2 <- subset(s_control_2, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 15)
s_cancer_1 <- subset(s_cancer_1, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 15)
s_cancer_2 <- subset(s_cancer_2, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 15)


#Combine without adjustment; Normalizarion, find highly variable genes, find PCA
m_harmony    <- merge(s_control_1, y = c(s_control_2, s_cancer_1,s_cancer_2),  add.cell.ids = c("control_1", "control_2", "cancer_1","cancer_2"), project = "Merged_noNoM")
#head(colnames(m_harmony))
#tail(colnames(m_harmony))
#unique(sapply(X = strsplit(colnames(m_harmony), split = "_"), FUN = "[", 1))
#table(m_harmony$orig.ident)
#table(m_harmony$batch)

m_harmony <- NormalizeData(m_harmony , normalization.method = "LogNormalize", scale.factor = 10000)
m_harmony <- FindVariableFeatures(m_harmony, selection.method = "vst", nfeatures = 5000)
all.genes <- rownames(m_harmony)
m_harmony<- ScaleData(m_harmony, features = all.genes)
m_harmony<- RunPCA(m_harmony, features = VariableFeatures(object = m_harmony))

#m_harmony <- FindNeighbors(m_harmony, dims = 1:20)
#m_harmony <- FindClusters(m_harmony, resolution = 0.5, n.iter = 200)

m_harmony  <- RunUMAP(m_harmony , dims = 1:20)
DimPlot(m_harmony, reduction = "umap")


######Try Hormony to adjust batch effect and sample differences
m_harmony <- m_harmony %>% RunHarmony(group.by.vars= c("orig.ident", "batch"), plot_convergence = T)

harmony_embeddings <- Embeddings(m_harmony, 'harmony')
harmony_embeddings[1:5, 1:5]

p1 <- DimPlot(object = m_harmony, reduction = "harmony", pt.size = .1, group.by = "orig.ident") + NoLegend()
p2 <- VlnPlot(object = m_harmony, features = "harmony_1", group.by = "orig.ident", pt.size = .1) + NoLegend()
plot_grid(p1,p2)

m_harmony <- FindNeighbors(m_harmony,reduction = "harmony", dims = 1:20)
m_harmony <- FindClusters(m_harmony, resolution = 0.5, n.iter = 200)

m_harmony  <- RunUMAP(m_harmony,reduction = "harmony", dims = 1:20)
DimPlot(m_harmony, reduction = "umap")


########Compared the results after Harmony adjustment
m_harmony <- SetIdent(m_harmony,value = "orig.ident")
DimPlot(m_harmony,reduction = "umap") + plot_annotation(title = "Control vs early endometrial cancer, after integration (Harmony)")
DimPlot(m_harmony, reduction = "umap", group.by = "orig.ident", pt.size = .1, split.by = 'orig.ident') + NoLegend()

m_harmony <- SetIdent(m_harmony,value = "seurat_clusters")
DimPlot(m_harmony, reduction = "umap")



######Calculate the cell numbers in each cluster
library(data.table)
## extract meta data
md <- m_harmony@meta.data %>% as.data.table
head(md)

library(plyr)
summary(md)
ddply(md,~seurat_clusters,summarise, freq=length(seurat_clusters) )
a1 <- ddply(md,~orig.ident + seurat_clusters,summarise, freq=length(seurat_clusters) )
a2 <- ddply(md,~seurat_clusters +orig.ident ,summarise, freq=length(orig.ident) )


########MANUALLY ANNOTATE ALL CLUSTERS BY USING CANONICAL MARKERS##############
current_cluster_ids <- c(0, 1, 2, 3, 4, 5, 6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)

# List of new cluster IDs
new_cluster_ids <- c("PMN group A", "PMN group B", "Prap1+ luminal epithelium",
                     "PMN group C", "LTF+ luminal epithelium", "Adgre1+ Cd11b+ monocytes",
                     "Pax8+ secretory epithelium","T cells", "Dendritic cells",
                     "Mcam- endothelial cells","Fibroblast group A","Fibroblast group B",
                     "Adgre- Cd11b+ monocytes","Cbr2+ luminal epithelium","M-MDSC",
                     "Fibroblast group C ","Mcam+ endothelial cells", "B cells",
                     "Smooth muscle cells","Fibroblast Group D", "Adgre1+ Cd11b- macrophages")
names(new_cluster_ids) <- levels(x= m_harmony)
m_harmony <- RenameIdents(m_harmony, new_cluster_ids)

b <- DimPlot(m_harmony, reduction = "umap", label = "TRUE", repel = TRUE)+plot_annotation(title = "ctrl_vs_cancer")
plot(b)



#####Reorder the cluster
m_harmony@active.ident <- factor(m_harmony@active.ident, 
                            levels=rev(c("Prap1+ luminal epithelium", "LTF+ luminal epithelium","Pax8+ secretory epithelium","Cbr2+ luminal epithelium",
                                     "Fibroblast group A","Fibroblast group B","Fibroblast group C ", "Fibroblast Group D", "Smooth muscle cells",
                                     "Mcam+ endothelial cells", "Mcam- endothelial cells",
                                     "PMN group A", "PMN group B", "PMN group C", "M-MDSC",
                                     "Adgre1+ Cd11b- macrophages","Adgre1+ Cd11b+ monocytes","Adgre- Cd11b+ monocytes", "Dendritic cells", 
                                      "B cells","T cells"
                                     )))


####### Dot plots - the size of the dot corresponds to the percentage of cells expressing the
########## feature in each cluster. The color represents the average expression level
ft_combined <- c("Sprr2f","Epcam",  "Prap1", "Acta2", "Pecam1", "Mcam", "Cd34", "Pdgfra", "Mfap5","Ptprc",
                 "Cd2", "Cd3g", "Cd19", "Itgam","Cd14", "Cd84", "Cxcr2","Ly6g", "S100a9",
                 "Ly6c1","Icam1","Lgals3","Tnf", "Csf1r", "Adgre1", "Itgax","Flt3", "H2-Ab1"
                 )
DotPlot(m_harmony, features = ft_combined) #+ RotatedAxis()


#######Create the proportion plot for each group###########################################
##################Calculate the proportion from different clusters per sample group########
pt <- table(Idents(m_harmony), m_harmony$orig.ident)
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)

pt$Cluster_name <- new_cluster_ids

ggplot(pt, aes(x = Var2, y = Freq, fill = Cluster_name)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  xlab("Sample") +
  ylab("Proportion") +
  plot_annotation(title = "Proprotions of different clusters among four sample groups: 2 controls vs 2 cancer")

#####Manually annotate the cell type of different clusters#######
for(i in 1:length(pt$Var1)){
  if (pt$Var1[i] %in% c(0,1,3)) {pt$Var3[i]  <- "PMN"}
  else if (pt$Var1[i] %in% c(5,8,12,14,20)) {pt$Var3[i]  <- "Monocytic lineage"}
  else if (pt$Var1[i] %in% c(7,17)) {pt$Var3[i]  <- "Lymphocytes"}
  else if (pt$Var1[i] %in% c(2,4,6,13)) {pt$Var3[i]  <- "Epithelium"}
  else if (pt$Var1[i] %in% c(9,16)) {pt$Var3[i]  <- "Endothelium"}
  else pt$Var3[i] <- "Stromal cells"
}

pt$Var3 <- as.character(pt$Var3)
pt$Var2 <- as.character(pt$Var2)

library(plyr)
c2 <- ddply(pt,~Var2 + Var3 ,summarise, Freq =sum(Freq) )
myresults_pct <- c2 %>% dplyr::group_by(Var2) %>% dplyr::mutate(pct = round(Freq/sum(Freq),4)*100)
myresults_pct$Var3 <- as.factor(myresults_pct$Var3)
myresults_pct$pct
myresults_pct$pct <- c(3,69,5,7,0,16,4,49,11,19,1,16,4,7,2,9,76,2,2,14,1,7,74,2)
myresults_pct <- myresults_pct[-5,]

ggplot()+
  theme_bw(base_size = 15) +
  geom_bar(aes(Var2, pct,  fill=Var3), data=myresults_pct, stat="identity", position = position_stack(reverse=TRUE))  +
  geom_text(data=myresults_pct, aes( x=Var2, y= pct, label = paste0(pct,"%")),
            position = position_stack(vjust = 0.5))+
  xlab("Sample") +
  ylab("Proportion") +
  labs(fill="Cell types")+
  plot_annotation(title = "Proprotions of components among four sample groups: 2 controls vs 2 cancer")+
  scale_fill_discrete(breaks = c('Stromal cells', 'Monocytic lineage', 'PMN', 'Lymphocytes','Epithelium', 'Endothelium'))

####################################################################################################################
#######Figure S1B Create heatmaps for PMNs######
###Creating a subset of only PMNs####
a1<- subset(m_harmony, idents = c("PMN group A","PMN group B","PMN group C") )


####PMN-MDSC genes
ft_combined <- c("Ngp","Camp",  "Ltf", "Chil3", "Lcn2", "Ifitm6", "Lyz2", "Cybb", "Serpinb1a","Cd177",
                 "Anxa1", "Aldh2", "Ly6c2", "Mmp8","Adpgk", "Dstn", "Arhgdib","Ly6g", "S100a8",
                 "S100a8","AA467197","Tkt","Wfdc21", "Capg", "Cebpe", "S100a9",
                 "Syne1","Pglyrp1", "Ltb4r1","Lgals3","Tmsb10", "Olfm4","Hmgn2","Ceacam10","H2afz"
)

DoHeatmap(object = a1, features = ft_combined,size =3, disp.min = -2.5)+  scale_fill_gradient2( low = rev(c('#d1e5f0','#67a9cf','#2166ac')), mid = "white", high = rev(c('#b2182b','#ef8a62','#fddbc7')), midpoint = 0, guide = "colourbar", aesthetics = "fill")


###Activated PMN-MDSC genes
ft_combined <- c("Hba-a1",  "Hba-a2", "Ccl4", "Ccl3", "Cxcl3", "Jun", "Ccrl2", "Saa3","Spp1",
                 "Gadd45b", "Il1b", "Ninj1", "Clec4n","Hcar2", "Basp1", "Nfkbia","Btg1", "Il1rn",
                 "Ifrd1","Txnip","Ccl9","Ier3", "Ier5", "Rgs1", "Thbs1",
                 "Cxcl1","Hilpda", "Hist1h1c","Srgn","Hspa5", "Csf1","Ptgs2","Xbp1"
)

DoHeatmap(object = a1, features = ft_combined,size =3, disp.min = -2.5)+  scale_fill_gradient2( low = rev(c('#d1e5f0','#67a9cf','#2166ac')), mid = "white", high = rev(c('#b2182b','#ef8a62','#fddbc7')), midpoint = 0, guide = "colourbar", aesthetics = "fill")


#######################################################################################################################

#saveRDS(m_harmony, file = "~/ctrl_vs_cancer.rds")
m_harmony <- readRDS(file = "~/ctrl_vs_cancer.rds")

######################################################################################################################
#########################Trajectory inference of PMN goup ABC#########################################################
######################################################################################################################
#BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
#                       'limma', 'S4Vectors', 'SingleCellExperiment',
#                       'SummarizedExperiment', 'batchelor', 'Matrix.utils'))
#install.packages("devtools")
#devtools::install_github("cole-trapnell-lab/leidenbase")
#devtools::install_github('cole-trapnell-lab/monocle3')

library(monocle3)
library(Seurat)
library(ggplot2)
library(SingleR)
library(dplyr)
library(celldex)
library(RColorBrewer)
library(SingleCellExperiment)
library(patchwork)
library(harmony)
library(dplyr)
library(SeuratDisk)
library(SeuratWrappers)
library(cowplot)

a1<- subset(m_harmony, idents = c("PMN group A","PMN group B","PMN group C") )
a1  <- RunUMAP(a1,reduction = "harmony", dims = 1:15)
DimPlot(a1)

cds <- as.cell_data_set(a1)
cds <- cluster_cells(cds, resolution=5*1e-4,num_iter = 100)


cds <- learn_graph(cds, use_partition = FALSE)
plot_cells(cds, label_groups_by_cluster = TRUE, label_leaves = FALSE, label_branch_points = FALSE)

#> sessionInfo()
#R version 4.1.1 (2021-08-10)
#Platform: x86_64-apple-darwin17.0 (64-bit)
#Running under: macOS Big Sur 11.3.1

#Matrix products: default
#LAPACK: /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRlapack.dylib

#locale:
#[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

#attached base packages:
#[1] stats4    stats     graphics  grDevices utils     datasets  methods   base

#other attached packages:
# [1] plyr_1.8.6                  data.table_1.14.2           cowplot_1.1.1
# [4] SeuratWrappers_0.3.0        SeuratDisk_0.0.0.9019       harmony_0.1.0
# [7] Rcpp_1.0.8                  patchwork_1.1.1             RColorBrewer_1.1-2
#[10] celldex_1.2.0               dplyr_1.0.7                 SingleR_1.6.1
#[13] ggplot2_3.3.5               SeuratObject_4.0.4          Seurat_4.0.6
#[16] monocle3_1.0.0              SingleCellExperiment_1.14.1 SummarizedExperiment_1.24.0
#[19] GenomicRanges_1.44.0        GenomeInfoDb_1.30.0         IRanges_2.28.0
#[22] S4Vectors_0.32.3            MatrixGenerics_1.6.0        matrixStats_0.61.0
#[25] Biobase_2.54.0              BiocGenerics_0.40.0

#loaded via a namespace (and not attached):
#  [1] AnnotationHub_3.0.2           BiocFileCache_2.0.0           igraph_1.2.11
#  [4] lazyeval_0.2.2                splines_4.1.1                 BiocParallel_1.28.2
#  [7] listenv_0.8.0                 scattermore_0.7               digest_0.6.29
# [10] htmltools_0.5.2               viridis_0.6.2                 fansi_1.0.2
# [13] magrittr_2.0.1                memoise_2.0.1                 ScaledMatrix_1.0.0
# [16] tensor_1.5                    cluster_2.1.2                 ROCR_1.0-11
# [19] remotes_2.4.2                 Biostrings_2.62.0             globals_0.14.0
# [22] R.utils_2.11.0                spatstat.sparse_2.1-0         colorspace_2.0-2
# [25] rappdirs_0.3.3                blob_1.2.2                    ggrepel_0.9.1
# [28] crayon_1.4.2                  RCurl_1.98-1.5                jsonlite_1.7.2
# [31] spatstat.data_2.1-2           survival_3.2-13               zoo_1.8-9
# [34] glue_1.6.0                    polyclip_1.10-0               gtable_0.3.0
# [37] zlibbioc_1.40.0               XVector_0.34.0                leiden_0.3.9
# [40] DelayedArray_0.20.0           BiocSingular_1.8.1            future.apply_1.8.1
# [43] abind_1.4-5                   scales_1.1.1                  pheatmap_1.0.12
# [46] DBI_1.1.2                     miniUI_0.1.1.1                viridisLite_0.4.0
# [49] xtable_1.8-4                  reticulate_1.22               spatstat.core_2.3-2
# [52] proxy_0.4-26                  bit_4.0.4                     rsvd_1.0.5
# [55] htmlwidgets_1.5.4             httr_1.4.2                    ellipsis_0.3.2
# [58] ica_1.0-2                     farver_2.1.0                  R.methodsS3_1.8.1
# [61] pkgconfig_2.0.3               dbplyr_2.1.1                  uwot_0.1.11
# [64] deldir_1.0-6                  utf8_1.2.2                    labeling_0.4.2
# [67] tidyselect_1.1.1              rlang_0.4.12                  reshape2_1.4.4
# [70] later_1.3.0                   AnnotationDbi_1.56.2          BiocVersion_3.13.1
# [73] munsell_0.5.0                 tools_4.1.1                   cachem_1.0.6
# [76] cli_3.1.1                     ExperimentHub_2.0.0           generics_0.1.1
# [79] RSQLite_2.2.9                 ggridges_0.5.3                stringr_1.4.0
# [82] fastmap_1.1.0                 yaml_2.2.1                    goftest_1.2-3
# [85] bit64_4.0.5                   fitdistrplus_1.1-6            purrr_0.3.4
# [88] RANN_2.6.1                    KEGGREST_1.34.0               pbapply_1.5-0
# [91] future_1.23.0                 nlme_3.1-153                  sparseMatrixStats_1.4.2
# [94] mime_0.12                     R.oo_1.24.0                   hdf5r_1.3.5
# [97] compiler_4.1.1                rstudioapi_0.13               interactiveDisplayBase_1.30.0
#[100] filelock_1.0.2                curl_4.3.2                    plotly_4.10.0
#[103] png_0.1-7                     spatstat.utils_2.3-0          tibble_3.1.6
#[106] stringi_1.7.6                 RSpectra_0.16-0               lattice_0.20-45
#[109] Matrix_1.4-0                  vctrs_0.3.8                   pillar_1.6.4
#[112] lifecycle_1.0.1               BiocManager_1.30.16           spatstat.geom_2.3-1
#[115] lmtest_0.9-39                 RcppAnnoy_0.0.19              BiocNeighbors_1.10.0
#[118] bitops_1.0-7                  irlba_2.3.5                   httpuv_1.6.4
#[121] R6_2.5.1                      promises_1.2.0.1              KernSmooth_2.23-20
#[124] gridExtra_2.3                 parallelly_1.30.0             codetools_0.2-18
#[127] MASS_7.3-54                   assertthat_0.2.1              leidenbase_0.1.3
#[130] withr_2.4.3                   sctransform_0.3.2             GenomeInfoDbData_1.2.7
#[133] mgcv_1.8-38                   parallel_4.1.1                grid_4.1.1
#[136] rpart_4.1-15                  beachmat_2.8.1                tidyr_1.1.4
#[139] DelayedMatrixStats_1.14.3     Rtsne_0.15                    shiny_1.7.1z
