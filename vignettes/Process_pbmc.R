library(Seurat)
library(pbmc3k.SeuratData)
library(dplyr)

counts <- pbmc3k@assays$RNA$counts
pbmc <- CreateSeuratObject(counts = counts, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes, vars.to.regress = "percent.mt")

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunUMAP(pbmc, dims = 1:10)

new.cluster.ids <- c("Naive_CD4_T", "CD14_Mono", "Memory_CD4_T", "B", "CD8_T", "FCGR3A_Mono",
                     "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
pbmc$Clusters <- Idents(pbmc) %>% as.character()
DimPlot(pbmc, reduction = "umap", label = T, pt.size = 1) 

saveRDS(pbmc, file = "./Desktop/SigXTalk_demo_data.rds")
