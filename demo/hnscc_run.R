rm(list = ls())
gc()
library(Seurat)
library(dplyr)
library(CellChat)
library(MASS)
library(tibble)
library(ggplot2)
# library(hdf5r)
library(future)
options(stringsAsFactors = FALSE)

source("./Rcodes/Communications.R")
source("./Rcodes/Utilities.R")
source("./Rcodes/Crosstalk_analysis.R")


# Load data
hnscc_expression = readRDS("./datasets/hnscc_expression.rds")
expression = hnscc_expression$expression
sample_info = hnscc_expression$sample_info
tumors_remove = c("HN10","HN","HN12", "HN13", "HN24", "HN7", "HN8","HN23")
# tumors_remove = c("HN")
sample_info <- sample_info %>% filter(`Lymph node` == 0 & !(tumor %in% tumors_remove))
expression <- expression[sample_info$cell,] %>% t()
cell_anno <- data.frame(cell = sample_info$cell, cluster = sample_info$`non-cancer cell type`)
cell_anno$cluster[cell_anno$cluster == '0'] <- "malignant"
cell_anno$cluster[cell_anno$cluster == '-Fibroblast'] <- "Fibroblast"
cell_anno$cluster[cell_anno$cluster == 'B cell'] <- "B_cell"
cell_anno$cluster[cell_anno$cluster == 'T cell'] <- "T_cell"

SeuratObj <- CreateSeuratObject(expression)
rm(expression,hnscc_expression,sample_info)
Idents(SeuratObj) <- cell_anno$cluster
gc()

SeuratObj <- subset(SeuratObj, nCount_RNA < 25000 & nFeature_RNA<10000)
SeuratObj <- SeuratObj %>% NormalizeData() %>% FindVariableFeatures() %>% ScaleData()
SeuratObj <- SCTransform(SeuratObj,vst.flavor = "v1")

cell_anno <- data.frame(cell = names(Idents(SeuratObj)), cluster = Idents(SeuratObj) %>% as.character())

# load("/home/jiawen/myMLnet/pathways/LR_layer1_human.rda")
LRDB <- readRDS("./datasets/LR_layer1.rds") %>% distinct(from, to, .keep_all = T)
RecTFDB <- readRDS("./datasets/RT_layer2.rds") %>% distinct(from, to, .keep_all = T)
TFTGDB <- readRDS("./datasets/TT_layer3.rds") %>% distinct(from, to, .keep_all = T)

allgenes <- rownames(SeuratObj@assays$RNA$data)
LRDB <- Filter_DB(LRDB,allgenes)
RecTFDB <- Filter_DB(RecTFDB,allgenes)
TFTGDB <- Filter_DB(TFTGDB, allgenes)

write.table(RecTFDB[,1:2], file = "./inputs/RecTFDB.txt",quote = F,sep = " ")
write.table(TFTGDB[,1:2], file = "./inputs/TFTGDB.txt",quote = F,sep = " ")

cellchat_use <- T
LigRec_original <- Infer_CCI(SeuratObj, LRDB = LRDB, cellchat_output = cellchat_use, db_use = "human")

all_types <- aggregate(cell_anno$cell,list(cell_anno$cluster),length)

types_used <- (all_types %>% filter(x > 200))$Group.1

Seurat_info <- list("Object" = SeuratObj,
                    "Anno" = cell_anno,
                    "Cellchat" = cellchat_use,
                    "CCC" = LigRec_original,
                    "Celltypes" = types_used)
filename <- "./results/seuratinfo.rds"
saveRDS(Seurat_info, file = filename)
rm(Seurat_info)

#Seurat_info <- readRDS(filename)
#SeuratObj <- Seurat_info$Object
#cell_anno <- Seurat_info$Anno
#cellchat_use <- Seurat_info$Cellchat
#LigRec_original <- Seurat_info$CCC
#types_used <- Seurat_info$Celltypes

High_TGs <- list()
for(i in 1:length(types_used)){
  target_type <- types_used[i]
  Exp_clu <- Get_Exp_Clu(SeuratObj, clusterID = target_type, assay = "SCT", datatype = "counts", cutoff = 0.15)
  allcells <- colnames(Exp_clu)
  TG_used <- SeuratObj@assays$RNA$counts[rowSums(SeuratObj@assays$RNA$counts[,allcells] > 0) > 0.7*length(allcells),] %>% rownames()
  TG_used <- intersect(TG_used,TFTGDB$To)
  High_TGs[[i]] <- TG_used
}
common_TGs <- Reduce(intersect, High_TGs)
write.table(common_TGs, file = "./inputs/TG.csv",quote = F)

conda_python <- "/home/jiawen/anaconda3/envs/SigXTalk/bin/python"

target_type <- types_used[2]
message(paste0("Analysing the cell type ",target_type,"...\n"))

LR_Pairprob <- Extract_LR_Prob(LigRec_original, target_type = target_type, cellchat_use = cellchat_use)
Rec_act <- aggregate(LR_Pairprob$Weight, list(LR_Pairprob$To), sum)
colnames(Rec_act) <- c("Rec", "Weight")
Rec_act <- Rec_act[Rec_act$Weight > 0.1*max(Rec_act$Weight),]
LR_Pairprob <- LR_Pairprob[LR_Pairprob$To %in% Rec_act$Rec,]

write.table(LR_Pairprob, file = "./inputs/LigRec.csv",quote = F, sep = " ")
Exp_clu <- Get_Exp_Clu(SeuratObj, clusterID = target_type, assay = "SCT", datatype = "counts")
write.table(Exp_clu, file = './inputs/ExpressionCount.csv',quote = F, sep = " ")

system2(conda_python, args = c("./pythoncodes/main_new.py"))
gc()



