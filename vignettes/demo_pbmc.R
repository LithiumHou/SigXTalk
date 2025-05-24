library(SigXTalkR)
library(Seurat)
library(dplyr)
library(CellChat)

work_dir <- "./vignettes"
setwd("work_dir")

SeuratObj <- readRDS("./SigXTalk_demo_data.rds") # as the seurat object
cell_anno <- data.frame(cell = names(Idents(SeuratObj)), cluster = Idents(SeuratObj) %>% as.character()) # The metadata of the dataset

DimPlot(SeuratObj, reduction = "umap", pt.size = 1)

data("RTF_human")
data("TFT_human")
RecTFDB <- RTF_human %>% distinct(from, to, .keep_all = T) %>% Filter_DB(rownames(SeuratObj@assays$RNA$data))
TFTGDB <- TFT_human %>% distinct(from, to, .keep_all = T) %>% Filter_DB(rownames(SeuratObj@assays$RNA$data))

LR_original <- Infer_CCI(SeuratObj, cell_anno, use_spatial = F, db_use = "human")

target_type <- "CD14_Mono"
TG_used <- FindMarkers(SeuratObj, target_type, min.pct = 0.25, only.pos = T, logfc.threshold = 0.25)
TG_used <- filter(TG_used, p_val_adj<1e-3) %>% rownames()

input_dir <- "./inputs"
Prepare_Input(SeuratObj, target_type, TGs = TG_used, CCC_results = LR_original, RecTFDB, TFTGDB, data_dir = input_dir,
              assay = "RNA", datatype = "scale.data", exp_threshold = 0.05, CCC_threshold = 0.05)

args.project <- "PBMC" # The name of the project, e.g., the dataset, or any other name you like
python_script <- "./main.py" # The realtive path of the python script
args <- c("--project", args.project, "--target_type",target_type)

conda_env <- "SigXTalk_py" # The conda environment that is previously installed to train the hypergraph neural network
Run_py_script(python_script, conda_env, args)

output_dir <- './outputs/'
filen <- paste0(output_dir, args.project,'/pathways_',target_type,'.csv')
RTFTG_results <- read.csv(filen, header = T)
RTFTG_results <- RTFTG_results[RTFTG_results$pred_label > 0.75, ]
RTFTG_results <- RTFTG_results[,1:3] # The activated pathways
Exp_clu <- Get_Exp_Clu(SeuratObj, clusterID = target_type, assay = "RNA", datatype = "data", cutoff = 0.1)
ress <- PRS_calc(Exp_clu, RTFTG_results, cutoff = 0.1)
results_filtered <- Filter_results(ress, PRS_thres = 0.01)

filen <- paste0(output_dir, args.project,'/results_',target_type,'.csv')
write.table(results_filtered,file = filen, quote = F, sep = ",",row.names = F)

library(ggplot2)
library(patchwork)
library(ggalluvial)

groupSize <- as.numeric(table(LR_original@idents))
CellChat::netVisual_circle(LR_original@net$weight, vertex.weight = groupSize, weight.scale = T,
                           label.edge= F, title.name = "Interaction weights/strength",
                           vertex.label.cex = 2)
CCC_threshold <- 0.1
LR_Pairprob <- Extract_LR_Prob(LR_original, target_type = target_type)
LR_Pairprob <- LR_Pairprob[which(LR_Pairprob$Weight >= CCC_threshold * max(LR_Pairprob$Weight)), ]
PlotCCI_CirclePlot(LR_Pairprob, topk = 30)

ps <- PlotXT_Counts(results_filtered, top_percent = 1)
ps$outer
ps$inner

CC_pair_results <- Aggregate_Causality(results_filtered, data_type = "Target")
Counts_pathway <- Count_Crosstalk(results_filtered, KeyGenes = NULL, data_type = "Target", verbose = F)
TG_used <- sort(Counts_pathway, decreasing = T) %>% names()
TG_used <- TG_used[1:15] # Select the top-15 target genes with most crosstalk pathways
PlotXT_RecTGHeatmap(CC_pair_results, Exp_clu, KeyTG = TG_used,topk = 100)

PlotXT_HeatMap(results_filtered, TG_used[1], "TG") # The target with most crosstalk pathways
PlotXT_Alluvial(results_filtered, TG_used[1:3], min_weight = 0.8) # You may adjust the min_weight parameter to display more or less pathways
PlotXT_FidSpe(results_filtered, TG_used[1], threshold = 0.01)
