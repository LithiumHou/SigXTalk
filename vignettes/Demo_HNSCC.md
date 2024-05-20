# SigXTalk: a step-by-step tutorial

## Load required R packages and source codes 
```
library(Seurat)
library(dplyr)
library(CellChat)
library(MASS)
library(tibble)
library(ggplot2)

setwd("/home/jiawen/myMLnet") # Set the working directory where you store the SigXTalk source codes and inputs
source("./Rcodes/Communications.R")
source("./Rcodes/Utilities.R")
source("./Rcodes/Crosstalk_analysis.R")
```

## Load and preprocess the scRNA-seq datasets (as a Seurat Object) 
We study the HNSCC dataset. The raw data could be found at: https://zenodo.org/records/3260758/files/hnscc_expression.rds?download=1
```{r}
SeuratObj <- readRDS("./datasets/demo_seurat_object.rds") # The HNSCC dataset
cell_anno <- read.table("./datasets/demo_cell_annotation.csv",header = True)
Idents(SeuratObj) <- cell_anno$cluster # The cells are labelled using pre-assigned annotations

SeuratObj <- subset(SeuratObj, nCount_RNA < 25000 & nFeature_RNA<10000) # Quality control
SeuratObj <- SeuratObj %>% NormalizeData() %>% FindVariableFeatures() %>% ScaleData()
SeuratObj <- SCTransform(SeuratObj,vst.flavor = "v1")

cell_anno <- data.frame(cell = names(Idents(SeuratObj)), cluster = Idents(SeuratObj) %>% as.character())
Exp_clu <- Get_Exp_Clu(SeuratObj, clusterID = target_type, assay = "SCT", datatype = "counts") # The formatted expression matrix
write.table(Exp_clu, file = './inputs/ExpressionCount.csv',quote = F, sep = " ")
```

## Load and preprocess the prior gene-gene interaction database s
```{r}
RecTFDB <- readRDS("./datasets/RT_layer2.rds") %>% distinct(from, to, .keep_all = T)
TFTGDB <- readRDS("./datasets/TT_layer3.rds") %>% distinct(from, to, .keep_all = T)
allgenes <- rownames(SeuratObj@assays$RNA$data)
LRDB <- Filter_DB(LRDB,allgenes)
RecTFDB <- Filter_DB(RecTFDB,allgenes)
TFTGDB <- Filter_DB(TFTGDB, allgenes)
# The filtered databases are prepared 
write.table(RecTFDB[,1:2], file = "./inputs/RecTFDB.txt",quote = F,sep = " ")
write.table(TFTGDB[,1:2], file = "./inputs/TFTGDB.txt",quote = F,sep = " ")
```

## Infer the cell-cell communication and identify the activated receptors (signals) of receiver cells
```
LigRec_original <- Infer_CCI(SeuratObj, cellchat_output = T, db_use = "human")
target_type <- "malignant"
LR_Pairprob <- Extract_LR_Prob(LigRec_original, target_type = target_type, cellchat_use = cellchat_use)
Rec_act <- aggregate(LR_Pairprob$Weight, list(LR_Pairprob$To), sum)
colnames(Rec_act) <- c("Rec", "Weight")
Rec_act <- Rec_act[Rec_act$Weight > 0.1*max(Rec_act$Weight),]
LR_Pairprob <- LR_Pairprob[LR_Pairprob$To %in% Rec_act$Rec,]
write.table(LR_Pairprob, file = "./inputs/LigRec.txt",quote = F, sep = " ")
```

## Visualize the cell-cell communication (optional)
```
# Visualize the CCC patterns
CellChat::netVisual_circle(LigRec_original@net$weight, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, title.name = "Interaction weights/strength",
                 vertex.label.cex = 2) 
# Visualize the LR pairs that targets the target type
PlotCCI_CirclePlot(LR_Pairprob, topk = 20)
```

## Input the target genes
Here, we use the genes that are highly expressed in all types of cells as an example. You can also input the group of genes that you're interested in.
```
all_types <- aggregate(cell_anno$cell,list(cell_anno$cluster),length)
types_used <- (all_types %>% filter(x > 200))$Group.1
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
write.table(common_TGs, file = "./inputs/TG.txt",quote = F,sep = " ")
```

## Run the SigXTalk inference
Setup your python environment
```{r}
conda_python <- "/home/jiawen/anaconda3/envs/SigXTalk/bin/python"
```
Before you proceed, please make sure you've successfully generated the correct inputs, including: the expression matrix ("ExpressionCount.csv"), the filtered databases ("RecTFDB.txt" and "TFTGDB.txt"), the results of CCC inference ("LigRec.txt"), and the list of target genes ("TG.txt"). All the files should be in the same directory. Then, run the following:
```{r}
system2(conda_python, args = c("./pythoncodes/main.py"))
```
SigXTalk offers tunable hyperparameters using Argparser. If you want to adjust them (for example, the learning rate and the batch size), you can run the following:
```{r}
system2(conda_python, args = c("./pythoncodes/main.py", paste("--lr",0.02), paste("--batch_size",128)))

```
The output of SigXTalk is saved in "./results/result_demo.csv". The well-trained deep learning model is saved as "./results/demo_model.pth". You may edit these paths and file names in "./pythoncodes/main.py"

## Visualize the results
Load and filter the results:
```{r}
CC_results <- read.csv("./results/result_demo.csv")
colnames(CC_results) <- c("Signal","SSC","TG","PRS")

# Filter out the mitochondrial genes and ribosome genes
CC_results <- CC_results[!grepl("^MT", CC_results$TG), ]
CC_results <- CC_results[!grepl("^RPL", CC_results$TG), ]
CC_results <- CC_results[!grepl("^RPS", CC_results$TG), ]

# Filter out lowly activated pathways
CC_results <- CC_results[CC_results$Weight_all > 0.01*max(CC_results$PRS),]

# Calculate the TRS
CC_pair_results <- Aggregate_Causality(CC_results, sum_type = "Copula",data_type = "TG")
CC_pair_results <- CC_pair_results[CC_pair_results$Weight > 0.05,]
```
Visualize the number of crosstalk pathways that regulates the selected target genes:
```{r}
ps <- PlotXT_Counts(CC_results, datatype = "TG", top_percent = 5)
ps$outer
print(ps$inner, vp = grid::viewport(x = 0.58, y = 0.66, width = 0.4, height = 0.3, just = c("left","bottom")))
```
Visualize the contribution of each signal to the expression of most-expressed target genes:
```{r}
PlotXT_RecTGHeatmap(CC_pair_results, Exp_clu, topk = 15 , ranktype = "expression")
```
Visualize the information flow to these target genes. By default, we use the most-expressed target genes. You can also input targets you're interested in, so long as they have activated pathways.
```{r}
PlotXT_Alluvial(CC_results, KeyTG = NULL, topk = 5,ranktype = "expression",min_weight = 0.7)
```
Visualize the fidelity and specificity of a target gene:
```{r}
KeyTG <- "MALAT1"
PlotXT_FidSpe(CC_results,KeyTG)
```
