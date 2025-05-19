# A quick start for SigXTalk

## Load required packages
```
library(SigXTalkR)
library(Seurat)
library(dplyr)
library(CellChat)
```
## Set the working directory
The working directory is where you store the datasets, pathway information and pythoncodes.
The structure of the working directory should be like:
- **my_dir/**
  - **pythoncodes/**
  - **pathways/**
  - **demo**
  - **dataset.rds**
```
work_dir <- ".../my_dir"
setwd(work_dir)
```
## Load the example dataset
The COVID dataset (seurat_covid.rds) is avaliable in the vignettes/ folder (the same folder as this demo).
```
SeuratObj <- readRDS("./seurat_covid.rds") # as the seurat object
cell_anno <- data.frame(cell = names(Idents(SeuratObj)), cluster = Idents(SeuratObj) %>% as.character()) # The metadata of the dataset
```
Note: if you want to use your own dataset, please make sure the dataset is stored as a Seurat Object. The data needs to be normalized, scaled and well-annotated.
Here is a simplified pipeline for the data preprocessing. For a full turotial, visit [here](https://satijalab.org/seurat/articles/pbmc3k_tutorial).
```
# DO NOT run for this tutorial
# Pre-process the data starting from the expression matrix
SeuratObj <- CreateSeuratObject(expression)
SeuratObj[["percent.mt"]] <- PercentageFeatureSet(SeuratObj, pattern = "^MT-")
SeuratObj <- SeuratObj %>% NormalizeData() %>% FindVariableFeatures() 
SeuratObj <- ScaleData(SeuratObj, features = rownames(SeuratObj), vars.to.regress = "percent.mt")
SeuratObj <- SCTransform(SeuratObj,vst.flavor = "v1")
SeuratObj <- RunPCA(SeuratObj) %>% RunUMAP(dims = 1:10)
```

## Load and filter the prior database
The databases could be accessed in the /pathways directory in the GitHub repo.
```
RecTFDB <- readRDS("./pathways/RTF_human.rds") %>% distinct(from, to, .keep_all = T) %>% Filter_DB(rownames(SeuratObj@assays$RNA$data))
TFTGDB <- readRDS("./pathways/TFT_human.rds") %>% distinct(from, to, .keep_all = T) %>% Filter_DB(rownames(SeuratObj@assays$RNA$data))
```

## Infer the cell-cell communication
```
LR_original <- Infer_CCI(SeuratObj, cell_anno, use_spatial = F, db_use = "human")
```

## Prepare the input files for the HGNN module
Note: Please carefully check the input_dir, ensuring it matches the directory that you store the python codes.
The structure of the python codes directory should be like:
- **pythoncodes/**
  - **inputs/**
  - **outputs/**
  - **preprocessing.py**
  - **predictor.py**
  - **training.py**
  - **main_new.py**
```
target_type <- "Fibroblasts"

# Use the differentially expressed genes as the target genes
TG_used <- FindMarkers(SeuratObj, target_type, min.pct = 0.25, only.pos = T, logfc.threshold = 0.25)
TG_used <- filter(TG_used, p_val_adj<1e-3) %>% rownames()

input_dir <- "./pythoncodes/inputs"
Prepare_Input_New(SeuratObj, target_type, TGs = TG_used, CCC_results = LR_original, RecTFDB, TFTGDB, data_dir = input_dir,
                  assay = "RNA", datatype = "scale.data", exp_threshold = 0.05, CCC_threshold = 0.05)
```

## Run the HGNN module to infer activated pathways
Note: Please check the name of the pre-installed python environment which contains the HGNN module. 
Run the HGNN module:
```
args.project <- "COVID" # The name of the project, e.g., the dataset, or any other name you prefer
conda_env <- "SigXTalk_py" # The conda environment that is previously installed to train the hypergraph neural network
python_script <- "./pythoncodes/main_new.py" # The realtive path of the python script
args <- c("--project", args.project, "--target_type",target_type)
system2("conda", args = c("run", "-n", conda_env, "python", python_script, args))
```
Note: You can try different values of hyperparameters using the command line. For example, you can set a higher learning rate by calling:
```
# Optional
args <- c("--project", args.project, "--target_type",target_type, "--lr", 0.05)
system2("conda", args = c("run", "-n", conda_env, "python", python_script, args))
```
The results will be saved in the pythoncodes/outputs directory. You can also save it to other directories by directly modifying the paths in main_new.py file.

## Calculate the PRS
```
output_dir <- './pythoncodes/outputs/'
filen <- paste0(output_dir, args.project,'/pathways_',target_type,'.csv')
RTFTG_results <- read.csv(filen, header = T)
RTFTG_results <- RTFTG_results[RTFTG_results$pred_label > 0.75, ]
RTFTG_results <- RTFTG_results[,1:3] # The activated pathways
Exp_clu <- Get_Exp_Clu(SeuratObj, clusterID = target_type, assay = "RNA", datatype = "data", cutoff = 0.1)
ress <- PRS_calc(Exp_clu, RTFTG_results, cutoff = 0.1)
ress <- ress[!grepl("^MT", ress$Target), ]
ress <- ress[!grepl("^RPL", ress$Target), ]
ress <- ress[!grepl("^RPS", ress$Target), ]
# Filter out low-PRS pathways
results_filtered <- filter(ress, Weight > 0.05*max(ress$Weight))
```

## Save the results (optional)
```
filen <- paste0(output_dir, args.project,'/results_',target_type,'.csv')
write.table(results_filtered,file = filen, quote = F, sep = ",")
```

## Visualize the results
The crosstalk patterns could be visualized in various ways!

### Load packages for visualization
```
library(ggplot2)
library(circlize)
library(ComplexHeatmap)
library(ggalluvial)
library(patchwork)
```
### Visualization of cell-cell communication
Overview of cell-cell communication
```
groupSize <- as.numeric(table(LR_original@idents))
CellChat::netVisual_circle(LR_original@net$weight, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, title.name = "Interaction weights/strength",
                 vertex.label.cex = 2)
```
Ligand-Receptor pairs targeting the receiver
```
CCC_threshold <- 0.1
LR_Pairprob <- Extract_LR_Prob(LR_original, target_type = target_type)
LR_Pairprob <- LR_Pairprob[which(LR_Pairprob$Weight >= CCC_threshold * max(LR_Pairprob$Weight)), ]
PlotCCI_CirclePlot(LR_Pairprob, topk = 30)
```

### Overview of the crosstalk patterns
```
ps <- PlotXT_Counts(results_filtered, top_percent = 1)
ps$outer
ps$inner
```

### Visualize the genes with most crosstalk pathways
```
CC_pair_results <- Aggregate_Causality(results_filtered, data_type = "Target")
Counts_pathway <- Count_Crosstalk(results_filtered, KeyGenes = NULL, data_type = "Target", verbose = F)
TG_used <- sort(Counts_pathway, decreasing = T) %>% names()
TG_used <- TG_used[1:15]
PlotXT_RecTGHeatmap(CC_pair_results, Exp_clu, KeyTG = TG_used,topk = 100)
```

### Detailed crosstalk pattern of a certain target gene
```
TG_used <- "CEBPB"

# Visualize using a single heatmap of PRS
PlotXT_HeatMap(results_filtered, TG_used, "TG")
# Visualize the flow of regulation
PlotXT_Alluvial(results_filtered, TG_used, min_weight = 0.8) # You may adjust the min_weight parameter to display more or less pathways
# Visualize the fidelity and specificity of pathways
PlotXT_FidSpe(results_filtered, TG_used, threshold = 0.01)

```

