# A quick start for SigXTalk

A clean R script (no markdown instructions) of this tutorial is available through [here](./demo_pbmc.R).

## Load required packages
```
library(SigXTalkR)
library(Seurat)
library(dplyr)
library(CellChat)
```
## Set the working directory
The working directory is where you store the dataset and python script `main.py`. It should at least contain the following:

- **work_dir/**
  - **main.py**
  - **dataset.rds**

In this tutorial, we simply use the `vignettes` as the directory.
```R
work_dir <- ".../vignettes"
setwd(work_dir)
```
## Load the example dataset
The PBMC dataset (SigXTalk_demo_data.rds) is available via the [Google Drive Link](https://drive.google.com/file/d/1e019IYCU_jV90FmCjJsPj0f1kvnzRf7u/view?usp=sharing) or [Zenodo Link](https://zenodo.org/records/15531285).
```R
SeuratObj <- readRDS("./SigXTalk_demo_data.rds") # as the seurat object
cell_anno <- data.frame(cell = names(Idents(SeuratObj)), cluster = Idents(SeuratObj) %>% as.character()) # The metadata of the dataset
```
Note: the example data has been processed using the R script [here](Process_pbmc.R). For a full turotial on how to process raw data with Seurat, visit [Seurat's tutorial for pbmc3k data](https://satijalab.org/seurat/articles/pbmc3k_tutorial). 

<details>
  <summary>Here is a quick guide for processing your data</summary>
  
```R
# DO NOT run for this tutorial
# Pre-process the data starting from the expression matrix
# No quality control is performed here. Please see Seurat's tutorial for details on filtering out low-quality cells. 
SeuratObj <- CreateSeuratObject(expression) # Create a Seurat object from the count matrix
SeuratObj[["percent.mt"]] <- PercentageFeatureSet(SeuratObj, pattern = "^MT-")
SeuratObj <- SeuratObj %>% NormalizeData() %>% FindVariableFeatures() 
SeuratObj <- ScaleData(SeuratObj, features = rownames(SeuratObj), vars.to.regress = "percent.mt")
SeuratObj <- RunPCA(SeuratObj) %>% RunUMAP(dims = 1:10)
```
</details>
If you want to use your own dataset, please make sure the dataset is stored as a Seurat Object. The data needs to be normalized, scaled and well-annotated (as above).

## Visualize the data (OPTIONAL)
You may visualize the data using UMAP:
```R
DimPlot(SeuratObj, reduction = "umap", pt.size = 1)
```

## Load and filter the prior database
The databases could be easily accessed using SigXTalk's built-in data.
```R
data("RTF_human")
data("TFT_human")
RecTFDB <- RTF_human %>% distinct(from, to, .keep_all = T) %>% Filter_DB(rownames(SeuratObj@assays$RNA$data))
TFTGDB <- TFT_human %>% distinct(from, to, .keep_all = T) %>% Filter_DB(rownames(SeuratObj@assays$RNA$data))
```

## Infer the cell-cell communication using CellChat
```R
LR_original <- Infer_CCI(SeuratObj, cell_anno, use_spatial = F, db_use = "human")
```

## Prepare the input files for the HGNN module

We use the differentially expressed genes as the target genes. You may also adjust the arguments of `FindMarkers` if you want a different list of target genes.
```R
target_type <- "CD14_Mono"
TG_used <- FindMarkers(SeuratObj, target_type, min.pct = 0.25, only.pos = T, logfc.threshold = 0.25)
TG_used <- filter(TG_used, p_val_adj<1e-3) %>% rownames()

input_dir <- "./inputs"
if(!dir.exists(input_dir)){
  dir.create(input_dir)
}


Prepare_Input(SeuratObj, target_type, TGs = TG_used, CCC_results = LR_original, RecTFDB, TFTGDB, data_dir = input_dir,
                  assay = "RNA", datatype = "scale.data", exp_threshold = 0.05, CCC_threshold = 0.05)
```

## Run the HGNN module to infer activated pathways
Note: Please check the name of the pre-installed python environment which contains the HGNN module. 
Set up the arguments passed to Python:
```R
args.project <- "PBMC" # The name of the project, e.g., the dataset, or any other name you like
python_script <- "./main.py" # The realtive path of the python script
args <- c("--project", args.project, "--target_type",target_type)
```
Note: You may try different values of hyperparameters using the command line. For example, you can set a higher learning rate by calling:
```R
# Optional
args <- c("--project", args.project, "--target_type",target_type, "--lr", 0.05)
```
Next, we need to run the python script `main.py`. If you're using a conda environment:
```R
conda_env <- "SigXTalk_py" # The conda environment that is previously installed to train the hypergraph neural network
Run_py_script(python_script, conda_env, args)
```

If you're using a mamba environment, you need to tell `reticulate` package to use a mamba environment as if it is a conda environment. Here, you could find the path to mamba by calling `where mamba` in the command line.
```R
conda_binary("/path/to/mamba")
mamba_env <- "SigXTalk_py" # The mamba environment that is previously installed to train the hypergraph neural network
Run_py_script(python_script, mamba_env, args)
```
<details>
<summary>TROUBLESHOOTING: If you encounter errors with reticulate::usecondaenv()</summary>

On a Windows terminal, sometimes you cannot locate the conda environment and run into error:

```R
Run_py_script(python_script, conda_env, args)
Error in reticulate::use_condaenv(conda_env, required = TRUE): 
  Unable to locate conda environment 'SigXTalk_py'.
```
In this case, you may use the absolute path of SigXTalk_py:
```R
conda_env <- "/path/to/SigXTalk_py"
Run_py_script(python_script, conda_env, args)
```
The path of the conda environment could be accessed by running `conda env list` in the command line (not in R)
Alternatively, you may use the system2() function:
```R
conda_env <- "SigXTalk_py"
system2("conda", args = c("run", "-n", conda_env, "python", python_script, args))
```

</details>

The results will be automatically saved in the `/outputs` directory. You can also save it to other directories by directly modifying the paths in `main.py` file.

## Calculate and filter the PRS
```R
output_dir <- './outputs/'
filen <- paste0(output_dir, args.project,'/pathways_',target_type,'.csv')
RTFTG_results <- read.csv(filen, header = T)
RTFTG_results <- RTFTG_results[RTFTG_results$pred_label > 0.75, ]
RTFTG_results <- RTFTG_results[,1:3] # The activated pathways
Exp_clu <- Get_Exp_Clu(SeuratObj, clusterID = target_type, assay = "RNA", datatype = "data", cutoff = 0.1)
ress <- PRS_calc(Exp_clu, RTFTG_results, cutoff = 0.1)
# Filter out the low-PRS pathways
results_filtered <- Filter_results(ress, PRS_threshold = 0.01)
```

## Save the results (optional)
```R
filen <- paste0(output_dir, args.project,'/results_',target_type,'.csv')
write.table(results_filtered,file = filen, quote = F, sep = ",",row.names = F)
```

## Visualize the results
The crosstalk patterns could be visualized in various ways!

### Load packages for visualization
```R
library(ggplot2)
library(patchwork)
library(ggalluvial)
```
### Visualization of cell-cell communication
Overview of cell-cell communication
```R
groupSize <- as.numeric(table(LR_original@idents))
CellChat::netVisual_circle(LR_original@net$weight, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, title.name = "Interaction weights/strength",
                 vertex.label.cex = 2)
```
For more visualization of the cell-cell communication network, please visit [CellChat's tutorial](https://github.com/jinworks/CellChat).
Ligand-Receptor pairs targeting the receiver
```R
CCC_threshold <- 0.1
LR_Pairprob <- Extract_LR_Prob(LR_original, target_type = target_type)
LR_Pairprob <- LR_Pairprob[which(LR_Pairprob$Weight >= CCC_threshold * max(LR_Pairprob$Weight)), ]
PlotCCI_CirclePlot(LR_Pairprob, topk = 30)
```

### Overview of the crosstalk patterns
```R
ps <- PlotXT_Counts(results_filtered, top_percent = 1)
ps$outer
ps$inner
```

### Visualize the genes with most crosstalk pathways
```R
CC_pair_results <- Aggregate_Causality(results_filtered, data_type = "Target")
Counts_pathway <- Count_Crosstalk(results_filtered, KeyGenes = NULL, data_type = "Target", verbose = F)
TG_used <- sort(Counts_pathway, decreasing = T) %>% names()
TG_used <- TG_used[1:15] # Select the top-15 target genes with most crosstalk pathways
PlotXT_RecTGHeatmap(CC_pair_results, Exp_clu, KeyTG = TG_used,topk = 100)
```

### Detailed crosstalk pattern of a certain target gene
```R
# Visualize using a single heatmap of PRS
PlotXT_HeatMap(results_filtered, TG_used[1], "TG") # The target with most crosstalk pathways
# Visualize the flow of regulation
PlotXT_Alluvial(results_filtered, TG_used[1:3], min_weight = 0.8) # You may adjust the min_weight parameter to display more or less pathways
# Visualize the fidelity and specificity of pathways
PlotXT_FidSpe(results_filtered, TG_used[1], threshold = 0.01)
```
