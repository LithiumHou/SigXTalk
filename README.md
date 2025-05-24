# SigXTalk
![Figure](/Fig1.jpg)
## About
This directory contains the code and resources of the following paper:
"Dissecting crosstalk induced by cell-cell communication using single-cell transcriptomic data"

## Introduction
SigXTalk is a deep learning-based computational method to analyze potential crosstalk between multiple regulatory pathways induced by cell-cell communication (CCC). Based on single-cell transcriptomic data and prior knowledge of gene-gene interaction, SigXTalk employs a specialized hypergraph learning framework to identify the crosstalk pathways and further measure their fidelity and specificity using tree-based machine learning approaches. Specifically, SigXTalk aims to:
-	identify the regulatory pathways that are activated by CCC signals, using a specialized hypergraph learning approach. Among the activated pathways, the pathways with shared signals, intermediate regulators or targets are identified as crosstalk pathways;
-	estimate the regulatory strength of each activated pathway using a tree-based machine learning method; and 
- measure the regulatory selectivity between signals and targets, by introducing the concept of pathway fidelity and specificity.

## Installation
SigXTalk is based on R+Python. The preprocessing (filtering, normalization, scaling and dim-reduction) of datasets, CCC analysis, result visualization are processed with R, while the hypergraph construction and representative learning are processed with Python. Therefore, both the R and Python environments need to be correctly setup before using SigXTalk. In short, the de-novo installation of SigXTalk contains the following steps:
- Installation of the dependencies of SigXTalk R package
- Installation of the SigXTalk R package
- Installation of the SigXTalk Python module
For users with experience on scRNA-seq data analysis using R (especially Seurat), the installation of the SigXTalk and other necessary dependencies could be finished in several minutes. However, the de-novo installation (starting from a vanilla R without any external libraries) of SigXTalk R package could be quite time-consuming as every dependency needs a series of Runtime libraries. It may take up to ~30 minutes to get everything done from the very beginning (especially for R within conda environments under Linux system).
The installation of the SigXTalk Python module, from creating the conda environment, could be finished in ~5 minutes.


### Installation of the dependencies of SigXTalk R package
R >= 4.1.1 is required to correctly install SigXTalk and other dependencies. 
Some dependencies are required for installing and running CellChat/SigXTalk but are not availiable on CRAN. We suggest that it be installed manually.
```
if (!require("devtools", quietly = TRUE))
    install.packages("devtools")  # If you haven't installed devtools before, it may take several minutes.

package_list <- c("Biobase","BiocNeighbors","ComplexHeatmap","BiocGenerics")
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager") 
BiocManager::install(package_list)
```
* If you are using R within a conda/mamba environment (especially if R is newly installed), the installation of dependencies (especially Seurat and CellChat) may become quite time-consuming and even annoying. You need to install additional libraries using command lines (not in R) before installing the dependencies:
```
conda install -c conda-forge \
  r-devtools r-ggplot2 r-svglite r-ggrepel \
  r-cowplot r-patchwork r-ggpubr r-ggnetwork r-plotly \
  r-mass r-lattice freetype libpng libxml2 libcurl openssl libuv cmake 
# If you are using mamba, simply replace 'conda' with 'mamba' (but keep 'conda-forge' unchanged)
``` 
### Installation of the SigXTalk R package
If you install SigXTalk on Windows, you will need to install Rtools. If you did not install Rtools while installing R/Rstudio, please see this [guide]
To install the SigXTalk R package, you may either install from remote or from local.

OPTION 1: remote installation. Run the following command in R:
```
devtools::install_github("LithiumHou/SigXTalk", dependencies = T, upgrade = "always")
```
Note: using `devtools::install_github` in Rstudio sometimes causes a github's token issue. In this case, you may need to generate a token. Please see [here](https://usethis.r-lib.org/articles/git-credentials.html). Alternatively, you may try local installation (see below).

OPTION 2: install from local. You may download or clone the SigXTalk repository to your device and run:
```
if (!require("devtools", quietly = TRUE))
    install.packages("devtools") 
devtools::install("/path/to/SigXTalk") # Replace it with the path where you store the SigXTalk repository
```
* CellChat is one of the non-CRAN dependencies required by SigXTalk. If you encounter any issue while installing CellChat, please visit the [CellChat homepage](https://github.com/jinworks/CellChat) for troubleshooting.
* If you encounter error when installing Seurat & CellChat, it's highly possible that the package `Matrix` is not installed correctly. Please try:
```
install.packages("https://cran.r-project.org/src/contrib/Archive/Matrix/Matrix_1.6-4.tar.gz")
```
and then try `devtools::install_github("LithiumHou/SigXTalk", dependencies = T, upgrade = "always")` again.

### Installation of the SigXTalk Python module 
SigXTalk requires a Python module to operate correctly. We strongly recommend that an independent python environment be created to run SigXTalk.
If you are using conda (Anaconda or Miniconda) environments:
```
conda create -n SigXTalk_py python=3.8
conda activate SigXTalk_py
```
Alternatively, if you prefer mamba environments:
```
mamba create --name SigXTalk_py python=3.8
mamba activate SigXTalk_py
```
If you want to use your own environment, please make sure the version of your Python version is 3.8.X. 
The Python library torch is necessary to perform the deep learning in SigXTalk, which could be run on both CUDA and CPU device. We strongly recommend using CUDA (for Linux and Windows systems only) to accelerate the training of neural network using torch.
The installation command of torch depends on your operating system and device and may cause compatibility issues (which is why we prefer a separated installation of torch instead of integrating it to the installation of other dependencies). If you encounter any issue, or want to know more details, please visit the [torch installation guide](https://pytorch.org/get-started/locally/).

```
# On Linux or Windows
pip install torch==1.13.1+cu117 --extra-index-url https://download.pytorch.org/whl/cu117
# On OSX
pip install torch==1.13.1
```
If you do not have a CUDA device, you may use the CPU version of torch. However, it could be quite time-consuming.
```
# On Linux or Windows
pip install torch==1.13.1+cpu --extra-index-url https://download.pytorch.org/whl/cpu
# On OSX
pip install torch==1.13.1
```
Then, install the SigXTalk's python dependency from GitHub: 

OPTION 1： If you have `git` installed on your device, you may install it remotely:
```
# Please make sure you are still in the SigXTalk_py environment
pip install git+https://github.com/LithiumHou/SigXTalk.git#subdirectory=pythoncodes
```
OPTION 2: If the above command does not work, you may manually clone the `pythoncodes` directory to your device and run the following command:
```
cd .../pythoncodes
pip install .
```
You are all set! You are now ready to perform the SigXTalk analysis.

## Usage

### Basic Usage
```
# Infer the cell-cell communication
LR_original <- Infer_CCI(SeuratObj, cell_anno, use_spatial = F, db_use = "human")
# Prepare the input for the HGNN+ module
Prepare_Input_New(SeuratObj, target_type, TGs = TG_used, CCC_results = LR_original, RecTFDB, TFTGDB, data_dir = input_dir,
                  assay = "RNA", datatype = "scale.data")
# Perform the HGNN
system2(conda_python, args = c("pythoncodes/main_new.py", paste("--project",shQuote(args.project)), paste("--target_type",args.target)))
# Calculate the PRS
ress <- PRS_calc(Exp_clu, RTFTG_results, cutoff = 0.1)
# Visualize the results
PlotXT_Alluvial(results_filtered, TG_used, min_weight = 0.8)
```
### Demo
A step-by-step tutorial to show the functionality of SigXTalk could be viewed [here](/vignettes/demo.md).

## License
SigXTalk is licensed under the MIT License.



