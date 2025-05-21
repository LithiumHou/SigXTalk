# SigXTalk
![Figure](/vignettes/Fig1.jpg)
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
- Installing the SigXTalk R package
- Create the Python environment for SigXTalk
- Installing the SigXTalk Python package
Usually, the installation of the SigXTalk and other necessary dependencies could be finished in several minutes.


### Installation of the SigXTalk R package (REQUIRED)
R >= 4.1.1 is required to correctly install SigXTalk and other dependencies. All necessary dependencies will be automatically installed when installing the SigXTalk R package. 
To install the SigXTalk R package, run the following command in R:
```
remotes::install_github("LithiumHou/SigXTalk")
```
* CellChat is one of the non-CRAN dependencies required by SigXTalk. If you encounter any issue while installing CellChat, please visit the [CellChat homepage](https://github.com/jinworks/CellChat) for troubleshooting.
  
### Suggested R packages for plotting (OPTIONAL)
These packages are not necessary for performing SigXTalk algorithm, but are useful when visualizing the results!
* ggridges 0.5.6
* patchwork 1.2.0
* scales 1.3.0

### Installation of the SigXTalk Python code dependencies (REQUIRED)
SigXTalk requires a Python module to operate correctly. We recommend that an independent python environment be created to run SigXTalk.
```
conda create -n SigXTalk_py python=3.8
conda activate SigXTalk_py
pip install pandas==2.0.3 scikit-learn==1.3.0 scipy==1.10.1 numpy==1.24.3 argparse==1.4.0
```
The Python package torch is necessary to perform the deep learning in SigXTalk, which could be run on both CUDA and CPU devices. We strongly recommend using CUDA to accelerate the training of neural network using torch.
The installation command of torch depends on your operating system and device. If you encounter any issue, or want to know more details, please visit the [torch installation guide](https://pytorch.org/get-started/locally/).

```
# On Linux or Windows
pip install torch==1.13.1+cu117 torchvision==0.14.1+cu117 torchaudio==0.13.1 --extra-index-url https://download.pytorch.org/whl/cu117
# On OSX
pip install torch==1.13.1 torchvision==0.14.1 torchaudio==0.13.1
```
If you do not have a CUDA device, you may use the CPU version of torch. However, it could be quite time-consuming.
```
# On Linux or Windows
pip install torch==1.13.1+cpu torchvision==0.14.1+cpu torchaudio==0.13.1 --extra-index-url https://download.pytorch.org/whl/cpu
# On OSX
pip install torch==1.13.1 torchvision==0.14.1 torchaudio==0.13.1
```
Then, the dhg package is required to handle the hypergraph object:
```
pip install dhg
```

That's it! You are now ready to perform the SigXTalk analysis.

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
A step-by-step tutorial to show the functionality of SigXTalk could be viewed [here](https://github.com/LithiumHou/SigXTalk/blob/master/vignettes/demo.md).

## License
SigXTalk is licensed under the MIT License.



