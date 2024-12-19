# SigXTalk
![Figure](/vignettes/Figure.jpg)
## About
This directory contains the code and resources of the following paper:
"Dissecting crosstalk induced by cell-cell communication using single-cell transcriptomic data"

## Introduction
SigXTalk is a deep learning-based computational method to analyze potential crosstalk between multiple regulatory pathways induced by cell-cell communication (CCC). Based on single-cell transcriptomic data and prior knowledge of gene-gene interaction, SigXTalk employs a specialized hypergraph learning framework to identify the crosstalk pathways and further measure their fidelity and specificity using tree-based machine learning approaches. Specifically, SigXTalk aims to:
-	identify the regulatory pathways that are activated by CCC signals, using a specialized hypergraph learning approach. Among the activated pathways, the pathways with shared signals, intermediate regulators or targets are identified as crosstalk pathways;
-	estimate the regulatory strength of each activated pathway using a tree-based machine learning method; and 
- measure the regulatory selectivity between signals and targets, by introducing the concept of pathway fidelity and specificity.

## Installation
SigXTalk is based on R+Python. The preprocessing (filtering, normalization, scaling and dim-reduction) of datasets, CCC analysis, result visualization are processed with R, while the hypergraph construction and representative learning are processed with Python. Therefore, both the R and Python environments need to be correctly setup before using SigXTalk.

### R code dependencies before installing SigXTalk
- R 4.3.1
* Seurat 5.1.0
* CellChat 2.1.1
You may install the Seurat R package using CRAN:
'''
install.packages("Seurat")
'''
Please visit the [CellChat homepage](https://github.com/jinworks/CellChat) for the installation of CellChat.

### Suggested R packages for plotting only
* ComplexHeatmap 2.16.0
* ggalluvial 0.12.5
* ggridges 0.5.6
* patchwork 1.2.0
* scales 1.3.0

### Installation of the SigXTalk R package
Now, you can run the following code in R to install the SigXTalk R package:
'''

'''
  
## Installation of the SigXTalk Python code dependencies
We recommend that an independent python environment be used to run SigXTalk.
```
conda create -f my_env.yml
conda activate SigXTalk_py
```
SigXTalk could be run on both CUDA and CPU devices. We recommend using CUDA to accelerate the training of neural network using torch:

```
pip install torch==1.13.1+cu117 torchvision==0.14.1+cu117 torchaudio==0.13.1 --extra-index-url https://download.pytorch.org/whl/cu117
```
If you do not have a CUDA device, you may use the CPU version of torch:
'''
pip install torch==1.13.1+cpu torchvision==0.14.1+cpu torchaudio==0.13.1 --extra-index-url https://download.pytorch.org/whl/cpu
'''

Then, the dhg package is required to handle the hypergraph object:
'''
pip install dhg
'''

That's it! You are all set for the SigXTalk installation.

## Usage

### Demo
A step-by-step tutorial to show the functionality of SigXTalk could be viewed [here](https://github.com/LithiumHou/SigXTalk/blob/master/vignettes/Demo_HNSCC.md).

## License
SigXTalk is licensed under the MIT License.



