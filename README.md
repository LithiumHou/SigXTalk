# SigXTalk
## About
This directory contains the code and resources of the following paper:
"Dissecting crosstalk induced by cell-cell communication using single-cell transcriptomic data"

## Introduction
SigXTalk is a deep learning-based computational method to analyze potential crosstalk between multiple regulatory pathways induced by cell-cell communication (CCC). Based on single-cell transcriptomic data and prior knowledge of gene-gene interaction, SigXTalk employs a specialized hypergraph learning framework to identify the crosstalk pathways and further measure their fidelity and specificity using tree-based machine learning approaches. Specifically, SigXTalk aims to:
-	identify the regulatory pathways that are activated by CCC signals, using a specialized hypergraph learning approach. Among the activated pathways, the pathways with shared signals, intermediate regulators or targets are identified as crosstalk pathways;
-	estimate the regulatory strength of each activated pathway using a tree-based machine learning method; and 
- measure the regulatory selectivity between signals and targets, by introducing the concept of pathway fidelity and specificity. 

## Dependency

### R code dependency
- R 4.3.1
* Seurat 5.1.0
* SeuratObject 5.0.2
* Cellchat 2.1.1
* dplyr 1.1.4
* ggplot2 3.5.1

### Python code dependency  
- python==3.8.12  
* numpy==1.24.3  
* torch==1.13.1 (CUDA version)  
* pandas==2.0.3  
* dhg==0.9.5  
* xgboost==2.0.3  
* scipy==1.10.1  
* scikit-learn==1.3.0  
* argparse==1.4.0
  
All dependencies can be installed within a few minutes.
## Usage

### Installation
No extra installation is needed. 
### Demo
A step-by-step tutorial to show the functionality of SigXTalk could be viewed [here](https://github.com/LithiumHou/SigXTalk/blob/master/vignettes/Demo_HNSCC.md).

## License
SigXTalk is licensed under the MIT License.



