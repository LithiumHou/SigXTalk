# SigXTalk
![Figure](/Fig1.jpg)
## About
This directory contains the code and resources of the following paper:
"Dissecting crosstalk induced by cell-cell communication using single-cell transcriptomic data"

## Introduction
SigXTalk is a deep learning-based computational method to analyze potential crosstalk between multiple regulatory pathways induced by cell-cell communication (CCC). Based on single-cell transcriptomic data and prior knowledge of gene-gene interaction, SigXTalk employs a specialized hypergraph learning framework to identify the crosstalk pathways and further measure their fidelity and specificity using tree-based machine learning approaches. Specifically, SigXTalk aims to:
- identify the regulatory pathways that are activated by CCC signals, using a specialized hypergraph learning approach. Among the activated pathways, the pathways with shared signals, intermediate regulators or targets are identified as crosstalk pathways;
- estimate the regulatory strength of each activated pathway using a tree-based machine learning method; and 
- measure the regulatory selectivity between signals and targets, by introducing the concept of pathway fidelity and specificity.

## Installation
SigXTalk is based on R+Python. The preprocessing (filtering, normalization, scaling and dim-reduction) of datasets, CCC analysis, result visualization are processed with R, while the hypergraph construction and representative learning are processed with Python. Therefore, both the R and Python environments need to be correctly setup before using SigXTalk. In short, the installation of SigXTalk contains the following steps:
- Installation of the dependencies of SigXTalk R package
- Installation of the SigXTalk R package
- Installation of the SigXTalk Python module

If Seurat and/or CellChat have already been installed, the installation of the SigXTalk and other necessary dependencies could be finished in several minutes. However, the de-novo installation (starting from a vanilla R without any external libraries) of SigXTalk R package could be quite time-consuming as every dependency needs a series of Runtime libraries. It may take up to ~30 minutes to get everything done from the very beginning (especially for R within conda environments under Linux system).
The installation of the SigXTalk Python module, from creating the conda environment, could be finished in ~5 minutes.


### Installation of the dependencies of SigXTalk R package
R >= 4.1.1 is required to correctly install SigXTalk and other dependencies. 
Some dependencies are required for installing and running CellChat/SigXTalk but are not availiable on CRAN. We suggest that it be installed manually.
The installation of R dependencies may vary across different operating systems. Please check the corresponding guide that matches your system.

<details>
<summary> Windows users </summary>
  
If you haven't installed Rtools on Windows (which is usually not automatically installed with R), please see [here](https://cran.r-project.org/bin/windows/Rtools).
  
```r
install.packages(c("httpuv","rlang"), type = "source")
install.packages(c("devtools","shiny"), type = "source")  # If you haven't installed devtools before, it may take several minutes.

package_list2 <- c("Biobase","BiocNeighbors","ComplexHeatmap","BiocGenerics")
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager") 
BiocManager::install(package_list2)
```

</details>


<details>
<summary> MacOS users </summary>
  
```r
if (!require("devtools", quietly = TRUE))
    install.packages("devtools") 

package_list <- c("Biobase","BiocNeighbors","ComplexHeatmap","BiocGenerics")
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager") 
BiocManager::install(package_list)
```
<details>
<summary> Troubleshooting </summary>
  
1. Sometimes, you may need XQuartz for the installation. If such error occurs, please visit [here](https://www.xquartz.org/) to install XQuartz. After that, restart R and try the above code again.

2. If you install the dependencies on a vanilla R environment, you may encounter the issue: `library 'gfortran' not found` or `'libintl.h' file not found`. In this case, you will need to add environment variables to ~/.R/Makevars as follows:
STEP 1: Install compilers via Homebrew (using the Terminal):

```bash
brew install gcc gettext
```

STEP 2: Create and open the `~/.R/Makevars` file (using the Terminal):

```bash
mkdir -p ~/.R
touch ~/.R/Makevars
```

STEP 3: Add the following to the `~/.R/Makevars` file:

```
FC = /opt/homebrew/bin/gfortran
F77 = /opt/homebrew/bin/gfortran
CPPFLAGS = -I/opt/homebrew/opt/gettext/include
LDFLAGS = -L/opt/homebrew/Cellar/gcc/15.1.0/lib/gcc/current
FLIBS = -lgfortran
```

Please ensure the file ends with a newline.

STEP 4: Install additional packages (using the Terminal):

```bash
brew install freetype harfbuzz pkg-config
```

STEP 5: Restart Rstudio and retry the installation of dependencies.
</details>

<details>
  <summary>Linux users WITH conda/mamba environment</summary>
If you use R inside a conda/mamba environment, you need to install additional libraries using command lines (not in R) before installing the dependencies:

```bash
conda install -c conda-forge \
  r-devtools r-ggplot2 r-svglite r-ggrepel r-ragg r-systemfonts\
  r-cowplot r-patchwork r-ggpubr r-ggnetwork r-plotly \
  r-mass r-lattice freetype libpng libxml2 libcurl openssl libuv cmake
# If you are using mamba, simply replace 'conda install' with 'mamba install' (but keep 'conda-forge' unchanged)
```

After that, enter R and run the following to install the dependencies:

```R
install.packages("htmltools")

package_list <- c("Biobase","BiocNeighbors","ComplexHeatmap","BiocGenerics")
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager") 
BiocManager::install(package_list)
```
</details>

<details>
  <summary>Linux users WITHOUT conda/mamba environment</summary>
We strongly suggest that you use R within a conda/mamba environment! You can easily run the following command to install R inside a new conda environment:

```bash
conda create -n my_r_env r-base r-devtools
conda activate my_r_env
```

Then, you may refer to the installation guide for Linux users WITH conda/mamba environment
If you are using a system R (not in a conda/mamba environment), it would be quite troublesome to install various libraries.
For Ubuntu/Debian users:

```bash
sudo apt update
sudo apt install -y \
  libfreetype6-dev \
  libpng-dev \
  libxml2-dev \
  libcurl4-openssl-dev \
  libssl-dev \
  libuv1-dev \
  cmake
```

For CentOS users:
```bash
sudo dnf install -y \
  freetype-devel \
  libpng-devel \
  libxml2-devel \
  libcurl-devel \
  openssl-devel \
  libuv-devel \
  cmake
```
After that, enter R and run the following to install the dependencies:

```R
if (!require("devtools", quietly = TRUE))
    install.packages("devtools")  # If you haven't installed devtools before, it may take several minutes.

package_list <- c("Biobase","BiocNeighbors","ComplexHeatmap","BiocGenerics")
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager") 
BiocManager::install(package_list)
install.packages("ragg")
install.packages("svglite")
```
However, there may still be libraries that you need to install manually :(.
</details>

### Installation of the SigXTalk R package
To install the SigXTalk R package, you may either install from remote or from local.
<details>
  <summary>OPTION 1: remote installation</summary>

Run the following command in R:

```R
devtools::install_github("LithiumHou/SigXTalk", dependencies = T, upgrade = "always")
```

Note: using `devtools::install_github` in Rstudio sometimes causes a github's token issue. In this case, you may need to generate a token. Please see [here](https://usethis.r-lib.org/articles/git-credentials.html). Alternatively, you may try local installation (see below).

</details>

<details>
  <summary>OPTION 2: install from local</summary>
You may download or clone the SigXTalk repository to your device and run:
  
```R
if (!require("devtools", quietly = TRUE))
    install.packages("devtools") 
devtools::install("/path/to/SigXTalk") # Replace it with the path where you store the SigXTalk repository
```

</details>

### Installation of the SigXTalk Python module 
SigXTalk requires a Python module to operate correctly. We strongly recommend that an independent python environment (either conda or mamba) be created to run SigXTalk.

<details>
  <summary>Conda users</summary>
  
```bash
conda create -n SigXTalk_py python=3.8
conda activate SigXTalk_py
```
  
</details>

<details>
  <summary>Mamba users</summary>
  
```bash
mamba create -n SigXTalk_py python=3.8
mamba activate SigXTalk_py
```
  
</details>
If you want to use your existing environment, please make sure the version of your Python version is 3.8.X. 

The Python library torch is necessary to perform the deep learning in SigXTalk, which could be run on both CUDA and CPU device. We strongly recommend using CUDA (for Linux and Windows systems only) to accelerate the training of neural network using torch.
The installation command of torch depends on your operating system and device:

<details>
<summary>If you have CUDA device</summary>
  
```bash
# On Linux or Windows only
pip install torch==1.13.1+cu117 --extra-index-url https://download.pytorch.org/whl/cu117
```

</details>

<details>
<summary>If you are using CPU</summary>
If you do not have a CUDA device (especially for MacOS users), you have to use the CPU version of torch. However, it could be a little bit more time-consuming.
  
```bash
# On Linux or Windows
pip install torch==1.13.1+cpu --extra-index-url https://download.pytorch.org/whl/cpu
# On OSX
pip install torch==1.13.1
```

</details>

Then, install the SigXTalk's python module from GitHub: 

<details>
<summary>OPTION 1: using Git</summary>
  
If you have `git` installed on your device, you may install it remotely:
```
# Please make sure you are still in the SigXTalk_py environment
pip install git+https://github.com/LithiumHou/SigXTalk.git#subdirectory=pythoncodes
```

</details>

<details>
<summary>OPTION 2: using local installation</summary>
If the above command does not work, you may manually clone the `pythoncodes` directory to your device and run the following command:
  
```
cd .../pythoncodes
pip install .
```

</details>

You are all set! You are now ready to perform the SigXTalk analysis.

## Usage

A step-by-step tutorial to show the functionality of SigXTalk could be viewed [here](/vignettes/demo.md).

## License
SigXTalk is licensed under the MIT License.
