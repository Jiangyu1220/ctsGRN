## :compass:  ctsGRN
**Comprehensive inference reveals a map of transcriptional regulation at cellular resolution in Arabidopsis root**

## :fireworks: Background
Most plant GRN studies have been performed at the tissue level resolution, making it difficult to disentangle regulatory programs unique to specific cell types or states. At present, there is a lack of cell-type-specific GRNs derived from single-cell data in plants, which limits our understanding of the organization of specific cell types in model species and crops. In this study, we constructed cell- type-specific GRNs for Arabidopsis roots at cell type resolution utilized three distinct algorithms. Machine learning algorithms demonstrates outstanding performance, and validation through Y1H assays and publicly available ChIP-seq data confirms the high predictive accuracy of the GRNs. Utilizing these high-quality networks, we identified numerous previously unknown regulatory relationships; We also pinpointed key TFs and functional modules specific to cell types which play important roles in different biological processes; Moreover, significant heterogeneity exists in TFs and their target genes across GRNs of different cell types.

Our study provides a comprehensive cell-type-resolved overview of transcription factor-mediated gene regulatory networks in Arabidopsis root, offering valuable insights into the molecular mechanisms governing root development and establishing a reference framework for investigating plant biological processes at cellular resolution. 

A schematic illustration of cteGRN is shown in the folowing.

<p align="center">
  <img src="https://github.com/Jiangyu1220/ctsGRN/blob/main/ctsGRN.png" alt="Schematic illustration of Scan" border="0.1">
</p>

## :book: Description of each file in R and src/DeepSEM folders
- **Case_study.R**: Case study for transcriptional regulation at cellular resolution in Arabidopsis root.
- **model.py**:Contains the deep generative model architecture (Variational Autoencoder, VAE) for GRN inference.Implements the encoder-decoder structure with gene expression as input and regulatory probabilities as output.
- **train.py**:Training pipeline for DeepSEM, including data loading, loss computation, and gradient updates.Supports early stopping and learning rate scheduling.
- **utils.py**:Preprocessing (normalization, batch correction) and postprocessing (network thresholding).Helper functions for visualization (e.g., AUPRC curves) and data splitting.
- **inference.py**:Infers gene regulatory networks from the trained VAE’s latent space.Computes edge confidence scores (e.g., gradient-based importance).

## :gear: The usage of ctsGRN
Paste all files into a single folder (set the folder as the directory of R environment). The users can simply run the scripts as follows.

```{r echo=FALSE, results='hide', message=FALSE}
source("R/Case_study.R")
```
## :zap: Quick example to use ctsGRN
```
# Load scRNA-seq dataset
Set the working Directory (Switch the working directory to the specified path)

# Load packages
library(Seurat)
library(ggplot2)
library(clustree)
library(cowplot)
library(dplyr)
library(R.utils)
library(pROC)
library(dplyr)
library(org.At.tair.db)
library(clusterProfiler)
library(AnnotationDbi)
library(ggplot2)
library(aPEAR)
library(DOSE)
library(cols4all)
library(UpSetR)
library(ggplot2)
library(dplyr)
library(purrr)
library(miRspongeR)
##=======================Create a Seurat object========================

seurat_data<- read.csv(gzfile("Data/GSE123818_matrix.csv.gz"), row.names = 1)
seurat_obj <- CreateSeuratObject(counts = seurat_data,
                                 min.features = 200,
                                 min.cells = 3, 
                                 project = "GSE123818")

##=======================Data quality control and standardization================================
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern='^ATMG')
violin <- VlnPlot(seurat_obj,
                  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                  pt.size = 0.1,
                  ncol = 3)
##Set quality control standards
seurat_obj<-subset(seurat_obj,subset=nFeature_RNA>500 & nFeature_RNA<5000 &percent.mt<0.5)
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)

#########################Three algorithms are used to construct GRN#############################
######################### WGCNA#############################
######################### GENIE3#############################
######################### DeepSEM#############################
load("Data/celltype.RData")
load("Data/celltype_TFs.RData")
``` 
