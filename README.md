## :hammer:  ctsGRN
**Comprehensive inference reveals a map of transcriptional regulation at cellular resolution in Arabidopsis root**

## :boom: Background
Most plant GRN studies have been performed at the tissue level resolution, making it difficult to disentangle regulatory programs unique to specific cell types or states. At present, there is a lack of cell-type-specific GRNs derived from single-cell data in plants, which limits our understanding of the organization of specific cell types in model species and crops. In this study, we constructed cell- type-specific GRNs for Arabidopsis roots at cell type resolution utilized three distinct algorithms. Machine learning algorithms demonstrates outstanding performance, and validation through Y1H assays and publicly available ChIP-seq data confirms the high predictive accuracy of the GRNs. Utilizing these high-quality networks, we identified numerous previously unknown regulatory relationships; We also pinpointed key TFs and functional modules specific to cell types which play important roles in different biological processes; Moreover, significant heterogeneity exists in TFs and their target genes across GRNs of different cell types.

Our study provides a comprehensive cell-type-resolved overview of transcription factor-mediated gene regulatory networks in Arabidopsis root, offering valuable insights into the molecular mechanisms governing root development and establishing a reference framework for investigating plant biological processes at cellular resolution. 

<p align="center">
  <img src="https://github.com/Jiangyu1220/ctsGRN/blob/main/ctsGRN.png" alt="Schematic illustration of Scan" border="0.1">
</p>

## :book: Description of each file in R and src/DeepSEM folders
- **Case_study.R**: Case study for transcriptional regulation at cellular resolution in Arabidopsis root
## :gear: The usage of ctsGRN
Paste all files into a single folder (set the folder as the directory of R environment). The users can simply run the scripts as follows.

```{r echo=FALSE, results='hide', message=FALSE}
source("R/Case_study.R")
```
## :zap: Quick example to use ctsGRN
