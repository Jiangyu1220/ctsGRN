#########################一.scRNA-seq data processing#############################
library(Seurat)
library(ggplot2)
library(clustree)
library(cowplot)
library(dplyr)
library(R.utils)

# Set the working Directory (Switch the working directory to the specified path)
setwd("E:/论文/AtScGRN")

##=======================1.Create a Seurat object========================

seurat_data<- read.csv(gzfile("Data/GSE123818_matrix.csv.gz"), row.names = 1)
seurat_obj <- CreateSeuratObject(counts = seurat_data,
                                 min.features = 200,
                                 min.cells = 3, 
                                 project = "GSE123818")

##=======================2.Data quality control and standardization================================
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern='^ATMG')
violin <- VlnPlot(seurat_obj,
                  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                  pt.size = 0.1,
                  ncol = 3)
##Set quality control standards
seurat_obj<-subset(seurat_obj,subset=nFeature_RNA>500 & nFeature_RNA<5000 &percent.mt<0.5)
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)

##=======================3.Data dimension reduction and clustering==================================
seurat_obj <- FindVariableFeatures(seurat_obj,mean.cutoff=c(0.0125,3),dispersion.cutoff =c(1.5,Inf) )
top10 <- head(VariableFeatures(seurat_obj), 10)
scale.genes <-  rownames(seurat_obj)
seurat_obj <- ScaleData(seurat_obj, features = scale.genes)
## PCA reduces dimension and extracts principal components
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj),npcs = 20) 
## Cell clustering
seurat_obj <- FindNeighbors(object = seurat_obj, dims = 1:20)
seurat_obj <- FindClusters(object = seurat_obj, resolution = 0.59)
table(seurat_obj@meta.data$seurat_clusters)
## tsne
seurat_obj <- RunTSNE(seurat_obj, dims =1:20)
embed_tsne <- Embeddings(seurat_obj, 'tsne')
## UMAP'
seurat_obj <- RunUMAP(seurat_obj,n.neighbors = 30,metric = 'correlation',min.dist = 0.3,dims = 1:20)
embed_umap <- Embeddings(seurat_obj, 'umap')
##=======================3.Marker gene screening==================================
seurat_obj <- JoinLayers(object = seurat_obj, suffix = "_joined")
markers <- FindAllMarkers(object = seurat_obj, test.use="wilcox" ,
                          only.pos = TRUE,
                          min.pct = 0.25,
                          logfc.threshold = 0.25)

all.markers =markers %>% dplyr::select(gene, everything()) %>% subset(p_val<0.05)
top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

View(seurat_obj@meta.data)
table(seurat_obj@meta.data$seurat_clusters)
##=======================4.Cell type annotation==================================
all_markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
all_markers %>%  group_by(cluster) %>%  slice_max(n = 2, order_by = avg_log2FC)
Arabi_complete_set=read.delim("Data/ara_mkgs.txt", header = FALSE,sep = "\t",comment.char = "#")
colnames(Arabi_complete_set)[1] = c("gene")
con <- inner_join(all_markers,Arabi_complete_set,by ="gene")
subset_root <- filter(con, V3 == "Root")
subset_root=subset_root[,-c(1:5)]
distinct_root <- distinct(subset_root, gene, .keep_all = TRUE)

result <- distinct_root %>%
  group_by(cluster, V4) %>%
  summarize(count = n())

sorted_result <- result %>%
  group_by(cluster) %>%
  arrange(cluster, desc(count))

#细胞类型注释
new.cluster.ids <- c("0" = "Columella", "1" = "Lateral root cap", "2" = "Lateral root cap", "3" = "Atrichoblast", "4"= "Metaphloem sieve element",  "5" = "Metaphloem sieve element", "6" = "Endodermis" , "7" = "Xylem", "8" = "Lateral root cap", "9" = "Cortex", "10" = "Endodermis", "11" = "Endodermis","12" = "Trichoblast", "13" = "Mature","14" = "Cortex")
names(new.cluster.ids) <- levels(seurat_obj)
seurat_obj <- RenameIdents(seurat_obj, new.cluster.ids)
DimPlot(seurat_obj, reduction = "umap", label = TRUE, pt.size = 1) + NoLegend()
levels(Idents(seurat_obj))

#Lateral root cap
Lateral_root_cap <- subset(seurat_obj, idents = "Lateral root cap")@assays$RNA$data
Lateral_root_cap_expr=data.frame(Lateral_root_cap)
dim(Lateral_root_cap_expr)

#Columella
Columella <- subset(seurat_obj, idents = "Columella")@assays$RNA$data
Columella_expr=data.frame(Columella)
dim(Columella_expr)

#Atrichoblast
Atrichoblast <- subset(seurat_obj, idents = "Atrichoblast")@assays$RNA$data
Atrichoblast_expr=data.frame(Atrichoblast)
dim(Atrichoblast_expr)

#Metaphloem sieve element
Metaphloem_sieve_element <- subset(seurat_obj, idents = "Metaphloem sieve element")@assays$RNA$data
Metaphloem_sieve_element_expr=data.frame(Metaphloem_sieve_element)
dim(Metaphloem_sieve_element_expr)

#Endodermis
Endodermis <- subset(seurat_obj, idents = "Endodermis")@assays$RNA$data
Endodermis_expr=data.frame(Endodermis)
dim(Endodermis_expr)

#Xylem
Xylem <- subset(seurat_obj, idents = "Xylem")@assays$RNA$data
Xylem_expr=data.frame(Xylem)
dim(Xylem_expr)

#Cortex
Cortex <- subset(seurat_obj, idents = "Cortex")@assays$RNA$data
Cortex_expr=data.frame(Cortex)
dim(Cortex_expr)

#Trichoblast
Trichoblast <- subset(seurat_obj, idents = "Trichoblast")@assays$RNA$data
Trichoblast_expr=data.frame(Trichoblast)
dim(Trichoblast_expr)

#Mature
Mature <- subset(seurat_obj, idents = "Mature")@assays$RNA$data
Mature_expr=data.frame(Mature)
dim(Mature_expr)

Atrichoblast <- as.matrix(Atrichoblast_expr)
Columella <- as.matrix(Columella_expr)
Cortex <- as.matrix(Cortex_expr)
Endodermis <- as.matrix(Endodermis_expr)
Lateral_root_cap <- as.matrix(Lateral_root_cap_expr)
Mature <- as.matrix(Mature_expr)
Metaphloem_sieve_element <- as.matrix(Metaphloem_sieve_element_expr)
Trichoblast <- as.matrix(Trichoblast_expr)
Xylem <- as.matrix(Xylem_expr)


Atrichoblast <- Atrichoblast[rowSums(Atrichoblast != 0) > 0, ]
Columella <- Columella[rowSums(Columella != 0) > 0, ]
Cortex <- Cortex[rowSums(Cortex != 0) > 0, ]
Endodermis <- Endodermis[rowSums(Endodermis != 0) > 0, ]
Lateral_root_cap <- Lateral_root_cap[rowSums(Lateral_root_cap != 0) > 0, ]
Mature <- Mature[rowSums(Mature != 0) > 0, ]
Metaphloem_sieve_element <- Metaphloem_sieve_element[rowSums(Metaphloem_sieve_element != 0) > 0, ]
Trichoblast <- Trichoblast[rowSums(Trichoblast != 0) > 0, ]
Xylem <- Xylem[rowSums(Xylem != 0) > 0, ]

save(Lateral_root_cap,Columella,Atrichoblast, Metaphloem_sieve_element,
     Endodermis,Xylem,Cortex,Trichoblast,Mature,
     file = "Data/celltype.RData")


tf_all <- read.table("Data/Ath_TF_list.txt", header = TRUE, sep = "\t")
tf_all <- intersect(tf_all$Gene_ID,scale.genes) # 

# Find TFs in each tissue type.
Columella_TFs <- intersect(rownames(Columella), tf_all)
Cortex_TFs <- intersect(rownames(Cortex), tf_all)
Lateral_root_cap_TFs <- intersect(rownames(Lateral_root_cap), tf_all)
Metaphloem_sieve_element_TFs <- intersect(rownames(Metaphloem_sieve_element), tf_all)
Endodermis_TFs <- intersect(rownames(Endodermis), tf_all)
Xylem_TFs <- intersect(rownames(Xylem), tf_all)
Atrichoblast_TFs <- intersect(rownames(Atrichoblast), tf_all)
Trichoblast_TFs <- intersect(rownames(Trichoblast), tf_all)
Mature_TFs <- intersect(rownames(Mature), tf_all)


save(Columella_TFs,Cortex_TFs,Lateral_root_cap_TFs, Metaphloem_sieve_element_TFs,
     Endodermis_TFs,Xylem_TFs,Atrichoblast_TFs,Mature_TFs,Trichoblast_TFs,
     file = "Data/celltype_TFs.RData")

#########################二.Three algorithms are used to construct GRN#############################

######################### 2-1 WGCNA#############################

##=======================2.1.1-WGCNA_matrix_celltype========================
load("Data/celltype.RData")
load("Data/celltype_TFs.RData")

Atrichoblast_TFs <- data.frame(TF = Atrichoblast_TFs)
Columella_TFs <- data.frame(TF = Columella_TFs)
Endodermis_TFs <- data.frame(TF = Endodermis_TFs)
Cortex_TFs <- data.frame(TF = Cortex_TFs)
Mature_TFs <- data.frame(TF = Mature_TFs)
Metaphloem_sieve_element_TFs <- data.frame(TF = Metaphloem_sieve_element_TFs)
Lateral_root_cap_TFs <- data.frame(TF = Lateral_root_cap_TFs)
Trichoblast_TFs <- data.frame(TF = Trichoblast_TFs)
Xylem_TFs <- data.frame(TF = Xylem_TFs)

Atrichoblast_expr =data.frame(Atrichoblast)
Columella_expr =data.frame(Columella)
Cortex_expr =data.frame(Cortex)
Endodermis_expr =data.frame(Endodermis)
Lateral_root_cap_expr =data.frame(Lateral_root_cap)
Mature_expr =data.frame(Mature)
Metaphloem_sieve_element_expr =data.frame(Metaphloem_sieve_element)
Trichoblast_expr =data.frame(Trichoblast)
Xylem_expr =data.frame(Xylem)

Atrichoblast_TFs <- Atrichoblast_TFs$TF
Atrichoblast_TFs_expr <- Atrichoblast_expr[row.names(Atrichoblast_expr) %in% Atrichoblast_TFs, ]
Atrichoblast_matrix <- cor(t(Atrichoblast_TFs_expr), t(Atrichoblast_expr), use = "complete", method = 'pearson')

Columella_TFs <- Columella_TFs$TF
Columella_TFs_expr <- Columella_expr[row.names(Columella_expr) %in% Columella_TFs, ]
Columella_matrix <- cor(t(Columella_TFs_expr), t(Columella_expr), use = "complete", method = 'pearson')

Endodermis_TFs <- Endodermis_TFs$TF
Endodermis_TFs_expr <- Endodermis_expr[row.names(Endodermis_expr) %in% Endodermis_TFs, ]
Endodermis_matrix <- cor(t(Endodermis_TFs_expr), t(Endodermis_expr), use = "complete", method = 'pearson')

Cortex_TFs <- Cortex_TFs$TF
Cortex_TFs_expr <- Cortex_expr[row.names(Cortex_expr) %in% Cortex_TFs, ]
Cortex_matrix <- cor(t(Cortex_TFs_expr), t(Cortex_expr), use = "complete", method = 'pearson')

Lateral_root_cap_TFs <- Lateral_root_cap_TFs$TF
Lateral_root_cap_TFs_expr <- Lateral_root_cap_expr[row.names(Lateral_root_cap_expr) %in% Lateral_root_cap_TFs, ]
Lateral_root_cap_matrix <- cor(t(Lateral_root_cap_TFs_expr), t(Lateral_root_cap_expr), use = "complete", method = 'pearson')

Mature_TFs <- Mature_TFs$TF
Mature_TFs_expr <- Mature_expr[row.names(Mature_expr) %in% Mature_TFs, ]
Mature_matrix <- cor(t(Mature_TFs_expr), t(Mature_expr), use = "complete", method = 'pearson')

Metaphloem_sieve_element_TFs <- Metaphloem_sieve_element_TFs$TF
Metaphloem_sieve_element_TFs_expr <- Metaphloem_sieve_element_expr[row.names(Metaphloem_sieve_element_expr) %in% Metaphloem_sieve_element_TFs, ]
Metaphloem_sieve_element_matrix <- cor(t(Metaphloem_sieve_element_TFs_expr), t(Metaphloem_sieve_element_expr), use = "complete", method = 'pearson')

Trichoblast_TFs <- Trichoblast_TFs$TF
Trichoblast_TFs_expr <- Trichoblast_expr[row.names(Trichoblast_expr) %in% Trichoblast_TFs, ]
Trichoblast_matrix <- cor(t(Trichoblast_TFs_expr), t(Trichoblast_expr), use = "complete", method = 'pearson')

Xylem_TFs <- Xylem_TFs$TF
Xylem_TFs_expr <- Xylem_expr[row.names(Xylem_expr) %in% Xylem_TFs, ]
Xylem_matrix <- cor(t(Xylem_TFs_expr), t(Xylem_expr), use = "complete", method = 'pearson')

##################################################
##=======================2.1.2-WGCNA_biadjacency_matrix_celltype========================
library(WGCNA)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(dbplyr)
load("Data/celltype.RData")
Atrichoblast_expr =data.frame(Atrichoblast)
Columella_expr =data.frame(Columella)
Cortex_expr =data.frame(Cortex)
Endodermis_expr =data.frame(Endodermis)
Lateral_root_cap_expr =data.frame(Lateral_root_cap)
Mature_expr =data.frame(Mature)
Metaphloem_sieve_element_expr =data.frame(Metaphloem_sieve_element)
Trichoblast_expr =data.frame(Trichoblast)
Xylem_expr =data.frame(Xylem)

nSamples_Atrichoblast = ncol(Atrichoblast_expr);
nSamples_Columella = ncol(Columella_expr);
nSamples_Cortex = ncol(Cortex_expr);
nSamples_Endodermis = ncol(Endodermis_expr);
nSamples_Lateral_root_cap = ncol(Lateral_root_cap_expr);
nSamples_Mature = ncol(Mature_expr);
nSamples_Metaphloem_sieve_element = ncol(Metaphloem_sieve_element_expr);
nSamples_Trichoblast = ncol(Trichoblast_expr);
nSamples_Xylem = ncol(Xylem_expr);

#moduleTraitCor = cor(MEs, BLCA_group, use = "p")
Atrichoblast_TraitPvalue = corPvalueStudent(Atrichoblast_matrix, nSamples_Atrichoblast);
Columella_TraitPvalue = corPvalueStudent(Columella_matrix, nSamples_Columella);
Cortex_TraitPvalue = corPvalueStudent(Cortex_matrix, nSamples_Cortex);
Endodermis_TraitPvalue = corPvalueStudent(Endodermis_matrix, nSamples_Endodermis);
Lateral_root_cap_TraitPvalue = corPvalueStudent(Lateral_root_cap_matrix, nSamples_Lateral_root_cap);
Mature_TraitPvalue = corPvalueStudent(Mature_matrix, nSamples_Mature);
Metaphloem_sieve_element_TraitPvalue = corPvalueStudent(Metaphloem_sieve_element_matrix, nSamples_Metaphloem_sieve_element);
Trichoblast_TraitPvalue = corPvalueStudent(Trichoblast_matrix, nSamples_Trichoblast);
Xylem_TraitPvalue = corPvalueStudent(Xylem_matrix, nSamples_Xylem);

## Specify the significance level
significance_level <- 0.05

# Convert the P-value matrix to a 0-1 matrix based on the significance level
biadjacency_Atrichoblast <- ifelse(Atrichoblast_TraitPvalue < significance_level, 1, 0)
biadjacency_Columella <- ifelse(Columella_TraitPvalue < significance_level, 1, 0)
biadjacency_Cortex <- ifelse(Cortex_TraitPvalue < significance_level, 1, 0)
biadjacency_Endodermis <- ifelse(Endodermis_TraitPvalue < significance_level, 1, 0)
biadjacency_Lateral_root_cap <- ifelse(Lateral_root_cap_TraitPvalue < significance_level, 1, 0)
biadjacency_Mature <- ifelse(Mature_TraitPvalue < significance_level, 1, 0)
biadjacency_Metaphloem_sieve_element <- ifelse(Metaphloem_sieve_element_TraitPvalue < significance_level, 1, 0)
biadjacency_Trichoblast <- ifelse(Trichoblast_TraitPvalue < significance_level, 1, 0)
biadjacency_Xylem <- ifelse(Xylem_TraitPvalue < significance_level, 1, 0)


##=======================2.1.3-WGCNA_pvalue_igraph_egdes_graph========================
library(igraph)
library(dplyr)

graph_Atrichoblast <- graph_from_biadjacency_matrix(
  biadjacency_Atrichoblast,
  directed = TRUE,
  mode = "out",
  multiple = FALSE,
  weighted = TRUE,
  add.names = NULL
)
edges_Atrichoblast <- as.data.frame(get.edgelist(graph_Atrichoblast), stringsAsFactors = FALSE)

graph_Columella <- graph_from_biadjacency_matrix(
  biadjacency_Columella,
  directed = TRUE,
  mode = "out",
  multiple = FALSE,
  weighted = TRUE,
  add.names = NULL
)
edges_Columella <- as.data.frame(get.edgelist(graph_Columella), stringsAsFactors = FALSE)

graph_Cortex <- graph_from_biadjacency_matrix(
  biadjacency_Cortex,
  directed = TRUE,
  mode = "out",
  multiple = FALSE,
  weighted = TRUE,
  add.names = NULL
)
edges_Cortex <- as.data.frame(get.edgelist(graph_Cortex), stringsAsFactors = FALSE)

graph_Endodermis <- graph_from_biadjacency_matrix(
  biadjacency_Endodermis,
  directed = TRUE,
  mode = "out",
  multiple = FALSE,
  weighted = TRUE,
  add.names = NULL
)
edges_Endodermis <- as.data.frame(get.edgelist(graph_Endodermis), stringsAsFactors = FALSE)

graph_Lateral_root_cap <- graph_from_biadjacency_matrix(
  biadjacency_Lateral_root_cap,
  directed = TRUE,
  mode = "out",
  multiple = FALSE,
  weighted = TRUE,
  add.names = NULL
)
edges_Lateral_root_cap <- as.data.frame(get.edgelist(graph_Lateral_root_cap), stringsAsFactors = FALSE)

graph_Mature <- graph_from_biadjacency_matrix(
  biadjacency_Mature,
  directed = TRUE,
  mode = "out",
  multiple = FALSE,
  weighted = TRUE,
  add.names = NULL
)
edges_Mature <- as.data.frame(get.edgelist(graph_Mature), stringsAsFactors = FALSE)

graph_Metaphloem_sieve_element <- graph_from_biadjacency_matrix(
  biadjacency_Metaphloem_sieve_element,
  directed = TRUE,
  mode = "out",
  multiple = FALSE,
  weighted = TRUE,
  add.names = NULL
)
edges_Metaphloem_sieve_element <- as.data.frame(get.edgelist(graph_Metaphloem_sieve_element), stringsAsFactors = FALSE)

graph_Trichoblast <- graph_from_biadjacency_matrix(
  biadjacency_Trichoblast,
  directed = TRUE,
  mode = "out",
  multiple = FALSE,
  weighted = TRUE,
  add.names = NULL
)
edges_Trichoblast <- as.data.frame(get.edgelist(graph_Trichoblast), stringsAsFactors = FALSE)

graph_Xylem <- graph_from_biadjacency_matrix(
  biadjacency_Xylem,
  directed = TRUE,
  mode = "out",
  multiple = FALSE,
  weighted = TRUE,
  add.names = NULL
)
edges_Xylem <- as.data.frame(get.edgelist(graph_Xylem), stringsAsFactors = FALSE)


######################### 2-2 GENIE3#############################
load("Data/celltype.RData")
load("Data/celltype_TFs.RData")

library(GENIE3)
# Calculate weigh matrix.
# weigh matrix (wm) for Atrichoblast
ptm <- proc.time()
wm_Atrichoblast <- GENIE3(Atrichoblast,regulators = Atrichoblast_TFs,nCores = 4, K = "sqrt")
proc.time() - ptm

# weigh matrix (wm) for Columella
ptm <- proc.time()
wm_Columella <- GENIE3(Columella,regulators = Columella_TFs,nCores = 4, K = "sqrt")
proc.time() - ptm

# weigh matrix (wm) for Endodermis
ptm <- proc.time()
wm_Endodermis <- GENIE3(Endodermis,regulators = Endodermis_TFs,nCores = 4, K = "sqrt")
proc.time() - ptm 

# weigh matrix (wm) for Lateral_root_cap
ptm <- proc.time()
wm_Lateral_root_cap <- GENIE3(Lateral_root_cap,regulators = Lateral_root_cap_TFs, nCores = 4,K = "sqrt")
proc.time() - ptm 

# weigh matrix (wm) for Mature
ptm <- proc.time()
wm_Metaphloem_sieve_element <- GENIE3(Metaphloem_sieve_element,regulators = Metaphloem_sieve_element_TFs,nCores = 4, K = "sqrt")
proc.time() - ptm 

# weigh matrix (wm) for Trichoblast
ptm <- proc.time()
wm_Trichoblast <- GENIE3(Trichoblast,regulators = Trichoblast_TFs,nCores = 4, K = "sqrt")
proc.time() - ptm 

# weigh matrix (wm) for Cortex
ptm <- proc.time()
wm_Cortex <- GENIE3(Cortex,regulators = Cortex_TFs,nCores = 4,  K = "sqrt")
proc.time() - ptm 

# weigh matrix (wm) for Mature
ptm <- proc.time()
Mature <- GENIE3(Mature,regulators = tf_Mature_TFs,nCores = 4, K = "sqrt")
proc.time() - ptm 

# weigh matrix (wm) for Xylem
ptm <- proc.time()
wm_Xylem <- GENIE3(Xylem,regulators = Xylem_TFs,nCores = 4, K = "sqrt")
proc.time() - ptm 



## Function for calculating z-score of a matrix
# mat: Input matrix
# Output: mat.zscore is the transformed zscore matrix
transform_zscore_Atrichoblast <- function(GRN_Atrichoblast){
  
  GRN_Atrichoblast.mean <- mean(GRN_Atrichoblast[!is.na(GRN_Atrichoblast)])
  GRN_Atrichoblast.sd <- sd(GRN_Atrichoblast[!is.na(GRN_Atrichoblast)])
  GRN_Atrichoblast.zscore <- (GRN_Atrichoblast - GRN_Atrichoblast.mean)/GRN_Atrichoblast.sd
  return(GRN_Atrichoblast.zscore)
}

## Function for obtaining biadjacency matrix using zscore matrix
# mat.pvalue: zscore matrix
# pvalue.cutoff: p-value cutoff, default is 0.05
# Output: biadjacency matrix
biadjacency_matrix_Atrichoblast <- function(GRN_Atrichoblast.zscore, pvalue.cutoff = 0.05){
  
  zscore.cutoff <- -qnorm(pvalue.cutoff)
  GRN_Atrichoblast.biadjacency <- as.matrix(GRN_Atrichoblast.zscore > zscore.cutoff) * 1
  return(GRN_Atrichoblast.biadjacency)
}


wm_Atrichoblast_z <- wm_Atrichoblast
zscore_Atrichoblast <- transform_zscore_Atrichoblast(wm_Atrichoblast_z)
biadjacency_Atrichoblast <- biadjacency_matrix_Atrichoblast(zscore_Atrichoblast, pvalue.cutoff = 0.05)



## Function for calculating z-score of a matrix
# mat: Input matrix
# Output: mat.zscore is the transformed zscore matrix
transform_zscore_Columella <- function(GRN_Columella){
  
  GRN_Columella.mean <- mean(GRN_Columella[!is.na(GRN_Columella)])
  GRN_Columella.sd <- sd(GRN_Columella[!is.na(GRN_Columella)])
  GRN_Columella.zscore <- (GRN_Columella - GRN_Columella.mean)/GRN_Columella.sd
  return(GRN_Columella.zscore)
}

## Function for obtaining biadjacency matrix using zscore matrix
# mat.pvalue: zscore matrix
# pvalue.cutoff: p-value cutoff, default is 0.05
# Output: biadjacency matrix
biadjacency_matrix_Columella <- function(GRN_Columella.zscore, pvalue.cutoff = 0.05){
  
  zscore.cutoff <- -qnorm(pvalue.cutoff)
  GRN_Columella.biadjacency <- as.matrix(GRN_Columella.zscore > zscore.cutoff) * 1
  return(GRN_Columella.biadjacency)
}

wm_Columella_z <- wm_Columella
zscore_Columella <- transform_zscore_Columella(wm_Columella_z)
biadjacency_Columella <- biadjacency_matrix_Columella(zscore_Columella, pvalue.cutoff = 0.05)


## Function for calculating z-score of a matrix
# mat: Input matrix
# Output: mat.zscore is the transformed zscore matrix
transform_zscore_Cortex <- function(GRN_Cortex){
  
  GRN_Cortex.mean <- mean(GRN_Cortex[!is.na(GRN_Cortex)])
  GRN_Cortex.sd <- sd(GRN_Cortex[!is.na(GRN_Cortex)])
  GRN_Cortex.zscore <- (GRN_Cortex - GRN_Cortex.mean)/GRN_Cortex.sd
  return(GRN_Cortex.zscore)
}

## Function for obtaining biadjacency matrix using zscore matrix
# mat.pvalue: zscore matrix
# pvalue.cutoff: p-value cutoff, default is 0.05
# Output: biadjacency matrix
biadjacency_matrix_Cortex <- function(GRN_Cortex.zscore, pvalue.cutoff = 0.05){
  
  zscore.cutoff <- -qnorm(pvalue.cutoff)
  GRN_Cortex.biadjacency <- as.matrix(GRN_Cortex.zscore > zscore.cutoff) * 1
  return(GRN_Cortex.biadjacency)
}

wm_Cortex_z <- wm_Cortex
zscore_Cortex <- transform_zscore_Cortex(wm_Cortex_z)
biadjacency_Cortex <- biadjacency_matrix_Cortex(zscore_Cortex, pvalue.cutoff = 0.05)


## Function for calculating z-score of a matrix
# mat: Input matrix
# Output: mat.zscore is the transformed zscore matrix
transform_zscore_Endodermis <- function(GRN_Endodermis){
  
  GRN_Endodermis.mean <- mean(GRN_Endodermis[!is.na(GRN_Endodermis)])
  GRN_Endodermis.sd <- sd(GRN_Endodermis[!is.na(GRN_Endodermis)])
  GRN_Endodermis.zscore <- (GRN_Endodermis - GRN_Endodermis.mean)/GRN_Endodermis.sd
  return(GRN_Endodermis.zscore)
}

## Function for obtaining biadjacency matrix using zscore matrix
# mat.pvalue: zscore matrix
# pvalue.cutoff: p-value cutoff, default is 0.05
# Output: biadjacency matrix
biadjacency_matrix_Endodermis <- function(GRN_Endodermis.zscore, pvalue.cutoff = 0.05){
  
  zscore.cutoff <- -qnorm(pvalue.cutoff)
  GRN_Endodermis.biadjacency <- as.matrix(GRN_Endodermis.zscore > zscore.cutoff) * 1
  return(GRN_Endodermis.biadjacency)
}

wm_Endodermis_z <- wm_Endodermis
zscore_Endodermis <- transform_zscore_Endodermis(wm_Endodermis_z)
biadjacency_Endodermis <- biadjacency_matrix_Endodermis(zscore_Endodermis, pvalue.cutoff = 0.05)



## Function for calculating z-score of a matrix
# mat: Input matrix
# Output: mat.zscore is the transformed zscore matrix
transform_zscore_Lateral_root_cap <- function(GRN_Lateral_root_cap){
  
  GRN_Lateral_root_cap.mean <- mean(GRN_Lateral_root_cap[!is.na(GRN_Lateral_root_cap)])
  GRN_Lateral_root_cap.sd <- sd(GRN_Lateral_root_cap[!is.na(GRN_Lateral_root_cap)])
  GRN_Lateral_root_cap.zscore <- (GRN_Lateral_root_cap - GRN_Lateral_root_cap.mean)/GRN_Lateral_root_cap.sd
  return(GRN_Lateral_root_cap.zscore)
}

## Function for obtaining biadjacency matrix using zscore matrix
# mat.pvalue: zscore matrix
# pvalue.cutoff: p-value cutoff, default is 0.05
# Output: biadjacency matrix
biadjacency_matrix_Lateral_root_cap <- function(GRN_Lateral_root_cap.zscore, pvalue.cutoff = 0.05){
  
  zscore.cutoff <- -qnorm(pvalue.cutoff)
  GRN_Lateral_root_cap.biadjacency <- as.matrix(GRN_Lateral_root_cap.zscore > zscore.cutoff) * 1
  return(GRN_Lateral_root_cap.biadjacency)
}

wm_Lateral_root_cap_z <- wm_Lateral_root_cap
zscore_Lateral_root_cap <- transform_zscore_Lateral_root_cap(wm_Lateral_root_cap_z)
biadjacency_Lateral_root_cap <- biadjacency_matrix_Lateral_root_cap(zscore_Lateral_root_cap, pvalue.cutoff = 0.05)




## Function for calculating z-score of a matrix
# mat: Input matrix
# Output: mat.zscore is the transformed zscore matrix
transform_zscore_Metaphloem_sieve_element <- function(GRN_Metaphloem_sieve_element){
  
  GRN_Metaphloem_sieve_element.mean <- mean(GRN_Metaphloem_sieve_element[!is.na(GRN_Metaphloem_sieve_element)])
  GRN_Metaphloem_sieve_element.sd <- sd(GRN_Metaphloem_sieve_element[!is.na(GRN_Metaphloem_sieve_element)])
  GRN_Metaphloem_sieve_element.zscore <- (GRN_Metaphloem_sieve_element - GRN_Metaphloem_sieve_element.mean)/GRN_Metaphloem_sieve_element.sd
  return(GRN_Metaphloem_sieve_element.zscore)
}

## Function for obtaining biadjacency matrix using zscore matrix
# mat.pvalue: zscore matrix
# pvalue.cutoff: p-value cutoff, default is 0.05
# Output: biadjacency matrix
biadjacency_matrix_Metaphloem_sieve_element <- function(GRN_Metaphloem_sieve_element.zscore, pvalue.cutoff = 0.05){
  
  zscore.cutoff <- -qnorm(pvalue.cutoff)
  GRN_Metaphloem_sieve_element.biadjacency <- as.matrix(GRN_Metaphloem_sieve_element.zscore > zscore.cutoff) * 1
  return(GRN_Metaphloem_sieve_element.biadjacency)
}

wm_Metaphloem_sieve_element_z <- wm_Metaphloem_sieve_element
zscore_Metaphloem_sieve_element <- transform_zscore_Metaphloem_sieve_element(wm_Metaphloem_sieve_element_z)
biadjacency_Metaphloem_sieve_element <- biadjacency_matrix_Metaphloem_sieve_element(zscore_Metaphloem_sieve_element, pvalue.cutoff = 0.05)


## Function for calculating z-score of a matrix
# mat: Input matrix
# Output: mat.zscore is the transformed zscore matrix
transform_zscore_Trichoblast <- function(GRN_Trichoblast){
  
  GRN_Trichoblast.mean <- mean(GRN_Trichoblast[!is.na(GRN_Trichoblast)])
  GRN_Trichoblast.sd <- sd(GRN_Trichoblast[!is.na(GRN_Trichoblast)])
  GRN_Trichoblast.zscore <- (GRN_Trichoblast - GRN_Trichoblast.mean)/GRN_Trichoblast.sd
  return(GRN_Trichoblast.zscore)
}

## Function for obtaining biadjacency matrix using zscore matrix
# mat.pvalue: zscore matrix
# pvalue.cutoff: p-value cutoff, default is 0.05
# Output: biadjacency matrix
biadjacency_matrix_Trichoblast <- function(GRN_Trichoblast.zscore, pvalue.cutoff = 0.05){
  
  zscore.cutoff <- -qnorm(pvalue.cutoff)
  GRN_Trichoblast.biadjacency <- as.matrix(GRN_Trichoblast.zscore > zscore.cutoff) * 1
  return(GRN_Trichoblast.biadjacency)
}

wm_Trichoblast_z <- wm_Trichoblast
zscore_Trichoblast <- transform_zscore_Trichoblast(wm_Trichoblast_z)
biadjacency_Trichoblast <- biadjacency_matrix_Trichoblast(zscore_Trichoblast, pvalue.cutoff = 0.05)



## Function for calculating z-score of a matrix
# mat: Input matrix
# Output: mat.zscore is the transformed zscore matrix
transform_zscore_Xylem <- function(GRN_Xylem){
  
  GRN_Xylem.mean <- mean(GRN_Xylem[!is.na(GRN_Xylem)])
  GRN_Xylem.sd <- sd(GRN_Xylem[!is.na(GRN_Xylem)])
  GRN_Xylem.zscore <- (GRN_Xylem - GRN_Xylem.mean)/GRN_Xylem.sd
  return(GRN_Xylem.zscore)
}

## Function for obtaining biadjacency matrix using zscore matrix
# mat.pvalue: zscore matrix
# pvalue.cutoff: p-value cutoff, default is 0.05
# Output: biadjacency matrix
biadjacency_matrix_Xylem <- function(GRN_Xylem.zscore, pvalue.cutoff = 0.05){
  
  zscore.cutoff <- -qnorm(pvalue.cutoff)
  GRN_Xylem.biadjacency <- as.matrix(GRN_Xylem.zscore > zscore.cutoff) * 1
  return(GRN_Xylem.biadjacency)
}


wm_Xylem_z <- wm_Xylem
zscore_Xylem <- transform_zscore_Xylem(wm_Xylem_z)
biadjacency_Xylem <- biadjacency_matrix_Xylem(zscore_Xylem, pvalue.cutoff = 0.05)


## Function for calculating z-score of a matrix
# mat: Input matrix
# Output: mat.zscore is the transformed zscore matrix
transform_zscore_Mature <- function(GRN_Mature){
  
  GRN_Mature.mean <- mean(GRN_Mature[!is.na(GRN_Mature)])
  GRN_Mature.sd <- sd(GRN_Mature[!is.na(GRN_Mature)])
  GRN_Mature.zscore <- (GRN_Mature - GRN_Mature.mean)/GRN_Mature.sd
  return(GRN_Mature.zscore)
}

## Function for obtaining biadjacency matrix using zscore matrix
# mat.pvalue: zscore matrix
# pvalue.cutoff: p-value cutoff, default is 0.05
# Output: biadjacency matrix
biadjacency_matrix_Mature <- function(GRN_Mature.zscore, pvalue.cutoff = 0.05){
  
  zscore.cutoff <- -qnorm(pvalue.cutoff)
  GRN_Mature.biadjacency <- as.matrix(GRN_Mature.zscore > zscore.cutoff) * 1
  return(GRN_Mature.biadjacency)
}

wm_Mature_z <- Mature
zscore_Mature <- transform_zscore_Mature(wm_Mature_z)
biadjacency_Mature <- biadjacency_matrix_Xylem(zscore_Mature, pvalue.cutoff = 0.05)


load("biadjacency_nine_celltype.RData")

library(igraph)
library(dplyr)

graph_Atrichoblast <- graph_from_biadjacency_matrix(
  biadjacency_Atrichoblast,
  directed = TRUE,
  mode = "out",
  multiple = F,
  weighted = TRUE,
  add.names = NULL
)
edges_Atrichoblast <- as.data.frame(get.edgelist(graph_Atrichoblast), stringsAsFactors = FALSE)


graph_Columella <- graph_from_biadjacency_matrix(
  biadjacency_Columella,
  directed = TRUE,
  mode = "out",
  multiple = F,
  weighted = TRUE,
  add.names = NULL
)
edges_Columella <- as.data.frame(get.edgelist(graph_Columella), stringsAsFactors = FALSE)


graph_Cortex <- graph_from_biadjacency_matrix(
  biadjacency_Cortex,
  directed = TRUE,
  mode = "out",
  multiple = F,
  weighted = TRUE,
  add.names = NULL
)
edges_Cortex <- as.data.frame(get.edgelist(graph_Cortex), stringsAsFactors = FALSE)


graph_Endodermis <- graph_from_biadjacency_matrix(
  biadjacency_Endodermis,
  directed = TRUE,
  mode = "out",
  multiple = F,
  weighted = TRUE,
  add.names = NULL
)
edges_Endodermis <- as.data.frame(get.edgelist(graph_Endodermis), stringsAsFactors = FALSE)


graph_Lateral_root_cap <- graph_from_biadjacency_matrix(
  biadjacency_Lateral_root_cap,
  directed = TRUE,
  mode = "out",
  multiple = F,
  weighted = TRUE,
  add.names = NULL
)
edges_Lateral_root_cap <- as.data.frame(get.edgelist(graph_Lateral_root_cap), stringsAsFactors = FALSE)


graph_Mature <- graph_from_biadjacency_matrix(
  biadjacency_Mature,
  directed = TRUE,
  mode = "out",
  multiple = F,
  weighted = TRUE,
  add.names = NULL
)
edges_Mature <- as.data.frame(get.edgelist(graph_Mature), stringsAsFactors = FALSE)


graph_Metaphloem_sieve_element <- graph_from_biadjacency_matrix(
  biadjacency_Metaphloem_sieve_element,
  directed = TRUE,
  mode = "out",
  multiple = F,
  weighted = TRUE,
  add.names = NULL
)
edges_Metaphloem_sieve_element <- as.data.frame(get.edgelist(graph_Metaphloem_sieve_element), stringsAsFactors = FALSE)


graph_Trichoblast <- graph_from_biadjacency_matrix(
  biadjacency_Trichoblast,
  directed = TRUE,
  mode = "out",
  multiple = F,
  weighted = TRUE,
  add.names = NULL
)
edges_Trichoblast <- as.data.frame(get.edgelist(graph_Trichoblast), stringsAsFactors = FALSE)


graph_Xylem <- graph_from_biadjacency_matrix(
  biadjacency_Xylem,
  directed = TRUE,
  mode = "out",
  multiple = F,
  weighted = TRUE,
  add.names = NULL
)
edges_Xylem <- as.data.frame(get.edgelist(graph_Xylem), stringsAsFactors = FALSE)


######################### 2-3 DeepSEM#############################
##########
# The DeepSEM algorithm used in this script is based on the following open-source project：
# Hantao Shu et al., Modeling gene regulatory networks using neural network architectures 
# GitHub: https://github.com/HantaoShu/DeepSEM
# Reference tutorial: https://github.com/HantaoShu/DeepSEM/blob/master/tutorial/GRN_inference_tutorial.ipynb
# Note: This script relies on the DeepSEM environment. For the installation method, please refer to README.md

#########################三.Three algorithms are used to construct GRN evaluations#############################
#########################3.1 Evaluation of Network indicators#############################
library(pROC)
library(dplyr)
load("Data/celltype_edges_TF.RData")

gold_network <- read.csv("Data/gold_network.csv")

#######WGCNA
WGCNA_GRN_edges_predictions <- WGCNA_GRN_edges_unique$weight
WGCNA_GRN_edges_TF<- unique(WGCNA_GRN_edges_TF$V1)
Gene1 <- unique(gold_network$Gene1)

WGCNA_intersected <- intersect(WGCNA_GRN_edges_TF, Gene1)
WGCNA_GRN_edges_TF_COM <- data.frame(SharedValues = WGCNA_intersected)

WGCNA_GRN_edges_TF_COM <- unique(unlist(WGCNA_GRN_edges_TF_COM))
WGCNA_gold_network_label <- as.numeric(apply(gold_network, 1, function(x) any(x %in% WGCNA_GRN_edges_TF_COM)))

if (length(WGCNA_gold_network_label) != length(WGCNA_GRN_edges_predictions)) {
  if (length(WGCNA_gold_network_label) > length(WGCNA_GRN_edges_predictions)) {
    WGCNA_gold_network_label <- WGCNA_gold_network_label[1:length(WGCNA_GRN_edges_predictions)]
  } else {
    WGCNA_GRN_edges_predictions <- WGCNA_GRN_edges_predictions[1:length(WGCNA_gold_network_label)]
  }
}


#######GENIE3
GENIE3_GRN_edges_predictions <- GENIE3_GRN_edges_unique$weight
GENIE3_GRN_edges_TF<- unique(GENIE3_GRN_edges_TF$regulatoryGene)

Gene1 <- unique(gold_network$Gene1)
GENIE3_intersected <- intersect(GENIE3_GRN_edges_TF, Gene1)
GENIE3_GRN_edges_TF_COM <- data.frame(SharedValues = GENIE3_intersected)

GENIE3_GRN_edges_TF_COM <- unique(unlist(GENIE3_GRN_edges_TF_COM))
GENIE3_gold_network_label <- as.numeric(apply(gold_network, 1, function(x) any(x %in% GENIE3_GRN_edges_TF_COM)))

if (length(GENIE3_gold_network_label) != length(GENIE3_GRN_edges_predictions)) {
  if (length(GENIE3_gold_network_label) > length(GENIE3_GRN_edges_predictions)) {
    GENIE3_gold_network_label <- GENIE3_gold_network_label[1:length(GENIE3_GRN_edges_predictions)]
  } else {
    GENIE3_GRN_edges_predictions <- GENIE3_GRN_edges_predictions[1:length(GENIE3_gold_network_label)]
  }
}


#######DeepSEM
DeepSEM_GRN_edges_predictions <- DeepSEM_GRN_edges_unique$EdgeWeight
DeepSEM_GRN_edges_TF<- unique(DeepSEM_GRN_edges_TF$TF)

Gene1 <- unique(gold_network$Gene1)

DeepSEM_intersected <- intersect(DeepSEM_GRN_edges_TF, Gene1)
DeepSEM_GRN_edges_TF_COM <- data.frame(SharedValues = DeepSEM_intersected)

DeepSEM_GRN_edges_TF_COM <- unique(unlist(DeepSEM_GRN_edges_TF_COM))
DeepSEM_gold_network_label <- as.numeric(apply(gold_network, 1, function(x) any(x %in% DeepSEM_GRN_edges_TF_COM)))

if (length(DeepSEM_gold_network_label) != length(DeepSEM_GRN_edges_predictions)) {
  if (length(DeepSEM_gold_network_label) > length(DeepSEM_GRN_edges_predictions)) {
    DeepSEM_gold_network_label <- DeepSEM_gold_network_label[1:length(DeepSEM_GRN_edges_predictions)]
  } else {
    DeepSEM_GRN_edges_predictions <- DeepSEM_GRN_edges_predictions[1:length(DeepSEM_gold_network_label)]
  }
}



library(PRROC)
library(ggplot2)
WGCNA_PR <- pr.curve(scores.class0 = WGCNA_GRN_edges_predictions, weights.class0 = WGCNA_gold_network_label, curve = TRUE)
WGCNA_aupr_result <- WGCNA_PR$auc.integral

GENIE3_PR <- pr.curve(scores.class0 = GENIE3_GRN_edges_predictions, weights.class0 = GENIE3_gold_network_label, curve = TRUE)
GENIE3_aupr_result <- GENIE3_PR$auc.integral

DeepSEM_PR <- pr.curve(scores.class0 = DeepSEM_GRN_edges_predictions, weights.class0 = DeepSEM_gold_network_label, curve = TRUE)
DeepSEM_aupr_result <- DeepSEM_PR$auc.integral


#########################3.2 Network topology analysis#############################
library(igraph)
library(foreach)
library(doParallel)
load("celltype_edges_TF.RData")
#WGCNA
g_WGCNA <- make_graph(edges = NULL, directed = FALSE)
nodes_WGCNA <- unique(c(as.character(WGCNA_GRN_edges_unique$V1),as.character(WGCNA_GRN_edges_unique$V2)))
g_WGCNA <- add_vertices(g_WGCNA, nv = length(nodes_WGCNA), name = nodes_WGCNA)
edges_WGCNA <- as.matrix(WGCNA_GRN_edges_unique[, 1:2])
edges_WGCNA <- apply(edges_WGCNA, 1, as.character)
edge_weights_WGCNA <- WGCNA_GRN_edges_unique$weight
g_WGCNA <- add_edges(g_WGCNA, edges_WGCNA, weight = edge_weights_WGCNA)


#GENIE3
g_GENIE3 <- make_graph(edges = NULL, directed = FALSE)
nodes_GENIE3 <- unique(c(as.character(GENIE3_GRN_edges_unique$regulatoryGene),as.character(GENIE3_GRN_edges_unique$targetGene)))
g_GENIE3 <- add_vertices(g_GENIE3, nv = length(nodes_GENIE3), name = nodes_GENIE3)
edges_GENIE3 <- as.matrix(GENIE3_GRN_edges_unique[, 1:2])
edges_GENIE3 <- apply(edges_GENIE3, 1, as.character)
edge_weights_GENIE3 <- GENIE3_GRN_edges_unique$weight
g_GENIE3 <- add_edges(g_GENIE3, edges_GENIE3, weight = edge_weights_GENIE3)


#DeepSEM
g_DeepSEM <- make_graph(edges = NULL, directed = FALSE)
nodes_DeepSEM <- unique(c(as.character(DeepSEM_GRN_edges_unique$TF),as.character(DeepSEM_GRN_edges_unique$Target)))
g_DeepSEM <- add_vertices(g_DeepSEM, nv = length(nodes_DeepSEM), name = nodes_DeepSEM)
edges_DeepSEM <- as.matrix(DeepSEM_GRN_edges_unique[, 1:2])
edges_DeepSEM <- apply(edges_DeepSEM, 1, as.character)
edge_weights_DeepSEM <- DeepSEM_GRN_edges_unique$EdgeWeight
g_DeepSEM <- add_edges(g_DeepSEM, edges_DeepSEM, weight = edge_weights_DeepSEM)

#########random network
density_WGCNA <- numeric(100)
avg_path_length_WGCNA <- numeric(100)
clustering_coefficient_WGCNA <- numeric(100)

for (i in 1:100) {
  random_graph <- sample_gnp(n = vcount(g_WGCNA), p = 0.05)
  density_random <- edge_density(random_graph, loops = FALSE)
  density_WGCNA[i] <- density_random
  avg_path_length_random <- mean_distance(random_graph)
  avg_path_length_WGCNA[i] <- avg_path_length_random
  clustering_coefficient_random <- transitivity(random_graph, type = "global")
  clustering_coefficient_WGCNA[i] <- clustering_coefficient_random
}


density_GENIE3 <- numeric(100)
avg_path_length_GENIE3 <- numeric(100)
clustering_coefficient_GENIE3 <- numeric(100)

for (i in 1:100) {
  random_graph <- sample_gnp(n = vcount(g_GENIE3), p = 0.05)
  density_random <- edge_density(random_graph, loops = FALSE)
  density_GENIE3[i] <- density_random
  avg_path_length_random <- mean_distance(random_graph)
  avg_path_length_GENIE3[i] <- avg_path_length_random
  clustering_coefficient_random <- transitivity(random_graph, type = "global")
  clustering_coefficient_GENIE3[i] <- clustering_coefficient_random
}

density_DeepSEM <- numeric(100)
avg_path_length_DeepSEM <- numeric(100)
clustering_coefficient_DeepSEM <- numeric(100)

for (i in 1:100) {
  random_graph <- sample_gnp(n = vcount(g_DeepSEM), p = 0.05)
  density_random <- edge_density(random_graph, loops = FALSE)
  density_DeepSEM[i] <- density_random
  avg_path_length_random <- mean_distance(random_graph)
  avg_path_length_DeepSEM[i] <- avg_path_length_random
  clustering_coefficient_random <- transitivity(random_graph, type = "global")
  clustering_coefficient_DeepSEM[i] <- clustering_coefficient_random
}


load("celltype_edges_TF.RData")
library(igraph)
WGCNA_GRN <- graph.data.frame(WGCNA_GRN_edges_unique, directed = TRUE)
WGCNA_in_degrees <- degree(WGCNA_GRN, mode = "in")
WGCNA_out_degrees <- degree(WGCNA_GRN, mode = "out")
library(poweRlaw)
WGCNA_fit <- displ$new(WGCNA_in_degrees)
WGCNA_fit$setXmin(1)  
WGCNA_alpha <- WGCNA_fit$pars$alpha
WGCNA_in_degrees_unique <- unique(WGCNA_in_degrees)
WGCNA_in_degrees_unique <- WGCNA_in_degrees_unique[WGCNA_in_degrees_unique != 0]
WGCNA_ks_result <- ks.test(WGCNA_in_degrees_unique, "pexp", WGCNA_alpha)


library(igraph)
GENIE3_GRN <- graph.data.frame(GENIE3_GRN_edges_unique, directed = TRUE)
GENIE3_in_degrees <- degree(GENIE3_GRN, mode = "in")
GENIE3_out_degrees <- degree(GENIE3_GRN, mode = "out")
library(poweRlaw)
GENIE3_fit <- displ$new(GENIE3_in_degrees)
GENIE3_fit$setXmin(1) 
GENIE3_alpha <- GENIE3_fit$pars$alpha
GENIE3_in_degrees_unique <- unique(GENIE3_in_degrees)
GENIE3_in_degrees_unique <- GENIE3_in_degrees_unique[GENIE3_in_degrees_unique != 0]
GENIE3_ks_result <- ks.test(GENIE3_in_degrees_unique, "pexp", GENIE3_alpha)


library(igraph)
DeepSEM_GRN <- graph.data.frame(DeepSEM_GRN_edges_unique, directed = TRUE)
DeepSEM_in_degrees <- degree(DeepSEM_GRN, mode = "in")
DeepSEM_out_degrees <- degree(DeepSEM_GRN, mode = "out")
library(poweRlaw)
DeepSEM_in_degrees_shifted <- DeepSEM_in_degrees + 1
DeepSEM_fit <- displ$new(DeepSEM_in_degrees_shifted)
DeepSEM_fit$setXmin(1) 
DeepSEM_alpha <- DeepSEM_fit$pars$alpha
DeepSEM_in_degrees_unique <- unique(DeepSEM_in_degrees)
DeepSEM_in_degrees_unique <- DeepSEM_in_degrees_unique[DeepSEM_in_degrees_unique != 0]
DeepSEM_ks_result <- ks.test(DeepSEM_in_degrees_unique, "pexp", DeepSEM_alpha)

load("celltype_edges_TF.RData")
WGCNA_p_values <- rep(FALSE, length(avg_path_length_WGCNA))
WGCNA_t_test_path_length <- t.test(avg_path_length_WGCNA, rep(WGCNA_avg_path_length, length(avg_path_length_WGCNA)))$p.value
WGCNA_t_test_density <- t.test(density_WGCNA, rep(WGCNA_density_network, length(density_WGCNA)))$p.value
WGCNA_log_p_value_path_length <- ifelse(WGCNA_t_test_path_length == 0, -log10(.Machine$double.eps), -log10(WGCNA_t_test_path_length))
WGCNA_log_p_value_density <- ifelse(WGCNA_t_test_density == 0, -log10(.Machine$double.eps), -log10(WGCNA_t_test_density))
WGCNA_alpha <- 2.2e-16
WGCNA_is_small_world <- ifelse(WGCNA_t_test_path_length < WGCNA_alpha & WGCNA_t_test_density < WGCNA_alpha, TRUE, FALSE)
WGCNA_result <- data.frame(WGCNA_GRN_edges_unique, WGCNA_is_small_world = WGCNA_is_small_world, WGCNA_log_p_value_path_length, WGCNA_log_p_value_density)

GENIE3_p_values <- rep(FALSE, length(avg_path_length_GENIE3))
GENIE3_t_test_path_length <- t.test(avg_path_length_GENIE3, rep(GENIE3_avg_path_length, length(avg_path_length_GENIE3)))$p.value
GENIE3_t_test_density <- t.test(density_GENIE3, rep(GENIE3_density_network, length(density_GENIE3)))$p.value
GENIE3_log_p_value_path_length <- ifelse(GENIE3_t_test_path_length == 0, -log10(.Machine$double.eps), -log10(GENIE3_t_test_path_length))
GENIE3_log_p_value_density <- ifelse(GENIE3_t_test_density == 0, -log10(.Machine$double.eps), -log10(GENIE3_t_test_density))
GENIE3_alpha <- 2.2e-16
GENIE3_is_small_world <- ifelse(GENIE3_t_test_path_length < GENIE3_alpha & GENIE3_t_test_density < GENIE3_alpha, TRUE, FALSE)
GENIE3_result <- data.frame(GENIE3_GRN_edges_unique, GENIE3_is_small_world = GENIE3_is_small_world, GENIE3_log_p_value_path_length, GENIE3_log_p_value_density)


DeepSEM_p_values <- rep(FALSE, length(avg_path_length_DeepSEM))
DeepSEM_t_test_path_length <- t.test(avg_path_length_DeepSEM, rep(DeepSEM_avg_path_length, length(avg_path_length_DeepSEM)))$p.value
DeepSEM_t_test_density <- t.test(density_DeepSEM, rep(DeepSEM_density_network, length(density_DeepSEM)))$p.value
DeepSEM_log_p_value_path_length <- ifelse(DeepSEM_t_test_path_length == 0, -log10(.Machine$double.eps), -log10(DeepSEM_t_test_path_length))
DeepSEM_log_p_value_density <- ifelse(DeepSEM_t_test_density == 0, -log10(.Machine$double.eps), -log10(DeepSEM_t_test_density))
DeepSEM_alpha <- 2.2e-16
DeepSEM_is_small_world <- ifelse(DeepSEM_t_test_path_length < DeepSEM_alpha & DeepSEM_t_test_density < DeepSEM_alpha, TRUE, FALSE)
DeepSEM_result <- data.frame(DeepSEM_GRN_edges_unique, DeepSEM_is_small_world = DeepSEM_is_small_world, DeepSEM_log_p_value_path_length, DeepSEM_log_p_value_density)


#########################四.Downstream analysis（GENIE3)#############################
######Network structure analysis reveals the complexity of gene expression regulation
load("Data/GENIE3_GRN_edges.RData")
library(igraph)
library(tidyverse)
library(gridExtra)
library(pheatmap)
library(RColorBrewer)
################################################################################
# calclate degree for four tissues

cal_deg <- function(link_list){
  link_list %>% 
    select(-weight) %>% 
    group_by(regulatoryGene) %>% 
    summarise(n = n()) %>% 
    arrange(desc(n))
}

deg_Atrichoblast <- cal_deg(Atrichoblast_edges)
deg_Columella <- cal_deg(Columella_edges)
deg_Cortex <- cal_deg(Cortex_edges)
deg_Endodermis <- cal_deg(Endodermis_edges)
deg_Lateral_root_cap <- cal_deg(Lateral_root_cap_edges)
deg_Mature <- cal_deg(Mature_edges)
deg_Metaphloem_sieve_element <- cal_deg(Metaphloem_sieve_element_edges)
deg_Trichoblast <- cal_deg(Trichoblast_edges)
deg_Xylem <- cal_deg(Xylem_edges)

# I choose TF connected with more than 2000 genes as key TFs.
key_n_Atrichoblast <- dim(filter(deg_Atrichoblast, n > 2000))[1] # 30
key_n_Columella <- dim(filter(deg_Columella, n > 2000))[1] # 38
key_n_Cortex <- dim(filter(deg_Cortex, n > 2000))[1] # 24
key_n_Endodermis <- dim(filter(deg_Endodermis, n > 2000))[1] # 41
key_n_Lateral_root_cap <- dim(filter(deg_Lateral_root_cap, n > 2000))[1] # 47
key_n_Mature <- dim(filter(deg_Mature, n > 2000))[1] # 7
key_n_Metaphloem_sieve_element <- dim(filter(deg_Metaphloem_sieve_element, n > 2000))[1] # 41
key_n_Trichoblast <- dim(filter(deg_Trichoblast, n > 2000))[1] # 4
key_n_Xylem <- dim(filter(deg_Xylem, n > 2000))[1] # 17

# plot barplot at nine celltype
p_Atrichoblast <- ggplot(data = deg_Atrichoblast) +
  geom_col(mapping = aes(x = reorder(regulatoryGene, -n), y = n)) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  labs(x = "TFs in Atrichoblast", y = "Number of targets", 
       title = "TF Degree Centrality in Atrichoblast") +
  geom_segment(aes(x = key_n_Atrichoblast, xend = key_n_Atrichoblast, y = 0, yend = 5200),
               color = "red", linetype = "dashed",size = 0.3)

p_Columella <- ggplot(data = deg_Columella) +
  geom_col(mapping = aes(x = reorder(regulatoryGene, -n), y = n)) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  labs(x = "TFs in Columella", y = "Number of targets", 
       title = "TF Degree Centrality in Columella") +
  geom_segment(aes(x = key_n_Columella, xend = key_n_Columella, y = 0, yend = 3000),
               color = "red", linetype = "dashed",size = 0.3)

p_Cortex <- ggplot(data = deg_Cortex) +
  geom_col(mapping = aes(x = reorder(regulatoryGene, -n), y = n)) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  labs(x = "TFs in Cortex", y = "Number of targets",
       title = "TF Degree Centrality in Cortex") +
  geom_segment(aes(x = key_n_Cortex, xend = key_n_Cortex, y = 0, yend = 5000),
               color = "red", linetype = "dashed",size = 0.3)

p_Endodermis <- ggplot(data = deg_Endodermis) +
  geom_col(mapping = aes(x = reorder(regulatoryGene, -n), y = n)) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  labs(x = "TFs in Endodermis", y = " Number of targets", 
       title = "TF Degree Centrality in Endodermis") +
  geom_segment(aes(x = key_n_Endodermis, xend = key_n_Endodermis, y = 0, yend = 4000),
               color = "red", linetype = "dashed",size = 0.3)


p_Lateral_root_cap <- ggplot(data = deg_Lateral_root_cap) +
  geom_col(mapping = aes(x = reorder(regulatoryGene, -n), y = n)) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  labs(x = "TFs in Lateral_root_cap", y = "Number of targets", 
       title = "TF Degree Centrality in Lateral_root_cap") +
  geom_segment(aes(x = key_n_Lateral_root_cap, xend = key_n_Lateral_root_cap, y = 0, yend = 5200),
               color = "red", linetype = "dashed",size = 0.3)

p_Mature <- ggplot(data = deg_Mature) +
  geom_col(mapping = aes(x = reorder(regulatoryGene, -n), y = n)) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  labs(x = "TFs in Mature", y = "Number of targets", 
       title = "TF Degree Centrality in Mature") +
  geom_segment(aes(x = key_n_Mature, xend = key_n_Mature, y = 0, yend = 3000),
               color = "red", linetype = "dashed",size = 0.3)

p_Metaphloem_sieve_element <- ggplot(data = deg_Metaphloem_sieve_element) +
  geom_col(mapping = aes(x = reorder(regulatoryGene, -n), y = n)) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  labs(x = "TFs in Metaphloem_sieve_element", y = "Number of targets",
       title = "TF Degree Centrality in Metaphloem_sieve_element") +
  geom_segment(aes(x = key_n_Metaphloem_sieve_element, xend = key_n_Metaphloem_sieve_element, y = 0, yend = 5000),
               color = "red", linetype = "dashed",size = 0.3)

p_Trichoblast <- ggplot(data = deg_Trichoblast) +
  geom_col(mapping = aes(x = reorder(regulatoryGene, -n), y = n)) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  labs(x = "TFs in Trichoblast", y = " Number of targets", 
       title = "TF Degree Centrality in Trichoblast") +
  geom_segment(aes(x = key_n_Trichoblast, xend = key_n_Trichoblast, y = 0, yend = 4000),
               color = "red", linetype = "dashed",size = 0.3)

p_Xylem <- ggplot(data = deg_Xylem) +
  geom_col(mapping = aes(x = reorder(regulatoryGene, -n), y = n)) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  labs(x = "TFs in Xylem", y = " Number of targets", 
       title = "TF Degree Centrality in Xylem") +
  geom_segment(aes(x = key_n_Xylem, xend = key_n_Xylem, y = 0, yend = 4000),
               color = "red", linetype = "dashed",size = 0.3)



# Put them together
p_total <- grid.arrange(p_Atrichoblast, p_Columella, p_Cortex, p_Endodermis,p_Lateral_root_cap,p_Mature,
                        p_Metaphloem_sieve_element,p_Trichoblast,p_Xylem, ncol = 3)

################################################################################
# PLOT QUANTILE HEATMAP
################################################################################

# Because the highest number for TF connection is about 4000-5000, but most of 
# connections are very low. So if I used even color break, most regions will be 
# blue and I do not want to use `scale row`.
# So I tried a script found online using quantile break. Works great.

# define a function to calculate quantile breaks.
quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

# Prepare the count number matrix for four tissues.
deg_nine <- full_join(deg_Atrichoblast, deg_Columella, by = "regulatoryGene") %>%
  full_join(., deg_Cortex, by = "regulatoryGene") %>%
  full_join(., deg_Endodermis, by = "regulatoryGene") %>%
  full_join(., deg_Lateral_root_cap, by = "regulatoryGene") %>%
  full_join(., deg_Mature, by = "regulatoryGene") %>%
  full_join(., deg_Metaphloem_sieve_element, by = "regulatoryGene") %>%
  full_join(., deg_Trichoblast, by = "regulatoryGene") %>%
  full_join(., deg_Xylem, by = "regulatoryGene") %>%
  rename(
    Atrichoblast = n.x,
    Columella = n.y,
    Cortex = n.x.x,
    Endodermis = n.y.y,
    Lateral_root_cap = n.x.x.x,
    Mature = n.y.y.y,
    Metaphloem_sieve_element = n.x.x.x.x,
    Trichoblast = n.y.y.y.y,
    Xylem = n
  ) %>%
  replace_na(.,
             replace = list(
               Atrichoblast = 0,
               Columella = 0,
               Cortex = 0,
               Endodermis = 0,
               Lateral_root_cap = 0,
               Mature = 0,
               Metaphloem_sieve_element = 0,
               Trichoblast = 0,
               Xylem = 0
             )
  )


library(pheatmap)
library(RColorBrewer)
mat <- as.matrix(deg_nine[, -1])
mat <- mat[rowSums(mat) > 0, ]
quantile_breaks <- function(x, n = 10) {
  breaks <- quantile(x, probs = seq(0, 1, length.out = n), na.rm = TRUE)
  breaks <- unique(breaks)  
  if(length(breaks) < 2) {  
    breaks <- seq(min(x), max(x), length.out = n)
  }
  return(breaks)
}
mat_breaks <- quantile_breaks(mat, n = 11)

pheatmap(
  mat = mat,
  color = colorRampPalette(rev(brewer.pal(n = 10, name = "RdYlBu")))(length(mat_breaks)-1),
  breaks = mat_breaks,
  border_color = NA,
  show_colnames = TRUE,
  show_rownames = FALSE,
  clustering_method = "complete", 
  clustering_distance_rows = "euclidean",
  drop_levels = TRUE,
  fontsize = 14,
  fontsize_col = 20,
  main = "TFs Degree Centrality Heatmap"
)

#############CV
library(tidyverse)
library(gridExtra)
library(broom)
library(stringr)
# combine four link_list into one.
ll_Atrichoblast <- mutate(Atrichoblast_edges, celltype = "Atrichoblast,")
ll_Columella <- mutate(Columella_edges, celltype = "Columella,")
ll_Cortex <- mutate(Cortex_edges, celltype = "Cortex,")
ll_Endodermis <- mutate(Endodermis_edges, celltype = "Endodermis,")
ll_Lateral_root_cap <- mutate(Lateral_root_cap_edges, celltype = "Lateral_root_cap,")
ll_Mature <- mutate(Mature_edges, celltype = "Mature,")
ll_Metaphloem_sieve_element <- mutate(Metaphloem_sieve_element_edges, celltype = "Metaphloem_sieve_element,")
ll_Trichoblast <- mutate(Trichoblast_edges, celltype = "Trichoblast,")
ll_Xylem <- mutate(Xylem_edges, celltype = "Xylem,")

ll_nine <- bind_rows(ll_Atrichoblast, ll_Columella, ll_Cortex, ll_Endodermis,ll_Lateral_root_cap,
                     ll_Mature,ll_Metaphloem_sieve_element,ll_Trichoblast,ll_Xylem )

# make sure the number is correct
ll_nine %>% 
  group_by(celltype) %>% 
  summarise(n = n())

# Find TFs that inlcude in all nine celltype. 672 TFs have at least one
# interactions in one tissue.
Atrichoblast <- ll_nine %>%
  filter(celltype == "Atrichoblast,") %>% 
  distinct(regulatoryGene)

Columella <- ll_nine %>%
  filter(celltype == "Columella,") %>% 
  distinct(regulatoryGene)

Cortex <- ll_nine %>%
  filter(celltype == "Cortex,") %>% 
  distinct(regulatoryGene)

Endodermis <- ll_nine %>%
  filter(celltype == "Endodermis,") %>% 
  distinct(regulatoryGene)

Lateral_root_cap <- ll_nine %>%
  filter(celltype == "Lateral_root_cap,") %>% 
  distinct(regulatoryGene)

Mature <- ll_nine %>%
  filter(celltype == "Mature,") %>% 
  distinct(regulatoryGene)

Metaphloem_sieve_element <- ll_nine %>%
  filter(celltype == "Metaphloem_sieve_element,") %>% 
  distinct(regulatoryGene)

Trichoblast <- ll_nine %>%
  filter(celltype == "Trichoblast,") %>% 
  distinct(regulatoryGene)

Xylem <- ll_nine %>%
  filter(celltype == "Xylem,") %>% 
  distinct(regulatoryGene)

nine_tf <- Reduce(f = intersect, 
                  list(Atrichoblast$regulatoryGene, Columella$regulatoryGene,Cortex$regulatoryGene,Endodermis$regulatoryGene, Lateral_root_cap$regulatoryGene,Mature$regulatoryGene,
                       Metaphloem_sieve_element$regulatoryGene, Trichoblast$regulatoryGene,Xylem$regulatoryGene) )

# Find TFs with number of connections in each tissue
count_Atrichoblast <- ll_nine %>% 
  filter(celltype == "Atrichoblast,") %>% 
  count(regulatoryGene) 

count_Columella <- ll_nine %>% 
  filter(celltype == "Columella,") %>% 
  count(regulatoryGene) 

count_Endodermis <- ll_nine %>% 
  filter(celltype == "Endodermis,") %>% 
  count(regulatoryGene) 

count_Cortex <- ll_nine %>% 
  filter(celltype == "Cortex,") %>% 
  count(regulatoryGene) 

count_Lateral_root_cap <- ll_nine %>% 
  filter(celltype == "Lateral_root_cap,") %>% 
  count(regulatoryGene) 

count_Mature <- ll_nine %>% 
  filter(celltype == "Mature,") %>% 
  count(regulatoryGene) 

count_Metaphloem_sieve_element <- ll_nine %>% 
  filter(celltype == "Metaphloem_sieve_element,") %>% 
  count(regulatoryGene) 

count_Trichoblast <- ll_nine %>% 
  filter(celltype == "Trichoblast,") %>% 
  count(regulatoryGene) 

count_Xylem <- ll_nine %>% 
  filter(celltype == "Xylem,") %>% 
  count(regulatoryGene) 

count_nine <- inner_join(count_Atrichoblast, count_Columella, by = "regulatoryGene", suffix = c("Atrichoblast", "Columella")) %>% 
  inner_join(count_Endodermis, by = "regulatoryGene") %>% 
  inner_join(count_Cortex, by = "regulatoryGene", suffix = c("Endodermis", "Cortex")) %>%
  inner_join(count_Lateral_root_cap, by = "regulatoryGene") %>% 
  inner_join(count_Mature, by = "regulatoryGene", suffix = c("Lateral_root_cap", "Mature")) %>%
  inner_join(count_Trichoblast, by = "regulatoryGene") %>%
  inner_join(count_Xylem, by = "regulatoryGene", suffix = c("Trichoblast", "Xylem")) %>%
  inner_join(count_Metaphloem_sieve_element, by = "regulatoryGene", suffix = c("Metaphloem_sieve_element", "Metaphloem_sieve_element"))

################################################################################
# Save a table
# Table contains TF expression, centrality, CV, 
# the tissue with higest centrality and the difference between 
# highest and lowest centrality. 
################################################################################
# calculate diff in TF interactions and write the table.
# diff in degreee centrality

count_diff <- count_nine %>% 
  mutate(ndiff = pmax(nAtrichoblast, nColumella, nEndodermis, nCortex,nLateral_root_cap, nMature, nTrichoblast, nXylem,nMetaphloem_sieve_element) - 
           pmin(nAtrichoblast, nColumella, nEndodermis, nCortex,nLateral_root_cap, nMature, nTrichoblast, nXylem,nMetaphloem_sieve_element)) %>% 
  arrange(desc(ndiff))

# find the tissue with max interactions.
max_tissue <- as_tibble(colnames(count_diff[,c(2:10)])[max.col(count_diff[,c(2:10)],
                                                               ties.method = "first")]) %>% 
  rename(max_tissue = value)

# delete "n" in the tissue name and change to factors.
count_diff <- bind_cols(count_diff, max_tissue) %>% 
  mutate(max_tissue = str_replace(max_tissue, "^n", "")) %>% 
  mutate(max_tissue = as.factor(max_tissue))

# add CV
CV <- function(x){
  (sd(x)/mean(x))*100
}

cv_table <- as_tibble(apply(count_diff[,c(2:10)], 1, CV)) %>% 
  rename(cv = value)


count_diff <- count_diff %>% 
  bind_cols(., cv_table) %>% 
  arrange(desc(cv)) %>% 
  filter(., ndiff > 500)


summary(count_diff[1:20,]$max_tissue)
summary(count_diff[1:50,]$max_tissue)
summary(count_diff[1:100,]$max_tissue)


# Sort in descending order by cv and select the Top100
top100 <- count_diff %>%
  arrange(desc(cv)) %>%
  head(100)

summary(top100[1:100,]$max_tissue)

tissue_counts <- top100 %>%
  count(max_tissue) %>%
  mutate(max_tissue = gsub("^n", "", max_tissue))  # 去除组织名前缀的"n"


ggplot(tissue_counts, aes(x = reorder(max_tissue, -n), y = n, fill = max_tissue)) +
  geom_bar(stat = "identity", width = 0.4) +  
  geom_text(aes(label = n), vjust = -0.5, size = 3.5) +  
  labs(
    title = "Top 100 TFs with Highest CV: Cell Type Distribution",
    x = NULL,  # 隐藏x轴标题
    y = "TFs Count"
  ) +
  theme_minimal(base_size = 12) +  # 调整基础字体大小
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, size = 10), 
    plot.title = element_text(hjust = 0.5, face = "bold"),  
    panel.grid.major.x = element_blank(),  
    legend.position = "none"
  ) +
  scale_fill_brewer(palette = "Set1") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) 


plot <- ggplot(tissue_counts, aes(x = reorder(max_tissue, -n), y = n, fill = max_tissue)) +
  geom_bar(stat = "identity", width = 0.4) +
  geom_text(aes(label = n), vjust = -0.5, size = 3.5) +
  labs(
    title = "Top 100 TFs with Highest CV: Cell Type Distribution",
    x = NULL,
    y = "TFs Count"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, size = 10),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid.major.x = element_blank(),
    legend.position = "none"
  ) +
  scale_fill_brewer(palette = "Set1") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

ggsave("TF_CellType_Distribution.pdf", 
       plot = plot,
       width = 8,   
       height = 6,  
       device = "pdf")

##########Correlation analysis of key transcription factors
load("Data/D_key_n_nine_celltype_TFs.RData")

data_frames <- list(D_key_n_Atrichoblast$regulatoryGene, D_key_n_Columella$regulatoryGene, D_key_n_Cortex$regulatoryGene, D_key_n_Endodermis$regulatoryGene, 
                    D_key_n_Lateral_root_cap$regulatoryGene, D_key_n_Mature$regulatoryGene, D_key_n_Metaphloem_sieve_element$regulatoryGene, 
                    D_key_n_Trichoblast$regulatoryGene, D_key_n_Xylem$regulatoryGene)

Sim.hub <- function(hub1, hub2){
  
  if(class(hub1)!="list" | class(hub2)!="list") {
    stop("Please check your input hub! The input hub should be list object! \n")
  }
  
  m <- length(hub1)
  n <- length(hub2)
  Sim <- matrix(NA, m, n)
  for (i in seq(m)){
    for (j in seq(n)){	    
      overlap_interin <- length(intersect(hub1[[i]], hub2[[j]]))
      Sim[i, j] <- overlap_interin/min(length(hub1[[i]]), length(hub2[[j]]))
    }
  }
  
  return(Sim)
}
similarity_matrix <- Sim.hub(data_frames, data_frames)
print(similarity_matrix)


key_TF_names <- c("Atrichoblast", "Columella", "Cortex", "Endodermis", "Metaphloem sieve element", "Trichoblast", "Lateral root cap", "Mature", "Xylem")
rownames(similarity_matrix) <- colnames(similarity_matrix) <- key_TF_names

Redice_res_Sim_1 <- similarity_matrix
Redice_res_Sim_2 <- similarity_matrix
diag(Redice_res_Sim_1) <- 0
diag(Redice_res_Sim_2) <- 1

## Hierarchical cluster analysis of cell-specific miRNA-mRNA causal regulatory network
hclust_res_network <- hclust(as.dist(1- similarity_matrix), "complete") 
clust_network_membership <- cutree(hclust_res_network, 3)

## Heatmap of network similarity matrix
library(corrplot)
corrplot.mixed(similarity_matrix, lower = "number", upper = "pie", 
               lower.col = "black", upper.col = "black", tl.cex = 1.0, 
               tl.pos = "lt", addCoefAsPercent = TRUE)


########Module identification and analysis
library(miRspongeR)
load("Data/GENIE3_GRN_edges.RData")
# Module identification from BRCA-related miRNA sponge interaction network identified by integrative method
GENIE3_Atrichoblast_Cluster_MCL <- netModule(Atrichoblast_edges, method = "MCL",directed = TRUE,modulesize = 10)
GENIE3_Columella_Cluster_MCL <- netModule(Columella_edges, method = "MCL",directed = TRUE,modulesize = 10)
GENIE3_Cortex_Cluster_MCL <- netModule(Cortex_edges, method = "MCL",directed = TRUE,modulesize = 10)
GENIE3_Endodermis_Cluster_MCL <- netModule(Endodermis_edges, method = "MCL",directed = TRUE,modulesize = 10)
GENIE3_Metaphloem_sieve_element_Cluster_MCL <- netModule(Metaphloem_sieve_element_edges, method = "MCL",directed = TRUE,modulesize = 10)
GENIE3_Trichoblast_Cluster_MCL <- netModule(Trichoblast_edges, method = "MCL",directed = TRUE,modulesize = 10)
GENIE3_Lateral_root_cap_Cluster_MCL <- netModule(Lateral_root_cap_edges, method = "MCL",directed = TRUE,modulesize = 10)
GENIE3_Mature_Cluster_MCL <- netModule(Mature_edges, method = "MCL",directed = TRUE,modulesize = 10)
GENIE3_Xylem_Cluster_MCL <- netModule(Xylem_edges, method = "MCL",directed = TRUE,modulesize = 10)

###########GO
library(org.At.tair.db)
library(clusterProfiler)
library(AnnotationDbi)
library(ggplot2)
library(aPEAR)
library(DOSE)
library(cols4all)

GO_BP_Atrichoblast <- list()
for (i in 1:21) {
  current_module <- GENIE3_Atrichoblast_Cluster_MCL[[i]]
  current_module_df <- data.frame(Gene = current_module)
  colnames(current_module_df) <- c("TAIR")
  rownames(current_module_df) <- current_module_df$TAIR
  ids_current_module <- bitr(rownames(current_module_df), 
                             fromType = "TAIR", 
                             toType = c("ENTREZID"), 
                             OrgDb = org.At.tair.db, 
                             drop = TRUE)
  
  
  GO_current_module <- enrichGO(ids_current_module$ENTREZID,
                                OrgDb = org.At.tair.db,
                                keyType = "ENTREZID",
                                ont = "BP", 
                                pvalueCutoff = 0.05,
                                qvalueCutoff = 0.05,
                                readable = TRUE)
  
  
  GO_BP_Atrichoblast[[i]] <- data.frame(GO_current_module)
}

library(dplyr)
for (i in 1:length(GO_BP_Atrichoblast)) {
  if (nrow(GO_BP_Atrichoblast[[i]]) > 0) { 
    GO_BP_Atrichoblast[[i]]$Module <- paste0("Module_", i) 
  }
}

GO_BP_Columella <- list()

for (i in 1:11) {
  current_module <- GENIE3_Columella_Cluster_MCL[[i]]
  current_module_df <- data.frame(Gene = current_module)
  colnames(current_module_df) <- c("TAIR")
  rownames(current_module_df) <- current_module_df$TAIR
  ids_current_module <- bitr(rownames(current_module_df), 
                             fromType = "TAIR", 
                             toType = c("ENTREZID"), 
                             OrgDb = org.At.tair.db, 
                             drop = TRUE)
  GO_current_module <- enrichGO(ids_current_module$ENTREZID,
                                OrgDb = org.At.tair.db,
                                keyType = "ENTREZID",
                                ont = "BP", 
                                pvalueCutoff = 0.05,
                                qvalueCutoff = 0.05,
                                readable = TRUE)
  GO_BP_Columella[[i]] <- data.frame(GO_current_module)
}

library(dplyr)

for (i in 1:length(GO_BP_Columella)) {
  if (nrow(GO_BP_Columella[[i]]) > 0) { 
    GO_BP_Columella[[i]]$Module <- paste0("Module_", i)  
  }
}

GO_BP_Cortex <- list()

for (i in 1:68) {
  current_module <- GENIE3_Cortex_Cluster_MCL[[i]]
  current_module_df <- data.frame(Gene = current_module)
  colnames(current_module_df) <- c("TAIR")
  rownames(current_module_df) <- current_module_df$TAIR
  ids_current_module <- bitr(rownames(current_module_df), 
                             fromType = "TAIR", 
                             toType = c("ENTREZID"), 
                             OrgDb = org.At.tair.db, 
                             drop = TRUE)
  GO_current_module <- enrichGO(ids_current_module$ENTREZID,
                                OrgDb = org.At.tair.db,
                                keyType = "ENTREZID",
                                ont = "BP",  
                                pvalueCutoff = 0.05,
                                qvalueCutoff = 0.05,
                                readable = TRUE)
  GO_BP_Cortex[[i]] <- data.frame(GO_current_module)
}

library(dplyr)
for (i in 1:length(GO_BP_Cortex)) {
  if (nrow(GO_BP_Cortex[[i]]) > 0) { 
    GO_BP_Cortex[[i]]$Module <- paste0("Module_", i)  
  }
}

GO_BP_Endodermis <- list()

for (i in 1:72) {
  current_module <- GENIE3_Endodermis_Cluster_MCL[[i]]
  current_module_df <- data.frame(Gene = current_module)
  colnames(current_module_df) <- c("TAIR")
  rownames(current_module_df) <- current_module_df$TAIR
  ids_current_module <- bitr(rownames(current_module_df), 
                             fromType = "TAIR", 
                             toType = c("ENTREZID"), 
                             OrgDb = org.At.tair.db, 
                             drop = TRUE)
  GO_current_module <- enrichGO(ids_current_module$ENTREZID,
                                OrgDb = org.At.tair.db,
                                keyType = "ENTREZID",
                                ont = "BP",  
                                pvalueCutoff = 0.05,
                                qvalueCutoff = 0.05,
                                readable = TRUE)
  GO_BP_Endodermis[[i]] <- data.frame(GO_current_module)
}

library(dplyr)
for (i in 1:length(GO_BP_Endodermis)) {
  if (nrow(GO_BP_Endodermis[[i]]) > 0) { 
    GO_BP_Endodermis[[i]]$Module <- paste0("Module_", i) 
  }
}

GO_BP_Lateral_root_cap <- list()

for (i in 1:8) {
  current_module <- GENIE3_Lateral_root_cap_Cluster_MCL[[i]]
  current_module_df <- data.frame(Gene = current_module)
  colnames(current_module_df) <- c("TAIR")
  rownames(current_module_df) <- current_module_df$TAIR
  ids_current_module <- bitr(rownames(current_module_df), 
                             fromType = "TAIR", 
                             toType = c("ENTREZID"), 
                             OrgDb = org.At.tair.db, 
                             drop = TRUE)
  GO_current_module <- enrichGO(ids_current_module$ENTREZID,
                                OrgDb = org.At.tair.db,
                                keyType = "ENTREZID",
                                ont = "BP",  
                                pvalueCutoff = 0.05,
                                qvalueCutoff = 0.05,
                                readable = TRUE)
  GO_BP_Lateral_root_cap[[i]] <- data.frame(GO_current_module)
}

library(dplyr)
for (i in 1:length(GO_BP_Lateral_root_cap)) {
  if (nrow(GO_BP_Lateral_root_cap[[i]]) > 0) { 
    GO_BP_Lateral_root_cap[[i]]$Module <- paste0("Module_", i) 
  }
}

GO_BP_Metaphloem_sieve_element <- list()

for (i in 1:2) {
  current_module <- GENIE3_Metaphloem_sieve_element_Cluster_MCL[[i]]
  current_module_df <- data.frame(Gene = current_module)
  colnames(current_module_df) <- c("TAIR")
  rownames(current_module_df) <- current_module_df$TAIR
  ids_current_module <- bitr(rownames(current_module_df), 
                             fromType = "TAIR", 
                             toType = c("ENTREZID"), 
                             OrgDb = org.At.tair.db, 
                             drop = TRUE)
  GO_current_module <- enrichGO(ids_current_module$ENTREZID,
                                OrgDb = org.At.tair.db,
                                keyType = "ENTREZID",
                                ont = "BP",  
                                pvalueCutoff = 0.05,
                                qvalueCutoff = 0.05,
                                readable = TRUE)
  GO_BP_Metaphloem_sieve_element[[i]] <- data.frame(GO_current_module)
}

library(dplyr)

for (i in 1:length(GO_BP_Metaphloem_sieve_element)) {
  if (nrow(GO_BP_Metaphloem_sieve_element[[i]]) > 0) { 
    GO_BP_Metaphloem_sieve_element[[i]]$Module <- paste0("Module_", i) 
  }
}

GO_BP_Trichoblast <- list()

for (i in 1:60) {
  current_module <- GENIE3_Trichoblast_Cluster_MCL[[i]]
  current_module_df <- data.frame(Gene = current_module)
  colnames(current_module_df) <- c("TAIR")
  rownames(current_module_df) <- current_module_df$TAIR
  ids_current_module <- bitr(rownames(current_module_df), 
                             fromType = "TAIR", 
                             toType = c("ENTREZID"), 
                             OrgDb = org.At.tair.db, 
                             drop = TRUE)
  GO_current_module <- enrichGO(ids_current_module$ENTREZID,
                                OrgDb = org.At.tair.db,
                                keyType = "ENTREZID",
                                ont = "BP",  
                                pvalueCutoff = 0.05,
                                qvalueCutoff = 0.05,
                                readable = TRUE)
  GO_BP_Trichoblast[[i]] <- data.frame(GO_current_module)
}

library(dplyr)
for (i in 1:length(GO_BP_Trichoblast)) {
  if (nrow(GO_BP_Trichoblast[[i]]) > 0) { 
    GO_BP_Trichoblast[[i]]$Module <- paste0("Module_", i) 
  }
}

GO_BP_Mature <- list()

for (i in 1:71) {
  current_module <- GENIE3_Mature_Cluster_MCL[[i]]
  current_module_df <- data.frame(Gene = current_module)
  colnames(current_module_df) <- c("TAIR")
  rownames(current_module_df) <- current_module_df$TAIR
  ids_current_module <- bitr(rownames(current_module_df), 
                             fromType = "TAIR", 
                             toType = c("ENTREZID"), 
                             OrgDb = org.At.tair.db, 
                             drop = TRUE)
  GO_current_module <- enrichGO(ids_current_module$ENTREZID,
                                OrgDb = org.At.tair.db,
                                keyType = "ENTREZID",
                                ont = "BP", 
                                pvalueCutoff = 0.05,
                                qvalueCutoff = 0.05,
                                readable = TRUE)
  GO_BP_Mature[[i]] <- data.frame(GO_current_module)
}

library(dplyr)
for (i in 1:length(GO_BP_Mature)) {
  if (nrow(GO_BP_Mature[[i]]) > 0) {  
    GO_BP_Mature[[i]]$Module <- paste0("Module_", i)  
  }
}

GO_BP_Xylem <- list()

for (i in 1:90) {
  current_module <- GENIE3_Xylem_Cluster_MCL[[i]]
  current_module_df <- data.frame(Gene = current_module)
  colnames(current_module_df) <- c("TAIR")
  rownames(current_module_df) <- current_module_df$TAIR
  ids_current_module <- bitr(rownames(current_module_df), 
                             fromType = "TAIR", 
                             toType = c("ENTREZID"), 
                             OrgDb = org.At.tair.db, 
                             drop = TRUE)
  GO_current_module <- enrichGO(ids_current_module$ENTREZID,
                                OrgDb = org.At.tair.db,
                                keyType = "ENTREZID",
                                ont = "BP",  
                                pvalueCutoff = 0.05,
                                qvalueCutoff = 0.05,
                                readable = TRUE)

  GO_BP_Xylem[[i]] <- data.frame(GO_current_module)
}

library(dplyr)
for (i in 1:length(GO_BP_Xylem)) {
  if (nrow(GO_BP_Xylem[[i]]) > 0) {  
    GO_BP_Xylem[[i]]$Module <- paste0("Module_", i)  
  }
}

GO_Atrichoblast <- bind_rows(GO_BP_Atrichoblast)
GO_Columella <- bind_rows(GO_BP_Columella)
GO_Cortex <- bind_rows(GO_BP_Cortex)
GO_Endodermis <- bind_rows(GO_BP_Endodermis)
GO_Lateral_root_cap <- bind_rows(GO_BP_Lateral_root_cap)
GO_Metaphloem_sieve_element <- bind_rows(GO_BP_Metaphloem_sieve_element)
GO_Trichoblast <- bind_rows(GO_BP_Trichoblast)
GO_Mature <- bind_rows(GO_BP_Mature)
GO_Xylem <- bind_rows(GO_BP_Xylem)


library(UpSetR)
library(ggplot2)
library(dplyr)
library(purrr)

go_dataframes <- list(
  Atrichoblast = GO_Atrichoblast,
  Columella = GO_Columella,
  Cortex = GO_Cortex,
  Endodermis = GO_Endodermis,
  Lateral_root_cap = GO_Lateral_root_cap,
  Metaphloem_sieve_element = GO_Metaphloem_sieve_element,
  Trichoblast = GO_Trichoblast,
  Mature = GO_Mature,
  Xylem = GO_Xylem
)

process_go_data <- function(go_df, df_name) {
  module_go <- split(go_df$Description, go_df$Module)
  
  all_go <- unique(unlist(module_go))
  binary_mat <- matrix(0, nrow = length(all_go), ncol = length(module_go),
                       dimnames = list(all_go, names(module_go)))
  
  for (mod in names(module_go)) {
    binary_mat[module_go[[mod]], mod] <- 1
  }
  
  list(
    module_go = module_go,
    binary_matrix = binary_mat,
    upset_plot = upset(as.data.frame(binary_mat),
                       nsets = ncol(binary_mat),
                       mainbar.y.label = paste("Share the number of GO Terms -", df_name),
                       sets.x.label = "The number of GO terms")
  )
}

all_results <- map2(go_dataframes, names(go_dataframes), process_go_data)
all_go_terms <- unique(unlist(lapply(go_dataframes, function(x) x$Description)))
cross_dataset_mat <- matrix(0, nrow = length(all_go_terms), 
                            ncol = length(go_dataframes),
                            dimnames = list(all_go_terms, names(go_dataframes)))

for (i in seq_along(go_dataframes)) {
  cross_dataset_mat[unique(go_dataframes[[i]]$Description), i] <- 1
}

cross_dataset_df <- as.data.frame(cross_dataset_mat)

# Generate a UpSet graph with numerical labels
UpSetR::upset(cross_dataset_df,
              nsets = ncol(cross_dataset_df),
              nintersects = 30,
              mb.ratio = c(0.5, 0.3),
              order.by = "freq",
              decreasing = TRUE,
              text.scale = c(
                intersize = 1.0,   
                intsersects = 1.0,  
                set = 1.1,          
                set.size = 1.0     
              ),
              point.size = 3.0,
              line.size = 1.5,
              main.bar.color = "#E74C3C",
              sets.bar.color = "#F1C40F",
              matrix.color = "#27AE60",
              show.numbers = "yes",        
              number.angles = 0,          
              empty.intersections = NULL,  
              
              set_size.show = TRUE,       
              set_size.scale_max = max(colSums(cross_dataset_df)) * 1.1, 
          
              mainbar.y.label = "The number of shared GO Terms (numerical display)",
              sets.x.label = "The total number of GO terms in each network"
)


# 7. Extract GO Terms shared in at least N datasets
shared_go_terms <- cross_dataset_df[rowSums(cross_dataset_df) >= 2, ]  # 修改3为您需要的阈值
print("GO terms shared in at least two datasets:")
print(rownames(shared_go_terms))


# Calculate the Jaccard similarity matrix (similarity between datasets)

library(pheatmap)
library(viridis)
library(ggcorrplot)

go_term_sets <- lapply(go_dataframes, function(df) unique(df$Description))
dataset_names <- names(go_dataframes)
n <- length(go_term_sets)


jaccard_sim <- matrix(0, nrow = n, ncol = n, 
                      dimnames = list(dataset_names, dataset_names))

for (i in 1:(n-1)) {
  set_i <- go_term_sets[[i]]
  for (j in (i+1):n) {
    set_j <- go_term_sets[[j]]
    intersection <- length(intersect(set_i, set_j))
    union_size <- length(set_i) + length(set_j) - intersection
    jaccard <- intersection / union_size
    jaccard_sim[i,j] <- jaccard
    jaccard_sim[j,i] <- jaccard  
  }
}
diag(jaccard_sim) <- 1 


library(ggplot2)
library(reshape2)  

jaccard_melt <- melt(jaccard_sim)
colnames(jaccard_melt) <- c("Dataset1", "Dataset2", "Similarity")

# Create a circular heat map and limit the value range to 0-1
ggplot(jaccard_melt, aes(x = Dataset1, y = Dataset2, fill = Similarity)) +
  geom_point(aes(size = Similarity), shape = 21, color = "white") +
  scale_size_continuous(range = c(3, 10), limits = c(0, 1)) +  
  scale_fill_gradient2(low = "white", mid = "#6D9EC1", high = "#E46726",
                       midpoint = 0.5,
                       limits = c(0, 1), 
                       na.value = "grey90") +
  geom_text(aes(label = round(Similarity, 2)), color = "black", size = 3) +
  labs(title = "数据集相似性矩阵 (Jaccard指数 0-1)",
       x = "", y = "", 
       fill = "Jaccard\nIndex",
       size = "Jaccard\nIndex") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        panel.grid.major = element_blank(),
        legend.position = "right") +
  coord_fixed() +
  guides(size = "none")  


#####
library(ggplot2)
library(dplyr)
library(purrr)

# 1. Combine all dataframes into one with celltype identifier
combined_go <- imap_dfr(go_dataframes, ~ {
  .x %>% 
    mutate(celltype = .y)  # Add celltype column with the list name
})

# 2. Filter top significant GO Terms for each celltype
top_go <- combined_go %>%
  group_by(celltype) %>%
  slice_min(p.adjust, n = 10) %>%  # Get top 5 most significant per celltype
  ungroup()

# 3. Create the visualization
ggplot(top_go, aes(x = celltype, y = Description)) +
  geom_point(aes(size = Count, color = -log10(p.adjust))) +
  scale_size(range = c(4, 10)) +
  scale_color_gradient(low = "blue", high = "red") +
  labs(title = "Significant GO Terms by Cell Type",
       x = NULL, y = NULL,
       size = "Gene Count", 
       color = "-log10(adj)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5))  # Center the title

### Shared path relationship
library(UpSetR)
library(ggplot2)
# 1. Prepare the data (ensure that cross_dataset_df is a binary matrix)
cross_dataset_mat <- as.matrix(cross_dataset_df)
###2. Sankey Diagram (Showing the Flow of Shared Pathways)
library(ggalluvial)
library(dplyr)
library(tidyr)
library(ggalluvial)
library(dplyr)
library(tidyr)

# Prepare the Sankey plot data
sankey_data <- as.data.frame(cross_dataset_mat) %>%
  tibble::rownames_to_column("GO_Term") %>%
  tidyr::gather("CellType", "Present", -GO_Term) %>%
  filter(Present == 1) %>%
  group_by(GO_Term) %>%
  mutate(Shared_In = n()) %>%
  filter(Shared_In > 1) %>% 
  ungroup()

# Draw a Sankey diagram (Select the 20 paths with the highest sharing frequency)
  count(GO_Term, sort = TRUE) %>%
  head(15) %>%
  pull(GO_Term)

ggplot(sankey_data %>% filter(GO_Term %in% top_shared),
       aes(axis1 = reorder(GO_Term, -Shared_In), 
           axis2 = CellType)) +
  geom_alluvium(aes(fill = CellType), width = 1/12) +
  geom_stratum(width = 1/12, fill = "gray90", color = "gray40") +
  
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), 
            size = 3, angle = 0, vjust = 0.5, hjust = 0.5,
            position = position_nudge(y = 0)) + 
  scale_x_discrete(limits = c("GO Term", "CellType"), expand = c(0.05, 0.05)) +
  scale_fill_brewer(palette = "Set3") +
  labs(title = "The relationship between high-frequency shared GO Terms and cell types",
       y = "Count") +
  theme_minimal() +
  theme(legend.position = "bottom",
        axis.text.y = element_blank(),
        panel.grid = element_blank())


### Visual identification of unique pathways
unique_go <- combined_go %>%
  group_by(Description) %>%
  mutate(n_celltypes = n_distinct(celltype)) %>%
  ungroup() %>%
  filter(n_celltypes == 1) %>%
  group_by(celltype) %>%
  slice_min(p.adjust, n = 15) %>%
  ungroup()


library(ggplot2)
library(ggforce)
library(ggrepel)
library(dplyr)
library(stringr)
library(ggstance)  

# -log10(p.adjust)
unique_go$log_p_adj <- -log10(unique_go$p.adjust)

# Screen the top 10 most significant pathways for each cell type
top_go <- unique_go %>%
  group_by(celltype) %>%
  top_n(5, log_p_adj) %>% 
  ungroup()


p <- ggplot(top_go, aes(x = celltype, y = log_p_adj)) +

  geom_point(aes(size = Count, color = log_p_adj), alpha = 0.7, stroke = 0.5) +
  geom_text_repel(
    aes(label = str_wrap(Description, 25)),
    size = 2.8, 
    max.overlaps = Inf, 
    box.padding = 0.8,   
    point.padding = 0.6,
    force = 2,          
    force_pull = 0.1,    
    segment.size = 0.3,  
    segment.alpha = 0.5, 
    direction = "both",  
    nudge_y = 0.1,       
    nudge_x = 0.05       
  ) +
  

  scale_size(range = c(2.5, 10)) +

  scale_color_viridis_c(
    option = "D", 
    direction = -1,
    limits = c(min(top_go$log_p_adj), max(top_go$log_p_adj)),
    breaks = seq(round(min(top_go$log_p_adj)), round(max(top_go$log_p_adj)), by = 1),
    labels = scales::number_format(accuracy = 0.1)
  ) +
  
  labs(
    title = "Cell type-specific pathway cluster diagram",
    subtitle = "Each cell type shows the 10 most significantly enriched pathways",
    x = "cell type", 
    y = "-log10(corrected p value)",
    size = "The number of enriched genes",
    color = "-log10(corrected p value)"
  ) +

  theme_minimal(base_size = 12) +

  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray50"),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.key.size = unit(0.8, "lines"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.title.y = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.margin = margin(1, 2, 1, 1, "cm")  
  ) +
  scale_x_discrete(expand = expansion(mult = c(0.1, 0.3)))


print(p)
