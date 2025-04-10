######################################################
# CREATE SEURAT AND ANNDATA FILES
#
# GOAL: Use expression matrix and bead coordinate data
# to create Seurat- and Scanpy-compatible versions of
# the data with a spatial embedding and basic
# pre-processing.
######################################################

# ================= IMPORTS ==========================

library(Seurat)
library(SeuratDisk)

# ================= PARAMS ===========================

# get manifest params
manifestFile <- commandArgs(trailingOnly=TRUE)[1]
true <- TRUE
false <- FALSE
source(manifestFile)
if (!endsWith(writeDir,'/')) { writeDir <- paste0(writeDir,'/') }
if (!endsWith(syrahDir,'/')) { syrahDir <- paste0(syrahDir,'/') }

if (!exists('minUMI')) { minUMI <- 10 }

# ================= METHODS ==========================

# ================= INIT DATA ========================

puck <- read.delim(puckFile,row.names=1,header=FALSE)
syrah <- read.delim(paste0(writeDir,batchName,"_counts.tsv.gz"),row.names=1)

if (doNonSyrah) { std <-  read.delim(paste0(writeDir,batchName,"_nonSyrah_counts.tsv.gz"),row.names=1) }

# ================= SYRAH ============================

# SEURAT OBJECT ==================
counts <- syrah[,colSums(syrah)>=minUMI]
seu <- CreateSeuratObject(counts=counts,project="curio_test_data")

# SPATIAL EMBEDDING ==============
coords <- as.matrix(puck[Cells(seu),])
colnames(coords) <- c("SPATIAL_1","SPATIAL_2")
seu@reductions[["SPATIAL"]] <- CreateDimReducObject(embeddings=coords,assay="RNA",key="SPATIAL_")

# STANDARD PROCESSING ============
seu <- SCTransform(seu,conserve.memory=TRUE)
seu <- RunPCA(seu)
seu <- FindNeighbors(seu)
seu <- FindClusters(seu)
seu <- RunUMAP(seu,dims=1:30)

# SAVE SEURAT ====================
saveRDS(seu,paste0(writeDir,batchName,"_Seurat.rds"),compress=FALSE)

# CONVERT TO ANNDATA =============
SaveH5Seurat(seu, filename=paste0(writeDir,batchName,"_AnnData.h5Seurat"))
Convert(paste0(writeDir,batchName,"_AnnData.h5Seurat"), dest="h5ad")
file.remove(paste0(writeDir,batchName,"_AnnData.h5Seurat"))

# ================= NON-SYRAH ========================

if (doNonSyrah) {
  # SEURAT OBJECT ==================
  counts <- std[,colSums(std)>=minUMI]
  seu <- CreateSeuratObject(counts=counts,project="curio_test_data_nonSyrah")
  
  # SPATIAL EMBEDDING ==============
  coords <- as.matrix(puck[Cells(seu),])
  colnames(coords) <- c("SPATIAL_1","SPATIAL_2")
  seu@reductions[["SPATIAL"]] <- CreateDimReducObject(embeddings=coords,assay="RNA",key="SPATIAL_")
  
  # STANDARD PROCESSING ============
  seu <- SCTransform(seu,conserve.memory=TRUE)
  seu <- RunPCA(seu)
  seu <- FindNeighbors(seu)
  seu <- FindClusters(seu)
  seu <- RunUMAP(seu,dims=1:30)
  
  # SAVE SEURAT ====================
  saveRDS(seu,paste0(writeDir,batchName,"_nonSyrah_Seurat.rds"),compress=FALSE)
  
  # CONVERT TO ANNDATA =============
  SaveH5Seurat(seu, filename=paste0(writeDir,batchName,"_nonSyrah_AnnData.h5Seurat"))
  Convert(paste0(writeDir,batchName,"_nonSyrah_AnnData.h5Seurat"), dest="h5ad")
  file.remove(paste0(writeDir,batchName,"_nonSyrah_AnnData.h5Seurat"))
}
