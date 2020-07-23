---
title: "Seurat Polyp Comparison Report"
author: Kyle Kimler
date: July 22, 2020
output:
  html_document:
    keep_md: true
    toc: true
    number_sections: false
    toc_float: true
---

Introduction:
------------

Comparing UMAP with t-SNE using many different numbers of dimensions and annotating clusters that don't contain canonical markers by their first rank Wilcoxon ranked gene

------------


```r
library(tidyverse)
library(Seurat)
library(Matrix)

umifile <- read.csv(file="/Users/kkimler/Projects/NasalPolyp/data/polypALL.csv",header=TRUE,stringsAsFactors=FALSE)
umidata <- umifile %>% remove_rownames %>% column_to_rownames(colnames(umifile)[1])
#subset for testing
#umidata <- umidata %>% select(starts_with(c("Polyp1TOT","Polyp3TOT")))
#umidata <- umidata %>% slice(1:500)

metadata <- read.csv("/Users/kkimler/Projects/NasalPolyp/data/polypMetadata.csv")
#reduce metadata table to data you actually have
metadata <- metadata[1:12,]

#Need a metadata row for each cell, not for each patient.
metadatalist = NULL
for (id in metadata$orig.ident){
	meta <- data.frame(patient = rep(id,length(umidata %>% select(starts_with(id)))), row.names=colnames(umidata %>% select(starts_with(id))))
	metadatalist = bind_rows(metadatalist,meta)
}
metadatalist <- metadatalist %>% left_join(metadata,by=c("patient"="orig.ident"))
row.names(metadatalist) <- colnames(umidata)

#sparsify the dataset
umisparse <- umidata %>% data.matrix %>% Matrix(sparse=TRUE)
rm(umidata)

polyp <- CreateSeuratObject(umisparse, meta.data = metadatalist, project = "polyp_scRNAseq")

polyp <- NormalizeData(polyp, normalization.method = "LogNormalize", scale.factor = 10000)

polyp.list <- SplitObject(polyp, split.by = "polyp")

polyp.list <- lapply(X = polyp.list, FUN = function(x) {
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

for (z in c(12,15,20,25,30)){
anchors <- FindIntegrationAnchors(object.list = polyp.list, dims=1:z)
integratedpolyp <- IntegrateData(anchorset=anchors, dims=1:z)
all.genes <- rownames(polyp)
integratedpolyp <- ScaleData(integratedpolyp, features=all.genes, vars.to.regress=NULL)
integratedpolyp <- RunPCA(integratedpolyp, npcs=z, verbose=FALSE)
integratedpolyp <- FindNeighbors(integratedpolyp, dims = 1:z)
integratedpolyp <- FindClusters(integratedpolyp, resolution = 0.25)
integratedpolyp <- RunUMAP(integratedpolyp, dims = 1:z)
DefaultAssay(integratedpolyp) <- "RNA"
integratedpolyp.markers <- FindAllMarkers(integratedpolyp, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
integratedpolyp.markers %>% group_by(cluster) %>% filter(p_val_adj<0.001) %>% top_n(n = 1, wt = avg_logFC) -> biomarkers
biomarkers <- unlist(as.list(biomarkers$gene))
#This doesn't work!!!
#integratedpolyp.markers %>% group_by(cluster) %>% top_n(n = -1, wt = p_val) -> biomarkers
#biomarkers <- unlist(as.list(biomarkers$gene))
classic.biomarkers <- data.frame(row.names=c("KRT5","KRT8","LTF","FOXJ1","COL1A2","DARC","CD79A","TRBC2","HLA-DRA","TPSAB1"),celltype=c("Basal","Apical","Glandular","Ciliated","Fibroblast","Endothelial","Plasma_Cell","T_Cell","Myeloid","Mast_Cell"))
classic.biomarkers$clusterID <- integratedpolyp.markers[rownames(classic.biomarkers),]$cluster
classic.biomarkers <- classic.biomarkers %>% group_by(clusterID) %>% summarise(celltype=toString(celltype))
classic.biomarkers <- classic.biomarkers %>% complete(clusterID) %>% data.frame 
classic.biomarkers$newbiomarkers <- paste(biomarkers,"characterized cells")
classic.biomarkers <- classic.biomarkers %>% mutate(celltype=coalesce(celltype,newbiomarkers))
IdentityRenamingList <- split(classic.biomarkers$celltype,classic.biomarkers$clusterID)
integratedpolyp <- RenameIdents(integratedpolyp, IdentityRenamingList)
markers.to.plot <- unique(c("KRT5","KRT8","LTF","FOXJ1","COL1A2","DARC","CD79A","TRBC2","HLA-DRA","TPSAB1",unlist(as.list((integratedpolyp.markers %>% group_by(cluster) %>% top_n(n=2,wt=avg_logFC))$gene))))
integratedpolyp <- RunTSNE(integratedpolyp, dims= 1:z)
plt1 <- DimPlot(integratedpolyp, reduction = "tsne", pt.size=0.2, split.by="polyp") + labs(title=paste("t-sne from",z,"PC's"))
plt2 <- DimPlot(integratedpolyp,reduction="umap", pt.size=0.2, split.by="polyp") + labs(title=paste("UMAP from",z,"PC's"))
print(plt1)
print(plt2)
}
```

```
## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
## 
## Number of nodes: 19196
## Number of edges: 658607
## 
## Running Louvain algorithm...
## Maximum modularity in 10 random starts: 0.9502
## Number of communities: 12
## Elapsed time: 2 seconds
```

![](UMAPTsnePolypDims5s_files/figure-html/unnamed-chunk-1-1.png)<!-- -->![](UMAPTsnePolypDims5s_files/figure-html/unnamed-chunk-1-2.png)<!-- -->

```
## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
## 
## Number of nodes: 19196
## Number of edges: 684126
## 
## Running Louvain algorithm...
## Maximum modularity in 10 random starts: 0.9518
## Number of communities: 12
## Elapsed time: 2 seconds
```

![](UMAPTsnePolypDims5s_files/figure-html/unnamed-chunk-1-3.png)<!-- -->![](UMAPTsnePolypDims5s_files/figure-html/unnamed-chunk-1-4.png)<!-- -->

```
## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
## 
## Number of nodes: 19196
## Number of edges: 716233
## 
## Running Louvain algorithm...
## Maximum modularity in 10 random starts: 0.9535
## Number of communities: 14
## Elapsed time: 2 seconds
```

![](UMAPTsnePolypDims5s_files/figure-html/unnamed-chunk-1-5.png)<!-- -->![](UMAPTsnePolypDims5s_files/figure-html/unnamed-chunk-1-6.png)<!-- -->

```
## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
## 
## Number of nodes: 19196
## Number of edges: 733446
## 
## Running Louvain algorithm...
## Maximum modularity in 10 random starts: 0.9519
## Number of communities: 13
## Elapsed time: 2 seconds
```

![](UMAPTsnePolypDims5s_files/figure-html/unnamed-chunk-1-7.png)<!-- -->![](UMAPTsnePolypDims5s_files/figure-html/unnamed-chunk-1-8.png)<!-- -->

```
## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
## 
## Number of nodes: 19196
## Number of edges: 760505
## 
## Running Louvain algorithm...
## Maximum modularity in 10 random starts: 0.9524
## Number of communities: 13
## Elapsed time: 3 seconds
```

![](UMAPTsnePolypDims5s_files/figure-html/unnamed-chunk-1-9.png)<!-- -->![](UMAPTsnePolypDims5s_files/figure-html/unnamed-chunk-1-10.png)<!-- -->
