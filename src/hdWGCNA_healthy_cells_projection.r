library(WGCNA)
library(hdWGCNA)
library(glue)
library(Seurat)

setwd("../Data/")

healthy_cells_path = ""  # Modify by the path of your healthy cells path
libs_info <- list(
  list(
    "path"="../data/hdwgcna_output/hdwgcna_processed_lib1.rds",
    "name"="M127",
    "diag_name"="M104",
    "relapse_name"="M127"
  ),
  list(
    "path"="../data/hdwgcna_output/hdwgcna_processed_lib2.rds",
    "name"="M148",
    "diag_name"="143",
    "relapse_name"="M148"
  ),
  list(
    "path"="../data/hdwgcna_output/hdwgcna_processed_lib3.rds",
    "name"="M187r",
    "diag_name"="M187",
    "relapse_name"="M187r"
  )
)

for (i in c(1,2,3)){
  # Load libraries -------------------------------------------------------------
  seurat_query <- readRDS(libs_info[i]["path"])
  seurat_ref <- readRDS(healthy_cells_path)
  name <- libs_info[i]["name"]
  
  ## Project modules from query to reference dataset ---------------------------
  conflict_prefer("cor", "WGCNA")
  seurat_query <- ProjectModules(
    seurat_obj = seurat_query,
    seurat_ref = seurat_ref,
    group.by.vars = "seurat_clusters",
    wgcna_name_proj=glue("healthy2lib{i}"),
    wgcna_name = name
  )
  
  ## Compute module hubGene expression scores for projected modules ------------
  seurat_query <- ModuleExprScore(
    seurat_query,
    n_genes = 25,
    method='Seurat'
  )
  
  ## Intramodular connectivity -------------------------------------------------
  seurat_query <- ModuleConnectivity(seurat_query, assay="RNA", slot="data")
  plot_list <- ModuleFeaturePlot(
    seurat_query,
    features='hMEs',
    order=TRUE
  )
  p <- wrap_plots(plot_list, ncol=5)
  plot(p)
  
  ## Get projected hMEs --------------------------------------------------------
  projected_hMEs <-  GetMEs(seurat_query, harmonized=TRUE)
  seurat_query@meta.data <- cbind(
    seurat_query@meta.data,
    projected_hMEs
  )
  p <- DotPlot(
    seurat_query,
    features = colnames(projected_hMEs),
    group.by = columnToGroup
  )
  
  p <- p +
    RotatedAxis() +
    scale_color_gradient2(high='red', mid='grey95', low='blue') +
    xlab('') + ylab('')
  plot(p)
  
}
