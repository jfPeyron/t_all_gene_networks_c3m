# Setup ------------------------------------------------------------------------

library(dplyr)
library(conflicted)
library(Seurat)
library(patchwork)
library(gridExtra)
library(ggplot2)
library(grid)
library(RColorBrewer)
library(igraph)
library(WGCNA)
library(org.Hs.eg.db)
library(GOplot)
library(hdWGCNA)
library(UCell)
library(ggrepel)
library(enrichR)

setwd("/src")

UCELL_METHOD <- TRUE
GO_ANALYSIS <- TRUE

libs_info <- list(
  list(
    "path"="../data/librairie1Clusterised.rds",
    "name"="M127",
    "diag_name"="M104",
    "relapse_name"="M127"
  ),
  list(
    "path"="../data/librairie2Clusterised.rds",
    "name"="M148",
    "diag_name"="143",
    "relapse_name"="M148"
  ),
  list(
    "path"="../data/librairie3Clusterised.rds",
    "name"="M187r",
    "diag_name"="M187",
    "relapse_name"="M187r"
  )
)

# Start analysis ---------------------------------------------------------------

for (i in c(1,2,3)) {
  singlecell_cluster <- readRDS(libs_info[[i]][["path"]])
  name <- libs_info[[i]][["name"]]
  
  ## Loading Seurat File clusterised --------------------------------------------
  Idents(singlecell_cluster) <- singlecell_cluster$hash.ID
  seurat_obj <- singlecell_cluster
  seurat_obj <- SetupForWGCNA(
    seurat_obj,
    gene_select = "fraction",
    fraction = 0.05,
    wgcna_name = name,
  )
  
  ## Construction metacells for each group --------------------------------------
  seurat_obj <- MetacellsByGroups(
    seurat_obj = seurat_obj,
    group.by = c("hash.ID"),
    k = 25,
    max_shared = 10,
    ident.group = 'hash.ID'
  )
  
  ## Normalize metacells expression matrix --------------------------------------
  seurat_obj <- NormalizeMetacells(seurat_obj)
  metacell_obj <- GetMetacellObject(seurat_obj)
  seurat_obj <- NormalizeMetacells(seurat_obj)
  seurat_obj <- ScaleMetacells(seurat_obj, features=VariableFeatures(seurat_obj))
  seurat_obj <- RunPCAMetacells(seurat_obj, features=VariableFeatures(seurat_obj))
  seurat_obj <- RunHarmonyMetacells(seurat_obj, group.by.vars='hash.ID')
  seurat_obj <- RunUMAPMetacells(seurat_obj, reduction='harmony', dims=1:15)
  
  p1 <- DimPlotMetacells(seurat_obj, group.by='hash.ID') + umap_theme() + ggtitle("Cell Type")
  plot(p1) 
  
  ## Select the cluster to analyze ----------------------------------------------
  seurat_obj <- SetDatExpr(
    seurat_obj,
    group_name = name,
    group.by='hash.ID',
    assay = 'RNA',
    slot = 'data'
  )
  
  ## Soft powers analysis -------------------------------------------------------
  seurat_obj <- TestSoftPowers(
    seurat_obj,
    networkType = 'signed'
  )
  plot_list <- PlotSoftPowers(seurat_obj)
  wrap_plots(plot_list, ncol=2)
  power_table <- GetPowerTable(seurat_obj)
  head(power_table)
  
  ## Co-expression network ------------------------------------------------------
  conflict_prefer("cor", "WGCNA")
  seurat_obj <- ConstructNetwork(
    seurat_obj, soft_power=4,
    setDatExpr=FALSE,
    tom_name = name,
    overwrite_tom = TRUE
  )
  PlotDendrogram(seurat_obj, main=paste0(name,'hdWGCNA Dendrogram'))
  seurat_obj <- ScaleData(seurat_obj, features=VariableFeatures(seurat_obj))
  
  ## Get ModuleEigengenes and harmonized ModuleEigengenes -----------------------
  seurat_obj <- ModuleEigengenes(
    seurat_obj,
    group.by.vars="hash.ID",
    wgcna_name = name
  )
  MEs <- GetMEs(seurat_obj, harmonized=FALSE)
  hMEs <- GetMEs(seurat_obj)
  
  ## Get Eigengene-based connectivity (kME): ------------------------------------
  seurat_obj <- ModuleConnectivity(
    seurat_obj,
    group.by = 'hash.ID',
    group_name = name
  )
  
  ## Rename modules -------------------------------------------------------------
  seurat_obj <- ResetModuleNames(
    seurat_obj,
    new_name = paste(name,"-M")
  )
  
  ## Plots ----------------------------------------------------------------------
  conflict_prefer("select", "dplyr")
  p <- PlotKMEs(seurat_obj, ncol=2)
  p
  
  ## Get the module assignment table --------------------------------------------
  modules <- GetModules(seurat_obj)
  head(modules[,1:6])
  
  ## Get hub genes --------------------------------------------------------------
  hub_df <- GetHubGenes(seurat_obj, n_hubs = 10)
  head(hub_df)
  
  ## Compute gene scoring for the top 25 hub genes by kME for each module -------
  seurat_obj <- ModuleExprScore(
    seurat_obj,
    n_genes = 25,
    method='Seurat'
  )
  
  if (UCELL_METHOD){
    seurat_obj <- ModuleExprScore(
      seurat_obj,
      n_genes = 25,
      method='UCell'
    )
  }
  
  ## Generate TOM ---------------------------------------------------------------
  TOM = GetTOM(seurat_obj)
  
  ## Plots ----------------------------------------------------------------------
  plot_list <- ModuleFeaturePlot(
    seurat_obj,
    features='hMEs',
    order=TRUE
  )
  wrap_plots(plot_list, ncol=3)
  
  
  plot_list <- ModuleFeaturePlot(
    seurat_obj,
    features='scores', 
    order= TRUE,
    ucell = TRUE
  )
  wrap_plots(plot_list, ncol=3)
  
  
  ModuleCorrelogram(seurat_obj)
  ## Filter modules -------------------------------------------------------------
  mods <- colnames(MEs); mods <- mods[mods != 'grey']
  ## Add hMEs to Seurat metadata ------------------------------------------------
  seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs)
  
  p <- DotPlot(seurat_obj, features=mods, group.by = 'hash.ID')
  p <- p +
    coord_flip() +
    RotatedAxis() +
    scale_color_gradient2(high='red', mid='grey95', low='blue')
  p
  
  
  ## Plot module ----------------------------------------------------------------
  ModuleNetworkPlot(seurat_obj, outdir = paste0("TOM/",name))

  graphics.off() # To solve the graphical interfering issue
  HubGeneNetworkPlot(
    seurat_obj,
    n_hubs = 3, n_other=5,
    edge_prop = 0.75,
    mods = 'all'
  )
  g <- HubGeneNetworkPlot(seurat_obj,  return_graph=TRUE)
  
  ## Plot Module UMAP -----------------------------------------------------------
  seurat_obj <- RunModuleUMAP(
    seurat_obj,
    n_hubs = 10,
    n_neighbors=15, 
    min_dist=0.1
  )
  # umap_df <- GetModuleUMAP(seurat_obj)
  # ggplot(umap_df, aes(x=UMAP1, y=UMAP2)) +
  #   geom_point(
  #     color=umap_df$color,
  #     size=umap_df$kME*2
  #   ) +
  #   umap_theme()
  
  ModuleUMAPPlot(
    seurat_obj,
    edge.alpha=0.25,
    sample_edges=TRUE,
    edge_prop=0.1,
    label_hubs=2 ,
    keep_grey_edges=FALSE
    
  )
  
  ## DME analysis comparing two groups (Relapse and Diag) -----------------------
  diag <- seurat_obj@meta.data %>% subset(hash.ID == libs_info[[i]][["diag_name"]]) %>% rownames()
  rechute <- seurat_obj@meta.data %>% subset(hash.ID == libs_info[[i]][["relapse_name"]]) %>% rownames()
  
  DMEs <- FindDMEs(
    seurat_obj,
    barcodes1 = rechute,
    barcodes2 = diag,
    test.use='wilcox',
    wgcna_name=name
  )
  
  PlotDMEsVolcano(
    seurat_obj,
    DMEs,
    wgcna_name = name
  )
  
  
  ## Gene Ontology --------------------------------------------------------------
  if (GO_ANALYSIS) {
    dbs <- c('GO_Biological_Process_2021','GO_Cellular_Component_2021','GO_Molecular_Function_2021')
    seurat_obj <- RunEnrichr(
      seurat_obj,
      dbs=dbs,
      max_genes = 100
    )
    enrich_df <- GetEnrichrTable(seurat_obj)
    conflict_prefer("select", "dplyr")
    for (db in dbs){
      p <- EnrichrDotPlot(
        seurat_obj,
        mods = "all",
        database = db,
        n_terms=3 
      )
      plot(p)
    }
  }

  saveRDS(seurat_obj, file = paste0("hdwgcna_processed_lib",i,".RDS"))
}
