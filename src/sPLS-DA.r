library(mixOmics)
library(gridExtra)
library(glue)
library(VennDiagram)
library(conflicted)

setwd("/src")

# Import the preselected hub genes ---------------------------------------------
hubGenesFiltered <- read.table("../data/hdwgcna/hubGenesFiltered.txt", quote="\"", comment.char="")
hubGenes.unique <- unique(hubGenesFiltered$V1)

# Import libraries -------------------------------------------------------------
lib1 <- readRDS("../data/hdwgcna/hdWCGNA_lib1.rds")
lib2 <- readRDS("../data/hdwgcna/hdWCGNA_lib2.rds")
lib3 <- readRDS("../data/hdwgcna/hdWCGNA_lib3.rds")
libs <- list(lib1, lib2, lib3)
lib1.genes <- names(as.data.frame(t(lib1@assays$RNA@scale.data)[,names(as.data.frame(t(lib1@assays$RNA@scale.data))) %in% hubGenes.unique]))
lib2.genes <- names(as.data.frame(t(lib2@assays$RNA@scale.data)[,names(as.data.frame(t(lib2@assays$RNA@scale.data))) %in% hubGenes.unique]))
lib3.genes <- names(as.data.frame(t(lib3@assays$RNA@scale.data)[,names(as.data.frame(t(lib3@assays$RNA@scale.data))) %in% hubGenes.unique]))

# Intersection of all the hub genes in our samples -----------------------------
intersect.genes <- Reduce(intersect, list(lib1.genes,lib2.genes,lib3.genes))


# Preservation / tune step -----------------------------------------------------
count = 1
for (raw in c(lib1, lib2, lib3)) {
  lib <- raw
  count <- count + 1
  ncomp <- 5
  
  lib.scaled.df.X <- as.data.frame(t(lib@assays$RNA@scale.data)[,names(as.data.frame(t(lib@assays$RNA@scale.data))) %in% intersect.genes])
  lib.scaled.df.Y <- lib$hash.ID
  
  result.splsda.lib <- splsda(lib.scaled.df.X,lib.scaled.df.Y,ncomp = ncomp)
  
  ## Find best components number -----------------------------------------------
  perf.splsda.lib <- perf(result.splsda.lib, validation = "Mfold", 
                          folds = 10, nrepeat = 100, 
                          progressBar = TRUE, auc = TRUE) 
  
  plot(perf.splsda.lib, col = color.mixo(5:7))
  perf.splsda.lib$choice.ncomp
  saveRDS(perf.splsda.lib,paste0("../data/spls_da_output/PLS-DA_perf/perf100repeat_sPLSDA_lib",count,".rds"))
  
  list.keepX <- c(1:5,  seq(1, 200, 1))
  tune.splsda.lib <- tune.splsda(lib.scaled.df.X, lib.scaled.df.Y, ncomp = 5,
                                 validation = 'Mfold',
                                 folds = 10, nrepeat = 100,
                                 dist = 'max.dist',
                                 measure = "BER",
                                 test.keepX = list.keepX,
                                 progressBar = TRUE)
  saveRDS(tune.splsda.lib,paste0("../data/spls_da_output/PLS-DA_tune/tune100repeat_sPLSDA_lib",count,".rds"))
}

# Plot sPLS-DA tuned biplot ----------------------------------------------------
perflib1 <- readRDS("../data/spls_da_output/PLS-DA_tune/perf100repeat_sPLSDA_lib1.rds")
perflib2 <- readRDS("../data/spls_da_output/PLS-DA_tune/perf100repeat_sPLSDA_lib2.rds")
perflib3 <- readRDS("../data/spls_da_output/PLS-DA_tune/perf100repeat_sPLSDA_lib3.rds")

tunelib1 <- readRDS("../data/spls_da_output/PLS-DA_tune/tune100repeat_sPLSDA_lib1.rds")
tunelib2 <- readRDS("../data/spls_da_output/PLS-DA_tune/tune100repeat_sPLSDA_lib2.rds")
tunelib3 <- readRDS("../data/spls_da_output/PLS-DA_tune/tune100repeat_sPLSDA_lib3.rds")


select.keepX.libs <- list("lib1" = c(),"lib2" = c(),"lib3" = c())
ncomp.libs <- list("lib1" = c(),"lib2" = c(),"lib3" = c())
for (i in c(1,2,3)) {
  perf.splsda.lib <- readRDS(glue("../data/spls_da_output/PLS-DA_tune/perf100repeat_sPLSDA_lib{i}.rds"))
  tune.splsda.lib <- readRDS(glue("../data/spls_da_output/PLS-DA_tune/tune100repeat_sPLSDA_lib{i}.rds"))
  lib <- libs[i]
  
  print(glue("Processing lib {i}"))
  print("N comp according to perf():", perf.splsda.lib$choice.ncomp)
  ## optimal number of components based on t-tests on the error rate -----------
  ncomp.lib <- tune.splsda.lib$choice.ncomp$ncomp 
  print("N comp according to t-test on the error rate:", ncomp.lib)
  ## optimal number of variables to select -------------------------------------
  select.keepX.lib <- na.omit(tune.splsda.lib$choice.keepX[1:ncomp.lib])
  print("Optimal number of variables to select:",select.keepX.lib)
  
  lib.scaled.df.X <- as.data.frame(t(lib@assays$RNA@scale.data)[,names(as.data.frame(t(lib@assays$RNA@scale.data))) %in% intersect.genes])
  lib.scaled.df.Y <- lib$hash.ID
  result.splsda.lib <- splsda(lib.scaled.df.X,lib.scaled.df.Y,ncomp = ncomp.lib, keepX = select.keepX.lib)
  png(glue("./data/spls_da_output/splsDABIPLOT_lib{i}.png"), width = 800, height = 600)
  plotIndiv(result.splsda.lib, ind.names = FALSE, legend=TRUE,ellipse = TRUE, ellipse.level = 0.95,title="sPLS-DA - M187/187r")
  dev.off()
  
  select.keepX.libs[glue("lib{i}")] <- select.keepX.lib
  ncomp.libs[glue("lib{i}")] <- ncomp.lib
}

# Récupération des top features ------------------------------------------------
lib1.scaled.df.X <- as.data.frame(t(lib1@assays$RNA@scale.data)[,names(as.data.frame(t(lib1@assays$RNA@scale.data))) %in% intersect.genes])
lib2.scaled.df.X <- as.data.frame(t(lib2@assays$RNA@scale.data)[,names(as.data.frame(t(lib2@assays$RNA@scale.data))) %in% intersect.genes])
lib3.scaled.df.X <- as.data.frame(t(lib3@assays$RNA@scale.data)[,names(as.data.frame(t(lib3@assays$RNA@scale.data))) %in% intersect.genes])

lib1.scaled.df.Y <- lib1$hash.ID
lib2.scaled.df.Y <- lib2$hash.ID
lib3.scaled.df.Y <- lib3$hash.ID

result.splsda.lib1 <- splsda(lib1.scaled.df.X,lib1.scaled.df.Y,ncomp = ncomp.libs["lib1"], keepX = select.keepX.libs["lib1"])
result.splsda.lib2 <- splsda(lib2.scaled.df.X,lib2.scaled.df.Y,ncomp = ncomp.libs["lib2"], keepX = select.keepX.libs["lib2"])
result.splsda.lib3 <- splsda(lib3.scaled.df.X,lib3.scaled.df.Y,ncomp = ncomp.libs["lib3"], keepX = select.keepX.libs["lib3"])

plotIndiv(result.splsda.lib1, comp = c(1,2), ind.names = FALSE, legend=TRUE, ellipse = TRUE, ellipse.level = 0.95,title="sPLS-DA")
plotIndiv(result.splsda.lib2, comp = c(1,2), ind.names = FALSE, legend=TRUE, ellipse = TRUE, ellipse.level = 0.95,title="sPLS-DA")
plotIndiv(result.splsda.lib3, comp = c(1,2), ind.names = FALSE, legend=TRUE, ellipse = TRUE, ellipse.level = 0.95,title="sPLS-DA")


lib1.vars <- unique(c(rownames(head(as.data.frame(result.splsda.lib1$loadings$X)[order(as.data.frame(result.splsda.lib1$loadings$X)$comp1, decreasing = TRUE),],10)),
                      rownames(head(as.data.frame(result.splsda.lib1$loadings$X)[order(as.data.frame(result.splsda.lib1$loadings$X)$comp2, decreasing = TRUE),],10)),
                      rownames(head(as.data.frame(result.splsda.lib1$loadings$X)[order(as.data.frame(result.splsda.lib1$loadings$X)$comp3, decreasing = TRUE),],10)),
                      rownames(head(as.data.frame(result.splsda.lib1$loadings$X)[order(as.data.frame(result.splsda.lib1$loadings$X)$comp4, decreasing = TRUE),],10)),
                      rownames(head(as.data.frame(result.splsda.lib1$loadings$X)[order(as.data.frame(result.splsda.lib1$loadings$X)$comp5, decreasing = TRUE),],10))))

lib2.vars <- unique(c(rownames(head(as.data.frame(result.splsda.lib2$loadings$X)[order(as.data.frame(result.splsda.lib2$loadings$X)$comp1, decreasing = TRUE),],10)),
                      rownames(head(as.data.frame(result.splsda.lib2$loadings$X)[order(as.data.frame(result.splsda.lib2$loadings$X)$comp2, decreasing = TRUE),],10)),
                      rownames(head(as.data.frame(result.splsda.lib2$loadings$X)[order(as.data.frame(result.splsda.lib2$loadings$X)$comp3, decreasing = TRUE),],10)),
                      rownames(head(as.data.frame(result.splsda.lib2$loadings$X)[order(as.data.frame(result.splsda.lib2$loadings$X)$comp4, decreasing = TRUE),],10)),
                      rownames(head(as.data.frame(result.splsda.lib2$loadings$X)[order(as.data.frame(result.splsda.lib2$loadings$X)$comp5, decreasing = TRUE),],10))))

lib3.vars <- unique(c(rownames(head(as.data.frame(result.splsda.lib3$loadings$X)[order(as.data.frame(result.splsda.lib3$loadings$X)$comp1, decreasing = TRUE),],10)),
                      rownames(head(as.data.frame(result.splsda.lib3$loadings$X)[order(as.data.frame(result.splsda.lib3$loadings$X)$comp2, decreasing = TRUE),],10)),
                      rownames(head(as.data.frame(result.splsda.lib3$loadings$X)[order(as.data.frame(result.splsda.lib3$loadings$X)$comp3, decreasing = TRUE),],10)),
                      rownames(head(as.data.frame(result.splsda.lib3$loadings$X)[order(as.data.frame(result.splsda.lib3$loadings$X)$comp4, decreasing = TRUE),],10)),
                      rownames(head(as.data.frame(result.splsda.lib3$loadings$X)[order(as.data.frame(result.splsda.lib3$loadings$X)$comp5, decreasing = TRUE),],10))))

vennProtocols <- venn.diagram(
  x = list(
    M104_M127= lib1.vars ,
    M143_M148= lib2.vars,
    M187_M187r = lib3.vars),
  filename = NULL,
  cex=1.5, cat.cex=1.5,
  fill = c('green', 'darkblue',  'yellow')
)
dev.off()
grid.draw(vennProtocols)

conflicts_prefer(dplyr::intersect)
intersect.keygenes <- Reduce(intersect, list(lib1.vars,lib2.vars,lib3.vars))
intersect.keygenes