# README

The data used in imput must be processed by Seurat in a previous step. You will therefore need a Seurat object as input, containing your processed data (according to the guidelines provided by https://satijalab.org/seurat/).

Scripts should be run in the following order:
1. hdWGCNA_script.r: Construction of the Gene Coexpression Modules (GCNs) through the package hdWGCNA (https://smorabit.github.io/hdWGCNA/articles/hdWGCNA.html)
2. hdWGCNA_healthy_cells_projection.r: Projection of the GCNs onto healthy cells through the package hdWGCNA (https://smorabit.github.io/hdWGCNA/articles/hdWGCNA.html)
3. sPLS-DA.r: Determination of the hub genes most involved in the differentiation of Diagnosis/Relapse states with the MixOmics Package (http://mixomics.org)
4. coxph_optimization.r: COXPH optimization using the MASS package (https://cran.r-project.org/web/packages/MASS/index.html)

**NOTE**: In each of these, you will need to modify the paths to the various input files, replacing them with those pointing to your data. 
