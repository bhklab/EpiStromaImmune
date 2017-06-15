# EpiStromaImmune

Abstract
--------
Engaging the immune system promises to be key for optimal cancer therapy, especially in hard-to-treat triple-negative breast cancer (TNBC). Using laser capture microdissection-coupled expression profiling, we identify three tumor-associated immune microenvironments with distinct CD8+ T cell localization and outcome. Approximately 25% of TNBCs possess an immunoreactive microenvironment, defined by enhanced infiltration of granzyme B+/CD8+ T cells into the tumor bed and a type I interferon signature. These display elevated expression of multiple immune checkpoint inhibitors but are associated with good outcome. In contrast, TNBCs with an “immune-cold” microenvironment restrict CD8+ T cells to tumor margins, possess elevated expression of the immunosuppressive marker B7-H4, and exhibit signatures of activated stroma and worse outcome. A third immunomodulatory microenvironment, also associated with worse outcome, is enriched for IL-17-producing cells and neutrophils and exhibits stromal localisation of CD8+ T cells and PD-L1. These distinct immune microenvironments have implications for TNBC patient stratification for immunotherapies.

Citation
--------

To be publsihed.

Full Reproducibility of the Analysis Results
--------------------------------------------

We describe below how to fully reproduce the figures and tables reported in the paper

1.  Set up the software environment

2.  Run the R scripts

3.  Generate figures

Set up the software environment
-------------------------------

We developed and tested our analysis pipeline using R running on linux and Mac OSX platforms. The following is a copy of `sessionInfo()` from the development environment in R

```
R version 3.1.3 (2015-03-09)
Platform: x86_64-apple-darwin13.4.0 (64-bit)
Running under: OS X 10.10.2 (Yosemite)

attached base packages:
[1] grid      parallel  stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] BiocInstaller_1.16.5 moments_0.14         GSA_1.03             VennDiagram_1.6.16   futile.logger_1.4.1  PharmacoGx_1.1.5     Hmisc_3.17-2         ggplot2_2.1.0       
 [9] Formula_1.2-1        lattice_0.20-33      gplots_2.17.0        xlsx_0.5.7           stringi_1.0-1        affy_1.44.1          Biobase_2.26.0       BiocGenerics_0.12.1 
[17] xlsxjars_0.6.1       rJava_0.9-8          pvclust_2.0-0        piano_1.6.2          genefu_1.16.0        biomaRt_2.22.0       mclust_5.1           survcomp_1.16.0     
[25] prodlim_1.5.7        survival_2.38-3      RCytoscape_1.16.0    XMLRPC_0.3-0         graph_1.44.1        

loaded via a namespace (and not attached):
 [1] acepack_1.3-3.3       affyio_1.34.0         amap_0.8-14           AnnotationDbi_1.28.2  bitops_1.0-6          bootstrap_2015.2      caTools_1.17.1        cluster_2.0.3        
 [9] colorspace_1.2-6      corrplot_0.73         DBI_0.3.1             digest_0.6.9          downloader_0.4        foreign_0.8-66        futile.options_1.0.0  gdata_2.17.0         
[17] GenomeInfoDb_1.2.5    gridExtra_2.2.1       gtable_0.2.0          gtools_3.4.2          igraph_1.0.1          IRanges_2.0.1         KernSmooth_2.23-15    lambda.r_1.1.7       
[25] latticeExtra_0.6-28   lava_1.4.1            limma_3.22.7          lsa_0.73.1            magicaxis_1.9.4       magrittr_1.5          marray_1.44.0         MASS_7.3-45          
[33] munsell_0.4.3         nnet_7.3-12           plotrix_3.6-1         plyr_1.8.3            preprocessCore_1.28.0 RColorBrewer_1.1-2    Rcpp_0.12.3           RCurl_1.95-4.8       
[41] relations_0.6-6       rmeta_2.16            rpart_4.1-10          RSQLite_1.0.0         S4Vectors_0.4.0       scales_0.4.0          sets_1.0-16           slam_0.1-32          
[49] sm_2.2-5.4            SnowballC_0.5.1       splines_3.1.3         stats4_3.1.3          SuppDists_1.1-9.2     survivalROC_1.0.3     tools_3.1.3           XML_3.98-1.4         
[57] zlibbioc_1.12.0
```

All these packages are available on [CRAN](http://cran.r-project.org) or [Bioconductor](http://www.bioconductor.org)

all necessary packages have `library(<package>)` calls within the R scripts themselves, or the script assumes a previous script has been run and thus should have loaded nessesary packages. 

Running R Scripts
-------------------------------
CoreDenScripts-

"CoreDen-GSEA-BT-IPA.R": Pathway significantly associated with core density phenotype on the bulktumor data, using Gene Set Enrichment Analysis (GSEA) as implemented in the Piano R package
"CoreDen-IPA-PathwayScores.R": To obtain Metagene Signatures for the Core Density phenotype (Supplementary Figure 9b in the manuscript)

SRIScripts-

"SRI-GSEA-BT-IPA.R": Pathway significantly associated with SRI phenotype on bulktumor data using Gene Set Enrichment Analysis (GSEA) as implemented in the Piano R package
"SRI-IPA-PathwayScores.R": To obtain Metagene Signatures for the SRI phenotype (Supplementary Figure 9c in the manuscript)

KM plot-

"KMPlot-Rody.R": Prognostic value of CoreDen-MetaSig1 on the Rody dataset (GSE31519) 

