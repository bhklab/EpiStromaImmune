# EpiStromaImmune

Abstract
--------


Citation
--------

To be publsihed.

Full Reproducibility of the Analysis Results
--------------------------------------------

We describe below how to fully reproduce the figures and tables reported in the paper

1.  Set up the software environment

2.  Run the R scripts

3.  Generate figures

Set up the software environment (needs to be updated)
-------------------------------

We developed and tested our analysis pipeline using R running on linux and Mac OSX platforms. The following is a copy of `sessionInfo()` from the development environment in R

```
R version 3.3.0 (2016-05-03)
Platform: x86_64-redhat-linux-gnu (64-bit)
Running under: CentOS Linux 7 (Core)

attached base packages:
[1] grid      stats     graphics  grDevices utils     datasets 
[7] methods   base  

other attached packages:
 [1] pvclust_2.0-0       vegan_2.4-0         lattice_0.20-33    
 [4] permute_0.9-0       survcomp_1.22.0     prodlim_1.5.7      
 [7] survival_2.39-4     piano_1.12.0        xlsx_0.5.7         
[10] xlsxjars_0.6.1      rJava_0.9-8         RColorBrewer_1.1-2 
[13] gplots_3.0.1        mgcv_1.8-12         nlme_3.1-128       
[16] ggplot2_2.1.0       reshape2_1.4.1      VennDiagram_1.6.17 
[19] futile.logger_1.4.1 PharmacoGx_1.1.6  

loaded via a namespace (and not attached):
 [1] gtools_3.5.0         lsa_0.73.1           slam_0.1-35         
 [4] sets_1.0-16          splines_3.3.0        colorspace_1.2-6    
 [7] SnowballC_0.5.1      marray_1.50.0        sm_2.2-5.4          
[10] magicaxis_1.9.4      BiocGenerics_0.18.0  lambda.r_1.1.7      
[13] plyr_1.8.4           lava_1.4.3           stringr_1.0.0       
[16] munsell_0.4.3        survivalROC_1.0.3    gtable_0.2.0        
[19] caTools_1.17.1       labeling_0.3         Biobase_2.32.0      
[22] parallel_3.3.0       Rcpp_0.12.5          KernSmooth_2.23-15  
[25] relations_0.6-6      scales_0.4.0         limma_3.28.11       
[28] gdata_2.17.0         rmeta_2.16           plotrix_3.6-2       
[31] bootstrap_2015.2     digest_0.6.9         stringi_1.1.1       
[34] SuppDists_1.1-9.2    tools_3.3.0          bitops_1.0-6        
[37] magrittr_1.5         cluster_2.0.4        futile.options_1.0.0
[40] MASS_7.3-45          Matrix_1.2-6         downloader_0.4      
[43] igraph_1.0.1 
```

All these packages are available on [CRAN](http://cran.r-project.org) or [Bioconductor](http://www.bioconductor.org)

all necessary packages have `library(<package>)` calls within the R scripts themselves, or the script assumes a previous script has been run and thus should have loaded nessesary packages. 

Running R Scripts
-------------------------------

