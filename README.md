# ThETA
ThETA: Transcriptome-driven Efficacy estimates for gene-based TArget discovery

This R package was built to provide the implementation of novel trasncriptome-based efficacy scores of target(gene)-disease associations, which have been recentely described in [Failli et al. 2019](https://www.nature.com/articles/s41598-019-46293-7). It provides utility functions to recompile the scores based on the selection of different disease-relevant gene sets and tissue-specific gene expresstion profiles. Morevoer, it provides the users with an easy access to the disease-gene association scores compiled by the [Open Target platform](https://www.targetvalidation.org/) and functions to merge the OT-based scores with our novel efficacy scores in order to provide a final prioritization of target(gene)-disease associations. Finally, the users can run basic visualization functions in order to visualize the tissue-specific gene networks and biological annotations associated to top drug targets (or genes) and closely related genes slected with a random walk algorithm. 

## Installation

As a pre-requisite, you need to have BioConductor and a set of BioConductor packages installed, you can do this with:

```r
# Install BioConductor
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()

# Install required BioConductor packages
BiocManager::install(c("MIGSA", "clusterProfiler", "enrichplot", "supraHex", "ReactomePA", "graph", "Rgraphviz", "MeSH.db", "BiocStyle"))
```

You can install ThETA from GitHub by using the [devtools](https://cran.r-project.org/web/packages/devtools/index.html) package. 

If you do not have devtools installed, you can do so by:
```r
install.packages("devtools")
```

and then to install ThETA:

```r
library(devtools)
install_github("vittoriofortino84/ThETA")
```
The above code will not build vignettes by default because they’re time consuming and may require additional packages.

To build the vignette use the following code. It will take approximately 13 minutes.  
```r
install_github("vittoriofortino84/ThETA", build_opts = c("--no-resave-data", "--no-manual"), build_vignettes = TRUE)
```
Note: the building process of the vignettes takes more than 15 minutes. This is why we provide a pre-compiled version accessible from this Git-hub page:

[Introduction](https://htmlpreview.github.io/?https://github.com/vittoriofortino84/ThETA/blob/master/Introduction.html)


```r
browseVignettes('ThETA')
vignette("Introduction", package = "ThETA")
vignette("DataProcessing", package = "ThETA")
```

### Contact Information
Mario Failli <m.failli@tigem.it>

Vittorio Fortino <vittorio.fortino@uef.fi>

