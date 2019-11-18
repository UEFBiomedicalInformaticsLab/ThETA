# ThETA
ThETA: Transcriptome-driven Efficacy estimates for gene-based TArget discovery

This R package was built to provide the implementation of novel trasncriptome-based efficacy scores of target(gene)-disease associations, which have been recentely described in [Failli et al. 2019](https://www.nature.com/articles/s41598-019-46293-7). It provides utility functions to recompile the scores based on the selection of different disease-relevant gene sets and tissue-specific gene expresstion profiles. Morevoer, it provides the users with an easy access to the disease-gene association scores compiled by the [Open Target platform](https://www.targetvalidation.org/) and functions to merge the OT-based scores with our novel efficacy scores in order to provide a final prioritization of target(gene)-disease associations. Finally, the users can run basic visualization functions in order to visualize the tissue-specific gene networks and biological annotations associated to top drug targets (or genes) and closely related genes slected with a random walk algorithm. 

## Pre-requisites

You must have the following packages installed:
* RCurl
* SPARQL
* XML
* MIGSA
* ReactomePA
* clusterProfiler
* enrichplot
* dnet
* mltools
* snow, doParallel

## Installation

You can install ThETA from GitHub by using the [devtools](https://cran.r-project.org/web/packages/devtools/index.html) package. 

```r
install_github("vittoriofortino84/ThETA/ThETA")
```

Please, open the available vignette files to learn how to use this R package.

```r
vignette("Introduction", package = "ThETA")
vignette("DataProcessing", package = "ThETA")
```

### Contact Information
Mario Failli <m.failli@tigem.it>

Vittorio Fortino <vittorio.fortino@uef.fi>

