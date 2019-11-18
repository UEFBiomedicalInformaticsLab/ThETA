#' Node centrality scores for the tissue-specific efficacy scores
#'
#' Three node centrality scores are measured (node degree, clustering coefficient and betweenness)
#' for each tissue-specific network. The networks are built by combining tissue-specific gene expression profiles
#' (which are obtained from GTEx; \url{https://gtexportal.org/}) and the human protein-protein interaction
#' (which is obtained from StringDB; \url{https://string-db.org/}). The Borda method is used to merge the three ndoe centrality scores.
#' These scores are finally discretized by replacing the rank values with discrete levels of node centrality: 0.25, 0.5, 0.75 and 1.
#'
#' It is a list of 4 items:
#'    - \code{scores}: three centrality scores compiled for each <gene, tissue> pair (list of matrices with genes on the rows and scores on the columns);
#'    - \code{ranks}: node rankings (list of matrices with genes on the rows and ranks on the columns);
#'    - \code{borda}: the borda count merging the node rankings (a matrix with genes on the rows and the borda-based rank on the columns/tissues);
#'    - \code{borda}: a list of discretized borda counts for each tissue.
#'
#' @docType data
#'
#' @usage data(centrality_score)
#'
#' @format A list containing node centrality scores.
#'
#' @keywords datasets
#'
#' @references \insertRef{Failli2019}{ThETA}
#'
#' @examples
#' data(centrality_score)
#' names(centrality_score)
#' head(centrality_score$scores[[1]])
#' table(centrality_score$borda.disc[[1]])
"centrality_score"

#' Disease-associated genes retrieved from DIsGeNET.
#'
#' Disease-associated genes retrieved from DIsGeNET (\url{www.disgenet.org}).
#'
#' It corresponds to a list of data frames. Each data frame contains three columns reporting entrez id,
#' symbol and  score of genes significantly associated with a disease. Each entry (or disease)
#' in the list is indexed by an EFO-ID.
#'
#' @docType data
#'
#' @usage data(dis_vrnts)
#'
#' @format A list of data frames.
#'
#' @keywords datasets
#'
#' @references \insertRef{Failli2019}{ThETA}
#'
#' @examples
#' data(dis_vrnts)
"dis_vrnts"

#' Disease-associated tissue scores.
#'
#' Significance of disease modules in tissue-specific interactomes.
#'
#' Significances were obtained by implementing the algorithm described in \insertRef{Kitsak2016}{ThETA},
#' which aims to assess tissue specificity of human disease modules.
#' It is a list of 5 matrices with EFO-id on the rows and GTEx-based tissues on the clomuns:
#'        - \code{maxc}: disease maximal connected components.
#'        - \code{mer}: averages of random maximal connected components.
#'        - \code{sdr}: standard deviations of random maximal connected components.
#'        - \code{z}: significances of disease-tissue associations in the form of Z-scores.
#'        - \code{genes}: cardinality of disease-relevant genes expressed in each tissue.
#'
#' @docType data
#'
#' @usage data(disease_tissue_zscores)
#'
#' @format A list of numerical matrices.
#'
#' @keywords datasets
#'
#' @references \insertRef{Failli2019}{ThETA}
#'
#' @examples
#' data(disease_tissue_zscores)
#' names(disease_tissue_zscores)
#' head(disease_tissue_zscores$z)
"disease_tissue_zscores"

#' Tissue-specific gene expression profiles retrieved from GTEx.
#'
#' Z-score of log(TPM+1) values.
#'
#' @docType data
#'
#' @usage data(gtexv7_zscore)
#'
#' @format A numerical matrix.
#'
#' @keywords datasets
#'
#' @references \insertRef{GTEx2013}{ThETA}
#'
#' @examples
#' data(gtexv7_zscore)
#' head(gtexv7_zscore)
"gtexv7_zscore"

#' Human protein-prtein interaction network retrieved from StringDB.
#'
#' This ppi was built by considering only connections with a combined score greater or equal than 700.
#'
#' @docType data
#'
#' @usage data(ppi_strdb_700)
#'
#' @format A numerical matrix.
#'
#' @keywords datasets
#'
#' @references \insertRef{Franceschini2013}{ThETA}
#'
#' @examples
#' data(ppi_strdb_700)
#' head(ppi_strdb_700)
"ppi_strdb_700"

#' Lists of gene sets induced by disease and gene perturbations.
#'
#' These lists were retrieved from EnrichR (\url{https://amp.pharm.mssm.edu/Enrichr/}) and formatted
#' in order to contain two items:
#'        - \code{geneIds}: a vector of gene IDs up- or down-regulated in diseases or gene perturbations
#'        - \code{Description}: a character specifying either the extended name of the disease or the gene modulation type.
#'
#' @docType data
#'
#' @usage data(geo_gene_sets)
#'
#' @format A numerical matrix.
#'
#' @keywords datasets
#'
#' @references \insertRef{Kuleshov2016}{ThETA}
#'
#' @examples
#' data(geo_gene_sets)
"geo_gene_sets"

#' Lists of top 50 genes (derived from the tissue-specific score) for T2D
#'
#' @docType data
#'
#' @usage data(tss_t2d_top50)
#'
#' @format A numerical matrix.
#'
#' @keywords datasets
#'
#' @examples
#' data(tss_t2d_top50)
"tss_t2d_top50"
