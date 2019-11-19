#'Calculate The L2 Norm
#'
#'Calculate the l2 norm.
#'
#'@keywords internal
#'
#'@param x numeric, a vector of values.
#'@param na.rm logical, whether or not to remove NA values from the calculation.
#'@return The L2 norm of \code{x}
l2norm <- function(x, na.rm = TRUE) {
  return(sqrt(mean(abs(x) ^ 2, na.rm = TRUE)))
}
#'Calculate The Geometric Mean
#'
#'Calculate the geometric mean.
#'
#'@keywords internal
#'
#'@param x numeric, a vector of values.
#'@param na.rm logical, whether or not to remove NA values from the calculation.
#'@return the geometric mean of \code{x}
geo.mean <- function(x, na.rm = TRUE) {
  return(exp(mean(log(x), na.rm = TRUE)))
}

#'Borda Based Rank Aggregation
#'
#'Computes Borda scores and rank based on four different aggregation functions.
#'
#'Computes Borda scores and rank based on four different aggregation functions. The four aggregation functions
#'are mean, median, geometric mean and L2 norm.
#'
#'@keywords internal
#'
#'@param input a list containing individual ranked lists.
#'@return a list with two components:\cr
#' -\strong{TopK} a matrix with 4 columns corresponding to rankings by each of the 4 aggregation functions;\cr
#' -\strong{Scores} a matrix with 4 columns corresponding to the Borda scores from each of the 4 aggregation functions.
Borda <- function(input) {
  aggreg.function = c("mean", "median", "geo.mean", "l2norm")
  allranks = data.frame(matrix(
    0,
    nrow = nrow(input),
    ncol = length(aggreg.function)
  ))
  names(allranks) = aggreg.function
  allfun = data.frame(matrix(
    0,
    nrow = nrow(input),
    ncol = length(aggreg.function)
  ))
  names(allfun) = aggreg.function
  for (i in 1:length(aggreg.function)) {
    allfun[, i] = apply(input, 1, aggreg.function[i],
                        na.rm = TRUE)
    allranks[, i] = rank(allfun[, i], ties.method = 'random')
  }
  results = list(allranks, allfun)
  names(results) = c("TopK", "Scores")
  return(results)
}

#'Compile Node Centrality Scores For A Given PPI Network.
#'
#'Determine three node centrality measures (node degree, clustering coefficient and betweenness) for a given ppi network.
#'
#'The three scores are merged with the Borda method. See \insertRef{Failli2019}{ThETA}, for more details.
#'
#'@param ppi_network a matrix or a data frame with at least two columns
#'reporting the ppi connections (or edges). Each line corresponds to a direct interaction.
#'Columns give the gene IDs of the two interacting proteins.
#'Only interactions whose destination nodes match with the IDs provided in \code{tissue_expr_data}
#'will be considered during the computation.
#'@param tissue_expr_data  numeric matrix or data frame indicating expression significances
#'in the form of Z-scores. Columns are tissues and rows are genes; colnames and rownames must be provided.
#'@param agg_function integer value between 1 and 4. Provides control over Borda aggregation functions.
#'1 for \code{mean}, 2 for \code{median}, 3 for \code{geo.mean} and 4 for \code{l2norm}. Defaults to 4.
#'@param parallel an integer indicating how many cores will be registered for parallel computation.
#'If NULL, the function will not use parallel computation.
#'@param verbose logical to indicate whether the messages will be displayed or not in the screen.
#'@return a list of 4 items:\cr
#'    - \strong{scores}: three centrality scores compiled for each <gene, tissue> pair (list of matrices with genes on the rows and scores on the columns);\cr
#'    - \strong{ranks}: node rankings (list of matrices with genes on the rows and ranks on the columns);\cr
#'    - \strong{borda}: the borda count merging the node rankings (a matrix with genes on the rows and the borda-based rank on the columns/tissues);\cr
#'    - \strong{borda}: a list of discretized borda counts for each tissue.
#'@export
#'@importFrom Rdpack reprompt
#'@importFrom igraph V E graph_from_edgelist as.undirected strength transitivity betweenness
#'@importFrom scales rescale
#'@importFrom snow makeCluster stopCluster
#'@importFrom doParallel registerDoParallel
#'@importFrom foreach foreach %dopar%
#'@importFrom mltools bin_data
node.centrality <-
  function(ppi_network,
           tissue_expr_data,
           agg_function = 4,
           directed_network = F,
           parallel = NULL,
           verbose = FALSE) {
    if (length(rownames(tissue_expr_data)) == 0 |
        length(colnames(tissue_expr_data)) == 0) {
      stop('Both colnames and rownames for tissue_expr_data must be provided')
    }
    for (i in 1:2)
      ppi_network[, i] <- as.character(ppi_network[, i])
    ppi_network <- ppi_network[!duplicated(ppi_network[, 1:2]), ]
    ppi_network_size <- nrow(ppi_network)
    if (directed_network)
      idx <- ppi_network[, 1] %in% rownames(tissue_expr_data)
    else
      idx <-
      ppi_network[, 1] %in% rownames(tissue_expr_data) &
      ppi_network[, 2] %in% rownames(tissue_expr_data)
    ppi_network <- ppi_network[idx, ]
    if (nrow(ppi_network) == 0)
      stop('No corresponding IDs between ppi_network and tissue_expr_data.')
    else if (ppi_network_size != nrow(ppi_network)) {
      if (verbose)
        print(paste(
          nrow(ppi_network),
          'out of',
          ppi_network_size,
          'network edges selected.',
          sep = ' '
        ))
    }
    if (!directed_network) {
      if (verbose)
        print('Undirect network. Duplicating network edges...')
      ppi_network_rev <- ppi_network[, c(2:1)]
      colnames(ppi_network_rev) <- colnames(ppi_network)[1:2]
      ppi_network <- rbind(ppi_network[, 1:2], ppi_network_rev)
    }
    # rescaling z-scores into positive values
    tissue_expr_scaled <-
      scales::rescale(tissue_expr_data, c(.Machine$double.eps, 1))
    tissue_expr_scaled_btw <-
      scales::rescale(tissue_expr_data, c(1, .Machine$double.eps))
    # building a graph from an edge list
    g = igraph::graph_from_edgelist(as.matrix(ppi_network[, 1:2]), directed =
                                      T)
    # computing node degree, clustering coefficient and betweenness per tissue
    node_degree <- c()
    node_cc <- c()
    node_btw <- c()
    if (verbose)
      print("Compiling node centrality scores...")
    if (is.null(parallel)) {
      node_degree <- lapply(colnames(tissue_expr_data), function(i) {
        igraph::E(g)$weight <- tissue_expr_scaled[ppi_network[, 1], i]
        igraph::strength(g, mode = 'in')
      })
      node_cc <- lapply(colnames(tissue_expr_data), function(i) {
        igraph::E(g)$weight <- tissue_expr_scaled[ppi_network[, 1], i]
        igraph::transitivity(igraph::as.undirected(g),
                             type = 'barrat',
                             isolates = 'zero')
      })
      node_btw <- lapply(colnames(tissue_expr_data), function(i) {
        igraph::E(g)$weight <- tissue_expr_scaled_btw[ppi_network[, 1], i]
        igraph::betweenness(g)
      })
    }
    else {
      cl <- snow::makeCluster(parallel)
      doParallel::registerDoParallel(cl)
      `%dopar%` <- foreach::`%dopar%`
      node_degree <-
        foreach::foreach(i = colnames(tissue_expr_data)) %dopar% {
          igraph::E(g)$weight <- tissue_expr_scaled[ppi_network[, 1], i]
          igraph::strength(g, mode = 'in')
        }
      names(node_degree) <- colnames(tissue_expr_data)
      node_cc <- foreach(i = colnames(tissue_expr_data)) %dopar% {
        igraph::E(g)$weight <- tissue_expr_scaled[ppi_network[, 1], i]
        igraph::transitivity(igraph::as.undirected(g),
                             type = 'barrat',
                             isolates = 'zero')
      }
      names(node_cc) <- colnames(tissue_expr_data)
      node_btw <- foreach(i = colnames(tissue_expr_data)) %dopar% {
        igraph::E(g)$weight <- tissue_expr_scaled_btw[ppi_network[, 1], i]
        igraph::betweenness(g)
      }
      names(node_btw) <- colnames(tissue_expr_data)
      snow::stopCluster(cl)
    }
    if (verbose)
      print("Merging different node centrality scores...")
    # aggregating individual centrality rankings with Borda (l2norm aggregation function)
    node_cent_scores <-
      mapply(function(x, y, z)
        cbind(degree = x, cc = y, btw = z),
        node_degree,
        node_cc,
        node_btw,
        SIMPLIFY = F)
    node_cent_ranks <-
      lapply(node_cent_scores, function(x)
        apply(x, 2, rank, ties.method = 'random'))
    borda_var <- lapply(node_cent_ranks, Borda)
    rank_agg_function <-
      sapply(borda_var, function(x)
        x$TopK[, agg_function])
    rownames(rank_agg_function) <- igraph::V(g)$name
    if (verbose)
      print("Discretizing the node centrality score...")
    # assigning each node a discrete value between 0.25 and 1 depending on its centrality ranking quartile
    chunk <- as.data.frame(apply(rank_agg_function, 2, function(x)
      mltools::bin_data(x, bins = quantile(x), boundaryType = "[lorc")))
    W <-
      lapply(chunk, function(x)
        as.numeric(factor(x, levels = rev(levels(
          x
        )))) / 4)
    for (i in 1:length(W))
      names(W[[i]]) <- rownames(node_cent_scores[[1]])
    return(
      list(
        scores = node_cent_scores,
        ranks = node_cent_ranks,
        borda = rank_agg_function,
        borda.disc = W
      )
    )
  }
