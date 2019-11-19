#'Compiling Z-scores From Gene-Level Expression Profiles
#'
#'Utility function to determine the Z-scores for tissue-specific gene expression profiles.
#'
#'@param tissue_expr_data a data matrix containing gene expression values.
#'@return a numeric matrix of gene expression values in the form of Z-scores. Columns are tissues and rows are diseases.
#'@export
tissue.expr <- function(tissue_expr_data) {
  if (is.null(rownames(tissue_expr_data)) |
      is.null(colnames(tissue_expr_data))) {
    stop('Both colnames and rownames for tissue_expr_data must be provided!')
  }
  tissue_expr_data <-
    tissue_expr_data[!rowSums(tissue_expr_data) == 0, ]
  tissue_expr_zscore <-
    t(apply(tissue_expr_data, 1, function(x)
      (x - mean(x)) / sd(x)))
  return(tissue_expr_zscore)
}

#'Retrieving Gene-Disease Assocations From DisGeNET
#'
#'Utility function to extract gene-disease association from DisGeNET \url{http://www.disgenet.org/}.
#'
#'@param x character indicating the disease ID (eg. efo ID).
#'@param id_type character indicating the type of disease ID. Possible values are:
#'\code{efo}, \code{mesh} and \code{omim}. Defaults to \code{efo}.
#'@param min_score numeric value from 0 to 1 indicating the cut-off for gene-disease scores. Defaults to 0.
#'@param curated a boolean value indicating if manually curated gene-disease pairs must be selected.
#'@return a three columns data frame:\cr
#'        - \strong{entrez}: a character variable containing entrez gene ids;\cr
#'        - \strong{symbol}: a character variable containing gene symbols;\cr
#'        - \strong{score}: a numeric variable containing gene-disease scores.
#'@export
#'@import RCurl
#'@import XML
#'@import SPARQL
disease.vrnts <- function(x,
                          id_type = 'efo',
                          min_score = 0,
                          curated = F) {
  endpoint <- "http://rdf.disgenet.org/sparql/"
  partI <- 'PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
  PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
  PREFIX owl: <http://www.w3.org/2002/07/owl#>
  PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>
  PREFIX dcterms: <http://purl.org/dc/terms/>
  PREFIX foaf: <http://xmlns.com/foaf/0.1/>
  PREFIX skos: <http://www.w3.org/2004/02/skos/core#>
  PREFIX void: <http://rdfs.org/ns/void#>
  PREFIX sio: <http://semanticscience.org/resource/>
  PREFIX ncit: <http://ncicb.nci.nih.gov/xml/owl/EVS/Thesaurus.owl#>
  PREFIX up: <http://purl.uniprot.org/core/>
  PREFIX dcat: <http://www.w3.org/ns/dcat#>
  PREFIX dctypes: <http://purl.org/dc/dcmitype/>
  PREFIX wi: <http://http://purl.org/ontology/wi/core#>
  PREFIX eco: <http://http://purl.obolibrary.org/obo/eco.owl#>
  PREFIX prov: <http://http://http://www.w3.org/ns/prov#>
  PREFIX pav: <http://http://http://purl.org/pav/>
  PREFIX obo: <http://purl.obolibrary.org/obo/>
  SELECT DISTINCT str(?DiseaseName) as ?DiseaseName
  ?gene
  ?GeneSymbol
  ?score
  WHERE {?disease skos:exactMatch '

  if (id_type == 'mesh')
    url <- paste('<http://id.nlm.nih.gov/mesh/', x, '> .', sep = '')
  else if (id_type == 'omim')
    url <- paste('<http://bio2rdf.org/omim:', x, '> .', sep = '')
  else if (id_type == 'efo')
    url <- paste('<http://www.ebi.ac.uk/efo/', x, '> .', sep = '')

  partII <- '?gene sio:SIO_000205 ?GeneSymbol .
  ?scoreIRI sio:SIO_000300 ?score .
  ?variant so:associated_with ?gene .
  ?vda sio:SIO_000628 ?variant,?disease ;
  sio:SIO_000253 ?source ;
  sio:SIO_000216 ?scoreIRI .
  ?gene rdf:type ncit:C16612 .
  ?GeneSymbol rdf:type ncit:C43568 .
  ?disease rdf:type ncit:C7057 ; dcterms:title ?DiseaseName .
  ?scoreIRI rdf:type ncit:C25338 .
  ?variant a so:0001060 .'

  filter <- paste('FILTER(?score >= ', min_score, ')', sep = '')

  if (curated == T) {
    filterII <- 'FILTER (!regex(?source,"BEFREE"))'
    filter <- paste(filter, filterII, sep = '\n      ')
  }

  partIII <- '} ORDER BY ?GeneSymbol'

  query <- paste(partI, url, partII, filter, partIII, sep = '\n      ')

  res <- SPARQL::SPARQL(endpoint, query)$results

  if (nrow(res != 0)) {
    etz <-
      sapply(strsplit(res$gene, '/'), function(x)
        gsub('>', '', tail(x, 1)))
    symb <-
      sapply(strsplit(res$GeneSymbol, '/'), function(x)
        gsub('>', '', tail(x, 1)))
    scr <- res$score
    df <-
      data.frame(
        entrez = etz,
        symbol = symb,
        score = scr,
        stringsAsFactors = F
      )
    df <- df[!duplicated(df), ]
    df <- df[order(df$score, decreasing = T), ]
  }

  else
    df <-
    data.frame(
      entrez = character(),
      symbol = character(),
      score = numeric(),
      stringsAsFactors = F
    )
  return(df)
}

#'Significance Of Disease Modules In Tissue-Specific Interactomes
#'
#'Compile significance of disease modules in tissue-specific interactomes.
#'
#'This function implements the method from \insertRef{Kitsak2016}{ThETA}.
#'The function tests the significance of disease modules in tissue-specific interactomes, consisting
#'of genes with expression significance >= \code{thr}. First, disease maximal connected
#'component is measured in tissue \emph{t}. Then, the average of \code{rand} maximal connected
#'components in \emph{t} is determined to compute the significance of real data.
#'As a major improvement, the radial distance between each gene and all disease-relevant genes is computed
#'to extend the disease gene set.
#'
#'@param disease_genes character vector containing the IDs of the genes related to a particular disease.
#'Gene IDs are expected to match with those provided in \code{ppi_network} and \code{tissue_expr_data}.
#'@param disease_gene_list a list of disease-associated genes. Each element of the list is a character vector
#'containing the IDs of the genes related to a particular disease.
#'@param ppi_network a matrix or a data frame with at least two columns
#'reporting the ppi connections (or edges). Each line corresponds to a direct interaction.
#'Columns give the gene IDs of the two interacting proteins.
#'@param weighted logical indicating whether the \code{ppi_network} contain edge weights? If \code{TRUE}, then the third column
#'is used to weight edges.
#'@param tissue_expr_data a numeric matrix or data frame indicating expression significances
#'in the form of Z-scores. Columns are tissues and rows are genes; colnames and rownames must be provided.
#'Gene IDs are expected to match with those provided in \code{ppi_network}.
#'@param thr an integer specifying the threshold for expression significance. Only genes with
#'expression significance >= \code{thr} will be considered expressed in a particular tissue.
#'Defaults to 1.
#'@param top an integer indicating the number of genes to be considered to extend the disease gene set.
#'Genes are selected based on the radial distance from the disease-relevant genes.
#'@param rand an integer indicating the number of random sets of genes needed to calculate the singificance of disease-tissue associations. Defaults to 100.
#'@param verbose logical indicating whether the messages will be displayed or not in the screen.
#'@param parallel an integer indicating how many cores will be registered for parallel computation.
#'If \code{NULL}, the function will not use parallel computation.
#'@return A list of five items:\cr
#'        - \strong{maxc}: a vector or a matrix specifying the disease maximal connected components by tissue;\cr
#'        - \strong{mer}: a vector or a matrix specifying the averages of \code{rand} maximal connected components by tissue;\cr
#'        - \strong{sdr}: a vector or a matrix specifying the standard deviations of \code{rand} maximal connected components by tissue;\cr
#'        - \strong{z}: a vector or a matrix specifying the significances of disease-tissue associations in the form of Z-scores;\cr
#'        - \strong{genes}: a vector or a matrix specifying the cardinality of disease-relevant genes expressed in each tissue.
#'@export
#'@importFrom Rdpack reprompt
#'@importFrom igraph V E graph_from_data_frame distances induced.subgraph clusters
#'@importFrom scales rescale
dis.rel.tissues <-
  function(disease_genes,
           ppi_network,
           weighted = FALSE,
           tissue_expr_data,
           thr = 1,
           top = 20,
           rand = 100,
           verbose = FALSE) {
    if (is.null(rownames(tissue_expr_data)) |
        is.null(colnames(tissue_expr_data))) {
      stop('Both colnames and rownames for tissue_expr_data must be provided!')
    }
    if (!is.character(ppi_network[, 1]))
      ppi_network[, 1] <- as.character(ppi_network[, 1])
    if (!is.character(ppi_network[, 2]))
      ppi_network[, 2] <- as.character(ppi_network[, 2])
    ppi_network <- ppi_network[!duplicated(ppi_network[, 1:2]), ]
    ppi_network_node <- unique(unlist(ppi_network[, 1:2]))
    universe <- intersect(ppi_network_node, rownames(tissue_expr_data))
    if (length(universe) == 0)
      stop('No corresponding IDs between ppi_network and tissue_expr_data!')
    else if (length(universe) != length(ppi_network_node) |
             length(universe) != nrow(tissue_expr_data)) {
      if (verbose)
        print(
          paste(
            length(universe),
            'IDs in common between ppi_network and tissue_expr_data will be considered.',
            sep = ' '
          )
        )
      tissue_expr_data <- tissue_expr_data[universe, ]
    }
    disease_genes_size <- length(disease_genes)
    disease_genes <- intersect(disease_genes, universe)
    if (length(disease_genes) == 0)
      stop('No disease-associated ID match with ppi_network and tissue_expr_data!')
    else if (disease_genes_size != length(disease_genes)) {
      if (verbose)
        print(
          paste(
            length(disease_genes),
            'disease-associated IDs match with ppi_network and tissue_expr_data. Computing z-scores for disease-tissue pairs...',
            sep = ' '
          )
        )
    }
    gene_sets <- list()
    max_comps <- c()
    sd_rand <- c()
    mean_rand <- c()
    zscores <- c()
    tissue_genes <-
      sapply(colnames(tissue_expr_data), function(i)
        rownames(tissue_expr_data)[tissue_expr_data[, i] >= thr])
    for (tgs in tissue_genes) {
      # intersect the ppi network with tissue-specific genes
      df <-
        ppi_network[ppi_network[, 1] %in% tgs & ppi_network[, 2] %in% tgs, ]
      g <- igraph::graph_from_data_frame(df[, 1:2], directed = FALSE)
      genes <- intersect(disease_genes, igraph::V(g)$name)
      if (length(genes) > 0) {
        if (top > 0) {
          # calculate the radial distance between each gene and all disease-relevant genes
          if (weighted)
            igraph::E(g)$weight <-
              scales::rescale(as.numeric(df[, 3]), c(1, .Machine$double.eps))
          rad_dist <-
            igraph::distances(
              g,
              to = intersect(disease_genes, igraph::V(g)$name),
              mode =  "out",
              weights = NULL,
              algorithm = "dijkstra"
            )
          rad_dist <- rowMeans(rad_dist)
          genes <-
            unique(c(genes, names(sort(rad_dist))[1:top])) # extended set of disease-relevant genes
        }
      }
      gene_sets[[length(gene_sets) + 1]] <- genes
      # calculate the maximal components by considering disease-relevant genes
      if (length(genes) > 1) {
        gs <- igraph::induced.subgraph(graph = g, vids = genes)
        max_comps <- c(max_comps, max(igraph::clusters(gs)$csize))
        # random
        rand_max_comps <- c()
        for (i in 1:rand) {
          gr <-
            igraph::induced.subgraph(graph = g,
                                     vids = sample(igraph::V(g)$name, length(genes), replace = T))
          rand_max_comps <-
            c(rand_max_comps, max(igraph::clusters(gr)$csize))
        }
        #hist(rand_max_comps)
        mean_rand <- c(mean_rand, mean(rand_max_comps))
        sd_rand <- c(sd_rand, sd(rand_max_comps))
        zscores <-
          c(zscores, (max_comps[length(max_comps)] - mean(rand_max_comps)) / sd(rand_max_comps))
      }
      else{
        max_comps <- c(max_comps, length(genes))
        mean_rand <- c(mean_rand, length(genes))
        sd_rand <- c(sd_rand, 0)
        zscores <- c(zscores, 0)
      }
    }
    final_df <-
      data.frame(
        maxc = max_comps,
        mer = mean_rand,
        sdr = sd_rand,
        z = zscores,
        genes = sapply(gene_sets, length)
      )
    rownames(final_df) <- colnames(tissue_expr_data)
    return(final_df)
  }

#'@rdname dis.rel.tissues
#'@export
#'@importFrom snow makeCluster stopCluster
#'@importFrom doParallel registerDoParallel
#'@importFrom foreach foreach %dopar%
list.dis.rel.tissues <-
  function(disease_gene_list,
           ppi_network,
           weighted = FALSE,
           tissue_expr_data,
           thr = 1,
           top = 20,
           rand = 100,
           parallel = NULL) {
    if (!is.list(disease_gene_list))
      stop('Argument disease_gene_list is not a list!')
    if (is.null(names(disease_gene_list)))
      stop('Names for disease_gene_list must be provided!')
    mat_set_tiss_zscore <- NULL
    if (!is.null(parallel)) {
      cl <- snow::makeCluster(parallel)
      doParallel::registerDoParallel(cl)
      `%dopar%` <- foreach::`%dopar%`
      dis_rlvnt_tiss <-
        foreach::foreach(i = disease_gene_list, .export = 'get.disease.relevant.tissues') %dopar%
        get.disease.relevant.tissues(i, ppi_network, weighted, tissue_expr_data, thr, top, rand)
      snow::stopCluster(cl)
    }
    else {
      warning("A parallel computation is highly recommended", immediate. = T)
      dis_rlvnt_tiss <- list()
      for (d in disease_gene_list) {
        dis_rlvnt_tiss[[length(dis_rlvnt_tiss) + 1]] <-
          get.disease.relevant.tissues(d, ppi_network, weighted, tissue_expr_data, thr, top, rand)
      }
    }
    names(dis_rlvnt_tiss) <- names(disease_gene_list)
    df <-
      do.call(function(...)
        mapply(rbind, ..., SIMPLIFY = F), dis_rlvnt_tiss)
    for (i in 1:length(df))
      colnames(df[[i]]) <- colnames(tissue_expr_data)
    return(df)
  }


#'Interactive visualization of tissue-specific networks.
#'
#'It uses visNetwork to build a neteworks showing all short pathways connecting disease genes with putative drug targets.
#'
#'@param tissue_scores a data.frame as the one compiled by \code{get.tissue.specific.scores}
#'@param disease_genes character vector containing the IDs of the genes related to a particular disease.
#'Gene IDs are expected to match with those provided in \code{ppi_network} and \code{tissue_expr_data}.
#'@param ppi_network a matrix or a data frame with at least two columns
#'reporting the ppi connections (or edges). Each line corresponds to a direct interaction.
#'Columns give the gene IDs of the two interacting proteins.
#'@param directed_network logical indicating whether the PPI is directed.
#'@param tissue_expr_data a numeric matrix or data frame indicating expression significances
#'in the form of Z-scores. Columns are tissues and rows are genes; colnames and rownames must be provided.
#'Gene IDs are expected to match with those provided in \code{ppi_network}.
#'@param top_targets character vector indicating a list of ENTREZ id to be used for the slection of the shortest paths.
#'@param db character indicating the database to consider for enrichment analysis.
#'Possible values are: \code{kegg}, \code{BP}, \code{MF} and \code{CC}. Defaults to \code{kegg}.
#'@param verbose logical indicating whether the messages will be displayed or not in the screen.
#'@export
#'@importFrom clusterProfiler enrichKEGG
#'@importFrom clusterProfiler enrichGO
visualize.graph <-
  function(tissue_scores,
           disease_genes,
           ppi_network,
           directed_network = FALSE,
           tissue_expr_data,
           top_targets = NULL,
           db = 'kegg',
           verbose = FALSE) {
    if (is.null(top_targets))
      stop('Please specifiy a set of targets (ENTREZ ids)!')
    if (is.null(rownames(tissue_expr_data)) |
        is.null(colnames(tissue_expr_data))) {
      stop('Both colnames and rownames for tissue_expr_data must be provided.')
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
        print('Undirect network. Converting to direct network...')
      ppi_network_rev <- ppi_network[, c(2:1)]
      colnames(ppi_network_rev) <- colnames(ppi_network)[1:2]
      ppi_network <- rbind(ppi_network[, 1:2], ppi_network_rev)
    }
    disease_genes_size <- length(disease_genes)
    disease_genes <- intersect(disease_genes, ppi_network[, 2])
    if (length(disease_genes) == 0)
      stop('No disease-associated ID match with ppi_network and tissue_expr_data!')
    else if (disease_genes_size != length(disease_genes)) {
      if (verbose)
        print(paste(
          length(disease_genes),
          'disease-associated IDs are reachable from the network.',
          sep = ' '
        ))
    }
    if (!(db %in% c('kegg', 'BP', 'MF', 'CC')))
      stop('Unused value for argument db')
    idx <- order(tissue_scores$avg_tissue_score, decreasing = T)
    tissue_scores <- tissue_scores[idx, ]
    g <-
      igraph::graph_from_edgelist(as.matrix(ppi_network[, 1:2]), directed = T)
    tissue_expr_data <-
      scales::rescale(tissue_expr_data, c(1, .Machine$double.eps))
    sign_tiss <- colnames(tissue_scores)[-ncol(tissue_scores)]
    sh.path <- lapply(sign_tiss, function(i) {
      igraph::E(g)$weight <- tissue_expr_data[ppi_network[, 1], i]
      lapply(top_targets, function(x)
        igraph::shortest_paths(
          g,
          from = x,
          to = disease_genes,
          mode =  "out",
          weights = NULL,
          output = 'epath'
        )$epath)
    })
    names(sh.path) <- sign_tiss
    for (i in names(sh.path))
      names(sh.path[[i]]) <- top_targets
    all_target_path <-
      lapply(sh.path, function(x)
        lapply(x, function(y)
          do.call(c, y)))
    all_tissue_path <- lapply(all_target_path, function(x)
      do.call(c, x))
    ids <- lapply(all_tissue_path, igraph::as_ids)
    shinyApp(
      ui = fluidPage(fluidRow(
        column(
          width = 2,
          checkboxGroupInput("tissue", "Select tissue/s", sign_tiss, selected = sign_tiss[1])
        ),
        column(width = 12,
               visNetworkOutput("network"))
      ),
      fluidRow(column(
        width = 12,
        tableOutput("table")
      ))),
      server = function(input, output) {
        output$network <- renderVisNetwork({
          validate(need(input$tissue != '', 'Please select one or more tissues'))
          if (length(input$tissue) > 1) {
            edge <-
              strsplit(Reduce(intersect, ids[input$tissue]),
                       split = '|',
                       fixed = T)
          }
          else if (length(input$tissue) == 1) {
            edge <- strsplit(ids[[input$tissue]], split = '|', fixed = T)
          }
          edge <- as.data.frame(do.call(rbind, edge))
          edge <- edge[!duplicated(edge), ]
          colnames(edge) <- c('from', 'to')
          edge$width <-
            scales::rescale(rowMeans(1 - tissue_expr_data[edge$from, input$tissue, drop =
                                                            F]), c(1, 5))
          nb <- unique(unlist(edge[, 1:2]))
          node <- data.frame(id = nb,
                             label = nb,
                             stringsAsFactors = F)
          node$value <-
            rowMeans(tissue_scores[nb, input$tissue, drop = F])
          gp <- rep('bridge gene', nrow(node))
          gp[node$label %in% top_targets] <- 'target gene'
          gp[node$label %in% disease_genes] <- 'disease gene'
          node$group <- gp
          node$shape <-
            c("circle", "star", "diamond")[as.numeric(as.factor(gp))]
          visNetwork::visNetwork(node, edge, height = "1500px", width = "500%") %>%
            visGroups(
              groupname = "target gene",
              color = list(background = "skyblue", border = "deepskyblue")
            ) %>%
            visGroups(
              groupname = "disease gene",
              color = list(background = "lightcoral", border = "red")
            ) %>%
            visGroups(groupname = "bridge gene", shape = "circle") %>%
            visLegend(
              addNodes = list(
                list(
                  label = "disease gene",
                  shape = "star",
                  color = list(background = "lightcoral", border = "red")
                ),
                list(
                  label = "target gene",
                  shape = "diamond",
                  color = list(background = "skyblue", border = "deepskyblue")
                ),
                list(label = "bridge gene", shape = "circle")
              ),
              useGroups = FALSE,
              width = 0.1
            ) %>%
            visLayout(hierarchical = TRUE) %>%  # visLayout(randomSeed = 123)
            visEvents(select = "function(nodes) {
                    Shiny.onInputChange('current_node_id', nodes.nodes);
                    ;}")
        })
        output$table <- renderTable({
          validate(
            need(input$tissue != '', NULL),
            need(
              input$current_node_id != '' &
                input$current_node_id %in% top_targets,
              'Please select a target gene'
            )
          )
          current_target_path <-
            sapply(all_target_path[input$tissue], function(x)
              x[input$current_node_id])
          current_target_path <-
            Reduce(igraph::intersection, current_target_path)
          target_interactors <-
            unique(c(igraph::ends(g, current_target_path)))
          if (db == 'kegg')
            res <-
            clusterProfiler::enrichKEGG(target_interactors, organism = 'hsa')@result
          else
            res <-
            clusterProfiler::enrichGO(target_interactors, 'org.Hs.eg.db', ont = db)@result
          res <- res[res$p.adjust < 0.05, ]
          if (nrow(res) > 25)
            res <- res[1:25, ]
          res
        })
      }
    )
  }
