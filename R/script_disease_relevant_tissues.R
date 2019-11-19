#'Compiling Z-scores From Gene-Level Expression Profiles
#'
#'Utility function to determine the Z-scores for tissue-specific gene expression profiles.
#'
#'@param tissue_expr_data a data matrix containing gene expression values.
#'@return a numeric matrix of gene expression values in the form of Z-scores. Columns are tissues and rows are diseases.
#'@export
tissue.expr <- function(tissue_expr_data){
  if(is.null(rownames(tissue_expr_data))|is.null(colnames(tissue_expr_data))){
    stop('Both colnames and rownames for tissue_expr_data must be provided!')
  }
  tissue_expr_data <- tissue_expr_data[!rowSums(tissue_expr_data)==0,]
  tissue_expr_zscore <- t(apply(tissue_expr_data, 1, function(x) (x-mean(x))/sd(x)))
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
disease.vrnts <- function(x, id_type='efo', min_score=0, curated=F){
  endpoint<-"http://rdf.disgenet.org/sparql/"
  partI<-'PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
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

  if(id_type=='mesh') url<-paste('<http://id.nlm.nih.gov/mesh/',x,'> .',sep='')
  else if(id_type=='omim') url<-paste('<http://bio2rdf.org/omim:',x,'> .',sep='')
  else if (id_type=='efo') url<-paste('<http://www.ebi.ac.uk/efo/',x,'> .',sep='')

  partII<-'?gene sio:SIO_000205 ?GeneSymbol .
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

  filter<-paste('FILTER(?score >= ',min_score,')',sep='')

  if(curated==T){
    filterII<-'FILTER (!regex(?source,"BEFREE"))'
    filter<-paste(filter,filterII,sep='\n      ')
  }

  partIII<-'} ORDER BY ?GeneSymbol'

  query<-paste(partI,url,partII,filter,partIII,sep='\n      ')

  res<-SPARQL::SPARQL(endpoint,query)$results

  if(nrow(res !=0)){
    etz <- sapply(strsplit(res$gene,'/'),function(x)gsub('>','',tail(x,1)))
    symb <- sapply(strsplit(res$GeneSymbol,'/'),function(x)gsub('>','',tail(x,1)))
    scr <- res$score
    df <- data.frame(entrez=etz,symbol=symb,score=scr,stringsAsFactors = F)
    df <- df[!duplicated(df),]
    df <- df[order(df$score,decreasing=T),]
  }

  else df <- data.frame(entrez=character(),symbol=character(),score=numeric(),stringsAsFactors = F)
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
dis.rel.tissues <- function(disease_genes, ppi_network, weighted = FALSE, tissue_expr_data, thr = 1, top = 20, rand = 100, verbose=FALSE){
  if(is.null(rownames(tissue_expr_data))|is.null(colnames(tissue_expr_data))){
    stop('Both colnames and rownames for tissue_expr_data must be provided!')
  }
  if(!is.character(ppi_network[,1])) ppi_network[,1] <- as.character(ppi_network[,1])
  if(!is.character(ppi_network[,2])) ppi_network[,2] <- as.character(ppi_network[,2])
  ppi_network <- ppi_network[!duplicated(ppi_network[,1:2]),]
  ppi_network_node<-unique(unlist(ppi_network[,1:2]))
  universe <- intersect(ppi_network_node,rownames(tissue_expr_data))
  if(length(universe)==0) stop('No corresponding IDs between ppi_network and tissue_expr_data!')
  else if(length(universe)!=length(ppi_network_node)|
          length(universe)!=nrow(tissue_expr_data)){
    if (verbose) print(paste(length(universe),'IDs in common between ppi_network and tissue_expr_data will be considered.', sep=' '))
    tissue_expr_data<-tissue_expr_data[universe,]
  }
  disease_genes_size<-length(disease_genes)
  disease_genes<-intersect(disease_genes,universe)
  if(length(disease_genes)==0) stop('No disease-associated ID match with ppi_network and tissue_expr_data!')
  else if( disease_genes_size!=length(disease_genes)){
    if(verbose) print(paste(length(disease_genes),'disease-associated IDs match with ppi_network and tissue_expr_data. Computing z-scores for disease-tissue pairs...',sep=' '))
  }
  gene_sets <- list()
  max_comps <- c()
  sd_rand <- c()
  mean_rand <- c()
  zscores <- c()
  tissue_genes <- sapply(colnames(tissue_expr_data), function(i) rownames(tissue_expr_data)[tissue_expr_data[,i]>=thr])
  for(tgs in tissue_genes) {
    # intersect the ppi network with tissue-specific genes
    df <- ppi_network[ppi_network[,1] %in% tgs & ppi_network[,2] %in% tgs,]
    g <- igraph::graph_from_data_frame(df[,1:2], directed = FALSE)
    genes <- intersect(disease_genes, igraph::V(g)$name)
    if(length(genes) > 0){
      if(top > 0) {
        # calculate the radial distance between each gene and all disease-relevant genes
        if(weighted) igraph::E(g)$weight <- scales::rescale(as.numeric(df[,3]),c(1,.Machine$double.eps))
        rad_dist <- igraph::distances(g, to=intersect(disease_genes,igraph::V(g)$name), mode =  "out", weights = NULL, algorithm = "dijkstra")
        rad_dist <- rowMeans(rad_dist)
        genes <- unique(c(genes, names(sort(rad_dist))[1:top])) # extended set of disease-relevant genes
      }
    }
    gene_sets[[length(gene_sets)+1]] <- genes
    # calculate the maximal components by considering disease-relevant genes
    if(length(genes) > 1){
      gs <- igraph::induced.subgraph(graph=g, vids=genes)
      max_comps <- c(max_comps, max(igraph::clusters(gs)$csize))
      # random
      rand_max_comps <- c()
      for(i in 1:rand) {
        gr <- igraph::induced.subgraph(graph=g, vids=sample(igraph::V(g)$name, length(genes), replace = T))
        rand_max_comps <- c(rand_max_comps, max(igraph::clusters(gr)$csize))
      }
      #hist(rand_max_comps)
      mean_rand <- c(mean_rand, mean(rand_max_comps))
      sd_rand <- c(sd_rand, sd(rand_max_comps))
      zscores <- c(zscores, (max_comps[length(max_comps)] - mean(rand_max_comps)) / sd(rand_max_comps))
    }
    else{
      max_comps <- c(max_comps, length(genes))
      mean_rand <- c(mean_rand, length(genes))
      sd_rand <- c(sd_rand, 0)
      zscores <- c(zscores, 0)
    }
  }
  final_df <- data.frame(maxc = max_comps, mer = mean_rand, sdr = sd_rand, z = zscores, genes = sapply(gene_sets, length))
  rownames(final_df) <- colnames(tissue_expr_data)
  return(final_df)
}

#'@rdname dis.rel.tissues
#'@export
#'@importFrom snow makeCluster stopCluster
#'@importFrom doParallel registerDoParallel
#'@importFrom foreach foreach %dopar%
list.dis.rel.tissues <- function(disease_gene_list, ppi_network, weighted = FALSE, tissue_expr_data, thr = 1, top = 20, rand = 100, parallel = NULL) {
  if(!is.list(disease_gene_list))stop('Argument disease_gene_list is not a list!')
  if(is.null(names(disease_gene_list)))stop('Names for disease_gene_list must be provided!')
  mat_set_tiss_zscore <- NULL
  if(!is.null(parallel)) {
    cl <- snow::makeCluster(parallel)
    doParallel::registerDoParallel(cl)
    `%dopar%` <- foreach::`%dopar%`
    dis_rlvnt_tiss <- foreach::foreach(i=disease_gene_list, .export = 'get.disease.relevant.tissues') %dopar%
      get.disease.relevant.tissues(i, ppi_network, weighted, tissue_expr_data, thr, top, rand)
    snow::stopCluster(cl)
  }
  else {
    warning("A parallel computation is highly recommended",immediate. = T)
    dis_rlvnt_tiss <- list()
    for(d in disease_gene_list) {
      dis_rlvnt_tiss[[length(dis_rlvnt_tiss)+1]] <- get.disease.relevant.tissues(d, ppi_network, weighted, tissue_expr_data, thr, top, rand)
    }
  }
  names(dis_rlvnt_tiss) <- names(disease_gene_list)
  df <- do.call(function(...)mapply(rbind,...,SIMPLIFY=F),dis_rlvnt_tiss)
  for (i in 1:length(df)) colnames(df[[i]]) <- colnames(tissue_expr_data)
  return(df)
}
