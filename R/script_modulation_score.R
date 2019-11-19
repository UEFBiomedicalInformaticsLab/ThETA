#'Modulation Scores Compiled Between Gene Modulations And Disease Perturbations.
#'
#'Determine if a gene modulation can revert disease-induced gene expression patterns.
#'
#'This function implements the core of the modulation (efficacy) score described in \insertRef{Failli2019}{ThETA}.
#'
#'@param geneSets a list of four lists including up or down gene sets identified in disease or gene perturbations
#' (e.g., gene knockout, gene knockdown, etc.). Lists of disease and gene perturbations shall not be shuffled and
#' the same ordering criterion shall be applied (e.g. disease-up-gene-sets, disease-down-gene-sets, gene-up-gene-sets and gene-down-gene-sets).
#' Each gene set object is a list of two items:\cr
#'        - \code{geneIds}: a vector of gene IDs up- or down-regulated in diseases or gene perturbations;\cr
#'        - \code{Description}: a character specifying either the extended name of the disease or the gene modulation type.
#' If the list object is empty, then the following list of genes sets will be downloaded from EnrichR
#' (\url{https://amp.pharm.mssm.edu/Enrichr/}) by using the R (Bioconductor) package \pkg{MIGSA} \insertRef{Rodriguez2019}{ThETA}.
#'@return a data frame with target-disease associations and related modulation scores. The gene modulation type and
#'the disease name are also reported.
#'@export
#'@importFrom Rdpack reprompt
#'@importFrom MIGSA downloadEnrichrGeneSets
#'@importFrom reshape2 melt
#'@importFrom data.table as.data.table
#'@importFrom scales rescale
modulation.score <- function(geneSets = NULL) {
  if (is.null(geneSets)) {
    diff_exp_gene_sets <-
      MIGSA::downloadEnrichrGeneSets(
        c(
          'Disease_Perturbations_from_GEO_up',
          'Disease_Perturbations_from_GEO_down',
          'Single_Gene_Perturbations_from_GEO_up',
          'Single_Gene_Perturbations_from_GEO_down'
        )
      )
    geneSets <- list()
    for (i in names(diff_exp_gene_sets)) {
      x <- diff_exp_gene_sets[[i]]
      geneSets[[i]] <- list()
      is_disease <- grepl('Disease_Perturbations', i)
      for (j in 1:length(x)) {
        geneSets[[i]][[j]] <-
          list(geneIds = x[[j]]@geneIds,
               Description = x[[j]]@setName)
        temp <- strsplit(geneSets[[i]][[j]]$Description, ' ')[[1]]
        idx <- length(temp) - 4
        if (is_disease) {
          geneSets[[i]][[j]]$Description <-
            tolower(paste(temp[1:(idx - 1)], collapse = ' '))
          names(geneSets[[i]])[j] <- temp[idx]
        }
        else{
          geneSets[[i]][[j]]$Description <- paste(temp[2:idx], collapse = ' ')
          names(geneSets[[i]])[j] <- toupper(temp[1])
        }
        if (is_disease) {
          names(geneSets[[i]])[grepl("^[0-9]+$", names(geneSets[[i]]))] <-
            paste('DOID:', names(geneSets[[i]])[grepl("^[0-9]+$", names(geneSets[[i]]))], sep =
                    '')
          names(geneSets[[i]])[grepl("^C[0-9]+$", names(geneSets[[i]]))] <-
            paste('CUI:', names(geneSets[[i]])[grepl("^C[0-9]+$", names(geneSets[[i]]))], sep =
                    '')
          names(geneSets[[i]]) <- gsub('-', ':', names(geneSets[[i]]))
        }
      }
    }
  }
  prtrb_specular <-
    lapply(seq(1, 4, by = 2), function(i)
      geneSets[c(i, i + 1)])
  prtrb_specular[[2]] <- rev(prtrb_specular[[2]])
  occ. <-
    mapply(
      function(x, y)
        sapply(x, function(i)
          sapply(y, function(j)
            length(
              intersect(i$geneIds, j$geneIds)
            ))),
      prtrb_specular[[1]],
      prtrb_specular[[2]],
      SIMPLIFY = F
    )
  names(occ.) <- c('up-down', 'down-up')
  ### Calculating the composite z-scores for each disease-gene perturbation interaction
  z1 <-
    lapply(occ., function(x)
      t(apply(x, 1, function(y)
        (y - mean(y))/sd(y))))
  z2 <-
    lapply(occ., function(x)
      apply(x, 2, function(y)
        (y - mean(y))/sd(y)))
  Z <- mapply('+', z1, z2, SIMPLIFY = F)
  ### Composite score
  Z[['both']] <- (Z[[1]] + Z[[2]]) / 2
  ### Removing indexes of perturbations whose gene signatures correspond to few human orthologs
  len <-
    lapply(geneSets, function(x)
      sapply(x, function(y)
        length(y$geneIds)))
  out <-
    sapply(len, function(x) {
      q <- quantile(x)
      q[2] - 1.5 * (q[4] - q[2])
    })
  idx <- mapply(function(x, i)
    x > i, len, out, SIMPLIFY = F)
  idx <- lapply(seq(1, 4, by = 2), function(i)
    do.call('&', idx[c(i, i + 1)]))
  mat <- Z$both[idx[[2]], idx[[1]]]
  ### Aggregating target-disease pairs by max score
  df <- reshape2::melt(mat, value.name = "modscore")
  colnames(df)[1:2] <- c('target.id', 'disease.id')
  annotation <-
    expand.grid(
      target.modulationType = unlist(lapply(geneSets[[3]][idx[[2]]], '[', 2), use.names = F),
      disease.name = unlist(lapply(geneSets[[1]][idx[[1]]], '[', 2), use.names = F),
      KEEP.OUT.ATTRS = F
    )
  df <- cbind(df, annotation)
  dt <- data.table::as.data.table(df)
  dt <-
    dt[, data.table::.SD[which.max(modscore)], by = list(target.id, disease.id)]
  ### Rescaling scores to [0-1] | Setting outliers above Q3 + 1.5 IQR equal to 1
  q <- quantile(dt$modscore)
  outliers <- q[4] + (1.5 * (q[4] - q[2]))
  dt$modscore <-
    scales::rescale(dt$modscore,
                    from = c(min(dt$modscore), outliers),
                    to = c(0, 1))
  dt$modscore[dt$modscore > 1] <- 1
  ### Generating final data frame
  df <- as.data.frame(dt)
  df[] <-
    lapply(df, function(x)
      if (is.factor(x))
        as.character(x)
      else
        x)
  return(df)
}
