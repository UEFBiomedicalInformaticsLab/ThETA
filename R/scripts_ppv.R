#' Plotting function showing positive predictive values. 
#'
#' Given a set of knonw disease-drug-target(gene) associations and efficacy estimates compiled for diseas-gene associations,
#' it determines an visualize the positive predictive values (ppv). PPVs are computed as the percentage of true positives 
#' normalized against the expected ppv for a random ordering of disease-gene pairs. 
#'
#' This function was used for the comparison of different efficacy scoring methods described in \insertRef{Failli2019}{ThETA}.
#'
#'@param gold_std a dataframe representing known drug target disease associations. 
#'It must contain at least two columns indicating the gene-entrez-id and disease-efo-id (EFO-ID). 
#'Ex:
#'target.entrez      disease.id                        description
#'         7039     EFO:0003884             chronic kidney disease
#'         1490     EFO:0009549 female reproductive system disease
#'         
#'@param test a dataframe containing efficacy estimates of disease-drug-target(gene) associations.
#'It must contain at least three columns indicating the gene-entrez-id, disease-efo-id (EFO-ID) and the efficacy scores. 
#'Ex:
#'target.entrez   disease.id  modscore
#'          5376 EFO:0000270 0.4614037
#'         28996 EFO:0000270 0.5404949
#'          7403 EFO:0000270 0.6425815
#'
#'@param entrez.col a character value indicating the column for the gene-entrez-id.  
#'@param disease.col a character value indicating the column for the disease-efo-id.   
#'@param score.col a character value indicating the column for the efficacy estimates.   
#'@param a.disease a character value to plot the ppv values for a given disease. 
#'@param plot a boolean value to plot the ppv curves. 
#'@return a ggplot showing curves of positive predictive values (PPVs).
#'@export
#'@import ggplot2
#'@importFrom reshape2 melt
ppvpercents <- function(gold_std, test, 
                        entrez.col = NULL, 
                        disease.col = NULL, 
                        score.col = NULL, 
                        a.disease = NULL,
                        plot = TRUE){
  if(is.null(entrez.col) | is.null(disease.col) | is.null(score.col))
    stop('Incomplete data input.')
  if(!is.null(a.disease)) {
    idd = which(test[,disease.col] == a.disease)
    if(length(idd) == 0) stop('Disease not included.')
    test = test[idd,]
  }
  n_score <- length(score.col)
  ppv <- vector(mode = 'list',length = n_score)
  names(ppv) <- score.col
  if(is.data.frame(gold_std)) gold_std <- list(gold_std)
  for(i in seq_len(n_score)){
    new_test <- test[!is.na(test[,score.col[i]]), 
                     c(entrez.col,disease.col,score.col)]
    new_test <- new_test[order(new_test[,score.col[i]], decreasing = T),
                         c(entrez.col,disease.col)]
    ppv [[i]] <- vector(mode = 'list',length = length(gold_std))
    if(!is.null(names(gold_std))) names(ppv[[i]]) <- names(gold_std)
    for(j in 1:length(gold_std)){
      if(i==1)
        gold_std[[j]] <- gold_std[[j]][!duplicated(gold_std[[j]]),]
      new_test <- new_test[new_test[,disease.col] %in% gold_std[[j]][,disease.col],]
      ind <- matching.rows(gold_std[[j]][,c(entrez.col,disease.col)], new_test, nomatch = NA)
      ind <- sort(ind[!is.na(ind)])
      nconn <- nrow(new_test)
      ppv[[i]][[j]] <- vector(mode='integer',length = 100)
      iters <- seq(1,100,1)
      for (k in 1:length(iters)){
        perc <- round(nconn/100 * iters[k])
        if(perc ==0) perc <- 1
        ppv[[i]][[j]][k] <- sum(ind<perc)/perc/(length(ind)/nconn)
      }
    }
  }
  if(plot){
    ppv <- reshape2::melt(ppv)
    p <- ggplot2::ggplot(data=ppv,
                         ggplot2::aes(x=rep(1:100,nrow(ppv)/100), y=value, colour= factor(L2))) +
                         ggplot2::geom_hline(yintercept = 1) + 
                         ggplot2::xlab('Percentile (%)') + 
                         ggplot2::ylab('Normalized PPV') +
                         ggplot2::geom_line() + 
                         ggplot2::facet_wrap(~L1, scales = "free") + 
                         ggplot2::theme(legend.title = element_blank())
  }
  else return(ppv)
} 

#'A function for matching rows between two matrices or data.frames. 
#'@keywords internal
matching.rows <- function (dat, tab, nomatch = NA)  {
  if (class(tab) == "matridat") 
    tab <- as.data.frame(tab)
  if (is.null(dim(dat))) 
    dat <- as.data.frame(matridat(dat, nrow = 1))
  cx <- do.call("paste", c(dat[, , drop = FALSE], sep = "\r"))
  ct <- do.call("paste", c(tab[, , drop = FALSE], sep = "\r"))
  match(cx, ct, nomatch = nomatch)
}
