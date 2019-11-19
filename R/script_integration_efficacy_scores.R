#'Computing Overall Efficacy Scores
#'
#'Utility function to compute overall efficacy scores by using max and harmonic sum functions.
#'
#'@param data data frame reporting different efficacy estimates for the same target-disease associations.
#'@param col.scores character (or numerical) vector indicating the columns that represent efficacy scores.
#' Please, specify, at least, two different efficacy scoring methods.
#'@return input data frame with two extra columns corresponding to overall efficacy scores calculated by
#'using max and harmonic sum functions (see \insertRef{Failli2019}{ThETA} for details).
#'@export
#'@importFrom Rdpack reprompt
integrate.scores <- function(data, col.scores = NULL){
  if(is.null(col.scores) | length(col.scores) < 2)
      stop("Please indicate, at least, two different column-scores.")
  data[is.na(data)] <- 0
  print(dim(data))
  data$MAX <- rep(0, nrow(data))
  data$HS <- rep(0, nrow(data))
  for(i in 1:nrow(data)){
    data$MAX[i] <- max(data[i, col.scores])
    data$HS[i] <- sum(sort(data[i,col.scores],decreasing=T)/(1:length(data[i,col.scores]))^2)
  }
  return(data)
}

# integrate.scores<-function(scores){
#   scores[is.na(scores)] <- 0
#   cols <- colnames(scores)[3:ncol(scores)]
#   scores <- as.data.table(scores)
#   scores <- scores[,  c('Max','HS'):= list(max(.SD),sum(sort(.SD,decreasing=T)/(1:length(.SD))^2)),
#                    by=1:nrow(scores),.SDcols=cols]
#   return(as.data.frame(scores))
# }
