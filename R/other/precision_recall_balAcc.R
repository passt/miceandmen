precall <- function(fmat) {
# #
# function to calculate balanced accuracy, precision and recall
# #
  # fmat: frequency table of predicted (rows) vs true (columns) identity
  
  nr <- nrow(fmat)
  nom <- colnames(fmat)
  
  # balanced accuracy: (tpr + tnr) / 2
  balacc    <- sapply(nom, FUN=function(i) { ni <- nom != i;
                                            (fmat[i,i]/sum(fmat[i,]) + sum(fmat[ni,ni])/sum(fmat[ni,])) / 2 } )
  
  # precision: fraction of true positives per sum of all positives: tp / (tp + fp)                                                                
  precision <- fmat %>% normalizeMat(1) %>% diag

  # recall: fraction of true positives per sum of condition positives: tp / (tp + fn)
  recall    <- fmat %>% normalizeMat(2) %>% diag
  
  # f1: combination of precision and recall. 
  f1        <- 2 * precision * recall / (precision + recall)
  
  return(list("balanced Accuracy" = balacc, "precision" = precision, "recall" = recall, 'F1'=f1))
}


normalizeMat <- function(x, dim) {
# # 
# Function to normalize matrix by row or column sum.
# #
  if (dim == 1) { # re-scale by sum over rows - diagonal entries are 'precision'
    x / rowSums(x) # divide by vector ()
  } else if (dim == 2) { # re-scale by sum over rows - diagonal entries are 'recall'
    x / colSums(x)[col(x)] # divide by 
  } else { stop('Invalid dimension') }
}