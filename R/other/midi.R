
midi <- function(x, p) {
  # discretize p (low:mid:high)
  p <- p %>% cut(breaks = (0:3)/3, include.lowest = TRUE)
  
  # joint probability mass tabulate instances and normalise by total number of elements
  crt <-table(x, p) / length(x)
  
  # calculate denominator for MI
  dno <- matrix(rowSums(crt), nrow=nrow(crt), ncol=1) %*%
         matrix(colSums(crt), nrow=1, ncol=ncol(crt))
  
  # calculate MI
  sum(crt * log(crt / dno), na.rm = T)
}