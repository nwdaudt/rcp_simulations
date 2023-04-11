##
## RCP simulations
## source
##
## Nicholas W Daudt
## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##

#### Helper functions for checking RCPs fit

# FUN: max value in row to 1 and others to 0
FUN_max1else0 <- function(x) {
  # All values other than the max to 0
  x <- replace(x, x < max(x), 0)
  # Then, 'max' values to 1
  x <- replace(x, x > 0, 1)
  return(x)
}

# FUN: any RCP assigned for less than "N" sites?
FUN_lessNsites <- function(x, n) {
  # Apply 'fun_max1_else0' by row, and transpose the matrix to correct format
  b <- t(apply(X = x, MARGIN = 1, FUN = FUN_max1else0))
  # colSum values of the matrix (cols are RCP groups), and if <=n replace with 1
  c <- replace((colSums(b) <= n), (colSums(b) <= n) == TRUE, 1)
  # sum to get the number of columns (RCP groups) with <=n sites allocated
  d <- sum(c)
  return(d)
}