#' Parameter estimates based on two-way approximation
#'
#' This function returns a vector consists of parameter estimates for overall offset,
#' cell effect, and gene effect.
#'
#' This is a function used to calculate parameter estimates based on
#' \eqn{\lambda_{gc}  = e^{\mu + \alpha_g + \beta_c}},
#' where \eqn{\mu} is the overall offset,
#' \eqn{\alpha} is a vector with the same length as the number of genes,
#' and \eqn{\beta} is a vector with the same length as the number of cells.
#' The order of elements in vectors \eqn{\alpha} or \eqn{\beta} is the same as rows (genes) or
#' cells (columns) from input data. Be sure to remove cells/genes with all zeros.
#'
#' @param test_set A UMI count data matrix with genes as rows and cells as columns
#'
#' @return A numeric vector containing parameter estimates from overall offset (first element), gene effect (same order as rows) and cell effect (same order as columns).
#'
#' @examples
#' # Matrix as input
#' test_set <- matrix(rpois(500, 0.5), nrow = 10)
#' para_est_new(test_set)
#'
#' @export
para_est_new <- function(test_set){

  stopifnot('Remove columns with only zero values' = min(colSums(test_set)) > 0)
  stopifnot('Remove rows with only zero values' = min(rowSums(test_set)) > 0)
  test_set <- test_set[, which(colSums(test_set) > 0)]
  n <- nrow(test_set)
  d <- ncol(test_set)
  N <- n+d+1

  mu <- log(sum(test_set) / (n*d))
  w <- log(rowMeans(test_set)) -mu
  r <- log(colMeans(test_set)) - mu

  para <- c(mu, w, r)

  return(para)
}

#' A novel data representation based on Poisson probability
#'
#' This function returns a matrix of a novel data representation with the same dimension as input data matrix.
#'
#' This is a function used to calculate model departure as a novel data representation.
#'
#' @param data A UMI count data matrix with genes as rows and cells as columns or an S3 object for class 'scppp'.
#' @param change A numeric value used to correct for exactly 0 and 1 before logit transformation.
#' Any values below \code{change} are set to be \code{change} and
#' any values above \eqn{1- change} are set to be \eqn{1- change}.
#' @param ... not used.
#'
#' @return A matrix of departure as a novel data representation (matrix as input) or an S3 object for class 'scppp' (scppp object as input; departure result will be stored in object scppp under "representation").
#'
#' @examples
#' # Matrix as input
#' test_set <- matrix(rpois(500, 0.5), nrow = 10)
#' adj_CDF_logit(test_set)
#' # scppp object as input
#' adj_CDF_logit(scppp(test_set))
#'
#' @export
adj_CDF_logit <- function(data, change = 1e-10, ...) UseMethod("adj_CDF_logit")

#' @export
#' @return scppp
adj_CDF_logit.scppp <- function(data, change = 1e-10, ...) {

  test_dat <- data[["data"]]
  data$representation[["departure"]] <- adj_CDF_logit.matrix(test_dat, change)
  return(data)
}

#' @export
#' @return scppp_departure
adj_CDF_logit.matrix <- function(data, change = 1e-10, ...){

  test_set <- data
  stopifnot('Require a matrix as input' = is.matrix(test_set))
  stopifnot('Remove columns with only zero values' = min(colSums(test_set)) > 0)
  stopifnot('Remove rows with only zero values' = min(rowSums(test_set)) > 0)
  n <- nrow(test_set)
  d <- ncol(test_set)
  para <- para_est_new(test_set)
  mu <- para[1]
  w <- para[2:(n+1)]
  r <- para[(n+2):(n+d+1)]
  log_Pe <- mu + w %*% matrix(rep(1, len = d), nrow = 1) + matrix(rep(1, len = n), ncol = 1) %*% r
  Pe <- exp(log_Pe)
  rownames(Pe) <- rownames(test_set)
  colnames(Pe) <- colnames(test_set)

  # adjust CDF
  F_adj <- function(a, b, .version = version){
    Fadj <- (ppois(a,b) + ppois(a-1,b)) / 2
    return(Fadj)
  }

  l_adj_test <- list(a = as.list(as.vector(as.matrix(test_set))), b = as.list(Pe))
  cdf_adj_test <- purrr::pmap(l_adj_test, function(a,b) F_adj(a,b)) %>%
    matrix(ncol = ncol(test_set))
  mcdf_adj_test <- unlist(cdf_adj_test) %>% matrix(ncol = ncol(test_set))
  rownames(mcdf_adj_test) <- rownames(test_set)
  colnames(mcdf_adj_test) <- colnames(test_set)
  mcdf_adj_test[mcdf_adj_test < change] <- change
  mcdf_adj_test[mcdf_adj_test > (1 - change)] <- 1 - change
  return(logit(mcdf_adj_test))
}

