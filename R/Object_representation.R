# Object representation

#' Generate New scppp object
#'
#' Define S3 class that stores scRNA-seq data and associated information
#' (e.g. model departure representation, cell clustering results) if
#' corresponding functions are called.
#'
#' @param data input data - Usually a matrix of counts
#' @param sample by rows or columns
#' @return S3 object for class 'scppp'.
#' @export
scppp <- function(data, sample = c("columns", "rows")) {
  new_scppp(data, sample)
}


#'
new_scppp <- function(data, sample = c("columns", "rows")) {
  UseMethod("new_scppp")
}

new_scppp.integer <- function(data, sample = NULL) {
  stopifnot(typeof(data)=="integer")
  structure(
    list(data = data),
    class = "scppp"
  )
}

new_scppp.matrix <- function(data,  sample = c("columns", "rows")) {
  stopifnot(typeof(data)=="integer")
  sample <- match.arg(sample[1], choices = c("columns", "rows"))
  if (sample == "rows") {
    data <- t(data)
  }
  structure(
    list(data = data),
    class = "scppp"
  )
}

new_scppp.data.frame <- function(data,  sample = c("columns", "rows")) {
  .types <- vapply(data, typeof, character(1))
  if (any(.types!="integer")) {
    stop(sprintf("columns %s are not integers", paste(which(.types!='integer'),collapse=', ')))
  }
  data <- as.matrix(data)
  sample <- match.arg(sample[1], choices = c("rows", "columns"))
  if (sample == "rows") {
    data <- t(data)
  }
  structure(
    list(data = data),
    class = "scppp"
  )
}


# print.scppp <- function() {}
