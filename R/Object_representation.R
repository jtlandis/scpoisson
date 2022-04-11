# Object representation

#' Generate New scppp object
#' @param data input data - Usually a matrix of counts
#' @param sample by rows or columns
#' @export
scppp <- function(data, sample = c("columns", "rows"), ...) {
  new_scppp(data, sample)
}


#'
new_scppp <- function(data, sample = c("columns", "rows"), ...) {
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
    stop(glue::glue("columns {paste(which(.types!='integer'),collapse=', ')} are not integers"))
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
