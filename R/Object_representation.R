# Object representation


#' @export
scppp <- function(data, sample = c("rows", "columns"), ...) {
  new_scppp(data, sample)
}

#' @param data input data - Usually a matrix of counts
#'
new_scppp <- function(data, sample = c("rows", "columns"), ...) {
  UseMethod("new_scppp")
}

new_scppp.integer <- function(data, sample = "rows") {
  stopifnot(typeof(data)=="integer")
  sample <- match.arg(sample, choices = c("rows", "columns"))
  structure(
    list(data = data),
    class = "scppp"
  )
}

new_scppp.matrix <- function(data,  sample = c("rows", "columns")) {
  stopifnot(typeof(data)=="integer")
  sample <- match.arg(sample, choices = c("rows", "columns"))
  if (sample == "columns") {
    data <- t(data)
  }
  structure(
    list(data = data),
    class = "scppp"
  )
}

new_scppp.data.frame <- function(data,  sample = c("rows", "columns")) {
  .types <- vapply(data, typeof, character(1))
  if (any(.types!="integer")) {
    stop(glue::glue("columns {paste(which(.types!='integer'),collapse=', ')} are not integers"))
  }
  data <- as.matrix(data)
  sample <- match.arg(sample, choices = c("rows", "columns"))
  if (sample == "columns") {
    data <- t(data)
  }
  structure(
    list(data = data),
    class = "scppp"
  )
}


# print.scppp <- function() {}
