
print.scppp <- function(x, ...) {
  cat("An object of class scppp\n")
  cat("  * data = [",paste0(dim(x$data), collapse = " x ") ,
      "] <",paste0(class(x$data),collapse = "/"),"> ", sep = "")
}

#' get example data
#' @param x data set to choose
#' @return A data set from example data
#' @export
get_example_data <- function(x = c("p5","p56")) {
  switch (x[1],
    p5 = p5,
    p56 = p56
  )
}
