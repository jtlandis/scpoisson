
print.scppp <- function(x, ...) {
  cat("An object of class scppp\n")
  cat("  * data = [",paste0(dim(x$data), collapse = " x ") ,
      "] <",paste0(class(x$data),collapse = "/"),"> ", sep = "")
}



