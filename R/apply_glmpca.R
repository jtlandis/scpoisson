

apply_glmpca <- function(rawdata, L = 10){

  set.seed(1234)
  test_df <- rawdata[which(rowSums(rawdata) > 0), ]
  ctl <- list(maxIter=500,eps=1e-4)
  res <- glmpca::glmpca(test_df, L=L, fam = "poi", sz=colSums(test_df),verbose=TRUE,ctl=ctl)
  factors <- res$factors

  U <- as.matrix(res$factors)
  V <- as.matrix(res$loadings)
  Ni <- colSums(test_df)
  Vj <- as.vector(as.matrix(res$coefX))
  Gfit_dat <- exp(t(t(U %*% t(V)) + Vj) + log(Ni))
  return((Gfit_dat))
}
