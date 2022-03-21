#' A more "continuous" approximation of quantile of samples with a few integer case
#'
#' new_quantile returns a data frame including data points and corresponding quantile.
#'
#' This is a function developed to get quantile for samples with only a few integer values.
#' Define both \eqn{p_{-1} = 0} and \eqn{q_{-1} = 0}.
#' Replace the point mass at each integer \eqn{z} by a bar on the interval \eqn{[z – \frac{1}{2}, z+ \frac{1}{2}]}
#' with height \eqn{P(X = z)}. This is a more "continuous" approximation of quantiles in this case.
#'
#' @param data a numeric vector of sampled data points.
#' @param sample a character string denotes which sample data points come from.
#'
#' @return a data frame. First column is the corresponding probability from cumulative distribution function (CDF),
#' second column is sample name, and third column is correponsding respective quantiles.
#'
#' @examples
#'
#' P <- sample(0:2,1000, replace=TRUE, prob=c(1/3, 1/2, 1/6))
#' new_quantile(P, "P")
#'
#' @import dplyr
#' @import stats
#'
new_quantile <- function(data, sample){

  if(!is.numeric(data)){
    warning("need numeric values for input data")
  }

  dfp <- data.frame(value = data)
  dfp <- dfp %>%
    do(data.frame(., fval = ecdf(.$value)(.$value)))
  dfp <- dfp %>% distinct()
  dfp <- dfp[with(dfp, order(fval)),]
  dfp$cat <- sample
  dfp <- rbind(data.frame(dfp), data.frame(value = c(min(data)-1), cat = sample, fval = c(0)))
  dfp$value_new <- dfp$value + 1/2
  dfp <- dfp[with(dfp, order(fval)),]
  return(dfp)
}

#' Linear interpolation for one sample given reference sample
#'
#' interpolate returns a data frame with interpolated data points
#'
#' This is a function developed to do linear interpolation for corresponding probability
#' from empirical cumulative distribution function (CDF) and corresponding quanitles.
#' Given a reference data frame and a data frame needed to do interpolation,
#' if there are any CDF values in reference but not in object data frame,
#' do the linear interpolation and insert both CDF values and respective quantiles
#' to the original object data frame.
#'
#' @param df The object data frame need to do interpolation
#' @param reference The reference data frame to make comparison
#' @param sample_id A character to denote object data frame
#'
#' @return A data frame. First column is CDF, second column is sample name,
#' and third column is correponsding quantiles.
#'
#' @examples
#'
#' P <- sample(0:2,1000, replace=TRUE, prob=c(1/3, 1/2, 1/6))
#' Q <- sample(0:1,1000, replace=TRUE, prob=c(2/3, 1/3))
#' dfp <- new_quantile(P, "P")
#' dfq <- new_quantile(Q, "Q")
#' df_tq <- interpolate(dfp, dfq, "P")
#'
#' @import purrr
#' @import stats
#'
interpolate <- function(df, reference, sample_id){
  list_ <- c(list(df), lapply(1:nrow(reference), function(x) x))
  df <- purrr::reduce(list_, function(df, i, reference){
    if (df$fval[1] == 0 &reference$fval[i] < 1 & reference$fval[i] %notin% df$fval){
      cat_add <- sample_id
      fval_add <- reference$fval[i]
      index <- max(which(df$fval < fval_add))
      if (index < nrow(df))
      {value_new_add <- approx(df$fval[index : (index+1)], df$value_new[index : (index+1)], fval_add)$y
      df <- rbind(data.frame(df), data.frame(value = (value_new_add-1/2), fval = fval_add, cat = cat_add, value_new = value_new_add))}
      df <- df[with(df, order(fval)),]
    }
    return(df)
  }, reference = reference)
  return(df)
}

#' Paired quantile after interpolation between two samples
#'
#' qq_interpolation returns a data frame with paired quantiles in two samples after interpolation.
#'
#' This is a function for quantile interpolation of two samples.
#' For each unique quantile value that has original data
#' point in one sample but no corresponding original data point in another sample,
#' apply a linear interpolation. So the common quantile values after interpolation
#' should have unique points the same as unique quantile points from either sample.
#'
#' @param dfp a data frame generated from function new_quantile() with first distribution.
#' @param dfq a data frame generated from function new_quantile() with second distribution.
#' @param sample1 a character to denote sample name of first distribution.
#' @param sample2 a character to denote sample name of second distribution.
#'
#' @return a data frame. First column is corresponding probability
#' from cumulative distribution function (CDF),
#' second column is corresponding quanitles from the first sample,
#' and third  column is correponsding quanitles from the second sample.
#'
#' @examples
#'
#' P <- sample(0:2,1000, replace=TRUE, prob=c(1/3, 1/2, 1/6))
#' Q <- sample(0:1,1000, replace=TRUE, prob=c(2/3, 1/3))
#' dfp <- new_quantile(P, "P")
#' dfq <- new_quantile(Q, "Q")
#' df_tq <- qq_interpolation(dfp, dfq, "P", "Q")
#'
qq_interpolation <- function(dfp, dfq, sample1, sample2){


  if (!all(colnames(dfp) == c("value", "fval", "cat", "value_new"))){
    warning("column names for first data not match, run new_quantile() function first")
  }

  if (!all(colnames(dfq) == c("value", "fval", "cat", "value_new"))){
    warning("column names for second data not match, run new_quantile() function first")
  }

  if (dfp$fval[1] != 0){
    warning("fvalue for first data not start from 0, run new_quantile() function first")
  }

  if (dfq$fval[1] != 0){
    warning("fvalue for second data not start from 0, run new_quantile() function first")
  }

  # interpolation on dfq
  dfq <- interpolate(df = dfq, reference = dfp, sample_id = sample2)

  # interpolation on dfp
  dfp <- interpolate(df = dfp, reference = dfq, sample_id = sample1)

  colnames(dfp) <- c("value", "fval", "cat", "value_new_p")
  dfp1 <- dfp[,c("fval", "value_new_p")]
  colnames(dfq) <- c("value", "fval", "cat", "value_new_q")
  dfq1 <- dfq[,c("fval", "value_new_q")]
  df_pq <- merge(dfp1, dfq1, by = "fval")
  #df_tq <- df_pq
  return(df_pq)
}

#' Q-Q plot comparing two samples with small discrete counts
#'
#' qqplot_small_test returns a ggplot object visualize quantiles comapring distributions of two samples
#'
#' This is a function for quantile-quantile plot comapring comparing samples from two discrete distributions
#' after \emph{continuity correction} and linear interpolation
#'
#' @param P a numeric vector from one sample
#' @param Q a numeric vector from the other sample
#' @param sample1 a character to denote sample name of first distribution.
#' @param sample2 a character to denote sample name of second distribution.
#'
#' @return a ggplot object with quantiles from one sample on the horizontal axis and corresponding qunatiles
#' from the other sample on the vertical axis.
#'
#' @examples
#'
#' P <- sample(0:2,1000, replace=TRUE, prob=c(1/3, 1/2, 1/6))
#' Q <- sample(0:1,1000, replace=TRUE, prob=c(2/3, 1/3))
#' qqplot_small_test(P, Q, "P", "Q")
#'
#' @export
qqplot_small_test <- function(P, Q, sample1, sample2){

  dfp <- new_quantile(P, sample1)
  dfq <- new_quantile(Q, sample2)

  df_tq <- qq_interpolation(dfp, dfq, sample1, sample2)
  ggplot(data = df_tq) +
    geom_line(data = df_tq, aes(x=value_new_q, y=value_new_p)) +
    geom_point(data = df_tq, aes(x=value_new_q, y=value_new_p)) +
    geom_abline(intercept = 0, slope = 1, color = "red") +
    xlab("Q") + ylab("P") +
    coord_fixed(ratio = 1)
}

#' A more "continuous" approximation of quantile of samples from theoretical Poisson distribution.
#'
#' new_quantile_pois returns a data frame including the probability
#' from cumulative distribution function (CDF) and corresponding quantiles.
#'
#' This is a function developed to get corresponding quantiles from theoretical Poisson distribution.
#' The data points ranges from 0 to maximum value of sampled data used to compare with theoretical Poisson.
#'
#' @param data a numeric vector of sampled data points to comapre with theoretical Poisson.
#' @param lambda a numeric value for theoretical Poisson distribution parameter (equal to mean).
#'
#' @return a data frame. First column is quantile, second column is "theoretical",
#' and third column is correponsding quantiles from sample data.
#'
#' @examples
#'
#' P <- rpois(100, 3)
#' new_quantile_pois(P, 3)
#'
new_quantile_pois <- function(data, lambda){

  if(!is.numeric(data)){
    warning("need numeric values for input data")
  }

  if(min(data < 0)){
    warning("Poisson random number is non negative")
  }

  dfp <- new_quantile(data, "sample")
  dfq <- data.frame(value = seq(0:max(dfp$value)), fval = ppois(q=seq(0:max(dfp$value)), lambda = lambda))
  dfq <- dfq[with(dfq, order(fval)),]
  dfq$cat <- "theoretical"
  dfq <- rbind(data.frame(dfq), data.frame(value = c(-1, 0), cat = c("theoretical"), fval = c(0, dpois(0, lambda))))
  dfq$value_new <- dfq$value + 1/2
  dfq <- dfq[with(dfq, order(fval)), ]
  return(dfq)
}

#' Random sample generation function to generate sets of samples from theoretical Poisson distribution.
#'
#' nboot_small returns a data frame with generate sets of samples and simulation index.
#'
#' This is a function used to simulate a given number sets of samples from a theoretical Poisson distribution
#' that match input samples on sample size and sample mean (or theoretical Poisson parameter).
#' Plotting these as envelopes in Q-Q plot shows the variability in shapes we can expect when
#' sampling from the theoretical Poisson distribution.
#'
#' @param x a numeric vector of sampled data points to comapre with theoretical Poisson.
#' @param lambda a numeric value for mean of theoretical Poisson.
#' @param R a numeric value for mean of theoretical Poisson.
#'
#' @return a numeric vector of number of simulation sets that match input samples on sample size
#' and sample mean (or theoretical Poisson parameter).
#'
#' @examples
#' nboot_small(rpois(100, 3), 3, 200)
#'
nboot_small <- function(x, lambda, R) {
  n <- length(x)
  do.call(rbind,
          lapply(1 : R,
                 function(i) {
                   xx <- sort(rpois(n, lambda))
                   data.frame(value = xx, sim = i)
                 }))
}

#' Q-Q plot comparing samples with a theoretical Poisson distribution
#'
#' qqplot_small_pois returns a Q-Q plot using a more "continuous" approximation of quantiles.
#'
#' This is a function for Q-Q envelope plot used to compare whether given sample data points come from the
#' theoretical Poisson distribution.  By simulating repeated samples of the same size from the candidate
#' theoretical distribution, and overlaying the envelope on the same figure, it provides a feeling of
#' understanding the natural variation from the theoretical distribution.
#'
#' @param sample_data a numeric vector of sample data points
#' @param lambda a numeric value for theoretical Poisson parameter
#' @param envelope_size a numeric value of size of envelope on Q-Q plot (default 100)
#'
#' @return a ggplot object
#'
#' @examples
#'
#' qqplot_env_pois(rpois(200, 3), 3, 100)
#'
#' @import ggplot2
#'
#' @export
qqplot_env_pois <- function(sample_data, lambda, envelope_size = 100, ...) {
  UseMethod("qqplot_env_pois")
}

#' @export
qqplot_env_pois.numeric <- function(sample_data, lambda, envelope_size = 100){

  test_raw <- sample_data

  dfp <- new_quantile(test_raw, "sample")
  dfq <- new_quantile_pois(test_raw, lambda)

  df_tq <- qq_interpolation(dfp, dfq, "sample", "theoretical")
  df_tq <- subset(df_tq, df_tq$value_new_p >= (dfp$value[2]-0.5) & df_tq$value_new_p <= max(dfp$value)
                  & df_tq$value_new_q >= (dfp$value[2]-0.5) & df_tq$value_new_q <= max(dfp$value))
  df_tq$source <- "sample"

  test_all_qq <- df_tq

  gb <- nboot_small(sample_data, lambda, envelope_size)

  df_all <- purrr::map_dfr(
    .x = 1:envelope_size,
    .f = function(j, gb, lambda){
      dfp <- new_quantile(gb$value[which(gb$sim == j)], j)
      dfq <- new_quantile_pois(dfp$value, lambda)
      df_pq <- qq_interpolation(dfp, dfq, j, "theoretical")
      df_pq$sim <- j
      return(df_pq)
    }, gb = gb, lambda = lambda)
  df_all <- subset(df_all, df_all$value_new_p >= min(test_all_qq$value_new_p-1) & df_all$value_new_p <= (max(test_all_qq$value_new_q, test_all_qq$value_new_p) + 1)
                   & df_all$value_new_q >= min(test_all_qq$value_new_p-1) & df_all$value_new_q <= (max(test_all_qq$value_new_q, test_all_qq$value_new_p) + 1))

  min <- min(df_all$value_new_q, df_all$value_new_p, test_all_qq$value_new_q, test_all_qq$value_new_p)
  max <- max(df_all$value_new_q, df_all$value_new_p, test_all_qq$value_new_q, test_all_qq$value_new_p)
  #size <- 18
  p <- ggplot(data = df_all) +
    geom_line(aes(x = value_new_q, y=value_new_p, group = sim),
              color = "gray", size = 1.2) +
    geom_line(data = test_all_qq, aes(x = value_new_q, y = value_new_p, color = source), size = 1.2) +
    geom_point(data = test_all_qq, aes(x = value_new_q, y = value_new_p, color = source), size = 1.5) +
    scale_color_manual(values=c("#E69F00")) +
    geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
    #coord_cartesian(xlim =c((dfp$value[2]-1), max(dfp$value + 1)), ylim = c((dfp$value[2]-1), max(dfp$value + 1))) +
    xlab("theoretical quantile") + ylab("sample quantile") + coord_fixed(ratio = 1) +
    labs(title = paste(expression("\u03BB"), "=", lambda)) +
    xlim(c(min, max)) + ylim(c(min, max))
  return(p)
}

select_data <- function(test_dat, Gfit_dat, lambda,
                        select_by = "entry", entry_size = 200){

  if(select_by == "entry"){
    test_raw <- as.vector(test_dat)
    fitdata <- as.vector(Gfit_dat)
    order <- order(abs(fitdata - lambda), decreasing = F)
    test_raw <- test_raw[order]
    return(test_raw[1:entry_size])
  }

  if(select_by == "cell"){
    cell_mean <- rowMeans(test_dat)
    measure <- abs(cell_mean - lambda)
    index <- which(measure == min(measure))

    cell <- rownames(test_dat)[index]
    cell_sd <- rowSds(test_dat)[index]
    select_cell <- cell[which(cell_sd == max(cell_sd))]
    return(test_dat[which(rownames(test_dat) == select_cell[1]), ])
  }

  if(select_by == "gene"){
    gene_mean <- colMeans(test_dat)
    measure <- abs(gene_mean - lambda)
    index <- which(measure == min(measure))

    gene <- colnames(test_dat)[index]
    gene_sd <- colSds(test_dat)[index]
    select_gene <- gene[which(gene_sd == max(gene_sd))]
    return(test_dat[, which(colnames(test_dat) == select_gene[1])])
  }

}

apply_glmpca <- function(rawdata, L = 10){

  set.seed(1234)
  test_df <- t(rawdata)
  test_df <- test_df[which(rowSums(test_df) > 0), ]
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

#' @export
qqplot_env_pois.scppp <- function(scppp_obj, lambda,
                                  L = 10,
                                  select_by = "entry",
                                  entry_size = 200, envelope_size = 100){
  test_data <- scppp_obj[["data"]]
  if(is.vector(test_data) && is.atomic(test_data)){
    qqplot_env_pois.numeric(test_data, lambda, envelope_size)
  }
  else {
    test <- as.matrix(test_data)
    Gfit_dat <- apply_glmpca(test_data, L)
    sample_data <- select_data(test_data, Gfit_dat, lambda, select_by, entry_size)
    qqplot_env_pois.numeric(sample_data, lambda, envelope_size)
  }
}









