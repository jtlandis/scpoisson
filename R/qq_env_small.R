#' A more "continuous" approximation of quantiles of samples with a few integer case
#'
#' This function returns a data frame including data points and corresponding quantile.
#'
#' This is a function developed to get quantile for samples with only a few integer values.
#' Define both \eqn{p_{-1} = 0} and \eqn{q_{-1} = 0}.
#' Replace the point mass at each integer \eqn{z} by a bar on the interval \eqn{[z – \frac{1}{2}, z+ \frac{1}{2}]}
#' with height \eqn{P(X = z)}. This is a more "continuous" approximation of quantiles in this case.
#'
#' @param data A numeric vector of sampled data points.
#' @param sample A character string denotes which sample data points come from.
#'
#' @return A data frame contains the corresponding probability from cumulative distribution function (CDF), sample name, and corresponding respective quantiles.
#'
#'
new_quantile <- function(data, sample){

  if(!is.numeric(data)){
    warning("need numeric values for input data")
  }

  dfp <- data.frame(value = data)
  dfp <- dfp %>%
    dplyr::do(data.frame(.data, fval = ecdf(.data$value)(.data$value)))
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
#' This function returns a data frame with interpolated data points.
#'
#' This is a function developed to do linear interpolation for corresponding probability
#' from empirical cumulative distribution function (CDF) and corresponding quantiles.
#' Given a reference data frame and a data frame needed to do interpolation,
#' if there are any CDF values in reference but not in object data frame,
#' do the linear interpolation and insert both CDF values and respective quantiles
#' to the original object data frame.
#'
#' @param df The object data frame requires interpolation.
#' @param reference The reference data frame to make comparison.
#' @param sample_id A character to denote the object data frame.
#'
#' @return A data frame contains CDF, the sample name, and the corresponding quantiles.
#'
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
#' This function returns a data frame with paired quantiles in two samples after interpolation.
#'
#' This is a function for quantile interpolation of two samples.
#' For each unique quantile value that has original data
#' point in one sample but no corresponding original data point in another sample,
#' apply a linear interpolation. So the common quantile values after interpolation
#' should have unique points the same as unique quantile points from either sample.
#'
#' @param dfp A data frame generated from function new_quantile() based on a specific distribution.
#' @param dfq Another data frame generated from function new_quantile() based on a specific distribution.
#' @param sample1 A character to denote sample name of distribution used to generate \code{dfp}.
#' @param sample2 A character to denote sample name of distribution used to generate \code{dfq}.
#'
#' @return A data frame contains corresponding probability from cumulative distribution function (CDF),
#' corresponding quantiles from the first sample (\code{dfp}),
#' and corresponding quantiles from the second sample (\code{dfq}).
#'
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
#' This function returns a ggplot object used to visualize quantiles comparing distributions of two samples.
#'
#' This is a function for quantile-quantile plot comparing comparing samples from two discrete distributions
#' after \emph{continuity correction} and linear interpolation
#'
#' @param P A numeric vector from one sample.
#' @param Q A numeric vector from the other sample.
#' @param sample1 A character to denote sample name of one distribution \code{P} generated from.
#' @param sample2 A character to denote sample name of the other distribution \code{Q} generated from.
#'
#' @return A ggplot object. Q-Q plot with continuity correction. Quantiles from one sample on the horizontal axis and corresponding quantiles
#' from the other sample on the vertical axis.
#'

qqplot_small_test <- function(P, Q, sample1, sample2){

  dfp <- new_quantile(P, sample1)
  dfq <- new_quantile(Q, sample2)

  df_tq <- qq_interpolation(dfp, dfq, sample1, sample2)
  ggplot(data = df_tq) +
    geom_line(data = df_tq, aes(x=.data$value_new_q, y=.data$value_new_p)) +
    geom_point(data = df_tq, aes(x=.data$value_new_q, y=.data$value_new_p)) +
    geom_abline(intercept = 0, slope = 1, color = "red") +
    xlab("Q") + ylab("P") +
    coord_fixed(ratio = 1)
}

#' A more "continuous" approximation of quantiles from the theoretical Poisson distribution.
#'
#' This function returns a data frame including the probability
#' from cumulative distribution function (CDF) and corresponding quantiles.
#'
#' This is a function developed to get corresponding quantiles from theoretical Poisson distribution.
#' The data points ranges from 0 to maximum value of sampled data used to compare with the theoretical Poisson distribution.
#'
#' @param data A numeric vector of sampled data points to compare with theoretical Poisson.
#' @param lambda A numeric value for theoretical Poisson distribution parameter (equal to mean).
#'
#' @return A data frame contains CDF probability and corresponding quantiles from the theoretical Poisson distribution.
#'
#'
new_quantile_pois <- function(data, lambda){

  if(!is.numeric(data)){
    warning("need numeric values for input data")
  }

  if(min(data < 0)){
    warning("Poisson random number should be non negative")
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
#' This function returns a data frame with generated sets of samples and simulation index.
#'
#' This is a function used to simulate a given number sets of samples from a theoretical Poisson distribution
#' that match input samples on sample size and sample mean (or theoretical Poisson parameter).
#' Plotting these as envelopes in Q-Q plot shows the variability in shapes we can expect when
#' sampling from the theoretical Poisson distribution.
#'
#' @param x A numeric vector of sampled data points to compare with theoretical Poisson.
#' @param lambda A numeric value specifying mean for theoretical Poisson distribution.
#' @param R A numeric value specifying the number of simulated sets.
#'
#' @return A data frame contains simulated data and corresponding simulation index.

#' Random sample generation function to generate sets of samples from theoretical Poisson distribution.
#'
#' nboot_small returns a data frame with generate sets of samples and simulation index.
#'
#' This is a function used to simulate a given number sets of samples from a theoretical Poisson distribution
#' that match input samples on sample size and sample mean (or theoretical Poisson parameter).
#' Plotting these as envelopes in Q-Q plot shows the variability in shapes we can expect when
#' sampling from the theoretical Poisson distribution.
#'
#' @param x a numeric vector of sampled data points to compare with theoretical Poisson.
#' @param lambda a numeric value for mean of theoretical Poisson.
#' @param R a numeric value for mean of theoretical Poisson.
#'
#' @return a numeric vector of number of simulation sets that match input samples on sample size
#' and sample mean (or theoretical Poisson parameter).
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
#' This function returns a Q-Q plot with envelope using a more "continuous" approximation of quantiles.
#'
#' This is a function for Q-Q envelope plot used to compare whether given sample data points come from the
#' theoretical Poisson distribution.  By simulating repeated samples of the same size from the candidate
#' theoretical distribution, and overlaying the envelope on the same figure, it provides a feeling of
#' understanding the natural variation from the theoretical distribution.
#'
#' If an S3 object for class 'scppp' is used as input and the stored result under "data" is a matrix,
#' The GLM-PCA algorithm will be applied to estimate the Poisson parameter for each matrix entry.
#' Then a specific number of entries will be selected as sample data points to compare with the theoretical Poisson distribution.
#'
#' @param sample_data A numeric vector of sample data points or an S3 object for class 'scppp'.
#' @param lambda A numeric value specifying the theoretical Poisson parameter.
#' @param envelope_size A numeric value specifying the size of envelope on Q-Q plot (default 100).
#' @param ... not used.
#'
#' @return A ggplot object.
#'
#' @references
#' \insertRef{glmpca}{scpoisson}
#'
#' @import ggplot2
#'
#' @export
#'
qqplot_env_pois <- function(sample_data, lambda, envelope_size = 100, ...) {
  UseMethod("qqplot_env_pois")
}

#' @export
qqplot_env_pois.numeric <- function(sample_data, lambda, envelope_size = 100, ...){

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
    geom_line(aes(x = .data$value_new_q, y=.data$value_new_p, group = .data$sim),
              color = "gray", size = 1.2) +
    geom_line(data = test_all_qq, aes(x = .data$value_new_q, y = .data$value_new_p, color = source), size = 1.2) +
    geom_point(data = test_all_qq, aes(x = .data$value_new_q, y = .data$value_new_p, color = source), size = 1.5) +
    scale_color_manual(values=c("#E69F00")) +
    geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
    #coord_cartesian(xlim =c((dfp$value[2]-1), max(dfp$value + 1)), ylim = c((dfp$value[2]-1), max(dfp$value + 1))) +
    xlab("theoretical quantile") + ylab("sample quantile") + coord_fixed(ratio = 1) +
    labs(title = paste(expression("\u03BB"), "=", lambda)) +
    xlim(c(min, max)) + ylim(c(min, max))
  return(p)
}

#' @param L A numeric value specifying the number of latent vectors included when estimate the Poisson parameter for each matrix entry. This is not useful if a numeric vector is used as input.
#' @param select_by A character indicating whether entries should be selected by
#' \itemize{
#' \item{entry}: {independent matrix entries with Poisson parameter estimates closest to \code{lambda}}
#' \item{\code{cell}}: {one particular cell with UMI count mean closest to \code{lambda}}
#' \item{\code{gene}}: {one particular gene with UMI count mean closest to \code{lambda}}
#' }
#' This is not useful if a numeric vector is used as input.
#' @param entry_size A numeric value specifying the number of entries used to compare with the theoretical Poisson distribution.
#' This is not useful if a numeric vector is used as input, or the entries are selected by cell or gene.
#' @export
qqplot_env_pois.scppp <- function(sample_data, lambda, envelope_size = 100,
                                  L = 10,
                                  select_by = "entry",
                                  entry_size = 200, ...){
  test_data <- sample_data[["data"]]
  if(is.vector(test_data) && is.atomic(test_data)){
    qqplot_env_pois.numeric(test_data, lambda, envelope_size)
  }
  else {

    select_data <- function(test_dat, Gfit_dat, lambda,
                            select_by = "entry", entry_size = 200){

      if(select_by == "entry"){
        test_raw <- as.vector(test_dat)
        fitdata <- as.vector(Gfit_dat)
        order <- order(abs(fitdata - lambda), decreasing = F)
        test_raw <- test_raw[order]
        return(test_raw[1:entry_size])
      }

      if(select_by == "gene"){
        cell_mean <- rowMeans(test_dat)
        measure <- abs(cell_mean - lambda)
        index <- which(measure == min(measure))

        cell <- rownames(test_dat)[index]
        cell_sd <- matrixStats::rowSds(test_dat)[index]
        select_cell <- cell[which(cell_sd == max(cell_sd))]
        return(test_dat[which(rownames(test_dat) == select_cell[1]), ])
      }

      if(select_by == "cell"){
        gene_mean <- colMeans(test_dat)
        measure <- abs(gene_mean - lambda)
        index <- which(measure == min(measure))

        gene <- colnames(test_dat)[index]
        gene_sd <- matrixStats::colSds(test_dat)[index]
        select_gene <- gene[which(gene_sd == max(gene_sd))]
        return(test_dat[, which(colnames(test_dat) == select_gene[1])])
      }

    }

    apply_glmpca <- function(rawdata, L = 10){

      #set.seed(1234)
      test_df <- rawdata
      test_df <- test_df[which(rowSums(test_df) > 0), ]
      ctl <- list(maxIter=500,eps=1e-4)
      res <- glmpca::glmpca(test_df, L=L, fam = "poi", sz=colSums(test_df),verbose=TRUE,ctl=ctl)
      factors <- res$factors

      U <- as.matrix(res$factors)
      V <- as.matrix(res$loadings)
      Ni <- colSums(test_df)
      Vj <- as.vector(as.matrix(res$coefX))
      Gfit_dat <- exp(t(t(U %*% t(V)) + Vj) + log(Ni))
      return(t(Gfit_dat))
    }

    test <- as.matrix(test_data)
    Gfit_dat <- apply_glmpca(test_data, L)
    sample_data <- select_data(test_data, Gfit_dat, lambda, select_by, entry_size)
    qqplot_env_pois.numeric(sample_data, lambda, envelope_size)
  }
}









