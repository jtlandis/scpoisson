

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

nboot_small <- function(n, lambda, R) {
  do.call(rbind,
          lapply(1 : R,
                 function(i) {
                   xx <- sort(nsim(n, lambda))
                   data.frame(value = xx, sim = i)
                 }))
}

