#' Differential expression analysis
#'
#' This function returns a data frame with differential expression analysis results.
#'
#' This is a function used to find deferentially expressed genes between two clusters.
#'
#' @param data A departure matrix generated from adj_CDF_logit() or an S3 object for class 'scppp'.
#' @param final_clust_res A data frame with clustering results generated from HclustDepart(). It contains two columns: names (cell names) and clusters (cluster label).
#' @param clust1 One of the cluster label used to make comparison, default "1".
#' @param clust2 The other cluster label used to make comparison, default "2".
#' @param t_test A logical value indicating whether the t-test should be used to make comparison. In general, for large cluster (\eqn{n \ge 30}), the t-test should be used. Otherwise, the Wilcoxon test might be more appropriate.
#' @param ... not used.
#'
#' @return A data frame contains genes (ranked by decreasing order of mean difference), and associated statistics (p-values, FDR adjusted p-values, etc.).
#' If the input is an S3 object for class 'scppp', differential expression analysis results will be stored in object scppp under "de_results".
#'
#' @export
diff_gene_list <- function(data, final_clust_res = NULL,
                           clust1 = "1", clust2 = "2",
                           t_test = FALSE, ...) UseMethod("diff_gene_list")

#' @export
#' @return scppp
diff_gene_list.scppp <- function(data, final_clust_res = NULL,
                                 clust1 = "1", clust2 = "2",
                                 t_test = FALSE, ...) {
  stopifnot('Require a clustering result, run HclustDepart() first' = "clust_results" %in% attributes(data)$names)
  final_clust_res <- data[["clust_results"]]$Hclust[[1]]
  stopifnot('Clust1 not match the cluster label from HclustDepart' = clust1 %in% final_clust_res$cluster)
  stopifnot('Clust2 not match the cluster label from HclustDepart' = clust2 %in% final_clust_res$cluster)
  test_dat2 <- data[["representation"]]$departure
  data$de_results[["Hclust"]] <- diff_gene_list.matrix(test_dat2, final_clust_res,
                                                          clust1, clust2, t_test)
  return(data)
}

#' @export
#' @return scppp_de_results
diff_gene_list.matrix <- function(data, final_clust_res = NULL,
                                  clust1 = "1", clust2 = "2",
                                  t_test = FALSE, ...){
  stopifnot('Clust1 not match the cluster label from HclustDepart' = clust1 %in% final_clust_res$cluster)
  stopifnot('Clust2 not match the cluster label from HclustDepart' = clust2 %in% final_clust_res$cluster)
  test_dat2 <- t(data)
  cell_in <- final_clust_res$names[which(final_clust_res$cluster %in%
                                           c(clust1, clust2))]
  test_dat2 <- test_dat2[which(rownames(test_dat2) %in% cell_in),]

  dat_long3 <- data.frame(names = rownames(test_dat2), test_dat2)
  #dat_long3 <- reshape2::melt(dat_long3, id.vars=c("names"))
  dat_long3 <- dat_long3 %>%
    tidyr::gather("variable", "value", -names)
  dat_long3 <- merge(dat_long3, final_clust_res, by = "names")
  dat_long3$clust_test <- ifelse(dat_long3$cluster %in% clust1, "A", "B")

  fc_df2 <- dat_long3 %>%
    dplyr::group_by(.data$variable, .data$clust_test) %>%
    dplyr::summarize(Mean = mean(.data$value, na.rm=TRUE))


  #fc_df_wide <- reshape2::dcast(fc_df2,  variable ~ clust_test, value.var="Mean")
  fc_df_wide <- fc_df2 %>%
    tidyr::spread(key = .data$clust_test, value = .data$Mean)
  fc_df_wide$mean_diff <- fc_df_wide$A - fc_df_wide$B

  if(!t_test){
    k <- dat_long3 %>%
      dplyr::group_by(.data$variable, .data$clust_test) %>%
      tidyr::nest() %>%
      tidyr::spread(key = .data$clust_test, value = data) %>%
      dplyr::mutate(
        test = purrr::map2(.data$A, .data$B, ~{wilcox.test(.x$value, .y$value) %>% broom::tidy()}),
        A = purrr::map(.data$A, nrow),
        B = purrr::map(.data$B, nrow)
      ) %>%
      tidyr::unnest(cols = c(.data$A, .data$B, .data$test))
  }

  if(t_test){
    k <- dat_long3 %>%
      dplyr::group_by(.data$variable, .data$clust_test) %>%
      tidyr::nest() %>%
      tidyr::spread(key = .data$clust_test, value = data) %>%
      dplyr::mutate(
        test = purrr::map2(.data$A, .data$B, ~{t.test(.x$value, .y$value) %>% broom::tidy()}),
        A = purrr::map(.data$A, nrow),
        B = purrr::map(.data$B, nrow)
      ) %>%
      tidyr::unnest(cols = c(.data$A, .data$B, .data$test))
  }

  fc_df_wide <- merge(fc_df_wide, k, by = "variable")
  fc_df_wide$abs_diff <- abs(fc_df_wide$mean_diff)
  fc_df_wide$FC <- fc_df_wide$A.x / fc_df_wide$B.x
  fc_df_wide$padj <- p.adjust(fc_df_wide$p.value, "fdr")
  fc_df_wide <- fc_df_wide %>%
    dplyr::rename(clust1_mean = "A.x",
                  clust2_mean = "B.x",
                  clust1_n = "A.y",
                  clust2_n = "B.y")
  return(fc_df_wide %>%
           dplyr::select(.data$variable, .data$clust1_mean, .data$clust2_mean, .data$clust1_n, .data$clust2_n,
                         .data$mean_diff, .data$statistic, .data$p.value, .data$padj, .data$abs_diff) %>%
           #filter(padj < 0.05) %>%
           arrange(desc(.data$mean_diff)))
}


