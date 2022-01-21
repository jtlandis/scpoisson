

logit <- function (x, min = 0, max = 1) {
  p <- (x - min)/(max - min)
  log(p/(1 - p))
}

`%notin%` <- Negate(`%in%`)

nsim <- rpois

clust_clean <- function(clust){
  if(grepl("NA", clust, fixed = TRUE)){
    gsub("(.*?)(-NA.*)", "\\1", clust)
  }
  else(gsub('.{2}$', '', clust))
}


theme_dirk <- function(base_size = 22,
                       base_family = "",
                       base_line_size = base_size/22,
                       base_rect_size = base_size/22, time_stamp = FALSE){
  
  obj <- theme_classic(base_size = base_size,
                       base_family = base_family,
                       base_line_size = base_line_size,
                       base_rect_size = base_rect_size) %+replace%
    theme(
      axis.text = element_text(size = base_size),
      axis.text.x = element_text(vjust = 0.5),
      axis.title.x = element_text(margin = margin(t = .8 * base_size), vjust = 0),
      axis.title.x.top = element_text(margin = margin(b = .8 * base_size), vjust = 1),
      axis.title.y = element_text(margin = margin(r = .8 * base_size), vjust = 1, angle = 90),
      axis.title.y.right = element_text(margin = margin(l = .8 * base_size), vjust = 0, angle = 90),
      strip.text = element_text(size = base_size),
      strip.background = element_rect(colour = "white", fill = "white"),
      panel.grid = element_blank(),
      panel.background = element_rect(colour = "black", fill = "white"),
      plot.caption = element_text(size = rel(.5), hjust = 1, vjust = 0),
      plot.caption.position = "plot",
      complete = T
    )
  
  if (time_stamp) {
    obj <- list(obj, labs(caption = Sys.time()))
  }
  
  return(obj)
  
}


# do_data <- function(x) {
#   env <- new.env()
#   load(x, env)
#   name <- paste0("Zhengmix4eq_craft_",gsub("\\.Rdata", replacement = "", x = basename(x)))
#   assign(name, value = env$res)
#   eval(substitute(usethis::use_data(.x, overwrite = T), list(.x = as.name(name))))
# }
# 
# do_data <- function(x) {
#   env <- new.env()
#   load(x, env)
#   name <- paste0("Siyao_",gsub("\\.Rdata", replacement = "", x = basename(x)))
#   assign(name, value = env$res)
#   eval(substitute(usethis::use_data(.x, overwrite = T), list(.x = as.name(name))))
# }
# 
# lapply(list.files(pattern = "artificial_try", full.names = T), do_data)
