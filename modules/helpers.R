myColors <- c(
  "NS" = "grey30", 
  "IL1B" = "#FC4E07", 
  "GC" = "#E7B800", 
  "combo" = "#2E9FDF"
)

cleanup_objects <- function(keep = readLines("cleanup_ignore.txt")) {
  all_objs <- ls(envir = .GlobalEnv)
  to_remove <- setdiff(all_objs, keep)
  rm(list = to_remove, envir = .GlobalEnv)
  # gc()
}

plot_upset <- function(logicalMatrix, minSize = 50, plotTitle = "title here"){
  
  # a wrapper for the ComplexUpset::upset function that applies my custom theme
  # to create an UpSet plot given a logical matrix specifying samples as colnames,
  # presence/absence as 1/0, and gene names as rownames.
  
  suppressWarnings(
    upset(
      logicalMatrix,
      intersect = colnames(logicalMatrix),
      min_size = minSize,
      set_sizes = F,
      height_ratio = 0.75,
      wrap = T,
      base_annotations = list(
        "intersection size" = intersection_size(
          text = list(
            color = "black",
            angle = 90, 
            vjust = 0.5,
            hjust = -0.3,
            size = 3
          )
        )
      )) +
      theme_minimal() +
      theme(
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 20, hjust = 0.05)
      ) +
      labs(title = plotTitle)
  )
}




