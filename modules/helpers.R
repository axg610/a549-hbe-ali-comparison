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
        # plot.title = element_text(size = 20, hjust = 0.05)
      ) +
      labs(title = plotTitle)
  )
}

plot_euler <- function(
    logicalMatrix, 
    plotType = "euler",
    shapeType = "circle", 
    plotTitle = "title here",
    quants = c("counts"),     # vector of metrics to show (e.g., counts, percent)
    cutoff = 0,               # intersections smaller than this number will be removed
    returnAsFunction = FALSE, # helps for patchworking multiple eulers together
    showLabels = T,
    fillColors = c("#FFFFFF", "#D9D9D9", "#A6CEE3"),
    aspectRatio = 1
){
  
  # a wrapper for the eulerr:euler function to create a Euler diagram given a
  # logical matrix specifying groups as colnames, presence/absence as 1/0, and
  # gene names as rownames.
  
  # prepare intersections
  obj <- logicalMatrix %>%
    mutate(
      across(
        everything(),
        ~ ifelse(. == 1, cur_column(), NA)
      )
    ) %>%
    rowwise() %>%
    mutate(group = paste(na.omit(c_across(cols = everything())), collapse = "&")) %>%
    count(group) %>%
    deframe()
  
  # impose cutoff
  obj <- obj[obj >= cutoff]
  
  # calculate and plot
  if(plotType == "euler"){
    obj <- euler(obj, shape = shapeType)
  }
  else if(plotType == "venn"){
    obj <- venn(obj)
  }
  
  if(returnAsFunction){
    function() {
      plot(
        obj,
        quantities = list(type = quants),
        main = plotTitle,
        asp = aspectRatio,
        labels = showLabels,
        fills = list(fill = fillColors)
      )
    }
  }
  
  else if(!returnAsFunction){
    plot(
      obj,
      quantities = list(type = quants),
      main = plotTitle,
      asp = aspectRatio,
      labels = showLabels,
      fills = list(fill = fillColors)
    )
    
  }
  
}

generic_l2f_heatmap_colors <- function(mat){
  
  # return a pretty log2fold heatmap color palette given a matrix
  
  colors = circlize::colorRamp2(
    c(
      seq(
        quantile(mat, 0.02),
        -0.05,
        length = 75
      ),
      seq(
        -0.49,
        0.49,
        length = 50
      ),
      seq(
        0.5,
        quantile(mat, 0.98),
        length = 75
      )
    ),
    colorRampPalette(
      c(
        "steelblue4", "steelblue2" , 
        "white", 
        "firebrick2", "firebrick4")
      )(200)
  )
  
  colors
}

