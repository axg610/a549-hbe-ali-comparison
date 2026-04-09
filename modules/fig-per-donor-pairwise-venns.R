# plots

plot_pairwise_venns <- function(
    inputOverlaps, 
    referenceGroup = per_donor_data$a549_genes,
    plotTitle = "title here"){
  
  # take an overlaps dataframe from above and create an array of pairwise
  # venn digrams using ggVennDiagram.
  
  # must specify the reference group (the bottom venn lobe)
  
  plots <- inputOverlaps %>%
    mutate(
      p = map2(
        genes,
        title,
        ~ ggVennDiagram(
          list(
            " " = referenceGroup,
            "  " = .x
          ),
          label_alpha = 0,
          label_size = 3
        ) +
          ggtitle(.y) +
          scale_fill_gradient(low = "white", high = "white") +
          theme(
            legend.position = "none",
            plot.title = element_text(size = 10, hjust = 0.5)
          )
      )
    )
  
  print(
    wrap_plots(plots$p, ncol = 6) +
      plot_annotation(title = plotTitle)
  )
  
}

plot_pairwise_venns(
  per_donor_data$overlaps,         
  per_donor_data$a549_genes, 
  "(A) Pairwise Venns for IL1B-reg genes: abs(log2fold) >= 1, donor (top) vs A549 (bottom)"
)

plot_pairwise_venns(
  per_donor_data$overlaps_upreg,   
  per_donor_data$a549_genes_upreg, 
  "(B) Pairwise Venns for IL1B-upreg genes: log2fold >= 1, donor (top) vs A549 (bottom)"
)

plot_pairwise_venns(
  per_donor_data$overlaps_downreg, 
  per_donor_data$a549_genes_downreg, 
  "(C) Pairwise Venns for IL1B-downreg genes: log2fold <= -1, donor (top) vs A549 (bottom)"
)