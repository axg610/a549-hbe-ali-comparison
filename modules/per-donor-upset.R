


# all regulated genes, regardless of direction

p1 <- plot_upset(
  logicalMatrix = per_donor_data[["logical_matrix"]], 
  plotTitle = "(A) IL1B-regulated",
  minSize = 35
  )

# separated by direction

p2 <- plot_upset(
  logicalMatrix = per_donor_data[["logical_matrix_upreg"]], 
  plotTitle = "(B) IL1B-upregulated",
  minSize = 50
)

p3 <- plot_upset(
  logicalMatrix = per_donor_data[["logical_matrix_downreg"]], 
  plotTitle = "(C) IL1B-downregulated",
  minSize = 50
)

print(
  wrap_plots(
    p1,
    wrap_plots(
      p2, 
      p3, 
      nrow = 1
    ),
    nrow = 2
  )
)

