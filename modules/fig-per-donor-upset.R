# all regulated genes, regardless of direction

p1 <- plot_upset(
  logicalMatrix = per_donor_data[["logical_matrix"]], 
  plotTitle = "(A) IL1B-regulated",
  minSize = 40
  )

p2 <- plot_upset(
  per_donor_data[["logical_matrix_a5Exclude"]], 
  plotTitle = "(B) IL1B-regulated, primary only",
  minSize = 40
  )

print(
  wrap_plots(p1, p2, nrow = 2)
)

# separated by direction

p1 <- plot_upset(
  logicalMatrix = per_donor_data[["logical_matrix_upreg"]], 
  plotTitle = "(C) IL1B-upregulated"
  )

p2 <- plot_upset(
  logicalMatrix = per_donor_data[["logical_matrix_downreg"]], 
  plotTitle = "(D) IL1B-downregulated"
  )

p3 <- plot_upset(
  per_donor_data[["logical_matrix_a5Exclude_upreg"]], 
  plotTitle = "(E) IL1B-upreg, primary only")

p4 <- plot_upset(
  per_donor_data[["logical_matrix_a5Exclude_downreg"]], 
  plotTitle = "(F) IL1B-downreg, primary only")

print(
  wrap_plots(p1, p2, p3, p4, nrow = 2)
)

