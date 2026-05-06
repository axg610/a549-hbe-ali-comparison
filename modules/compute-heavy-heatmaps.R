
# === only initialize container if it doesn't exist ============================

if (!exists("clustered_matrices")) {
  clustered_matrices <- list()
}


# === only compute matrix if missing ===========================================

if (!"all_donors_union_A5mean" %in% names(clustered_matrices)) {
  
  # the union of all 14621 genes regulated in any donor, log2fold, hclust by row
  # A549s collapsed into average of 4 reps.
  
  union_genes <- per_donor_data$thresholded_log2fold %>% 
    filter(meets_threshold) %>% 
    pull(Gene)
  
  all_donors_union_A5mean <- per_donor_data$thresholded_log2fold %>%
    filter(
      Gene %in% union_genes
    ) %>%
    mutate(name = paste(celltype, rep)) %>%
    select(Gene, name, log2fold) %>%
    pivot_wider(names_from = "name", values_from = "log2fold") %>%
    column_to_rownames("Gene") %>%
    as.matrix() %>%
    .[hclust(dist(.), method = "ward.D2")$order, ]
  
  clustered_matrices[["all_donors_union_A5mean"]] <- all_donors_union_A5mean
}


if (!"all_donors_union" %in% names(clustered_matrices)) {
  
  # the union of all 14621 genes regulated in any donor, log2fold
  # hclust rows and columns
  # A549s shown as separate reps
  
  union_genes <- per_donor_data$thresholded_log2fold_A549separate %>% 
    filter(meets_threshold) %>% 
    pull(Gene)
  
  all_donors_union <- per_donor_data$thresholded_log2fold_A549separate %>%
    filter(
      Gene %in% union_genes
    ) %>%
    mutate(name = paste(celltype, rep)) %>%
    select(Gene, name, log2fold) %>%
    pivot_wider(names_from = "name", values_from = "log2fold") %>%
    column_to_rownames("Gene") %>%
    as.matrix() %>%
    .[
      hclust(dist(.), method = "ward.D2")$order,
    ]
  
  # rowhc.all_donors_union <- hclust(dist(all_donors_union), method = "ward.D2")
  # colhc.all_donors_union <- hclust(dist(t(all_donors_union)), method = "ward.D2")
  
  # all_donors_union <- all_donors_union %>%
  #   .[
  #     rowhc.all_donors_union$order,
  #     colhc.all_donors_union$order
  #   ]
  
  clustered_matrices[["all_donors_union"]] <- all_donors_union
  # clustered_matrices[["rowhc.all_donors_union"]] <- rowhc.all_donors_union
  # clustered_matrices[["colhc.all_donors_union"]] <- colhc.all_donors_union
}

if (!"all_donors_union_comparativeTPMs" %in% names(clustered_matrices)) {

  # the union of all 14621 genes regulated in any donor, log2fold
  # hclust rows and columns
  # A549s shown as separate reps
  # BUT, we display the TPMs instead

  union_genes <- per_donor_data$thresholded_log2fold_A549separate %>%
    filter(meets_threshold) %>%
    pull(Gene)

  all_donors_union_comparativeTPMs <- tpm %>%
    filter(
      Gene %in% union_genes,
      treatment == "NS"
    ) %>%
    mutate(name = paste(celltype, rep)) %>%
    select(Gene, name, log2tpm) %>%
    pivot_wider(names_from = "name", values_from = "log2tpm") %>%
    column_to_rownames("Gene") %>%
    as.matrix() %>%
    .[
      rownames(clustered_matrices[["all_donors_union"]]),
      colnames(clustered_matrices[["all_donors_union"]])
    ]

  clustered_matrices[["all_donors_union_comparativeTPMs"]] <- all_donors_union_comparativeTPMs
}