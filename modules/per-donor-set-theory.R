
# === initialize per donor data object =========================================
per_donor_data <- list()

per_donor_data[["tpm_by_donor"]] <- tpm %>%
  filter(treatment == "IL1B") %>%
  mutate(rep = if_else(celltype == "A549", "mean", rep)) %>%
  mutate(name = paste(celltype, rep)) %>%
  select(Gene, name, rep, celltype, treatment, tpm, log2tpm, fold, log2fold)

per_donor_data[["thresholded_log2fold"]] <- per_donor_data[["tpm_by_donor"]] %>%
  group_by(Gene, name, rep, celltype) %>%
  summarize(log2fold = mean(log2fold), .groups = "drop") %>%
  mutate(meets_threshold = abs(log2fold) >= 1)

# === basic counts and spread on a per-donor basis =============================

# counts: how many threshold-meeting genes in each donor?
counts <- per_donor_data[["thresholded_log2fold"]] %>%
  filter(meets_threshold) %>%
  group_by(name, rep, celltype) %>%
  summarize(n = length(unique(Gene)), .groups = "drop") %>%
  arrange(name) %>%
  mutate(name = factor(name, levels = unique(name)))

# spread: sd and %meetingThreshold for each gene, in each celltype
spread <- per_donor_data[["thresholded_log2fold"]] %>%
  group_by(Gene, celltype) %>%
  filter(any(meets_threshold)) %>%
  ungroup() %>%
  select(-log2fold) %>%
  distinct() %>%
  left_join(
    per_donor_data[["tpm_by_donor"]]
  ) %>%
  group_by(Gene, celltype) %>%
  summarize(
    sd = sd(log2fold),
    n_meets_threshold = length(Gene[abs(log2fold) >= 1]),
    nreps = length(Gene),
    percent_meets_threshold = 100 * n_meets_threshold / nreps,
    .groups = "drop"
  )

# === calculate overlaps for Venns and UpSets ==================================

# sets: which genes meet the log2fold threshold in each donor?
# for two-tailed, upregulated, and downregulated cases

sets <- per_donor_data[["thresholded_log2fold"]] %>%
  filter(abs(log2fold) >= 1) %>%
  group_by(name, celltype, rep) %>%
  summarize(genes = list(unique(Gene)), .groups = "drop")

sets_upreg <- per_donor_data[["thresholded_log2fold"]] %>%
  filter(log2fold >= 1) %>%
  group_by(name, celltype, rep) %>%
  summarize(genes = list(unique(Gene)), .groups = "drop")

sets_downreg <- per_donor_data[["thresholded_log2fold"]] %>%
  filter(log2fold <= -1) %>%
  group_by(name, celltype, rep) %>%
  summarize(genes = list(unique(Gene)), .groups = "drop")


# a549 genes: the set against which all other sets will be compared
# for two-tailed, upregulated, and downregulated cases

a549_genes <- sets %>%
  filter(name == "A549 mean") %>%
  pull(genes) %>%
  .[[1]]

a549_genes_upreg <- sets_upreg %>%
  filter(name == "A549 mean") %>%
  pull(genes) %>%
  .[[1]]

a549_genes_downreg <- sets_downreg %>%
  filter(name == "A549 mean") %>%
  pull(genes) %>%
  .[[1]]

# overlaps: how does each donor interact with the A549 consensus?
# for two-tailed, upregulated, and downregulated cases

overlaps <- sets %>%
  mutate(
    n_genes = map_int(genes, length),
    overlap_genes = map(genes, ~ intersect(a549_genes, .x)),
    overlap = map_int(genes, ~ length(intersect(a549_genes, .x))),
    union = map_int(genes, ~ length(union(a549_genes, .x))),
    jaccard = overlap / union,
    title = paste0(
      name,
      "\nJ=", round(jaccard, 2),
      " | n=", union
    )
  )

overlaps_upreg <- sets_upreg %>%
  mutate(
    n_genes = map_int(genes, length),
    overlap_genes = map(genes, ~ intersect(a549_genes_upreg, .x)),
    overlap = map_int(genes, ~ length(intersect(a549_genes_upreg, .x))),
    union = map_int(genes, ~ length(union(a549_genes_upreg, .x))),
    jaccard = overlap / union,
    title = paste0(
      name,
      "\nJ=", round(jaccard, 2),
      " | n=", union
    )
  )

overlaps_downreg <- sets_downreg %>%
  mutate(
    n_genes = map_int(genes, length),
    overlap_genes = map(genes, ~ intersect(a549_genes_downreg, .x)),
    overlap = map_int(genes, ~ length(intersect(a549_genes_downreg, .x))),
    union = map_int(genes, ~ length(union(a549_genes_downreg, .x))),
    jaccard = overlap / union,
    title = paste0(
      name,
      "\nJ=", round(jaccard, 2),
      " | n=", union
    )
  )

# === express overlaps as logical matrices as input for plotting functions =====

# logical matrix: is each gene present (1/0) in each donor?
# for two-tailed, upregulated, and downregulated cases

logical_matrix <- overlaps %>%
  select(name, genes) %>%
  unnest(genes) %>%
  mutate(value = 1) %>%
  pivot_wider(
    names_from = "name",
    values_from = "value",
    values_fill = 0
  ) %>%
  column_to_rownames("genes")

logical_matrix_upreg <- overlaps_upreg %>%
  select(name, genes) %>%
  unnest(genes) %>%
  mutate(value = 1) %>%
  pivot_wider(
    names_from = "name",
    values_from = "value",
    values_fill = 0
  ) %>%
  column_to_rownames("genes")

logical_matrix_downreg <- overlaps_downreg %>%
  select(name, genes) %>%
  unnest(genes) %>%
  mutate(value = 1) %>%
  pivot_wider(
    names_from = "name",
    values_from = "value",
    values_fill = 0
  ) %>%
  column_to_rownames("genes")

# an extra matrix where A549 and related genes are excluded

a5_exclude <- logical_matrix %>%
  select(-`A549 mean`) %>%
  filter(!if_all(everything(), ~ .x == 0))

a5_exclude_upreg <- logical_matrix_upreg %>%
  select(-`A549 mean`) %>%
  filter(!if_all(everything(), ~ .x == 0))

a5_exclude_downreg <- logical_matrix_downreg %>%
  select(-`A549 mean`) %>%
  filter(!if_all(everything(), ~ .x == 0))


# === write calculations to main per-donor data object

per_donor_data[["counts"]] <- counts
per_donor_data[["spread"]] <- spread

per_donor_data[["overlaps"]] <- overlaps
per_donor_data[["overlaps_upreg"]] <- overlaps_upreg
per_donor_data[["overlaps_downreg"]] <- overlaps_downreg

per_donor_data[["overlaps_summary"]] <- rbind(
  per_donor_data[["overlaps"]] %>%
    mutate(criterion = "reg", .before = 1),
  per_donor_data[["overlaps_upreg"]] %>%
    mutate(criterion = "upreg", .before = 1),
  per_donor_data[["overlaps_downreg"]] %>%
    mutate(criterion = "downreg", .before = 1)
)

per_donor_data[["logical_matrix"]] <- logical_matrix
per_donor_data[["logical_matrix_upreg"]] <- logical_matrix_upreg
per_donor_data[["logical_matrix_downreg"]] <- logical_matrix_downreg

per_donor_data[["logical_matrix_a5Exclude"]] <- a5_exclude
per_donor_data[["logical_matrix_a5Exclude_upreg"]] <- a5_exclude_upreg
per_donor_data[["logical_matrix_a5Exclude_downreg"]] <- a5_exclude_downreg





