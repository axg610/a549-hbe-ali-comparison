myColors <- c(
  "NS" = "grey30", 
  "IL1B" = "#FC4E07", 
  "Bud" = "#E7B800", 
  "I+B" = "#2E9FDF",
  "IB" = "#2E9FDF"
)

cleanup_objects <- function(keep = readLines("cleanup_ignore.txt")) {
  all_objs <- ls(envir = .GlobalEnv)
  to_remove <- setdiff(all_objs, keep)
  rm(list = to_remove, envir = .GlobalEnv)
  # gc()
}