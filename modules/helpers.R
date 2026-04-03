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