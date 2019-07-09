file_find_replace <- function(filepath, pattern, replacement) {
  file_contents <- readLines(filepath)
  updated_contents <- gsub(x = file_contents, pattern = pattern, replacement = replacement)
  cat(updated_contents, file = filepath, sep = "\n")
}


file_find_replace('~/R/Clonal sizes/stop3.R',"_5","_6")