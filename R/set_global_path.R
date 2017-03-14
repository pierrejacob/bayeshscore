#'@export
set_global_path <- function(){
  if (Sys.getenv("USER") == "pierre"){
    rdatapath <<- "~/Dropbox/HyvarinenScore/RData/"
  } else {
    rdatapath <<- getwd()
  }
  cat("rdatapath set to", rdatapath, "\n")
}
