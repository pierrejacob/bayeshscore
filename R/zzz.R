#------------------------------------------------------------------------------#
#-- This automatically calls the (RCpp) Tree module upon package loading  -----#
#------------------------------------------------------------------------------#
.onLoad <- function(libname, pkgname) {
  cat("module_tree loaded: TreeClass is now available\n")
  module_tree <<- Module("module_tree", PACKAGE = "bayeshscore")
  TreeClass <<- module_tree$Tree
}

