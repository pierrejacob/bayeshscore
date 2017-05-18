#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
#----------------     Some useful functions to manipulate (RCpp) trees       -----------#
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
#'@rdname tree_getattributes
#'@title tree_getattributes
#'@description Get all the attributes (as a \code{list}) of an RCpp tree object. The attributes are respectively:
#'number of x-particles (\code{N}); upper bound on total size of tree (\code{M});
#'dimension of each x-particle (\code{dimx}); 'number of time "insert" has been called (\code{nsteps});
#'vector of ancestor indices (\code{a_star}); vector of offspring counts (\code{o_star});
#'vector of particles (\code{x_star}); vector of leaf indices (\code{l_star}).
#'@export
tree_getattributes = function(tree) {
  return (list(N = tree$N, M = tree$M, dimx = tree$dimx, nsteps = tree$nsteps, a_star = tree$a_star,
               o_star = tree$o_star, x_star = tree$x_star, l_star = tree$l_star))
}

#'@rdname tree_reconstruct
#'@title tree_reconstruct
#'@description Reconstruct an RCpp tree object from a \code{list} of attributes. The attributes are respectively:
#'number of x-particles (\code{N}); upper bound on total size of tree (\code{M});
#'dimension of each x-particle (\code{dimx}); 'number of time "insert" has been called (\code{nsteps});
#'vector of ancestor indices (\code{a_star}); vector of offspring counts (\code{o_star});
#'vector of particles (\code{x_star}); vector of leaf indices (\code{l_star}).
#'@export
tree_reconstruct = function(tree_attributes_list) {
  tree <- new(TreeClass, tree_attributes_list$N, tree_attributes_list$M, tree_attributes_list$dimx)
  tree$reconstruct(tree_attributes_list$N, tree_attributes_list$M, tree_attributes_list$dimx,
                   tree_attributes_list$nsteps, tree_attributes_list$a_star, tree_attributes_list$o_star,
                   tree_attributes_list$x_star, tree_attributes_list$l_star)
  return (tree)
}
#----------------------------------------------------------------------------------------
#-----------------------------  "vectoried" wrappers ------------------------------------
#----------------------------------------------------------------------------------------
#'@rdname trees_getattributes
#'@title trees_getattributes
#'@description This is a wrapper that takes a \code{list} of trees as input
#'and outputs a \code{list} of tree attributes \code{list} via the function \code{tree_getattributes}
#'@export
trees_getattributes = function(trees) {
  return (lapply(1:length(trees),function(i)tree_getattributes(trees[[i]])))
}

#'@rdname trees_reconstruct
#'@title trees_reconstruct
#'@description This is a wrapper that takes a \code{list} of tree attributes \code{list}
#'and outputs a \code{list} of reconstructed trees via the function \code{tree_reconstruct}
#'@export
trees_reconstruct = function(trees_attributes_list) {
  return (lapply(1:length(trees_attributes_list),function(i)tree_reconstruct(trees_attributes_list[[i]])))
}
