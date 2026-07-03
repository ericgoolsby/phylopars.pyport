#' Example 25-species phylogenetic tree
#'
#' A phylogenetic tree with 25 species for demonstrating phylopars functions.
#'
#' @format A phylo object with 25 tips
#' @examples
#' data(example_tree)
#' plot(example_tree)
"example_tree"


#' Example 25-species trait data
#'
#' Trait data for 25 species with 4 continuous traits. Includes missing values
#' and multiple observations per species (3 replicates each).
#'
#' @format A data frame with 75 rows (25 species x 3 replicates) and 5 columns:
#' \describe{
#'   \item{species}{Species name}
#'   \item{V1}{First trait}
#'   \item{V2}{Second trait}
#'   \item{V3}{Third trait}
#'   \item{V4}{Fourth trait}
#' }
#'
#' @examples
#' data(example_traits)
#' head(example_traits)
"example_traits"

#' Example 50-species phylogenetic tree
#'
#' A phylogenetic tree with 50 species for demonstrating phylopars functions.
#'
#' @format A phylo object with 50 tips
#' @examples
#' data(tree)
#' plot(tree)
"tree"

#' 50 LOO-CV Trees
#'
#' A named list of phylogenetic trees, each corresponding to a tree where a single species was dropped from `tree`.
#' The order of dropped species corresponds to the species vector `tree$tip.label`.
#'
#' @format A list of phylo objects
#' @examples
#' data(loo_trees)
#' plot(loo_trees$t1)
"loo_trees"

#' Example trait data
#'
#' Trait data for 50 species with 4 continuous traits. Includes missing values
#' and multiple observations per species (3 replicates each). This dataset has a 
#' true interspecific trait covariance of diag(c(1,10,10,1)), a true intraspecific
#' trait covariance of diag(c(1,1,1,1)), and a true root value of c(0,0,0,0). Useful
#' for diagnosing convergence on the correct values (e.g., trait covariance estimates
#' should differ from the true variances by no more than 2 in either direction).
#'
#' @format A data frame with 75 rows (25 species x 3 replicates) and 5 columns:
#' \describe{
#'   \item{species}{Species name}
#'   \item{V1}{First trait}
#'   \item{V2}{Second trait}
#'   \item{V3}{Third trait}
#'   \item{V4}{Fourth trait}
#' }
#'
#' @examples
#' data(Y)
#' head(Y)
"Y"