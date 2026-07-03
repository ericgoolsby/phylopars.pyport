#' Example phylogenetic tree
#'
#' A phylogenetic tree with 25 species for demonstrating phylopars functions.
#'
#' @format A phylo object with 25 tips
#' @examples
#' data(example_tree)
#' plot(example_tree)
"example_tree"


#' Example trait data
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
