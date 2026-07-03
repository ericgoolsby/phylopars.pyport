#' Ancestral State Reconstruction
#'
#' Compute ancestral state estimates, standard errors, and confidence intervals
#' for all internal nodes of the phylogeny using the Bastide et al. (2021)
#' tree traversal algorithm.
#'
#' @param fit A fitted phylopars object (output of \code{\link{phylopars}})
#' @param tree A phylo object (the same tree used to fit the model)
#' @param trait_data Data frame with species in first column, traits in remaining
#'   columns (the same data used to fit the model)
#' @param alpha Significance level for confidence intervals (default 0.05 for
#'   95\% CI)
#'
#' @return A list with class "anc_recon" containing:
#'   \item{estimates}{Matrix of trait estimates at internal nodes (nodes x traits)}
#'   \item{se}{Matrix of standard errors (nodes x traits)}
#'   \item{lower}{Matrix of CI lower bounds (nodes x traits)}
#'   \item{upper}{Matrix of CI upper bounds (nodes x traits)}
#'   \item{node_ids}{Integer vector of internal node IDs (ape convention:
#'     tips are 1:n, internal nodes are (n+1):(2n-1))}
#'   \item{variances}{List of full variance-covariance matrices at each internal
#'     node (for multivariate inference)}
#'
#' @details
#' This function re-runs the Bastide tree traversal (forward + backward pass)
#' using the fitted model parameters to compute conditional means and variances
#' at all internal nodes.
#'
#' For each internal node, the conditional mean M gives the ancestral state
#' estimate, and the diagonal of the conditional variance V gives the per-trait
#' variances. Standard errors are sqrt(diag(V)), and confidence intervals are
#' computed as estimate +/- qnorm(1 - alpha/2) * SE.
#'
#' @examples
#' \dontrun{
#' library(ape)
#' tree <- rtree(20)
#' Y <- data.frame(species = tree$tip.label, trait1 = rnorm(20))
#' fit <- phylopars(Y, tree)
#' asr <- anc_recon(fit, tree, Y)
#' print(asr)
#' }
#'
#' @export
anc_recon <- function(fit, tree, trait_data, alpha = 0.05) {

  # --- Input validation ---
  if (!is.list(fit) || is.null(fit$Sigma) || is.null(fit$B)) {
    stop("'fit' must be a fitted phylopars object (output of phylopars())")
  }
  if (!inherits(tree, "phylo")) {
    stop("'tree' must be a phylo object")
  }
  if (!is.data.frame(trait_data)) {
    stop("'trait_data' must be a data frame with species in the first column")
  }
  if (alpha <= 0 || alpha >= 1) {
    stop("'alpha' must be between 0 and 1 (exclusive)")
  }

  # --- Prepare tree for the model ---
  working_tree <- reorder(tree, "postorder")

  # Apply evolutionary model transform if needed
  if (fit$model != "BM" && !is.null(fit$evo_model_par)) {
    working_tree <- dist_from_root(working_tree)
    working_tree <- transf_branch_lengths(working_tree, fit$model, fit$evo_model_par)
  }

  # --- Run Bastide traversal with backward pass ---
  trav <- bastide_traversal(
    tree = working_tree,
    trait_data = trait_data,
    R = fit$Sigma,
    S = fit$B,
    mu = as.numeric(fit$mu),
    compute_backward = TRUE
  )

  # --- Extract results for internal nodes ---
  n_tips <- trav$n
  n_nodes <- working_tree$Nnode
  N_obs <- trav$N
  p <- trav$p  # number of traits
  trait_names <- colnames(trait_data)[-1]

  # Internal node IDs in ape convention: (n+1) to (n + Nnode)
  node_ids <- (n_tips + 1):(n_tips + n_nodes)

  # In the traversal indexing, internal node i corresponds to index (i + N_obs)
  # in the M and V arrays (observations 1:N_obs, then nodes shifted by N_obs)
  estimates <- matrix(NA_real_, nrow = n_nodes, ncol = p)
  se_mat <- matrix(NA_real_, nrow = n_nodes, ncol = p)
  var_list <- vector("list", n_nodes)

  for (i in seq_len(n_nodes)) {
    node_id <- node_ids[i]
    # Index into traversal arrays: node_id + N_obs
    trav_idx <- node_id + N_obs

    estimates[i, ] <- trav$M[trav_idx, ]
    V_node <- trav$V[[trav_idx]]
    var_list[[i]] <- V_node
    se_mat[i, ] <- sqrt(pmax(diag(V_node), 0))
  }

  # --- Compute confidence intervals ---
  z <- qnorm(1 - alpha / 2)
  lower <- estimates - z * se_mat
  upper <- estimates + z * se_mat

  # --- Add names ---
  rownames(estimates) <- rownames(se_mat) <- rownames(lower) <- rownames(upper) <- node_ids
  if (length(trait_names) == p) {
    colnames(estimates) <- colnames(se_mat) <- colnames(lower) <- colnames(upper) <- trait_names
  }
  names(var_list) <- node_ids

  # --- Return result ---
  result <- list(
    estimates = estimates,
    se = se_mat,
    lower = lower,
    upper = upper,
    node_ids = node_ids,
    variances = var_list
  )
  class(result) <- "anc_recon"
  result
}


#' Print method for anc_recon objects
#'
#' @param x An anc_recon object
#' @param ... Additional arguments (ignored)
#'
#' @export
print.anc_recon <- function(x, ...) {
  cat("Ancestral State Reconstruction\n")
  cat("==============================\n")
  cat("Internal nodes:", length(x$node_ids), "\n")
  cat("Traits:", ncol(x$estimates), "\n\n")
  cat("Estimates:\n")
  print(round(x$estimates, 4))
  cat("\nStandard Errors:\n")
  print(round(x$se, 4))
  cat("\nConfidence Intervals (lower):\n")
  print(round(x$lower, 4))
  cat("\nConfidence Intervals (upper):\n")
  print(round(x$upper, 4))
  invisible(x)
}
