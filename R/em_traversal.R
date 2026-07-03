#' EM-Specific Tree Traversal
#'
#' Implements the E-step traversal for the EM algorithm. This computes:
#' - exps: Conditional expectations at all nodes
#' - vars: Conditional variances at all nodes
#' - covars: Covariances between parent-child pairs (crucial for M-step)
#'
#' Based on Bastide et al. (2018) PhylogeneticEM approach.
#'
#' @name em-traversal
NULL


#' EM E-Step: Compute Conditional Expectations and Covariances
#'
#' Runs the forward and backward passes to compute quantities needed
#' for the EM M-step. Unlike the gradient-focused bastide_traversal,
#' this computes parent-child covariances needed for variance estimation.
#'
#' @param tree A phylo object (should be postorder)
#' @param trait_data Data frame with species in first column, traits in rest
#' @param R Phylogenetic rate matrix
#' @param S Residual covariance matrix
#' @param mu Current mean estimate
#'
#' @return List with:
#'   \item{exps}{Conditional expectations at all nodes (total x p matrix)}
#'   \item{vars}{Conditional variances at all nodes (list of p x p matrices)}
#'   \item{covars}{Covariances between child and parent (list of p x p matrices)}
#'   \item{logL}{Log-likelihood}
#'   \item{mu}{Mean estimate}
#'
#' @export
em_estep <- function(tree, trait_data, R, S, mu) {

  # Extract data
  Y <- as.matrix(trait_data[, -1, drop = FALSE])
  species <- trait_data[, 1]

  n <- length(tree$tip.label)      # number of species
  Nnode <- tree$Nnode              # number of internal nodes
  m_nodes <- n + Nnode             # total tree nodes
  N <- nrow(Y)                     # number of observations
  total <- N + m_nodes             # total entities
  p <- ncol(Y)                     # number of traits
  Ip <- diag(p)

  # Map observations to species
  node_obs <- match(species, tree$tip.label)
  parent_nodes_of_observations <- node_obs

  # Edge info
  nedge <- nrow(tree$edge)
  anc <- tree$edge[, 1]
  des <- tree$edge[, 2]
  lens <- tree$edge.length

  # Which edges lead to each species (for observations)
  edges_leading_to_obs_k <- match(node_obs, tree$edge[, 2])

  # Initialize storage for forward pass
  # pA_orig: original accumulated precision (preserved for backward pass)
  # pA_prop: propagated precision (used for tree propagation)
  pA_orig <- rep(list(matrix(0, p, p)), total)  # Original precision
  pA_prop <- rep(list(matrix(0, p, p)), total)  # Propagated precision
  pY_orig <- matrix(0, total, p)                 # Original precision-weighted sum
  pY_prop <- matrix(0, total, p)                 # Propagated precision-weighted sum
  cond_exp <- matrix(0, total, p)                # Conditional expectations from forward

  # Initialize storage for backward pass
  exps <- matrix(0, total, p)              # Conditional expectations
  vars <- rep(list(matrix(0, p, p)), total) # Conditional variances
  covars <- rep(list(matrix(0, p, p)), total) # Parent-child covariances

  #---------------------------------------------------------------------------
  # FORWARD PASS: Process observations, then edges (tips to root)
  #---------------------------------------------------------------------------

  # Process observations first
  for (i in 1:N) {
    anc_i <- des[edges_leading_to_obs_k[i]] + N  # Species node (offset by N)
    des_i <- i  # Observation index

    # Handle missing data
    obs_traits <- which(!is.na(Y[i, ]))
    if (length(obs_traits) == 0) next

    # For observed traits, compute precision contribution
    if (length(obs_traits) == p) {
      # Fully observed
      S_inv <- solve(S)
      pA_orig[[des_i]] <- S_inv
      pA_prop[[des_i]] <- S_inv
      cond_exp[des_i, ] <- Y[i, ]
      pY_orig[des_i, ] <- S_inv %*% Y[i, ]
      pY_prop[des_i, ] <- pY_orig[des_i, ]
    } else {
      # Partially observed - use submatrix
      S_sub <- S[obs_traits, obs_traits, drop = FALSE]
      S_sub_inv <- solve(S_sub)
      pA_orig[[des_i]][obs_traits, obs_traits] <- S_sub_inv
      pA_prop[[des_i]] <- pA_orig[[des_i]]
      cond_exp[des_i, obs_traits] <- Y[i, obs_traits]
      pY_orig[des_i, obs_traits] <- S_sub_inv %*% Y[i, obs_traits]
      pY_prop[des_i, ] <- pY_orig[des_i, ]
    }

    # Propagate to species node
    pA_orig[[anc_i]] <- pA_orig[[anc_i]] + pA_prop[[des_i]]
    pA_prop[[anc_i]] <- pA_orig[[anc_i]]
    pY_orig[anc_i, ] <- pY_orig[anc_i, ] + pY_prop[des_i, ]
    pY_prop[anc_i, ] <- pY_orig[anc_i, ]
  }

  # Process tree edges (tips to root, i.e., postorder)
  for (i in 1:nedge) {
    anc_i <- anc[i] + N  # Parent node (offset)
    des_i <- des[i] + N  # Child node (offset)
    t_e <- lens[i]

    # Build covariance matrix for this edge
    Sigma_e <- t_e * R

    # Compute conditional expectation at child from forward pass data
    if (any(diag(pA_orig[[des_i]]) != 0)) {
      nonzero <- which(diag(pA_orig[[des_i]]) != 0)
      cond_var_sub <- solve(pA_orig[[des_i]][nonzero, nonzero, drop = FALSE])
      cond_exp[des_i, nonzero] <- cond_var_sub %*% pY_orig[des_i, nonzero]
    }

    if (any(diag(pA_orig[[des_i]]) != 0)) {
      # Compute (I + Sigma * P)^{-1}
      itpa <- Ip + Sigma_e %*% pA_orig[[des_i]]
      itpainv <- solve(itpa)

      # Propagated precision (for tree traversal), stored separately
      pA_prop[[des_i]] <- pA_orig[[des_i]] %*% itpainv
      pY_prop[des_i, ] <- pY_orig[des_i, ] %*% itpainv
    }

    # Propagate to parent
    pA_orig[[anc_i]] <- pA_orig[[anc_i]] + pA_prop[[des_i]]
    pA_prop[[anc_i]] <- pA_orig[[anc_i]]
    pY_orig[anc_i, ] <- pY_orig[anc_i, ] + pY_prop[des_i, ]
    pY_prop[anc_i, ] <- pY_orig[anc_i, ]
  }

  # Compute root expectation
  root_idx <- (n + 1) + N  # Root is node n+1, offset by N
  if (any(diag(pA_orig[[root_idx]]) != 0)) {
    root_var <- solve(pA_orig[[root_idx]])
    cond_exp[root_idx, ] <- root_var %*% pY_orig[root_idx, ]
    vars[[root_idx]] <- root_var
  }
  exps[root_idx, ] <- mu  # Use provided mu

  #---------------------------------------------------------------------------
  # BACKWARD PASS: Compute expectations and covariances (root to tips)
  #---------------------------------------------------------------------------

  # Process tree edges in reverse (root to tips)
  for (i in nedge:1) {
    anc_i <- anc[i] + N
    des_i <- des[i] + N
    t_e <- lens[i]

    Sigma_e <- t_e * R

    # Recompute itpainv using ORIGINAL pA (not propagated)
    if (any(diag(pA_orig[[des_i]]) != 0)) {
      itpa <- Ip + Sigma_e %*% pA_orig[[des_i]]
      itpainv <- solve(itpa)

      # Covariance between child and parent: Cov(Y_des, Y_anc | data)
      covars[[des_i]] <- t(itpainv %*% vars[[anc_i]])

      # Conditional expectation at child
      new_p <- pA_orig[[des_i]] %*% itpainv
      exps[des_i, ] <- itpainv %*% exps[anc_i, ] +
                       Sigma_e %*% new_p %*% cond_exp[des_i, ]

      # Conditional variance at child
      vars[[des_i]] <- itpainv %*% Sigma_e + itpainv %*% covars[[des_i]]
    } else {
      exps[des_i, ] <- exps[anc_i, ]
      covars[[des_i]] <- vars[[anc_i]]
      vars[[des_i]] <- vars[[anc_i]] + Sigma_e
    }
  }

  # Process observations (species to obs)
  for (i in 1:N) {
    anc_i <- des[edges_leading_to_obs_k[i]] + N  # Species node
    des_i <- i  # Observation

    obs_traits <- which(!is.na(Y[i, ]))
    miss_traits <- which(is.na(Y[i, ]))

    if (length(obs_traits) == p) {
      # Fully observed - expectation equals observation
      exps[des_i, ] <- Y[i, ]
      vars[[des_i]] <- matrix(0, p, p)
      covars[[des_i]] <- matrix(0, p, p)
    } else if (length(obs_traits) > 0) {
      # Partially observed
      # Observed traits are fixed
      exps[des_i, obs_traits] <- Y[i, obs_traits]

      # Missing traits: conditional on observed
      if (length(miss_traits) > 0) {
        S_mm <- S[miss_traits, miss_traits, drop = FALSE]
        S_mo <- S[miss_traits, obs_traits, drop = FALSE]
        S_oo <- S[obs_traits, obs_traits, drop = FALSE]
        S_oo_inv <- solve(S_oo)

        # Regression of missing on observed
        resid <- Y[i, obs_traits] - exps[anc_i, obs_traits]
        exps[des_i, miss_traits] <- exps[anc_i, miss_traits] +
                                     S_mo %*% S_oo_inv %*% resid

        # Conditional variance
        vars[[des_i]][miss_traits, miss_traits] <-
          S_mm - S_mo %*% S_oo_inv %*% t(S_mo)
      }

      # Covariance with parent
      covars[[des_i]] <- S %*% solve(S + vars[[anc_i]]) %*% vars[[anc_i]]
    } else {
      # Fully missing - inherit from parent
      exps[des_i, ] <- exps[anc_i, ]
      vars[[des_i]] <- S + vars[[anc_i]]
      covars[[des_i]] <- vars[[anc_i]]
    }
  }

  list(
    exps = exps,
    vars = vars,
    covars = covars,
    pA_orig = pA_orig,
    pY_orig = pY_orig,
    cond_exp = cond_exp,
    tree = tree,
    trait_data = trait_data,
    N = N,
    n = n,
    m_nodes = m_nodes,
    total = total,
    p = p,
    node_obs = node_obs,
    parent_nodes_of_observations = parent_nodes_of_observations,
    edges_leading_to_obs_k = edges_leading_to_obs_k
  )
}


#' EM M-Step: Update Parameters
#'
#' Given the E-step results, update R, S, and mu using closed-form solutions.
#'
#' @param estep_result Result from em_estep()
#' @param tree A phylo object
#'
#' @return List with updated R, S, mu
#'
#' @export
em_mstep <- function(estep_result, tree) {

  exps <- estep_result$exps
  vars <- estep_result$vars
  covars <- estep_result$covars
  N <- estep_result$N
  n <- estep_result$n
  m_nodes <- estep_result$m_nodes
  p <- estep_result$p
  node_obs <- estep_result$node_obs
  parent_nodes_of_observations <- estep_result$parent_nodes_of_observations

  nedge <- nrow(tree$edge)
  anc <- tree$edge[, 1]
  des <- tree$edge[, 2]

  #---------------------------------------------------------------------------
  # Update mu: weighted average of root children expectations
  #---------------------------------------------------------------------------

  root_edges <- which(anc == (n + 1))
  if (length(root_edges) > 0) {
    weights <- 1 / tree$edge.length[root_edges]
    weights <- weights / sum(weights)

    mu_new <- numeric(p)
    for (i in seq_along(root_edges)) {
      edge_idx <- root_edges[i]
      child_node <- des[edge_idx]
      child_idx <- child_node + N
      mu_new <- mu_new + weights[i] * exps[child_idx, ]
    }
  } else {
    root_idx <- (n + 1) + N
    mu_new <- exps[root_idx, ]
  }

  #---------------------------------------------------------------------------
  # Update R: phylogenetic rate matrix
  #---------------------------------------------------------------------------

  # sum_exp: outer product of differences
  # sum_var: variance of differences

  sum_exp <- matrix(0, p, p)
  sum_var <- matrix(0, p, p)

  for (e in 1:nedge) {
    parent_node <- anc[e]
    child_node <- des[e]
    t_e <- tree$edge.length[e]

    parent_idx <- parent_node + N
    child_idx <- child_node + N

    # Difference in expectations
    diff_exp <- exps[child_idx, ] - exps[parent_idx, ]
    sum_exp <- sum_exp + outer(diff_exp, diff_exp) / t_e

    # Variance of difference
    var_child <- vars[[child_idx]]
    var_parent <- vars[[parent_idx]]
    cov_cp <- covars[[child_idx]]

    var_diff <- var_child + var_parent - cov_cp - t(cov_cp)
    sum_var <- sum_var + var_diff / t_e
  }

  R_new <- (sum_exp + sum_var) / (m_nodes - 1)
  R_new <- (R_new + t(R_new)) / 2  # Symmetrize

  #---------------------------------------------------------------------------
  # Update S: residual covariance matrix
  #---------------------------------------------------------------------------

  sum_exp_ind <- matrix(0, p, p)
  sum_var_ind <- matrix(0, p, p)

  for (i in 1:N) {
    species_node <- parent_nodes_of_observations[i]
    species_idx <- species_node + N
    obs_idx <- i

    # Difference: species - observation
    diff_exp <- exps[species_idx, ] - exps[obs_idx, ]
    sum_exp_ind <- sum_exp_ind + outer(diff_exp, diff_exp)

    # Variance of difference
    var_species <- vars[[species_idx]]
    var_obs <- vars[[obs_idx]]
    cov_so <- covars[[obs_idx]]

    var_diff <- var_species + var_obs - cov_so - t(cov_so)
    sum_var_ind <- sum_var_ind + var_diff
  }

  S_new <- (sum_exp_ind + sum_var_ind) / N
  S_new <- (S_new + t(S_new)) / 2  # Symmetrize

  #---------------------------------------------------------------------------
  # Ensure positive definiteness
  #---------------------------------------------------------------------------

  eig_R <- eigen(R_new, symmetric = TRUE)
  if (any(eig_R$values <= 1e-10)) {
    eig_R$values[eig_R$values <= 1e-10] <- 1e-10
    R_new <- eig_R$vectors %*% diag(eig_R$values) %*% t(eig_R$vectors)
  }

  eig_S <- eigen(S_new, symmetric = TRUE)
  if (any(eig_S$values <= 1e-10)) {
    eig_S$values[eig_S$values <= 1e-10] <- 1e-10
    S_new <- eig_S$vectors %*% diag(eig_S$values) %*% t(eig_S$vectors)
  }

  list(R = R_new, S = S_new, mu = mu_new)
}


#' Run EM Algorithm
#'
#' Iteratively runs E-step and M-step until convergence.
#'
#' @param tree A phylo object (should be postorder)
#' @param trait_data Data frame with species and traits
#' @param R_init Initial R matrix (default: identity)
#' @param S_init Initial S matrix (default: identity)
#' @param mu_init Initial mean (default: column means)
#' @param max_iter Maximum iterations
#' @param tol Convergence tolerance (on log-likelihood change)
#' @param verbose Print progress
#'
#' @return List with R, S, mu, logL, converged, iterations
#'
#' @export
em_fit <- function(tree, trait_data,
                    R_init = NULL, S_init = NULL, mu_init = NULL,
                    max_iter = 100, tol = 1e-6, verbose = FALSE) {

  # Dimensions
  Y <- as.matrix(trait_data[, -1, drop = FALSE])
  p <- ncol(Y)

  # Defaults
  if (is.null(R_init)) R_init <- diag(p)
  if (is.null(S_init)) S_init <- diag(p)
  if (is.null(mu_init)) mu_init <- colMeans(Y, na.rm = TRUE)

  R <- R_init
  S <- S_init
  mu <- mu_init

  logL_prev <- -Inf
  logL_history <- numeric(max_iter)
  converged <- FALSE

  for (iter in 1:max_iter) {
    # E-step
    estep <- em_estep(tree, trait_data, R, S, mu)

    # M-step
    mstep <- em_mstep(estep, tree)
    R <- mstep$R
    S <- mstep$S
    mu <- mstep$mu

    # Compute log-likelihood
    logL <- loglik_bastide(R, S, tree, trait_data)
    logL_history[iter] <- logL

    if (verbose) {
      cat("EM iter", iter, ": logL =", round(logL, 4), "\n")
    }

    # Check convergence
    if (abs(logL - logL_prev) < tol) {
      converged <- TRUE
      break
    }

    logL_prev <- logL
  }

  list(
    R = R,
    S = S,
    mu = mu,
    logL = logL,
    logL_history = logL_history[1:iter],
    converged = converged,
    iterations = iter
  )
}
