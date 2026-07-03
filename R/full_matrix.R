#' Build Full Covariance Matrix W, Response y, and Design Matrix X
#'
#' Constructs the full n_obs x n_obs covariance matrix W directly using the
#' Kronecker structure: W = (Sigma x C) + (B x I_ind), where C is the
#' phylogenetic covariance matrix from vcv().
#'
#' This is O(n^2) in the number of observations, so it's slow for large
#' datasets, but conceptually clearer than the tree traversal algorithm.
#' Useful for validation and teaching.
#'
#' @param Sigma Phylogenetic covariance matrix (m x m)
#' @param B Residual/within-species covariance matrix (m x m)
#' @param tree A phylo object
#' @param trait_data Data frame with species in first column, traits in rest
#'
#' @return List with components:
#'   \item{W}{Full covariance matrix (n_obs x n_obs)}
#'   \item{y}{Response vector (observed trait values)}
#'   \item{X}{Design matrix for mean (indicator for each trait)}
#'   \item{obs_info}{Data frame with individual, species, trait indices}
#'   \item{n_obs}{Number of non-missing observations}
#'
#' @export
build_W_y_X <- function(Sigma, B, tree, trait_data) {

  # Get dimensions
  tree_species <- tree$tip.label
  n_species <- length(tree_species)
  trait_matrix <- trait_data[, -1, drop = FALSE]
  m <- ncol(trait_matrix)  # number of traits
  n_ind <- nrow(trait_data)  # number of individuals

  # Ensure species is a factor in tree order
  species_vec <- factor(trait_data[, 1], levels = tree_species)
  species_idx <- as.integer(species_vec)

  # Compute phylogenetic covariance matrix
  C <- ape::vcv(tree, corr = FALSE)

  # Build observation index list (which observations are non-missing)
  obs_info <- data.frame(
    ind = integer(0),
    species = integer(0),
    trait = integer(0)
  )

  y_vals <- numeric(0)
  X_rows <- list()

  for (i in seq_len(n_ind)) {
    for (j in seq_len(m)) {
      val <- trait_matrix[i, j]
      if (!is.na(val)) {
        obs_info <- rbind(obs_info, data.frame(
          ind = i,
          species = species_idx[i],
          trait = j
        ))
        y_vals <- c(y_vals, val)

        # Row of X: indicator for which trait mean
        x_row <- rep(0, m)
        x_row[j] <- 1
        X_rows <- c(X_rows, list(x_row))
      }
    }
  }

  n_obs <- length(y_vals)
  y <- y_vals
  X <- do.call(rbind, X_rows)

  # Build full W matrix
  W <- matrix(0, n_obs, n_obs)

  for (i in seq_len(n_obs)) {
    for (j in i:n_obs) {
      idx_i <- obs_info[i, ]
      idx_j <- obs_info[j, ]

      # Phylogenetic covariance: Sigma[trait_i, trait_j] * C[species_i, species_j]
      cov_phylo <- Sigma[idx_i$trait, idx_j$trait] *
                   C[idx_i$species, idx_j$species]

      # Within-individual covariance: B[trait_i, trait_j] if same individual
      cov_ind <- if (idx_i$ind == idx_j$ind) B[idx_i$trait, idx_j$trait] else 0

      W[i, j] <- cov_phylo + cov_ind
      if (i != j) W[j, i] <- W[i, j]  # Symmetry
    }
  }

  list(W = W, y = y, X = X, obs_info = obs_info, n_obs = n_obs)
}


#' Compute REML Log-Likelihood Using Full Matrix Approach
#'
#' Computes the restricted log-likelihood directly using matrix operations.
#' This is the "naive" O(n^3) approach that inverts W directly.
#'
#' The REML log-likelihood is:
#' logL = -0.5 * [(n-p)*log(2*pi) + log|W| + log|X'W^{-1}X| + (y-X*mu)'W^{-1}(y-X*mu)]
#'
#' where mu = (X'W^{-1}X)^{-1} X'W^{-1}y
#'
#' @param Sigma Phylogenetic covariance matrix (m x m)
#' @param B Residual covariance matrix (m x m)
#' @param tree A phylo object
#' @param trait_data Data frame with species in first column, traits in rest
#' @param ret_mu If TRUE, also return the estimated mean
#'
#' @return Log-likelihood value, or list with logL and mu if ret_mu=TRUE
#'
#' @export
loglik_full_matrix <- function(Sigma, B, tree, trait_data, ret_mu = FALSE) {

  # Build full matrices
  WyX <- build_W_y_X(Sigma, B, tree, trait_data)
  W <- WyX$W
  y <- WyX$y
  X <- WyX$X
  n_obs <- WyX$n_obs
  m <- ncol(X)  # number of traits = number of means to estimate

  # Invert W (use Cholesky for numerical stability)
  W_chol <- tryCatch(
    chol(W),
    error = function(e) NULL
  )

  if (is.null(W_chol)) {
    if (ret_mu) return(list(logL = -1e100, mu = NULL))
    return(-1e100)
  }

  W_inv <- chol2inv(W_chol)
  log_det_W <- 2 * sum(log(diag(W_chol)))

  # Compute X'W^{-1}X and X'W^{-1}y
  XtWinv <- t(X) %*% W_inv
  XtWinvX <- XtWinv %*% X
  XtWinvy <- XtWinv %*% y

  # Invert X'W^{-1}X
  XtWinvX_chol <- tryCatch(
    chol(XtWinvX),
    error = function(e) NULL
  )

  if (is.null(XtWinvX_chol)) {
    if (ret_mu) return(list(logL = -1e100, mu = NULL))
    return(-1e100)
  }

  XtWinvX_inv <- chol2inv(XtWinvX_chol)
  log_det_XtWinvX <- 2 * sum(log(diag(XtWinvX_chol)))

  # REML mean estimate
  mu <- XtWinvX_inv %*% XtWinvy

  # Residuals
  resid <- y - X %*% mu

  # REML log-likelihood
  # n_obs - m because REML adjusts for estimating m fixed effects
  logL <- -0.5 * (
    (n_obs - m) * log(2 * pi) +
    log_det_W +
    log_det_XtWinvX +
    as.numeric(t(resid) %*% W_inv %*% resid)
  )

  if (!is.finite(logL)) {
    if (ret_mu) return(list(logL = -1e100, mu = NULL))
    return(-1e100)
  }

  if (ret_mu) {
    return(list(logL = logL, mu = mu))
  }

  logL
}


#' Fit Model Using Full Matrix Approach
#'
#' Wrapper around loglik_full_matrix for use with optim().
#' Transforms theta to Sigma and B, then computes likelihood.
#'
#' @param theta Parameter vector
#' @param tree A phylo object
#' @param trait_data Data frame with species and traits
#' @param m Number of traits
#' @param Sigma_diag Constrain Sigma to diagonal
#' @param B_diag Constrain B to diagonal
#' @param Sigma_fixed Fixed value for Sigma
#' @param B_fixed Fixed value for B
#' @param model Evolutionary model
#' @param evo_model_par_fixed Fixed evolutionary parameter
#'
#' @return Log-likelihood value
#'
#' @export
restricted_loglik_full <- function(theta, tree, trait_data, m = NULL,
                                    Sigma_diag = FALSE, B_diag = FALSE,
                                    Sigma_fixed = NULL, B_fixed = NULL,
                                    model = "BM", evo_model_par_fixed = NULL) {

  # Infer m from data if not provided
  if (is.null(m)) {
    m <- ncol(trait_data) - 1
  }

  # Convert theta to matrices
  pars <- theta2vals(
    theta = theta, m = m,
    Sigma_diag = Sigma_diag, B_diag = B_diag,
    Sigma_fixed = Sigma_fixed, B_fixed = B_fixed,
    model = model, evo_model_par_fixed = evo_model_par_fixed
  )

  Sigma <- pars$Sigma
  B <- pars$B

  # Handle evolutionary model transforms
  if (model != "BM") {
    tree <- dist_from_root(tree)
    evo_par <- if (!is.null(evo_model_par_fixed)) evo_model_par_fixed else pars$evo_model_par
    tree <- transf_branch_lengths(tree, model, evo_par)
  }

  # Compute likelihood
  loglik_full_matrix(Sigma, B, tree, trait_data, ret_mu = FALSE)
}
