# =============================================================================
# Parameter Constraint System
# =============================================================================

#' Create a Constraint Specification
#'
#' Builds an internal representation of parameter constraints for the R
#' (phylogenetic rate) and S (within-species covariance) matrices.
#'
#' @param p Integer, dimension of the covariance matrices (number of traits)
#' @param phylo_correlated Logical, if FALSE constrains R to be diagonal (default TRUE)
#' @param pheno_correlated Logical, if FALSE constrains S to be diagonal (default TRUE)
#' @param R_fixed List of lists, each with elements i, j, value. Fixes R[i,i] to
#'   a known value by constraining the corresponding Cholesky elements. Only diagonal
#'   entries (i == j) are supported. When R[i,i] is fixed to value v, L[i,i] = sqrt(v)
#'   and L[i,j] = 0 for j < i, which also forces R[i,j] = 0 for j < i.
#' @param S_fixed List of lists, each with elements i, j, value. Same as R_fixed but
#'   for the S matrix.
#' @param R_shared_var List of integer vectors. Each vector specifies a group of
#'   diagonal indices that share the same L diagonal element (exp(theta)).
#'   E.g., list(c(1,2,3)) means L[1,1] = L[2,2] = L[3,3] = exp(theta_shared).
#'   Note: this shares the Cholesky diagonal, so R[i,i] values may differ slightly
#'   due to off-diagonal L elements (R[i,i] = sum(L[i,:]^2)).
#' @param S_shared_var List of integer vectors, same as R_shared_var but for S.
#' @param R_blocks List of integer vectors specifying block-diagonal structure.
#'   Each vector lists the trait indices in one block. Traits in different blocks
#'   have zero covariance. E.g., list(1:3, 4:5) creates two independent blocks.
#' @param S_blocks List of integer vectors, same as R_blocks but for S.
#'
#' @return A list with class "ConstraintSpec" containing:
#'   \item{p}{Dimension of matrices}
#'   \item{R}{List with constraint details for R matrix}
#'   \item{S}{List with constraint details for S matrix}
#'   \item{n_free_R}{Number of free parameters for R}
#'   \item{n_free_S}{Number of free parameters for S}
#'   \item{n_free}{Total number of free matrix parameters}
#'
#' @export
make_constraints <- function(p,
                              phylo_correlated = TRUE,
                              pheno_correlated = TRUE,
                              R_fixed = NULL,
                              S_fixed = NULL,
                              R_shared_var = NULL,
                              S_shared_var = NULL,
                              R_blocks = NULL,
                              S_blocks = NULL) {

  R_spec <- build_matrix_constraint(
    p, is_correlated = phylo_correlated,
    fixed = R_fixed, shared_var = R_shared_var, blocks = R_blocks,
    matrix_name = "R"
  )

  S_spec <- build_matrix_constraint(
    p, is_correlated = pheno_correlated,
    fixed = S_fixed, shared_var = S_shared_var, blocks = S_blocks,
    matrix_name = "S"
  )

  result <- list(
    p = p,
    R = R_spec,
    S = S_spec,
    n_free_R = R_spec$n_free,
    n_free_S = S_spec$n_free,
    n_free = R_spec$n_free + S_spec$n_free
  )
  class(result) <- "ConstraintSpec"
  result
}


#' Build constraint specification for one matrix
#'
#' @param p Dimension
#' @param is_correlated If FALSE, diagonal-only (no off-diagonal)
#' @param fixed List of fixed entry specs
#' @param shared_var List of shared variance groups
#' @param blocks List of block-diagonal groups
#' @param matrix_name "R" or "S" (for error messages)
#'
#' @return List describing the constraint
#' @keywords internal
build_matrix_constraint <- function(p, is_correlated = TRUE,
                                     fixed = NULL, shared_var = NULL,
                                     blocks = NULL, matrix_name = "R") {

  has_blocks <- !is.null(blocks) && length(blocks) > 0
  has_fixed <- !is.null(fixed) && length(fixed) > 0
  has_shared <- !is.null(shared_var) && length(shared_var) > 0

  # Validate blocks
  if (has_blocks) {
    all_indices <- sort(unlist(blocks))
    if (length(all_indices) != p || !all(all_indices == 1:p)) {
      stop(matrix_name, "_blocks must partition all traits 1:", p,
           ". Got: ", paste(all_indices, collapse = ", "))
    }
    if (any(duplicated(all_indices))) {
      stop(matrix_name, "_blocks contains duplicate trait indices")
    }
  }

  # Validate fixed entries
  if (has_fixed) {
    for (f in fixed) {
      if (is.null(f$i) || is.null(f$j) || is.null(f$value)) {
        stop(matrix_name, "_fixed entries must have i, j, and value")
      }
      if (f$i != f$j) {
        stop(matrix_name, "_fixed: only diagonal entries (i == j) are supported. ",
             "Got i=", f$i, ", j=", f$j)
      }
      if (f$i < 1 || f$i > p) {
        stop(matrix_name, "_fixed: index ", f$i, " out of range [1, ", p, "]")
      }
      if (f$value <= 0) {
        stop(matrix_name, "_fixed: variance must be positive, got ", f$value)
      }
    }
  }

  # Validate shared variance groups
  if (has_shared) {
    all_shared <- unlist(shared_var)
    if (any(duplicated(all_shared))) {
      stop(matrix_name, "_shared_var: trait index appears in multiple groups")
    }
    if (any(all_shared < 1 | all_shared > p)) {
      stop(matrix_name, "_shared_var: indices must be in [1, ", p, "]")
    }
  }

  # Cannot fix and share the same diagonal entry
  if (has_fixed && has_shared) {
    fixed_indices <- sapply(fixed, function(f) f$i)
    shared_indices <- unlist(shared_var)
    overlap <- intersect(fixed_indices, shared_indices)
    if (length(overlap) > 0) {
      stop(matrix_name, ": cannot both fix and share variance for trait(s) ",
           paste(overlap, collapse = ", "))
    }
  }

  if (has_blocks) {
    spec <- build_block_constraint(p, blocks, is_correlated, fixed, shared_var, matrix_name)
  } else if (!is_correlated) {
    spec <- build_diagonal_constraint(p, fixed, shared_var, matrix_name)
  } else {
    spec <- build_full_constraint(p, fixed, shared_var, matrix_name)
  }

  spec
}


#' Build block-diagonal constraint
#' @keywords internal
build_block_constraint <- function(p, blocks, is_correlated, fixed, shared_var, matrix_name) {
  block_specs <- list()
  n_free <- 0
  theta_offset <- 0

  for (b_idx in seq_along(blocks)) {
    block_indices <- blocks[[b_idx]]
    bp <- length(block_indices)

    # Find fixed entries within this block
    block_fixed <- NULL
    if (!is.null(fixed)) {
      for (f in fixed) {
        if (f$i %in% block_indices) {
          local_i <- match(f$i, block_indices)
          block_fixed <- c(block_fixed, list(list(i = local_i, j = local_i, value = f$value)))
        }
      }
    }

    # Find shared var groups within this block
    block_shared <- NULL
    if (!is.null(shared_var)) {
      for (sg in shared_var) {
        in_block <- sg[sg %in% block_indices]
        if (length(in_block) > 1) {
          local_indices <- match(in_block, block_indices)
          block_shared <- c(block_shared, list(local_indices))
        }
      }
    }

    if (is_correlated) {
      sub_spec <- build_full_constraint(bp, block_fixed, block_shared,
                                         paste0(matrix_name, "_block", b_idx))
    } else {
      sub_spec <- build_diagonal_constraint(bp, block_fixed, block_shared,
                                              paste0(matrix_name, "_block", b_idx))
    }
    sub_spec$global_indices <- block_indices
    sub_spec$theta_offset <- theta_offset
    theta_offset <- theta_offset + sub_spec$n_free
    n_free <- n_free + sub_spec$n_free

    block_specs[[b_idx]] <- sub_spec
  }

  list(
    type = "block",
    p = p,
    blocks = blocks,
    block_specs = block_specs,
    n_free = n_free,
    is_correlated = is_correlated,
    is_diagonal = !is_correlated
  )
}


#' Build diagonal-only constraint
#' @keywords internal
build_diagonal_constraint <- function(p, fixed, shared_var, matrix_name) {
  fixed_indices <- if (!is.null(fixed)) sapply(fixed, function(f) f$i) else integer(0)
  fixed_values <- if (!is.null(fixed)) sapply(fixed, function(f) f$value) else numeric(0)
  if (length(fixed_indices) > 0) names(fixed_values) <- fixed_indices

  theta_map <- integer(p)
  fixed_map <- logical(p)
  fixed_val_map <- numeric(p)

  for (i in seq_len(p)) {
    if (as.character(i) %in% names(fixed_values)) {
      fixed_map[i] <- TRUE
      fixed_val_map[i] <- fixed_values[as.character(i)]
    }
  }

  n_free <- 0
  if (!is.null(shared_var)) {
    for (sg in shared_var) {
      n_free <- n_free + 1
      for (idx in sg) {
        if (fixed_map[idx]) next
        theta_map[idx] <- n_free
      }
    }
  }

  for (i in seq_len(p)) {
    if (!fixed_map[i] && theta_map[i] == 0) {
      n_free <- n_free + 1
      theta_map[i] <- n_free
    }
  }

  list(
    type = "diagonal",
    p = p,
    n_free = n_free,
    is_diagonal = TRUE,
    is_correlated = FALSE,
    theta_map = theta_map,
    fixed_map = fixed_map,
    fixed_val_map = fixed_val_map,
    shared_var = shared_var
  )
}


#' Build full (correlated) constraint
#' @keywords internal
build_full_constraint <- function(p, fixed, shared_var, matrix_name) {
  fixed_diag <- if (!is.null(fixed)) sapply(fixed, function(f) f$i) else integer(0)
  fixed_vals <- if (!is.null(fixed)) sapply(fixed, function(f) f$value) else numeric(0)
  if (length(fixed_diag) > 0) names(fixed_vals) <- fixed_diag

  n_chol <- p * (p + 1) / 2
  chol_positions <- vector("list", n_chol)
  idx <- 1
  for (i in 1:p) {
    for (j in 1:i) {
      chol_positions[[idx]] <- c(i, j)
      idx <- idx + 1
    }
  }

  theta_idx <- integer(n_chol)
  is_fixed_pos <- logical(n_chol)
  fixed_chol_values <- numeric(n_chol)
  is_diagonal_pos <- logical(n_chol)

  for (k in seq_len(n_chol)) {
    pos <- chol_positions[[k]]
    i <- pos[1]; j <- pos[2]
    is_diagonal_pos[k] <- (i == j)

    if (as.character(i) %in% names(fixed_vals)) {
      if (i == j) {
        is_fixed_pos[k] <- TRUE
        fixed_chol_values[k] <- sqrt(fixed_vals[as.character(i)])
      } else {
        is_fixed_pos[k] <- TRUE
        fixed_chol_values[k] <- 0
      }
    }
  }

  n_free <- 0

  shared_theta_map <- list()
  if (!is.null(shared_var)) {
    for (sg_idx in seq_along(shared_var)) {
      sg <- shared_var[[sg_idx]]
      n_free <- n_free + 1
      shared_theta_idx <- n_free
      shared_theta_map[[sg_idx]] <- shared_theta_idx

      for (trait_i in sg) {
        for (k in seq_len(n_chol)) {
          pos <- chol_positions[[k]]
          if (pos[1] == trait_i && pos[2] == trait_i && !is_fixed_pos[k]) {
            theta_idx[k] <- shared_theta_idx
          }
        }
      }
    }
  }

  for (k in seq_len(n_chol)) {
    if (!is_fixed_pos[k] && theta_idx[k] == 0) {
      n_free <- n_free + 1
      theta_idx[k] <- n_free
    }
  }

  list(
    type = "full",
    p = p,
    n_free = n_free,
    is_diagonal = FALSE,
    is_correlated = TRUE,
    n_chol = n_chol,
    chol_positions = chol_positions,
    theta_idx = theta_idx,
    is_fixed_pos = is_fixed_pos,
    fixed_chol_values = fixed_chol_values,
    is_diagonal_pos = is_diagonal_pos,
    shared_var = shared_var,
    fixed = fixed
  )
}


#' Convert constrained theta to a covariance matrix
#'
#' @param theta Numeric vector of free parameters for this matrix
#' @param spec Constraint spec for one matrix (from build_*_constraint)
#'
#' @return p x p positive-definite covariance matrix
#' @keywords internal
constrained_theta_to_matrix <- function(theta, spec) {
  if (spec$type == "block") {
    return(constrained_theta_to_matrix_block(theta, spec))
  }

  p <- spec$p

  if (spec$type == "diagonal") {
    mat <- matrix(0, p, p)
    for (i in 1:p) {
      if (spec$fixed_map[i]) {
        mat[i, i] <- spec$fixed_val_map[i]
      } else {
        mat[i, i] <- exp(theta[spec$theta_map[i]])
      }
    }
    return(mat)
  }

  # type == "full": build L from theta, then compute L %*% t(L)
  L <- matrix(0, p, p)
  for (k in seq_len(spec$n_chol)) {
    pos <- spec$chol_positions[[k]]
    i <- pos[1]; j <- pos[2]

    if (spec$is_fixed_pos[k]) {
      L[i, j] <- spec$fixed_chol_values[k]
    } else {
      tidx <- spec$theta_idx[k]
      if (spec$is_diagonal_pos[k]) {
        L[i, j] <- exp(theta[tidx])
      } else {
        L[i, j] <- theta[tidx]
      }
    }
  }

  L %*% t(L)
}


#' Convert constrained theta to block-diagonal matrix
#' @keywords internal
constrained_theta_to_matrix_block <- function(theta, spec) {
  p <- spec$p
  mat <- matrix(0, p, p)

  for (b_idx in seq_along(spec$block_specs)) {
    bspec <- spec$block_specs[[b_idx]]
    gi <- bspec$global_indices
    if (bspec$n_free > 0) {
      block_theta <- theta[(bspec$theta_offset + 1):(bspec$theta_offset + bspec$n_free)]
    } else {
      block_theta <- numeric(0)
    }
    mat[gi, gi] <- constrained_theta_to_matrix(block_theta, bspec)
  }

  mat
}


#' Convert covariance matrix to constrained theta
#'
#' @param mat p x p covariance matrix
#' @param spec Constraint spec for one matrix
#'
#' @return Numeric vector of free parameters
#' @keywords internal
constrained_matrix_to_theta <- function(mat, spec) {
  if (spec$type == "block") {
    return(constrained_matrix_to_theta_block(mat, spec))
  }

  p <- spec$p

  if (spec$type == "diagonal") {
    theta <- numeric(spec$n_free)
    assigned <- logical(spec$n_free)
    for (i in 1:p) {
      if (!spec$fixed_map[i]) {
        tidx <- spec$theta_map[i]
        if (!assigned[tidx]) {
          theta[tidx] <- log(mat[i, i])
          assigned[tidx] <- TRUE
        }
      }
    }
    return(theta)
  }

  # type == "full"
  L <- tryCatch(
    t(chol(mat)),
    error = function(e) t(chol(mat + 1e-8 * diag(p)))
  )

  theta <- numeric(spec$n_free)
  assigned <- logical(spec$n_free)

  for (k in seq_len(spec$n_chol)) {
    if (spec$is_fixed_pos[k]) next

    pos <- spec$chol_positions[[k]]
    i <- pos[1]; j <- pos[2]
    tidx <- spec$theta_idx[k]

    if (!assigned[tidx]) {
      if (spec$is_diagonal_pos[k]) {
        theta[tidx] <- log(L[i, j])
      } else {
        theta[tidx] <- L[i, j]
      }
      assigned[tidx] <- TRUE
    }
  }

  theta
}


#' Convert block-diagonal matrix to constrained theta
#' @keywords internal
constrained_matrix_to_theta_block <- function(mat, spec) {
  theta <- numeric(spec$n_free)

  for (b_idx in seq_along(spec$block_specs)) {
    bspec <- spec$block_specs[[b_idx]]
    gi <- bspec$global_indices
    block_theta <- constrained_matrix_to_theta(mat[gi, gi, drop = FALSE], bspec)
    if (length(block_theta) > 0) {
      theta[(bspec$theta_offset + 1):(bspec$theta_offset + bspec$n_free)] <- block_theta
    }
  }

  theta
}


#' Convert constrained gradient matrices to theta-space gradient
#'
#' @param grad_R Gradient w.r.t. R matrix (p x p)
#' @param grad_S Gradient w.r.t. S matrix (p x p)
#' @param R Current R matrix
#' @param S Current S matrix
#' @param constraints ConstraintSpec object
#'
#' @return Numeric vector of gradient in constrained theta space
#' @keywords internal
constrained_gradient_to_theta <- function(grad_R, grad_S, R, S, constraints) {
  grad_theta_R <- constrained_grad_one_matrix(grad_R, R, constraints$R)
  grad_theta_S <- constrained_grad_one_matrix(grad_S, S, constraints$S)
  c(grad_theta_R, grad_theta_S)
}


#' Compute constrained gradient for one matrix
#' @keywords internal
constrained_grad_one_matrix <- function(grad_mat, mat, spec) {
  if (spec$type == "block") {
    return(constrained_grad_block(grad_mat, mat, spec))
  }

  p <- spec$p

  if (spec$type == "diagonal") {
    theta_grad <- numeric(spec$n_free)
    for (i in 1:p) {
      if (!spec$fixed_map[i]) {
        tidx <- spec$theta_map[i]
        theta_grad[tidx] <- theta_grad[tidx] + grad_mat[i, i] * mat[i, i]
      }
    }
    return(theta_grad)
  }

  # type == "full"
  L <- tryCatch(t(chol(mat)), error = function(e) t(chol(mat + 1e-8 * diag(p))))
  grad_L <- 2 * grad_mat %*% L

  theta_grad <- numeric(spec$n_free)
  for (k in seq_len(spec$n_chol)) {
    if (spec$is_fixed_pos[k]) next

    pos <- spec$chol_positions[[k]]
    i <- pos[1]; j <- pos[2]
    tidx <- spec$theta_idx[k]

    if (spec$is_diagonal_pos[k]) {
      theta_grad[tidx] <- theta_grad[tidx] + L[i, j] * grad_L[i, j]
    } else {
      theta_grad[tidx] <- theta_grad[tidx] + grad_L[i, j]
    }
  }

  theta_grad
}


#' Compute constrained gradient for block-diagonal matrix
#' @keywords internal
constrained_grad_block <- function(grad_mat, mat, spec) {
  theta_grad <- numeric(spec$n_free)

  for (b_idx in seq_along(spec$block_specs)) {
    bspec <- spec$block_specs[[b_idx]]
    gi <- bspec$global_indices
    block_theta_grad <- constrained_grad_one_matrix(
      grad_mat[gi, gi, drop = FALSE],
      mat[gi, gi, drop = FALSE],
      bspec
    )
    if (length(block_theta_grad) > 0) {
      theta_grad[(bspec$theta_offset + 1):(bspec$theta_offset + bspec$n_free)] <- block_theta_grad
    }
  }

  theta_grad
}


#' Convert constrained theta (for both R and S) to matrices
#'
#' Convenience function: given a full constrained theta vector containing
#' both R and S free parameters, split and reconstruct both matrices.
#'
#' @param theta Numeric vector of free parameters (R params then S params)
#' @param constraints ConstraintSpec object
#'
#' @return List with components R and S
#' @keywords internal
constrained_theta_to_matrices <- function(theta, constraints) {
  n_R <- constraints$n_free_R
  n_S <- constraints$n_free_S

  theta_R <- if (n_R > 0) theta[1:n_R] else numeric(0)
  theta_S <- if (n_S > 0) theta[(n_R + 1):(n_R + n_S)] else numeric(0)

  R <- constrained_theta_to_matrix(theta_R, constraints$R)
  S <- constrained_theta_to_matrix(theta_S, constraints$S)

  list(R = R, S = S)
}


#' Convert R and S matrices to constrained theta vector
#'
#' Convenience function: given R and S matrices, produce the combined
#' constrained theta vector.
#'
#' @param R Phylogenetic rate matrix
#' @param S Residual covariance matrix
#' @param constraints ConstraintSpec object
#'
#' @return Numeric vector of free parameters
#' @keywords internal
constrained_matrices_to_theta <- function(R, S, constraints) {
  theta_R <- constrained_matrix_to_theta(R, constraints$R)
  theta_S <- constrained_matrix_to_theta(S, constraints$S)
  c(theta_R, theta_S)
}


# =============================================================================
# Original parameterization functions (unchanged for backwards compatibility)
# =============================================================================

#' Convert parameter vector theta to covariance matrices
#'
#' Converts an unconstrained parameter vector to positive-definite covariance
#' matrices using Cholesky decomposition parameterization.
#'
#' @param theta Numeric vector of parameters
#' @param m Integer, dimension of the covariance matrices
#' @param Sigma_diag Logical, if TRUE Sigma is constrained to be diagonal
#' @param B_diag Logical, if TRUE B is constrained to be diagonal
#' @param Sigma_fixed Matrix or NULL, fixed value for Sigma (not optimized)
#' @param B_fixed Matrix or NULL, fixed value for B (not optimized)
#' @param model Character, evolutionary model ("BM" or "lambda")
#' @param evo_model_par_fixed Numeric or NULL, fixed evolutionary parameter
#'
#' @return List with components: Sigma, B, evo_model_par
#' @export
theta2vals <- function(theta, m = NULL, Sigma_diag = FALSE, B_diag = FALSE,
                       Sigma_fixed = NULL, B_fixed = NULL,
                       model = "BM", evo_model_par_fixed = NULL) {

  # Infer m if not provided
  if (is.null(m)) {
    theta_len <- length(theta)
    if (model != "BM" && is.null(evo_model_par_fixed)) {
      theta_len <- theta_len - 1
    }
    if (!is.null(Sigma_fixed)) {
      m <- nrow(Sigma_fixed)
    } else if (!is.null(B_fixed)) {
      m <- nrow(B_fixed)
    } else {
      if (Sigma_diag && B_diag) {
        m <- theta_len %/% 2
      } else if (Sigma_diag || B_diag) {
        m <- as.integer((-3 + sqrt(9 + 8 * theta_len)) / 2)
      } else {
        m <- as.integer((-1 + sqrt(1 + 4 * theta_len)) / 2)
      }
    }
  }

  next_idx <- 1
  evo_model_par <- NULL

  # Extract Sigma
  if (!is.null(Sigma_fixed)) {
    Sigma <- Sigma_fixed
  } else {
    if (Sigma_diag) {
      Sigma_params <- theta[1:m]
      Sigma <- diag(exp(Sigma_params), nrow = m)
      next_idx <- m + 1
    } else {
      L_sigma <- matrix(0, m, m)
      idx <- 1
      for (i in 1:m) {
        for (j in 1:i) {
          if (i == j) {
            # Diagonal elements are exponentiated
            L_sigma[i, j] <- exp(theta[idx])
          } else {
            # Off-diagonal elements used directly
            L_sigma[i, j] <- theta[idx]
          }
          idx <- idx + 1
        }
      }
      Sigma <- L_sigma %*% t(L_sigma)
      next_idx <- idx
    }
  }

  # Extract B
  if (!is.null(B_fixed)) {
    B <- B_fixed
  } else {
    if (B_diag) {
      B_params <- theta[next_idx:(next_idx + m - 1)]
      B <- diag(exp(B_params), nrow = m)
      next_idx <- next_idx + m
    } else {
      L_b <- matrix(0, m, m)
      idx <- next_idx
      for (i in 1:m) {
        for (j in 1:i) {
          if (i == j) {
            L_b[i, j] <- exp(theta[idx])
          } else {
            L_b[i, j] <- theta[idx]
          }
          idx <- idx + 1
        }
      }
      B <- L_b %*% t(L_b)
      next_idx <- idx
    }
  }

  # Extract evolutionary model parameter
  if (model != "BM") {
    if (!is.null(evo_model_par_fixed)) {
      evo_model_par <- evo_model_par_fixed
    } else {
      raw <- theta[next_idx]
      evo_model_par <- switch(model,
        "lambda" = plogis(raw),              # [0, 1] via inverse logit
        "kappa"  = exp(raw),                 # (0, Inf) via exp
        "delta"  = exp(raw),                 # (0, Inf) via exp
        "EB"     = -exp(raw),                # (-Inf, 0] via -exp
        "OU"     = exp(raw),                 # (0, Inf) via exp
        stop("Unknown model for theta extraction: ", model)
      )
    }
  }

  list(Sigma = Sigma, B = B, evo_model_par = evo_model_par)
}


#' Convert covariance matrices to parameter vector
#'
#' Converts positive-definite covariance matrices to an unconstrained parameter
#' vector using Cholesky decomposition.
#'
#' @param Sigma Phylogenetic covariance matrix
#' @param B Residual/within-species covariance matrix
#' @param evo_model_par Evolutionary model parameter (e.g., lambda)
#' @param ntheta Length of parameter vector (computed if NULL)
#' @param Sigma_diag Logical, if TRUE Sigma is diagonal
#' @param B_diag Logical, if TRUE B is diagonal
#' @param Sigma_fixed Matrix or NULL, fixed value for Sigma
#' @param B_fixed Matrix or NULL, fixed value for B
#' @param model Character, evolutionary model
#' @param evo_model_par_fixed Numeric or NULL, fixed evolutionary parameter
#'
#' @return Numeric vector of parameters
#' @export
vals2theta <- function(Sigma = NULL, B = NULL, evo_model_par = NULL,
                       ntheta = NULL, Sigma_diag = FALSE, B_diag = FALSE,
                       Sigma_fixed = NULL, B_fixed = NULL,
                       model = "BM", evo_model_par_fixed = NULL) {

  # Determine m
  if (!is.null(Sigma_fixed)) {
    m <- nrow(Sigma_fixed)
  } else if (!is.null(B_fixed)) {
    m <- nrow(B_fixed)
  } else if (!is.null(Sigma)) {
    m <- nrow(Sigma)
  } else if (!is.null(B)) {
    m <- nrow(B)
  } else {
    stop("Must provide at least one of Sigma, B, Sigma_fixed, or B_fixed")
  }

  # Compute ntheta if not provided
  if (is.null(ntheta)) {
    ntheta <- 0
    optimize_Sigma <- is.null(Sigma_fixed)
    optimize_B <- is.null(B_fixed)
    optimize_evo <- model != "BM" && is.null(evo_model_par_fixed)

    if (optimize_Sigma) {
      ntheta <- ntheta + if (Sigma_diag) m else m * (m + 1) / 2
    }
    if (optimize_B) {
      ntheta <- ntheta + if (B_diag) m else m * (m + 1) / 2
    }
    if (optimize_evo) {
      ntheta <- ntheta + 1
    }
  }

  if (ntheta == 0) {
    return(numeric(0))
  }

  theta <- numeric(ntheta)
  optimize_Sigma <- is.null(Sigma_fixed)
  optimize_B <- is.null(B_fixed)
  optimize_evo <- model != "BM" && is.null(evo_model_par_fixed)

  next_idx <- 1

  # Convert Sigma to parameters
  if (optimize_Sigma && !is.null(Sigma)) {
    if (Sigma_diag) {
      theta[1:m] <- log(diag(Sigma))
      next_idx <- m + 1
    } else {
      # Cholesky decomposition
      L_sigma <- tryCatch(
        chol(Sigma),
        error = function(e) chol(Sigma + 1e-8 * diag(m))
      )
      L_sigma <- t(L_sigma)  # Lower triangular

      idx <- next_idx
      for (i in 1:m) {
        for (j in 1:i) {
          if (i == j) {
            theta[idx] <- log(L_sigma[i, j])
          } else {
            theta[idx] <- L_sigma[i, j]
          }
          idx <- idx + 1
        }
      }
      next_idx <- idx
    }
  }

  # Convert B to parameters
  if (optimize_B && !is.null(B)) {
    if (B_diag) {
      theta[next_idx:(next_idx + m - 1)] <- log(diag(B))
      next_idx <- next_idx + m
    } else {
      L_b <- tryCatch(
        chol(B),
        error = function(e) chol(B + 1e-8 * diag(m))
      )
      L_b <- t(L_b)  # Lower triangular

      idx <- next_idx
      for (i in 1:m) {
        for (j in 1:i) {
          if (i == j) {
            theta[idx] <- log(L_b[i, j])
          } else {
            theta[idx] <- L_b[i, j]
          }
          idx <- idx + 1
        }
      }
      next_idx <- idx
    }
  }

  # Convert evolutionary parameter
  if (optimize_evo && !is.null(evo_model_par)) {
    theta[next_idx] <- switch(model,
      "lambda" = qlogis(evo_model_par),       # [0, 1] -> R via logit
      "kappa"  = log(evo_model_par),           # (0, Inf) -> R via log
      "delta"  = log(evo_model_par),           # (0, Inf) -> R via log
      "EB"     = log(-evo_model_par),           # (-Inf, 0] -> R via log(-x)
      "OU"     = log(evo_model_par),           # (0, Inf) -> R via log
      stop("Unknown model for theta encoding: ", model)
    )
  }

  theta
}


#' Generate default starting values for optimization
#'
#' @param m Dimension of covariance matrices
#' @param Sigma_diag Logical, constrain Sigma to diagonal
#' @param B_diag Logical, constrain B to diagonal
#' @param Sigma_fixed Fixed Sigma matrix
#' @param B_fixed Fixed B matrix
#' @param Sigma_start Starting value for Sigma
#' @param B_start Starting value for B
#' @param model Evolutionary model
#' @param evo_model_par_start Starting value for evolutionary parameter
#' @param evo_model_par_fixed Fixed evolutionary parameter
#'
#' @return Numeric vector of starting parameter values
#' @export
get_default_starting_theta <- function(m, Sigma_diag = FALSE, B_diag = FALSE,
                                        Sigma_fixed = NULL, B_fixed = NULL,
                                        Sigma_start = NULL, B_start = NULL,
                                        model = "BM", evo_model_par_start = NULL,
                                        evo_model_par_fixed = NULL) {

  optimize_Sigma <- is.null(Sigma_fixed)
  optimize_B <- is.null(B_fixed)
  optimize_evo <- model != "BM" && is.null(evo_model_par_fixed)

  # Set defaults
  if (optimize_Sigma && is.null(Sigma_start)) {
    Sigma_start <- diag(m)
  }
  if (optimize_B && is.null(B_start)) {
    B_start <- diag(m)
  }
  if (optimize_evo && is.null(evo_model_par_start)) {
    evo_model_par_start <- switch(model,
      "lambda" = 0.5,
      "kappa"  = 1.0,
      "delta"  = 1.0,
      "EB"     = -0.1,
      "OU"     = 0.1,
      stop("Unknown model for default starting value: ", model)
    )
  }

  # Compute ntheta
  ntheta <- 0
  if (optimize_Sigma) {
    ntheta <- ntheta + if (Sigma_diag) m else m * (m + 1) / 2
  }
  if (optimize_B) {
    ntheta <- ntheta + if (B_diag) m else m * (m + 1) / 2
  }
  if (optimize_evo) {
    ntheta <- ntheta + 1
  }

  vals2theta(
    Sigma = Sigma_start, B = B_start, evo_model_par = evo_model_par_start,
    ntheta = ntheta, Sigma_diag = Sigma_diag, B_diag = B_diag,
    Sigma_fixed = Sigma_fixed, B_fixed = B_fixed,
    model = model, evo_model_par_fixed = evo_model_par_fixed
  )
}


#' Convert R and S matrices to theta vector (convenience wrapper)
#'
#' Simple interface for converting two covariance matrices (R and S) to a theta
#' parameter vector. When constraints are provided, uses the constrained
#' parameterization; otherwise uses the standard log-Cholesky parameterization.
#'
#' @param R Phylogenetic rate matrix
#' @param S Residual covariance matrix
#' @param constraints Optional ConstraintSpec object (from make_constraints)
#'
#' @return Numeric vector of parameters
#'
#' @keywords internal
matrices_to_theta <- function(R, S, constraints = NULL) {
  if (!is.null(constraints)) {
    return(constrained_matrices_to_theta(R, S, constraints))
  }
  vals2theta(Sigma = R, B = S)
}


#' Convert theta vector back to R and S matrices (convenience wrapper)
#'
#' When constraints are provided, uses the constrained parameterization;
#' otherwise uses the standard log-Cholesky parameterization.
#'
#' @param theta Numeric vector of parameters
#' @param p Dimension of the covariance matrices
#' @param constraints Optional ConstraintSpec object (from make_constraints)
#'
#' @return List with components R and S
#'
#' @keywords internal
theta_to_matrices <- function(theta, p, constraints = NULL) {
  if (!is.null(constraints)) {
    return(constrained_theta_to_matrices(theta, constraints))
  }
  result <- theta2vals(theta, m = p)
  list(R = result$Sigma, S = result$B)
}


#' Convert gradient matrices to theta-space gradient
#'
#' Uses chain rule: d/d(theta) = d/d(L) * d(L)/d(theta)
#' where L is the lower-triangular Cholesky factor with log-diagonal.
#'
#' Optionally includes the evolutionary parameter gradient, chaining through
#' the parameter-to-theta transform (logit for lambda, log for others).
#'
#' When constraints are provided, uses the constrained gradient mapping instead.
#'
#' @param grad_R Gradient w.r.t. R matrix
#' @param grad_S Gradient w.r.t. S matrix
#' @param R Current R matrix
#' @param S Current S matrix
#' @param grad_evo Gradient w.r.t. evolutionary parameter (scalar, or NULL)
#' @param evo_model_par Current evolutionary parameter value (needed for chain rule)
#' @param model Evolutionary model name (needed for chain rule)
#' @param constraints Optional ConstraintSpec object (from make_constraints)
#'
#' @return Numeric vector of gradient in theta space
#'
#' @keywords internal
gradient_to_theta <- function(grad_R, grad_S, R, S,
                              grad_evo = NULL, evo_model_par = NULL,
                              model = "BM", constraints = NULL) {

  # Use constrained gradient mapping if constraints are provided
  if (!is.null(constraints)) {
    theta_grad <- constrained_gradient_to_theta(grad_R, grad_S, R, S, constraints)
  } else {
    p <- nrow(R)

    # Cholesky factors
    L_R <- tryCatch(t(chol(R)), error = function(e) t(chol(R + 1e-8 * diag(p))))
    L_S <- tryCatch(t(chol(S)), error = function(e) t(chol(S + 1e-8 * diag(p))))

    # Gradient w.r.t. L: grad_L = 2 * grad_Sigma * L (for lower triangular)
    grad_L_R <- 2 * grad_R %*% L_R
    grad_L_S <- 2 * grad_S %*% L_S

    # Convert to theta-space
    theta_grad_R <- numeric(p * (p + 1) / 2)
    theta_grad_S <- numeric(p * (p + 1) / 2)

    idx <- 1
    for (i in 1:p) {
      for (j in 1:i) {
        if (i == j) {
          # d/d(log L_ii) = L_ii * d/d(L_ii)
          theta_grad_R[idx] <- L_R[i, j] * grad_L_R[i, j]
          theta_grad_S[idx] <- L_S[i, j] * grad_L_S[i, j]
        } else {
          theta_grad_R[idx] <- grad_L_R[i, j]
          theta_grad_S[idx] <- grad_L_S[i, j]
        }
        idx <- idx + 1
      }
    }

    theta_grad <- c(theta_grad_R, theta_grad_S)
  }

  # Append evolutionary parameter gradient if present
  if (!is.null(grad_evo) && !is.null(evo_model_par) && model != "BM") {
    dparam_dtheta <- switch(model,
      "lambda" = evo_model_par * (1 - evo_model_par),  # logistic derivative
      "kappa"  = evo_model_par,                         # exp derivative
      "delta"  = evo_model_par,                         # exp derivative
      "EB"     = evo_model_par,                         # r = -exp(theta), dr/dtheta = -exp(theta) = r
      "OU"     = evo_model_par,                         # exp derivative
      stop("Unknown model for evo gradient chain rule: ", model)
    )
    theta_grad <- c(theta_grad, grad_evo * dparam_dtheta)
  }

  theta_grad
}


#' Calculate distance from root to each node
#'
#' @param tree A phylo object
#'
#' @return The tree with added 'root_dist' element
#' @export
dist_from_root <- function(tree) {
  nspecies <- length(tree$tip.label)
  Nnode <- tree$Nnode
  total_nodes <- nspecies + Nnode

  root_dist <- numeric(total_nodes)

  # Process edges in reverse order (root to tips)
  for (i in rev(seq_len(nrow(tree$edge)))) {
    anc <- tree$edge[i, 1]
    des <- tree$edge[i, 2]
    root_dist[des] <- root_dist[anc] + tree$edge.length[i]
  }

  tree$root_dist <- root_dist
  tree
}


#' Transform branch lengths for evolutionary models
#'
#' Delegates to the model-specific functions in tree_transforms.R.
#'
#' @param tree A phylo object with root_dist computed (via dist_from_root)
#' @param model Evolutionary model: "BM", "lambda", "kappa", "delta", "EB", "OU"
#' @param evo_model_par Parameter for the model
#'
#' @return Tree with transformed branch lengths
#' @export
transf_branch_lengths <- function(tree, model, evo_model_par) {
  transform_tree(tree, model, evo_model_par)
}
