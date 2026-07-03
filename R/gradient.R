#' Analytic Gradient for Phylogenetic Mixed Model
#'
#' Computes the analytic gradient of the log-likelihood with respect to the
#' phylogenetic rate matrix (R/Sigma) and residual covariance matrix (S/B)
#' using the Bastide 2021 algorithm.
#'
#' The gradient formulas (from Bastide et al. 2021) are:
#'
#' For the phylogenetic rate R (summed over internal edges k):
#'   dlogL/dR = sum_k -0.5 * t_k * Q_k %*% (Q_k^{-1} - (M_k - qn_k)(M_k - qn_k)' - V_k) %*% Q_k
#'
#' For the residual covariance S (summed over observations k):
#'   dlogL/dS = sum_k -0.5 * Q_k %*% (Q_k^{-1} - (M_k - qn_k)(M_k - qn_k)' - V_k) %*% Q_k
#'
#' Where:
#' - Q_k is the backward precision matrix at node k
#' - M_k is the conditional mean (ASR estimate)
#' - qn_k is the backward conditional mean
#' - V_k is the conditional variance
#' - t_k is the edge length
#'
#' @name gradient
NULL


#' Compute Analytic Gradient from Bastide Traversal Results
#'
#' Given the results from bastide_traversal(), computes the gradient of the
#' log-likelihood with respect to R and S.
#'
#' @param b Result from bastide_traversal() (must include backward pass quantities)
#'
#' @return List with components:
#'   \item{grad_R}{Gradient with respect to phylogenetic rate matrix (p x p)}
#'   \item{grad_S}{Gradient with respect to residual covariance matrix (p x p)}
#'
#' @export
gradient_from_bastide <- function(b) {

  p <- b$p
  N <- b$N
  total <- b$total
  edge <- b$edge

  # Initialize gradients
  grad_R <- matrix(0, p, p)
  grad_S <- matrix(0, p, p)

  # Gradient with respect to R (phylogenetic rate)
  # Sum over internal edges (k = N+1 to total)
  for (k in (N + 1):total) {
    des <- edge[k - N, 2] + N
    t_k <- b$len[k]

    if (t_k > 0) {  # Skip zero-length edges
      Q_des <- b$Q[[des]]
      M_des <- b$M[des, ]
      qn_des <- b$qn[des, ]
      V_des <- b$V[[des]]

      # Compute (M - qn)(M - qn)'
      diff <- M_des - qn_des
      outer_diff <- outer(diff, diff)

      # Gradient contribution: -0.5 * t_k * Q %*% (Q^{-1} - diff*diff' - V) %*% Q
      Q_inv <- mp_inv(Q_des)$inv
      inner <- Q_inv - outer_diff - V_des
      grad_R <- grad_R - 0.5 * t_k * Q_des %*% inner %*% Q_des
    }
  }

  # Gradient with respect to S (residual covariance)
  # Sum over observations (k = 1 to N)
  for (k in 1:N) {
    des <- k  # For observations, des = k

    Q_des <- b$Q[[des]]
    M_des <- b$M[des, ]
    qn_des <- b$qn[des, ]
    V_des <- b$V[[des]]

    # Compute (M - qn)(M - qn)'
    diff <- M_des - qn_des
    outer_diff <- outer(diff, diff)

    # Gradient contribution: -0.5 * Q %*% (Q^{-1} - diff*diff' - V) %*% Q
    Q_inv <- mp_inv(Q_des)$inv
    inner <- Q_inv - outer_diff - V_des
    grad_S <- grad_S - 0.5 * Q_des %*% inner %*% Q_des
  }

  list(grad_R = grad_R, grad_S = grad_S)
}


#' Compute Gradient of Log-Likelihood w.r.t. Evolutionary Parameter
#'
#' Uses the chain rule: dlogL/d(evo_param) = sum_k (dlogL/dt'_k) * (dt'_k/d(evo_param))
#'
#' The per-edge dlogL/dt'_k is computed from the backward pass quantities:
#'   A_k = Q_k %*% (Q_k^{-1} - (M_k - qn_k)(M_k - qn_k)' - V_k) %*% Q_k
#'   dlogL/dt'_k = -0.5 * tr(A_k %*% R)
#'
#' Only tree edges contribute (not observation edges or the fake root edge).
#'
#' @param b Result from bastide_traversal() (must include backward pass quantities)
#' @param R Phylogenetic rate matrix (p x p)
#' @param tree The original (pre-transform) phylo object with root_dist computed
#' @param model Character string: evolutionary model
#' @param evo_par Current evolutionary parameter value
#'
#' @return Scalar gradient dlogL/d(evo_param)
#'
#' @keywords internal
evo_param_gradient <- function(b, R, tree, model, evo_par) {

  p <- b$p
  N <- b$N
  total <- b$total
  edge <- b$edge

  # Compute dt'_k/d(evo_param) for each tree edge
  jac <- transform_jacobian(tree, model, evo_par)

  # Sum dlogL/dt'_k * dt'_k/d(evo_param) over tree edges
  grad_evo <- 0

  for (k in (N + 1):total) {
    # k indexes augmented edges: first N are obs, then tree edges, last is fake root
    tree_edge_idx <- k - N  # index into b$edge / tree$edge

    # Skip the fake root edge (last augmented edge)
    if (tree_edge_idx > nrow(tree$edge)) next

    des <- edge[tree_edge_idx, 2] + N
    t_k <- b$len[k]

    if (t_k > 0) {
      Q_des <- b$Q[[des]]
      M_des <- b$M[des, ]
      qn_des <- b$qn[des, ]
      V_des <- b$V[[des]]

      # Compute A_k = Q %*% (Q^{-1} - (M-qn)(M-qn)' - V) %*% Q
      diff <- M_des - qn_des
      outer_diff <- outer(diff, diff)
      Q_inv <- mp_inv(Q_des)$inv
      inner <- Q_inv - outer_diff - V_des
      A_k <- Q_des %*% inner %*% Q_des

      # dlogL/dt'_k = -0.5 * tr(A_k %*% R) = -0.5 * sum(A_k * R)
      dlogL_dt_k <- -0.5 * sum(A_k * R)

      # Chain rule with Jacobian
      grad_evo <- grad_evo + dlogL_dt_k * jac[tree_edge_idx]
    }
  }

  grad_evo
}


#' Compute Full Gradient (Bastide Traversal + Gradient Computation)
#'
#' Convenience function that runs the Bastide traversal and computes the
#' gradient in one call. Supports both BM and non-BM models.
#'
#' @param tree A phylo object (should be postorder). For non-BM models, this
#'   should be the ORIGINAL tree (pre-transform) with root_dist computed.
#' @param trait_data Data frame with species in first column, traits in rest
#' @param R Phylogenetic rate matrix (Sigma)
#' @param S Residual covariance matrix (B)
#' @param model Evolutionary model: "BM", "lambda", "kappa", "delta", "EB", "OU"
#' @param evo_model_par Evolutionary parameter value (ignored for BM)
#'
#' @return List with components:
#'   \item{logL}{Log-likelihood}
#'   \item{mu}{Mean estimate}
#'   \item{grad_R}{Gradient with respect to R}
#'   \item{grad_S}{Gradient with respect to S}
#'   \item{grad_evo}{Gradient with respect to evolutionary parameter (NULL for BM)}
#'
#' @export
bastide_gradient <- function(tree, trait_data, R, S,
                             model = "BM", evo_model_par = NULL) {

  # For non-BM models, transform the tree before running traversal
  if (model != "BM" && !is.null(evo_model_par)) {
    transformed_tree <- transform_tree(tree, model, evo_model_par)
  } else {
    transformed_tree <- tree
  }

  # Run traversal with backward pass on transformed tree
  b <- bastide_traversal(transformed_tree, trait_data, R, S,
                         compute_backward = TRUE)

  # Compute R and S gradients
  grad <- gradient_from_bastide(b)

  # Compute evolutionary parameter gradient
  grad_evo <- NULL
  if (model != "BM" && !is.null(evo_model_par)) {
    grad_evo <- evo_param_gradient(b, R, tree, model, evo_model_par)
  }

  list(
    logL = b$logL,
    mu = b$mu,
    grad_R = grad$grad_R,
    grad_S = grad$grad_S,
    grad_evo = grad_evo
  )
}


#' Numerical Gradient for Comparison
#'
#' Computes the gradient numerically using finite differences.
#' Useful for validating the analytic gradient. Supports non-BM models.
#'
#' @param tree A phylo object (should be postorder). For non-BM models,
#'   should have root_dist computed.
#' @param trait_data Data frame with species in first column, traits in rest
#' @param R Phylogenetic rate matrix
#' @param S Residual covariance matrix
#' @param model Evolutionary model: "BM", "lambda", "kappa", "delta", "EB", "OU"
#' @param evo_model_par Evolutionary parameter value (ignored for BM)
#' @param eps Step size for finite differences (default 1e-6)
#'
#' @return List with grad_R, grad_S, and optionally grad_evo (numerical gradients)
#'
#' @export
numerical_gradient <- function(tree, trait_data, R, S,
                               model = "BM", evo_model_par = NULL,
                               eps = 1e-6) {

  p <- nrow(R)

  # Numerical gradient for R
  grad_R <- matrix(0, p, p)
  for (i in 1:p) {
    for (j in 1:p) {
      # Symmetric perturbation
      R_plus <- R
      R_plus[i, j] <- R_plus[i, j] + eps
      if (i != j) R_plus[j, i] <- R_plus[j, i] + eps

      R_minus <- R
      R_minus[i, j] <- R_minus[i, j] - eps
      if (i != j) R_minus[j, i] <- R_minus[j, i] - eps

      logL_plus <- loglik_bastide(R_plus, S, tree, trait_data,
                                  model = model, evo_model_par = evo_model_par)
      logL_minus <- loglik_bastide(R_minus, S, tree, trait_data,
                                   model = model, evo_model_par = evo_model_par)

      grad_R[i, j] <- (logL_plus - logL_minus) / (2 * eps)
    }
  }

  # Numerical gradient for S
  grad_S <- matrix(0, p, p)
  for (i in 1:p) {
    for (j in 1:p) {
      S_plus <- S
      S_plus[i, j] <- S_plus[i, j] + eps
      if (i != j) S_plus[j, i] <- S_plus[j, i] + eps

      S_minus <- S
      S_minus[i, j] <- S_minus[i, j] - eps
      if (i != j) S_minus[j, i] <- S_minus[j, i] - eps

      logL_plus <- loglik_bastide(R, S_plus, tree, trait_data,
                                  model = model, evo_model_par = evo_model_par)
      logL_minus <- loglik_bastide(R, S_minus, tree, trait_data,
                                   model = model, evo_model_par = evo_model_par)

      grad_S[i, j] <- (logL_plus - logL_minus) / (2 * eps)
    }
  }

  # Symmetrize (for off-diagonal, we computed the sum of partials)
  for (i in 1:p) {
    for (j in 1:p) {
      if (i != j) {
        grad_R[i, j] <- grad_R[i, j] / 2
        grad_S[i, j] <- grad_S[i, j] / 2
      }
    }
  }

  result <- list(grad_R = grad_R, grad_S = grad_S)

  # Numerical gradient for evolutionary parameter
  if (model != "BM" && !is.null(evo_model_par)) {
    logL_plus <- loglik_bastide(R, S, tree, trait_data,
                                model = model, evo_model_par = evo_model_par + eps)
    logL_minus <- loglik_bastide(R, S, tree, trait_data,
                                 model = model, evo_model_par = evo_model_par - eps)
    result$grad_evo <- (logL_plus - logL_minus) / (2 * eps)
  }

  result
}
