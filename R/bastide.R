#' Bastide 2021 Tree Traversal Algorithm
#'
#' Implements the efficient tree traversal algorithm from Bastide et al. (2021)
#' for computing the likelihood and quantities needed for analytic gradients.
#'
#' The algorithm has two passes:
#' - Forward (tips to root): accumulates precision P, mean m, likelihood r
#' - Backward (root to tips): computes Q, qn, M, V for gradients and ASR
#'
#' References:
#' Bastide P, Ho LST, Baele G, Lemey P, Suchard MA (2021). Efficient Bayesian
#' inference of general Gaussian models on large phylogenetic trees.
#' Annals of Applied Statistics.
#'
#' @name bastide-algorithm
NULL


#' Moore-Penrose Pseudoinverse with Log Determinant
#'
#' Computes the pseudoinverse and log determinant simultaneously using SVD.
#' Handles rank-deficient matrices gracefully.
#'
#' @param x A matrix
#' @param tol Tolerance for determining rank
#'
#' @return List with components:
#'   \item{inv}{Moore-Penrose pseudoinverse}
#'   \item{logd}{Log of product of positive singular values}
#'   \item{rank}{Numerical rank}
#'
#' @keywords internal
mp_inv <- function(x, tol = sqrt(.Machine$double.eps)) {
  xsvd <- svd(x)
  positive <- xsvd$d > max(tol * xsvd$d[1], 0)

  if (all(positive)) {
    inv <- xsvd$v %*% (1 / xsvd$d * t(xsvd$u))
    logd <- sum(log(xsvd$d))
    rank <- length(xsvd$d)
  } else if (!any(positive)) {
    inv <- matrix(0, ncol(x), nrow(x))
    logd <- -Inf
    rank <- 0
  } else {
    inv <- xsvd$v[, positive, drop = FALSE] %*%
           ((1 / xsvd$d[positive]) * t(xsvd$u[, positive, drop = FALSE]))
    logd <- sum(log(xsvd$d[positive]))
    rank <- sum(positive)
  }

  list(inv = inv, logd = logd, rank = rank)
}


#' Log Pseudo-Determinant
#'
#' Computes the log of the product of positive singular values.
#'
#' @param x A matrix
#' @param tol Tolerance for determining positive eigenvalues
#'
#' @return Log pseudo-determinant
#'
#' @keywords internal
plogdet <- function(x, tol = sqrt(.Machine$double.eps)) {
  xsvd <- svd(x)
  positive <- xsvd$d > max(tol * xsvd$d[1], 0)

  if (all(positive)) {
    sum(log(xsvd$d))
  } else if (!any(positive)) {
    -Inf
  } else {
    sum(log(xsvd$d[positive]))
  }
}


#' Bastide Tree Traversal (Forward and Backward Pass)
#'
#' Main function implementing the Bastide 2021 algorithm. Performs forward pass
#' (tips to root) to compute likelihood, then backward pass (root to tips) to
#' compute quantities needed for gradients and ancestral state reconstruction.
#'
#' IMPORTANT: Tree must be in postorder for this to work correctly.
#'
#' @param tree A phylo object (MUST be in postorder via ape::reorder(tree, "postorder"))
#' @param trait_data Data frame with species in first column, traits in rest
#' @param R Phylogenetic rate matrix (m x m), called Sigma in our other code
#' @param S Residual covariance matrix (m x m), called B in our other code
#' @param mu Optional fixed mean vector. If NULL, ML mean is computed.
#' @param compute_backward Logical, whether to do backward pass (needed for gradients)
#'
#' @return List with components:
#'   \item{logL}{ML log-likelihood (at root node)}
#'   \item{mu}{ML mean estimate}
#'   \item{P}{List of precision matrices at each node}
#'   \item{Pstar}{List of "starred" precision matrices}
#'   \item{m_mat}{Matrix of conditional means (forward pass)}
#'   \item{Q}{List of backward precision matrices}
#'   \item{qn}{Matrix of backward conditional means}
#'   \item{M}{Matrix of final conditional means (ASR)}
#'   \item{V}{List of final conditional variances}
#'   \item{r}{Log-likelihood contributions}
#'   \item{len}{Edge lengths vector}
#'   \item{edge}{Augmented edge matrix}
#'   \item{n}{Number of species}
#'   \item{m}{Number of nodes (n + Nnode)}
#'   \item{N}{Number of observations}
#'   \item{total}{Total = N + m}
#'   \item{p}{Number of traits}
#'   \item{a, b, c}{Linear Gaussian model parameters per node}
#'   \item{delta}{Missing data indicator matrices}
#'
#' @export
bastide_traversal <- function(tree, trait_data, R, S,
                               mu = NULL, compute_backward = TRUE) {

  # Extract data
  data_species_by_row <- trait_data[, 1]
  parent_nodes_of_observations <- match(data_species_by_row, tree$tip.label)

  Y <- trait_data
  Yspecies <- Y[, 1]
  Y <- as.matrix(Y[, -1, drop = FALSE])

  n <- length(tree$tip.label)   # number of species
  m <- n + tree$Nnode           # number of nodes
  N <- nrow(Y)                  # number of observations
  total <- N + m
  p <- ncol(Y)                  # number of traits
  Ip <- diag(p)

  # Species node corresponding to each observation
  node_obs <- match(Yspecies, tree$tip.label)

  # Augmented edge matrix: tree edges + fake root edge
  edge <- rbind(tree$edge, c(m + 1, n + 1))

  # Edge lengths: 0 for obs, tree lengths, 0 for fake root
  len <- c(rep(0, N), tree$edge.length, 0)

  # Initialize storage
  m_mat <- matrix(0, total + 1, p)
  P <- Pstar <- rep(list(matrix(0, p, p)), total + 1)
  delta <- a <- b <- c_mat <- q <- d_vec <- Sigma <- vector("list", total + 1)
  r <- r1 <- r2 <- numeric(total + 1)

  # For backward pass
  Pstar_mkmk <- qn <- qnstar <- M <- matrix(0, total + 1, p)
  Pstar_mk <- Q <- Qstar <- V <- C_back <- rep(list(matrix(0, p, p)), total + 1)

  #---------------------------------------------------------------------------
  # FORWARD PASS (tips to root)
  #---------------------------------------------------------------------------

  for (k in 1:total) {
    if (k <= N) {
      # Processing observation k
      des <- k
      anc <- node_obs[k] + N

      # For observations: a = I, b = 0, c = S
      a[[des]] <- Ip
      b[[des]] <- numeric(p)
      c_mat[[des]] <- S

      # Handle missing data with delta matrix
      not_NA <- which(!is.na(Y[k, ]))
      if (all(is.na(Y[k, ]))) {
        is_NA <- 1:p
      } else {
        is_NA <- (1:p)[-not_NA]
      }

      m_mat[des, not_NA] <- Y[k, not_NA]

      delta[[des]] <- Ip
      if (length(is_NA) > 0) {
        diag(delta[[des]])[is_NA] <- 0
      }

      # Pstar = ginv(delta %*% S %*% delta)
      Pstar[[des]] <- mp_inv(delta[[des]] %*% S %*% delta[[des]])$inv

      # Likelihood contribution: -0.5 * rank * log(2pi) + 0.5 * log|Pstar|
      rank_Pstar <- sum(diag(delta[[des]]))  # = number of observed traits
      r[des] <- -0.5 * rank_Pstar * log(2 * pi) + 0.5 * plogdet(Pstar[[des]])

    } else {
      # Processing internal branch k-N
      des <- edge[k - N, 2] + N
      anc <- edge[k - N, 1] + N

      # For BM: a = I, b = 0, c = len * R
      a[[des]] <- Ip
      b[[des]] <- numeric(p)
      Sigma[[des]] <- len[k] * R
      c_mat[[des]] <- Sigma[[des]]

      # Update m_mat: m = ginv(P) %*% m
      m_mat[des, ] <- mp_inv(P[[des]])$inv %*% m_mat[des, ]

      if (k < (m + N)) {
        # Not at root yet - compute Pstar
        Sigma_inv <- mp_inv(Sigma[[des]])$inv
        invPSigma_P <- solve(P[[des]] + Sigma_inv) %*% P[[des]]
        Pstar[[des]] <- P[[des]] - P[[des]] %*% invPSigma_P

        # Likelihood contribution
        r[des] <- r1[des] +
                  0.5 * as.numeric(t(m_mat[des, ]) %*% P[[des]] %*% m_mat[des, ]) -
                  0.5 * r2[des] +
                  0.5 * plogdet(Ip - invPSigma_P)
      }
    }

    # Propagate to ancestor
    if (k < (m + N)) {
      ta_Pstar <- t(a[[des]]) %*% Pstar[[des]]
      P[[anc]] <- P[[anc]] + ta_Pstar %*% a[[des]]
      m_mat[anc, ] <- m_mat[anc, ] + ta_Pstar %*% (m_mat[des, ] - b[[des]])
      r1[anc] <- r1[anc] + r[des]
      r2[anc] <- r2[anc] + as.numeric(m_mat[des, ] %*% Pstar[[des]] %*% m_mat[des, ])
    } else {
      # At root: finalize likelihood
      r[des] <- r1[des] +
                0.5 * as.numeric(t(m_mat[des, ]) %*% P[[des]] %*% m_mat[des, ]) -
                0.5 * r2[des]
    }
  }

  # Get root index
  root_idx <- edge[total - N, 2] + N

  # ML mean at root
  if (is.null(mu)) {
    mu <- m_mat[root_idx, ]
  }

  # Final log-likelihood is stored at root node
  logL <- r[root_idx]

  #---------------------------------------------------------------------------
  # BACKWARD PASS (root to tips) - for gradients and ASR
  #---------------------------------------------------------------------------

  if (compute_backward) {
    for (k in total:1) {
      if (k > N) {
        des <- edge[k - N, 2] + N
        anc <- edge[k - N, 1] + N
      } else {
        des <- k
        anc <- node_obs[k] + N
      }

      if ((anc - N) == (n + 1)) {
        # Direct descendant of root (anc-N = n+1 is the root node index)
        Q[[des]] <- mp_inv(c_mat[[des]])$inv
        qn[des, ] <- a[[des]] %*% mu + b[[des]]

      } else if (k != total) {
        # Find sibling nodes
        if (des > N) {
          # Internal node: siblings are other tree edges from same parent
          not_des <- edge[-(k - N), ][(edge[-(k - N), 1] == (anc - N)), 2] + N
        } else {
          # Observation: siblings are other obs with same species
          pk <- parent_nodes_of_observations[k]
          ls <- which(parent_nodes_of_observations == pk)
          not_des <- ls[ls != k]
        }

        # Accumulate sibling contributions
        if (length(not_des) > 0) {
          for (l in not_des) {
            Pstar_mk[[des]] <- Pstar_mk[[des]] +
                               t(a[[l]]) %*% Pstar[[l]] %*% a[[l]]
            Pstar_mkmk[des, ] <- Pstar_mkmk[des, ] +
                                 t(a[[l]]) %*% Pstar[[l]] %*% (m_mat[l, ] - b[[l]])
          }
        }

        # Compute Q and qn
        Qstar[[des]] <- Pstar_mk[[des]] + Q[[anc]]
        qnstar[des, ] <- mp_inv(Qstar[[des]])$inv %*%
                         (Pstar_mkmk[des, ] + Q[[anc]] %*% qn[anc, ])
        Q[[des]] <- mp_inv(a[[des]] %*% mp_inv(Qstar[[des]])$inv %*% t(a[[des]]) +
                           c_mat[[des]])$inv
        qn[des, ] <- a[[des]] %*% qnstar[des, ] + b[[des]]
      }

      # Compute M and V (conditional mean and variance)
      if (k == total) {
        M[des, ] <- m_mat[des, ]
      } else if ((k > N) || (sum(diag(delta[[k]])) == 0)) {
        # Internal node or fully missing observation
        V[[des]] <- mp_inv(P[[des]] + Q[[des]])$inv
        M[des, ] <- V[[des]] %*% (P[[des]] %*% m_mat[des, ] + Q[[des]] %*% qn[des, ])
      } else {
        # Observation with some observed traits
        n_obs_traits <- sum(diag(delta[[k]]))
        if (n_obs_traits == p) {
          # Fully observed
          M[k, ] <- Y[k, ]
        } else {
          # Partially observed
          is_na <- is.na(Y[k, ])
          PIkm <- Ip[is_na, , drop = FALSE]   # Selector for missing
          PIko <- Ip[!is_na, , drop = FALSE]  # Selector for observed

          V[[k]][is_na, is_na] <- mp_inv(PIkm %*% Q[[k]] %*% t(PIkm))$inv
          M[k, !is_na] <- Y[k, !is_na]
          M[k, is_na] <- PIkm %*% qn[k, ] -
                         V[[k]][is_na, is_na] %*% PIkm %*% Q[[k]] %*% t(PIko) %*%
                         (Y[k, !is_na] - PIko %*% qn[k, ])
        }
      }
    }
  }

  #---------------------------------------------------------------------------
  # Return results
  #---------------------------------------------------------------------------

  list(
    logL = logL,
    mu = mu,
    P = P,
    Pstar = Pstar,
    m_mat = m_mat,
    r = r,
    r1 = r1,
    r2 = r2,
    Q = Q,
    Qstar = Qstar,
    qn = qn,
    qnstar = qnstar,
    M = M,
    V = V,
    a = a,
    b = b,
    c = c_mat,
    Sigma = Sigma,
    delta = delta,
    len = len,
    edge = edge,
    n = n,
    m = m,
    N = N,
    total = total,
    p = p,
    node_obs = node_obs,
    parent_nodes_of_observations = parent_nodes_of_observations
  )
}


#' Compute ML Log-Likelihood Using Bastide Traversal
#'
#' Wrapper function that computes only the likelihood (no backward pass).
#' Note: This computes ML likelihood, not REML.
#'
#' Optionally applies a tree transformation for non-BM models before
#' running the traversal. The tree should already have root_dist computed
#' (via dist_from_root) if using a non-BM model.
#'
#' @param R Phylogenetic rate matrix (Sigma)
#' @param S Residual covariance matrix (B)
#' @param tree A phylo object (should be postorder)
#' @param trait_data Data frame with species and traits
#' @param model Evolutionary model: "BM", "lambda", "kappa", "delta", "EB", "OU"
#' @param evo_model_par Model parameter value (ignored for BM)
#'
#' @return Log-likelihood value
#'
#' @export
loglik_bastide <- function(R, S, tree, trait_data,
                           model = "BM", evo_model_par = NULL) {
  # Apply tree transform for non-BM models
  if (model != "BM" && !is.null(evo_model_par)) {
    tree <- transform_tree(tree, model, evo_model_par)
  }
  result <- bastide_traversal(tree, trait_data, R, S, compute_backward = FALSE)
  result$logL
}
