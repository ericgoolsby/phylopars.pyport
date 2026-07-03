#' Advanced Backends for Bastide Traversal and EM
#'
#' R wrappers for C++ implementations of the Bastide algorithm,
#' gradient computation, and EM algorithm.
#'
#' @name backends-advanced
NULL


#' Compute Log-Likelihood Using C++ Bastide Traversal
#'
#' Fast C++ implementation of the Bastide 2021 ML likelihood computation.
#'
#' @param R Phylogenetic rate matrix (p x p)
#' @param S Residual covariance matrix (p x p)
#' @param tree A phylo object (should be postorder)
#' @param trait_data Data frame with species in first column, traits in rest
#' @param backend Which backend to use: "armadillo" or "eigen"
#'
#' @return ML log-likelihood value
#'
#' @export
loglik_bastide_cpp <- function(R, S, tree, trait_data, backend = "armadillo") {

  # Prepare data
  Y <- as.matrix(trait_data[, -1, drop = FALSE])
  species <- as.character(trait_data[, 1])
  species_idx <- match(species, tree$tip.label)

  n_tips <- length(tree$tip.label)
  n_nodes <- n_tips + tree$Nnode

  if (backend == "armadillo") {
    loglik_bastide_arma(
      R = R,
      S = S,
      edge = tree$edge,
      edge_length = tree$edge.length,
      traits = Y,
      species_idx = species_idx,
      n_tips = n_tips,
      n_nodes = n_nodes
    )
  } else if (backend == "eigen") {
    loglik_bastide_eigen(
      R = R,
      S = S,
      edge = tree$edge,
      edge_length = tree$edge.length,
      traits = Y,
      species_idx = species_idx,
      n_tips = n_tips,
      n_nodes = n_nodes
    )
  } else {
    stop("Unknown backend: ", backend)
  }
}


#' Compute Gradient Using C++ Bastide Traversal
#'
#' Fast C++ implementation of the Bastide 2021 gradient computation.
#'
#' @param R Phylogenetic rate matrix (p x p)
#' @param S Residual covariance matrix (p x p)
#' @param tree A phylo object (should be postorder)
#' @param trait_data Data frame with species in first column, traits in rest
#' @param backend Which backend to use: "armadillo" or "eigen"
#'
#' @return List with logL, mu, grad_R, grad_S
#'
#' @export
bastide_gradient_cpp <- function(R, S, tree, trait_data, backend = "armadillo") {

  # Prepare data
  Y <- as.matrix(trait_data[, -1, drop = FALSE])
  species <- as.character(trait_data[, 1])
  species_idx <- match(species, tree$tip.label)

  n_tips <- length(tree$tip.label)
  n_nodes <- n_tips + tree$Nnode

  if (backend == "armadillo") {
    bastide_gradient_arma(
      R = R,
      S = S,
      edge = tree$edge,
      edge_length = tree$edge.length,
      traits = Y,
      species_idx = species_idx,
      n_tips = n_tips,
      n_nodes = n_nodes
    )
  } else if (backend == "eigen") {
    bastide_gradient_eigen(
      R = R,
      S = S,
      edge = tree$edge,
      edge_length = tree$edge.length,
      traits = Y,
      species_idx = species_idx,
      n_tips = n_tips,
      n_nodes = n_nodes
    )
  } else {
    stop("Unknown backend: ", backend)
  }
}


#' Run EM Algorithm Using C++ Backend
#'
#' Fast C++ implementation of the EM algorithm for phylogenetic mixed models.
#'
#' @param tree A phylo object (should be postorder)
#' @param trait_data Data frame with species in first column, traits in rest
#' @param R_init Initial R matrix (default: identity)
#' @param S_init Initial S matrix (default: identity)
#' @param mu_init Initial mean (default: column means)
#' @param max_iter Maximum iterations
#' @param tol Convergence tolerance
#' @param verbose Print progress
#' @param backend Which backend to use: "armadillo" or "eigen"
#'
#' @return List with R, S, mu, logL, converged, iterations
#'
#' @export
em_fit_cpp <- function(tree, trait_data,
                       R_init = NULL, S_init = NULL, mu_init = NULL,
                       max_iter = 100, tol = 1e-6, verbose = FALSE,
                       backend = "armadillo") {

  # Prepare data
  Y <- as.matrix(trait_data[, -1, drop = FALSE])
  p <- ncol(Y)
  species <- as.character(trait_data[, 1])
  species_idx <- match(species, tree$tip.label)

  n_tips <- length(tree$tip.label)
  n_nodes <- n_tips + tree$Nnode

  # Defaults
  if (is.null(R_init)) R_init <- diag(p)
  if (is.null(S_init)) S_init <- diag(p)
  if (is.null(mu_init)) mu_init <- colMeans(Y, na.rm = TRUE)

  if (backend == "armadillo") {
    em_fit_arma(
      R_init = R_init,
      S_init = S_init,
      mu_init = mu_init,
      edge = tree$edge,
      edge_length = tree$edge.length,
      traits = Y,
      species_idx = species_idx,
      n_tips = n_tips,
      n_nodes = n_nodes,
      max_iter = max_iter,
      tol = tol,
      verbose = verbose
    )
  } else if (backend == "eigen") {
    em_fit_eigen(
      R_init = R_init,
      S_init = S_init,
      mu_init = mu_init,
      edge = tree$edge,
      edge_length = tree$edge.length,
      traits = Y,
      species_idx = species_idx,
      n_tips = n_tips,
      n_nodes = n_nodes,
      max_iter = max_iter,
      tol = tol,
      verbose = verbose
    )
  } else {
    stop("Unknown backend: ", backend)
  }
}


#' Unified Optimizer with C++ Backend
#'
#' @param tree A phylo object (should be postorder)
#' @param trait_data Data frame with species and traits
#' @param method Character vector of methods to run
#' @param likelihood Which likelihood to use: "ML" (default) or "REML"
#' @param control Named list with per-method settings
#' @param R_init Initial R matrix
#' @param S_init Initial S matrix
#' @param mu_init Initial mean
#' @param verbose Print progress
#' @param backend "R", "armadillo", or "eigen"
#'
#' @return List with R, S, mu, logL, likelihood, history, converged
#'
#' @details
#' When likelihood="REML", the EM method is not available (skipped automatically).
#' REML uses numerical gradients (no analytic gradients implemented yet).
#'
#' @export
optimize_phylopars_cpp <- function(tree, trait_data,
                                   method = c("EM", "BFGS"),
                                   likelihood = c("ML", "REML"),
                                   control = list(),
                                   R_init = NULL, S_init = NULL, mu_init = NULL,
                                   verbose = FALSE,
                                   backend = "armadillo") {

  # Validate likelihood
  likelihood <- match.arg(likelihood)

  # Validate methods
  valid_methods <- c("EM", "BFGS", "Nelder-Mead", "L-BFGS-B")
  for (m in method) {
    if (!(m %in% valid_methods)) {
      stop("Unknown method: ", m, ". Valid methods: ", paste(valid_methods, collapse = ", "))
    }
  }

  # EM not available for REML
  if (likelihood == "REML" && "EM" %in% method) {
    if (verbose) message("Note: EM not available for REML, skipping EM steps")
    method <- method[method != "EM"]
    if (length(method) == 0) {
      method <- "BFGS"
    }
  }

  # Prepare data
  Y <- as.matrix(trait_data[, -1, drop = FALSE])
  p <- ncol(Y)
  species <- as.character(trait_data[, 1])
  species_idx <- match(species, tree$tip.label)

  n_tips <- length(tree$tip.label)
  n_nodes <- n_tips + tree$Nnode

  # Initialize parameters
  if (is.null(R_init)) R_init <- diag(p)
  if (is.null(S_init)) S_init <- diag(p)
  if (is.null(mu_init)) mu_init <- colMeans(Y, na.rm = TRUE)

  R <- R_init
  S <- S_init
  mu <- mu_init

  # Default control settings
  default_control <- list(
    EM = list(max_iter = 50, tol = 1e-4),
    BFGS = list(max_iter = 100, tol = 1e-6, use_gradient = TRUE),
    "Nelder-Mead" = list(max_iter = 500, tol = 1e-8),
    "L-BFGS-B" = list(max_iter = 100, tol = 1e-6)
  )

  # Merge with user control
  for (m in names(default_control)) {
    if (is.null(control[[m]])) {
      control[[m]] <- default_control[[m]]
    } else {
      for (nm in names(default_control[[m]])) {
        if (is.null(control[[m]][[nm]])) {
          control[[m]][[nm]] <- default_control[[m]][[nm]]
        }
      }
    }
  }

  # Likelihood function based on backend and likelihood type
  if (likelihood == "ML") {
    loglik_fn <- if (backend == "R") {
      function(R, S) loglik_bastide(R, S, tree, trait_data)
    } else {
      function(R, S) loglik_bastide_cpp(R, S, tree, trait_data, backend)
    }
  } else {
    # REML - use C++ backend if armadillo, else R
    loglik_fn <- function(R, S) {
      theta <- matrices_to_theta(R, S)
      if (backend == "armadillo") {
        restricted_loglik(theta, tree, trait_data, backend = "armadillo")
      } else if (backend == "eigen") {
        restricted_loglik(theta, tree, trait_data, backend = "eigen")
      } else {
        restricted_loglik(theta, tree, trait_data, backend = "R")
      }
    }
  }

  # Track history
  history <- list()
  overall_converged <- FALSE

  # Run methods in sequence
  for (i in seq_along(method)) {
    m <- method[i]
    ctrl <- control[[m]]

    if (verbose) {
      cat("\n=== Running", m, "(", likelihood, ", step", i, "of", length(method), ") ===\n")
    }

    start_logL <- loglik_fn(R, S)

    if (m == "EM") {
      # Use C++ EM
      if (backend == "R") {
        result <- em_fit(tree, trait_data, R, S, mu, ctrl$max_iter, ctrl$tol, verbose)
      } else {
        result <- em_fit_cpp(tree, trait_data, R, S, mu,
                             ctrl$max_iter, ctrl$tol, verbose, backend)
      }
      R <- result$R
      S <- result$S
      if (!is.null(result$mu)) mu <- as.numeric(result$mu)
      end_logL <- result$logL
      converged <- result$converged
      iterations <- result$iterations

    } else {
      # BFGS/Nelder-Mead/L-BFGS-B
      theta <- matrices_to_theta(R, S)

      # Objective function
      fn <- function(theta) {
        mats <- tryCatch(theta_to_matrices(theta, p), error = function(e) NULL)
        if (is.null(mats)) return(1e100)
        logL <- tryCatch(loglik_fn(mats$R, mats$S), error = function(e) -1e100)
        if (!is.finite(logL) || logL < -1e99) return(1e100)
        -logL
      }

      # Gradient function (if using C++ and BFGS) - only for ML
      # REML doesn't have analytic gradients, so use numerical
      gr <- NULL
      if (isTRUE(ctrl$use_gradient) && m %in% c("BFGS", "L-BFGS-B") && likelihood == "ML") {
        gr <- function(theta) {
          mats <- tryCatch(theta_to_matrices(theta, p), error = function(e) NULL)
          if (is.null(mats)) return(rep(0, length(theta)))
          grad <- tryCatch({
            if (backend == "R") {
              g <- bastide_gradient(tree, trait_data, mats$R, mats$S)
            } else {
              g <- bastide_gradient_cpp(mats$R, mats$S, tree, trait_data, backend)
            }
            -gradient_to_theta(g$grad_R, g$grad_S, mats$R, mats$S)
          }, error = function(e) rep(0, length(theta)))
          grad[!is.finite(grad)] <- 0
          grad
        }
      }

      # Run optim
      opt_result <- tryCatch({
        optim(
          par = theta,
          fn = fn,
          gr = gr,
          method = if (m == "Nelder-Mead") "Nelder-Mead" else m,
          control = list(maxit = ctrl$max_iter, reltol = ctrl$tol)
        )
      }, error = function(e) {
        if (verbose) cat("  Optimizer failed:", e$message, "\n")
        list(par = theta, value = fn(theta), convergence = 1)
      })

      mats <- tryCatch(theta_to_matrices(opt_result$par, p),
                       error = function(e) list(R = R, S = S))
      R <- mats$R
      S <- mats$S
      end_logL <- -opt_result$value
      converged <- (opt_result$convergence == 0)
      iterations <- NA
    }

    # Store in history
    history[[i]] <- list(
      method = m,
      start_logL = start_logL,
      end_logL = end_logL,
      converged = converged,
      iterations = iterations
    )

    if (verbose) {
      cat("  logL:", round(start_logL, 4), "->", round(end_logL, 4),
          "(", ifelse(converged, "converged", "not converged"), ")\n")
    }

    if (converged) {
      overall_converged <- TRUE
    }
  }

  list(
    R = R,
    S = S,
    mu = mu,
    logL = loglik_fn(R, S),
    likelihood = likelihood,
    history = history,
    converged = overall_converged
  )
}
