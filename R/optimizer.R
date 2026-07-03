#' Unified Optimizer for Phylogenetic Mixed Models
#'
#' Provides a flexible optimization interface that supports:
#' - Multiple optimization methods: EM, BFGS (with analytic gradients), Nelder-Mead
#' - Method sequencing: run methods in user-specified order
#' - Per-method control parameters
#'
#' @name optimizer
NULL


#' Optimize Phylogenetic Model Parameters
#'
#' Main optimization function that supports running multiple methods in sequence.
#'
#' @param tree A phylo object (should be postorder)
#' @param trait_data Data frame with species in first column, traits in rest
#' @param method Character vector of methods to run in order.
#'   Options: "EM", "BFGS", "Nelder-Mead", "L-BFGS-B"
#'   Default: c("EM", "BFGS")
#' @param likelihood Which likelihood to use: "ML" (default) or "REML"
#' @param model Evolutionary model: "BM", "lambda", "kappa", "delta", "EB", "OU"
#' @param evo_model_par_start Starting value for evolutionary parameter (NULL = default)
#' @param evo_model_par_fixed Fixed value for evolutionary parameter (NULL = optimize)
#' @param control Named list with per-method settings. Each element should be
#'   a list with max_iter, tol, etc. See details.
#' @param R_init Initial phylogenetic rate matrix
#' @param S_init Initial residual covariance matrix
#' @param mu_init Initial mean vector
#' @param verbose Print progress
#'
#' @return List with:
#'   \item{R}{Final phylogenetic rate matrix}
#'   \item{S}{Final residual covariance matrix}
#'   \item{mu}{Final mean vector}
#'   \item{logL}{Final log-likelihood}
#'   \item{likelihood}{Which likelihood was used}
#'   \item{model}{Evolutionary model used}
#'   \item{evo_model_par}{Estimated evolutionary model parameter}
#'   \item{history}{List with results from each method step}
#'   \item{converged}{Overall convergence status}
#'
#' @details
#' ## Likelihood types
#'
#' - **ML**: Maximum likelihood (Bastide 2021 algorithm with analytic gradients)
#' - **REML**: Restricted maximum likelihood (integrates out mean parameter)
#'
#' Note: REML does not have analytic gradients implemented, so gradient-based
#' methods (BFGS, L-BFGS-B) will use numerical gradients when likelihood="REML".
#' EM is not available for REML.
#'
#' ## Control parameters
#'
#' The control list can have entries for each method type:
#'
#' ```
#' control = list(
#'   EM = list(max_iter = 50, tol = 1e-4),
#'   BFGS = list(max_iter = 100, tol = 1e-6, use_gradient = TRUE),
#'   "Nelder-Mead" = list(max_iter = 500, tol = 1e-8),
#'   "L-BFGS-B" = list(max_iter = 100, tol = 1e-6)
#' )
#' ```
#'
#' ## Method characteristics
#'
#' - **EM**: Most stable, good for finding general region. Slow near optimum. (ML only)
#' - **BFGS**: Fast with analytic gradients (ML) or numerical (REML).
#' - **Nelder-Mead**: No gradients needed. Robust but slow.
#' - **L-BFGS-B**: Limited memory BFGS. Good for high dimensions.
#'
#' ## Example sequences
#'
#' - `c("EM")`: Pure EM (stable but slow, ML only)
#' - `c("BFGS")`: Pure gradient descent (fast but may fail)
#' - `c("EM", "BFGS")`: EM to get close, BFGS to polish (recommended for ML)
#' - `c("BFGS", "Nelder-Mead")`: Recommended for REML
#' - `c("EM", "BFGS", "Nelder-Mead")`: EM, then BFGS, then polish with NM
#'
#' @export
optimize_phylopars <- function(tree, trait_data,
                                method = c("EM", "BFGS"),
                                likelihood = c("ML", "REML"),
                                model = "BM",
                                evo_model_par_start = NULL,
                                evo_model_par_fixed = NULL,
                                control = list(),
                                R_init = NULL, S_init = NULL, mu_init = NULL,
                                verbose = FALSE,
                                constraints = NULL) {

  # Validate likelihood
  likelihood <- match.arg(likelihood)

  # Validate methods
  valid_methods <- c("EM", "BFGS", "Nelder-Mead", "L-BFGS-B")
  for (m in method) {
    if (!(m %in% valid_methods)) {
      stop("Unknown method: ", m, ". Valid methods: ", paste(valid_methods, collapse = ", "))
    }
  }

  # EM not available for REML, non-BM models, or constrained optimization
  em_excluded <- (likelihood == "REML" || model != "BM" || !is.null(constraints))
  if (em_excluded && "EM" %in% method) {
    if (verbose) {
      if (likelihood == "REML") message("Note: EM not available for REML, skipping EM steps")
      if (model != "BM") message("Note: EM not available for non-BM models, skipping EM steps")
      if (!is.null(constraints)) message("Note: EM not available with constraints, skipping EM steps")
    }
    method <- method[method != "EM"]
    if (length(method) == 0) {
      method <- "BFGS"  # Fall back to BFGS
    }
  }

  # Precompute root distances for non-BM models
  if (model != "BM") {
    tree <- dist_from_root(tree)
  }

  # If evo_model_par is fixed, apply the transform to the tree once
  # and treat it as BM from that point on
  effective_model <- model
  evo_model_par <- evo_model_par_fixed
  if (model != "BM" && !is.null(evo_model_par_fixed)) {
    tree <- transform_tree(tree, model, evo_model_par_fixed)
    effective_model <- "BM"  # After transform, it's just BM
  }

  # Set default starting value for evo_model_par if needed
  optimize_evo <- (effective_model != "BM" && is.null(evo_model_par_fixed))
  if (optimize_evo && is.null(evo_model_par_start)) {
    evo_model_par_start <- switch(model,
      "lambda" = 0.5,
      "kappa"  = 1.0,
      "delta"  = 1.0,
      "EB"     = -0.1,
      "OU"     = 0.1
    )
    evo_model_par <- evo_model_par_start
  } else if (optimize_evo) {
    evo_model_par <- evo_model_par_start
  }

  # Get dimensions
  Y <- as.matrix(trait_data[, -1, drop = FALSE])
  p <- ncol(Y)

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
      # Merge: user settings override defaults
      for (nm in names(default_control[[m]])) {
        if (is.null(control[[m]][[nm]])) {
          control[[m]][[nm]] <- default_control[[m]][[nm]]
        }
      }
    }
  }

  # Track history
  history <- list()
  overall_converged <- FALSE

  # Helper to compute likelihood based on type
  compute_logL <- function(R, S, evo_par = NULL) {
    if (likelihood == "ML") {
      loglik_bastide(R, S, tree, trait_data,
                     model = effective_model, evo_model_par = evo_par)
    } else {
      # REML: need theta parameterization
      theta <- matrices_to_theta(R, S)
      restricted_loglik(theta, tree, trait_data, backend = "R",
                        model = effective_model,
                        evo_model_par_fixed = evo_par)
    }
  }

  # Run methods in sequence
  for (i in seq_along(method)) {
    m <- method[i]
    ctrl <- control[[m]]

    if (verbose) {
      cat("\n=== Running", m, "(", likelihood, ", step", i, "of", length(method), ") ===\n")
    }

    start_logL <- compute_logL(R, S, evo_model_par)

    result <- switch(m,
      "EM" = run_em_step(tree, trait_data, R, S, mu, ctrl, verbose),
      "BFGS" = run_bfgs_step(tree, trait_data, R, S, ctrl, verbose,
                             optim_method = "BFGS", likelihood = likelihood,
                             model = effective_model,
                             evo_model_par = evo_model_par,
                             optimize_evo = optimize_evo,
                             constraints = constraints),
      "Nelder-Mead" = run_nm_step(tree, trait_data, R, S, ctrl, verbose,
                                   likelihood = likelihood,
                                   model = effective_model,
                                   evo_model_par = evo_model_par,
                                   optimize_evo = optimize_evo,
                                   constraints = constraints),
      "L-BFGS-B" = run_bfgs_step(tree, trait_data, R, S, ctrl, verbose,
                                  optim_method = "L-BFGS-B", likelihood = likelihood,
                                  model = effective_model,
                                  evo_model_par = evo_model_par,
                                  optimize_evo = optimize_evo,
                                  constraints = constraints)
    )

    # Update parameters
    R <- result$R
    S <- result$S
    if (!is.null(result$mu)) mu <- result$mu
    if (!is.null(result$evo_model_par)) evo_model_par <- result$evo_model_par

    # Store in history
    history[[i]] <- list(
      method = m,
      start_logL = start_logL,
      end_logL = result$logL,
      converged = result$converged,
      iterations = result$iterations
    )

    if (verbose) {
      cat("  logL:", round(start_logL, 4), "->", round(result$logL, 4),
          "(", ifelse(result$converged, "converged", "not converged"), ")\n")
    }

    # Check if we've converged to the overall tolerance
    if (result$converged) {
      overall_converged <- TRUE
    }
  }

  list(
    R = R,
    S = S,
    mu = mu,
    logL = compute_logL(R, S, evo_model_par),
    likelihood = likelihood,
    model = model,
    evo_model_par = evo_model_par,
    history = history,
    converged = overall_converged
  )
}


#' Run EM Step
#'
#' @keywords internal
run_em_step <- function(tree, trait_data, R, S, mu, ctrl, verbose) {

  result <- em_fit(
    tree = tree,
    trait_data = trait_data,
    R_init = R,
    S_init = S,
    mu_init = mu,
    max_iter = ctrl$max_iter,
    tol = ctrl$tol,
    verbose = verbose
  )

  list(
    R = result$R,
    S = result$S,
    mu = result$mu,
    logL = result$logL,
    converged = result$converged,
    iterations = result$iterations
  )
}


#' Run BFGS Step with Analytic Gradients
#'
#' @keywords internal
run_bfgs_step <- function(tree, trait_data, R, S, ctrl, verbose,
                          optim_method = "BFGS", likelihood = "ML",
                          model = "BM", evo_model_par = NULL,
                          optimize_evo = FALSE, constraints = NULL) {

  p <- ncol(as.matrix(trait_data[, -1]))

  # Convert to theta (R and S only), using constraints if provided
  theta <- matrices_to_theta(R, S, constraints = constraints)

  # If optimizing evo param, append it to theta
  if (optimize_evo && !is.null(evo_model_par)) {
    evo_raw <- switch(model,
      "lambda" = qlogis(evo_model_par),
      "kappa"  = log(evo_model_par),
      "delta"  = log(evo_model_par),
      "EB"     = log(-evo_model_par),
      "OU"     = log(evo_model_par)
    )
    theta <- c(theta, evo_raw)
  }

  # Number of matrix params depends on constraints
  if (!is.null(constraints)) {
    n_mat_params <- constraints$n_free
  } else {
    n_mat_params <- p * (p + 1)  # number of R + S params (unconstrained)
  }

  # Helper to extract evo param from theta
  extract_evo <- function(theta) {
    if (!optimize_evo) return(evo_model_par)
    raw <- theta[length(theta)]
    switch(model,
      "lambda" = plogis(raw),
      "kappa"  = exp(raw),
      "delta"  = exp(raw),
      "EB"     = -exp(raw),
      "OU"     = exp(raw)
    )
  }

  # Negative log-likelihood function with singularity guard
  fn <- function(theta) {
    mat_theta <- theta[1:n_mat_params]
    mats <- tryCatch(
      theta_to_matrices(mat_theta, p, constraints = constraints),
      error = function(e) NULL
    )
    if (is.null(mats)) return(1e100)  # Penalty for invalid matrices

    cur_evo <- extract_evo(theta)

    logL <- tryCatch({
      if (likelihood == "ML") {
        loglik_bastide(mats$R, mats$S, tree, trait_data,
                       model = model, evo_model_par = cur_evo)
      } else {
        # REML with constraints: need to convert to unconstrained theta for REML
        reml_theta <- matrices_to_theta(mats$R, mats$S)
        restricted_loglik(reml_theta, tree, trait_data, backend = "R",
                          model = model, evo_model_par_fixed = cur_evo)
      }
    }, error = function(e) -1e100)
    if (!is.finite(logL) || logL < -1e99) return(1e100)
    -logL
  }

  # Gradient function (if using) with singularity guard
  # Use analytic gradients for ML (BM and non-BM); REML uses numerical
  gr <- NULL
  if (isTRUE(ctrl$use_gradient) && likelihood == "ML") {
    gr <- function(theta) {
      mat_theta <- theta[1:n_mat_params]
      mats <- tryCatch(
        theta_to_matrices(mat_theta, p, constraints = constraints),
        error = function(e) NULL
      )
      if (is.null(mats)) return(rep(0, length(theta)))  # Zero gradient for invalid

      cur_evo <- extract_evo(theta)

      grad <- tryCatch({
        g <- bastide_gradient(tree, trait_data, mats$R, mats$S,
                              model = model, evo_model_par = cur_evo)
        -gradient_to_theta(g$grad_R, g$grad_S, mats$R, mats$S,
                           grad_evo = g$grad_evo,
                           evo_model_par = cur_evo,
                           model = model,
                           constraints = constraints)
      }, error = function(e) {
        rep(0, length(theta))  # Zero gradient on error
      })

      # Replace any NA/Inf with 0
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
      method = optim_method,
      control = list(
        maxit = ctrl$max_iter,
        reltol = ctrl$tol
      )
    )
  }, error = function(e) {
    if (verbose) cat("  BFGS failed:", e$message, "\n")
    # Try to compute fn(theta) safely
    init_val <- tryCatch(fn(theta), error = function(e) 1e100)
    list(par = theta, value = init_val, convergence = 1)
  })

  # Convert back safely
  mat_theta <- opt_result$par[1:n_mat_params]
  mats <- tryCatch(
    theta_to_matrices(mat_theta, p, constraints = constraints),
    error = function(e) list(R = R, S = S)
  )

  final_evo <- extract_evo(opt_result$par)

  list(
    R = mats$R,
    S = mats$S,
    mu = NULL,  # BFGS doesn't update mu directly
    evo_model_par = final_evo,
    logL = -opt_result$value,
    converged = (opt_result$convergence == 0),
    iterations = NA  # optim doesn't report this directly
  )
}


#' Run Nelder-Mead Step
#'
#' @keywords internal
run_nm_step <- function(tree, trait_data, R, S, ctrl, verbose, likelihood = "ML",
                        model = "BM", evo_model_par = NULL, optimize_evo = FALSE,
                        constraints = NULL) {
  # Use BFGS function with Nelder-Mead method
  run_bfgs_step(tree, trait_data, R, S,
                list(max_iter = ctrl$max_iter, tol = ctrl$tol, use_gradient = FALSE),
                verbose, optim_method = "Nelder-Mead", likelihood = likelihood,
                model = model, evo_model_par = evo_model_par,
                optimize_evo = optimize_evo, constraints = constraints)
}


# NOTE: matrices_to_theta, theta_to_matrices, and gradient_to_theta
# are now defined in utils.R as the single source of truth for
# theta parameterization. optimizer.R calls them from there.
