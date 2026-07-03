#' Fit a phylogenetic mixed model using REML
#'
#' Estimates phylogenetic (Sigma) and residual (B) covariance matrices using
#' restricted maximum likelihood. Handles missing data and within-species
#' variation (multiple observations per species).
#'
#' @param Y Data frame with species names in first column, trait values in
#'   remaining columns. Missing values (NA) are allowed.
#' @param tree A phylo object (from ape package)
#' @param Sigma_start Starting value for Sigma matrix (default: identity)
#' @param B_start Starting value for B matrix (default: identity)
#' @param Sigma_fixed Fixed value for Sigma (not optimized if provided)
#' @param B_fixed Fixed value for B (not optimized if provided)
#' @param Sigma_diag If TRUE, constrain Sigma to be diagonal
#' @param B_diag If TRUE, constrain B to be diagonal
#' @param method Optimization method for optim() (default: "BFGS")
#' @param control List of control parameters for optim()
#' @param model Evolutionary model: "BM", "lambda", "kappa", "delta", "EB", or "OU"
#' @param evo_model_par_start Starting value for evolutionary model parameter
#' @param evo_model_par_fixed Fixed value for evolutionary model parameter
#'
#' @return List containing:
#'   \item{mu}{Estimated mean vector}
#'   \item{Sigma}{Estimated phylogenetic covariance matrix}
#'   \item{B}{Estimated residual covariance matrix}
#'   \item{evo_model_par}{Estimated evolutionary model parameter (if applicable)}
#'   \item{logL}{Restricted log-likelihood at optimum}
#'   \item{AIC}{Akaike Information Criterion}
#'   \item{convergence}{Convergence information from optim()}
#'   \item{model}{Evolutionary model used}
#'   \item{npars}{Number of parameters estimated}
#'
#' @examples
#' \dontrun{
#' library(ape)
#'
#' # Simulate a tree
#' tree <- rtree(20)
#'
#' # Create trait data with some missing values
#' Y <- data.frame(
#'   species = tree$tip.label,
#'   trait1 = rnorm(20),
#'   trait2 = rnorm(20)
#' )
#' Y$trait1[c(3, 7)] <- NA
#'
#' # Fit model
#' result <- phylopars(Y, tree)
#' print(result$Sigma)
#' print(result$B)
#' }
#'
#' @export
phylopars <- function(Y, tree, Sigma_start = NULL, B_start = NULL,
                      Sigma_fixed = NULL, B_fixed = NULL,
                      Sigma_diag = FALSE, B_diag = FALSE,
                      method = "BFGS", control = list(),
                      model = "BM", evo_model_par_start = NULL,
                      evo_model_par_fixed = NULL, backend = c("armadillo", "R", "eigen")) {

  # Validate inputs
  if (!inherits(tree, "phylo")) {
    stop("tree must be a phylo object")
  }
  backend <- match.arg(backend)

  if (!is.data.frame(Y)) {
    stop("Y must be a data frame")
  }

  # Validate model
  valid_models <- c("BM", "lambda", "kappa", "delta", "EB", "OU")
  if (!(model %in% valid_models)) {
    stop("Unknown model: '", model, "'. Valid models: ",
         paste(valid_models, collapse = ", "))
  }

  # Handle evolutionary models
  if (model != "BM") {
    tree <- dist_from_root(tree)

    # Validate fixed parameter if provided
    if (!is.null(evo_model_par_fixed)) {
      validate_evo_par(model, evo_model_par_fixed, "evo_model_par_fixed")
      # Transform branch lengths using fixed parameter
      tree <- transf_branch_lengths(tree, model, evo_model_par_fixed)
    }

    # Validate starting parameter if provided
    if (!is.null(evo_model_par_start)) {
      validate_evo_par(model, evo_model_par_start, "evo_model_par_start")
    }
  }

  # Precompute values for efficiency
  traits <- Y[, -1, drop = FALSE]
  m <- ncol(traits)
  species <- Y[, 1]
  observations <- split(seq_len(nrow(Y)), species)
  nind <- nrow(Y)
  Nnode <- tree$Nnode + length(tree$tip.label)
  n_obs <- sum(!is.na(traits))

  precompute_list <- list(
    observations = observations,
    nind = nind,
    Nnode = Nnode,
    species = species,
    traits = traits,
    m = m,
    n_obs = n_obs,
    original_tree = if (model != "BM" && is.null(evo_model_par_fixed)) tree else NULL
  )

  # Get starting parameter values
  theta0 <- get_default_starting_theta(
    m = m, Sigma_diag = Sigma_diag, B_diag = B_diag,
    Sigma_fixed = Sigma_fixed, B_fixed = B_fixed,
    Sigma_start = Sigma_start, B_start = B_start,
    model = model, evo_model_par_start = evo_model_par_start,
    evo_model_par_fixed = evo_model_par_fixed
  )

  ntheta <- length(theta0)

  if (ntheta > 0) {
    # Define negative log-likelihood function
    nlogl <- function(theta) {
      -restricted_loglik(
        theta = theta, tree = tree, Y = Y, backend = backend,
        Sigma_diag = Sigma_diag, B_diag = B_diag,
        Sigma_fixed = Sigma_fixed, B_fixed = B_fixed,
        precompute_list = precompute_list,
        model = model, evo_model_par_fixed = evo_model_par_fixed
      )
    }

    # Optimize
    res <- optim(theta0, nlogl, method = method, control = control)
    opt_theta <- res$par

  } else {
    # All parameters fixed, no optimization needed
    opt_theta <- vals2theta(
      Sigma = Sigma_fixed, B = B_fixed, evo_model_par = evo_model_par_fixed,
      ntheta = ntheta, Sigma_diag = Sigma_diag, B_diag = B_diag,
      Sigma_fixed = Sigma_fixed, B_fixed = B_fixed,
      model = model, evo_model_par_fixed = evo_model_par_fixed
    )
    res <- list(convergence = 0, message = "No optimization (all parameters fixed)")
  }

  # Extract final matrices
  pars <- theta2vals(
    theta = opt_theta, m = m, Sigma_diag = Sigma_diag, B_diag = B_diag,
    Sigma_fixed = Sigma_fixed, B_fixed = B_fixed,
    model = model, evo_model_par_fixed = evo_model_par_fixed
  )
  Sigma <- pars$Sigma
  B <- pars$B
  evo_model_par <- pars$evo_model_par

  # Get final results
  result_dict <- restricted_log_likelihood(
    theta = opt_theta, tree = tree, Y = Y,
    ret_logL = FALSE, Sigma_diag = Sigma_diag, B_diag = B_diag,
    Sigma_fixed = Sigma_fixed, B_fixed = B_fixed,
    precompute_list = precompute_list,
    model = model, evo_model_par_fixed = evo_model_par_fixed
  )

  # Compute AIC: AIC = 2k - 2ln(L)
  k <- ntheta
  AIC <- 2 * k - 2 * result_dict$logL

  # Add row/column names to matrices
  trait_names <- names(traits)
  rownames(Sigma) <- colnames(Sigma) <- trait_names
  rownames(B) <- colnames(B) <- trait_names
  rownames(result_dict$mu) <- trait_names

  # Return results
  list(
    mu = result_dict$mu,
    Sigma = Sigma,
    B = B,
    evo_model_par = evo_model_par,
    logL = result_dict$logL,
    AIC = AIC,
    convergence = list(
      code = if (exists("res") && !is.null(res$convergence)) res$convergence else 0,
      message = if (exists("res") && !is.null(res$message)) res$message else NULL,
      counts = if (exists("res") && !is.null(res$counts)) res$counts else NULL
    ),
    model = model,
    npars = k
  )
}


#' Print method for phylopars results
#'
#' @param x Result from phylopars()
#' @param ... Additional arguments (ignored)
#'
#' @export
print.phylopars <- function(x, ...) {
  cat("Phylogenetic Mixed Model Results\n")
  cat("================================\n")
  cat("Model:", x$model, "\n")
  cat("Log-likelihood:", round(x$logL, 4), "\n")
  cat("AIC:", round(x$AIC, 4), "\n")
  cat("Parameters estimated:", x$npars, "\n")

  if (!is.null(x$evo_model_par)) {
    cat("\nEvolutionary parameter:", round(x$evo_model_par, 4), "\n")
  }

  cat("\nMean (mu):\n")
  print(round(x$mu, 4))

  cat("\nPhylogenetic covariance (Sigma):\n")
  print(round(x$Sigma, 4))

  cat("\nResidual covariance (B):\n")
  print(round(x$B, 4))

  invisible(x)
}
