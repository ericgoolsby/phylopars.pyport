#' Backend-aware log-likelihood computation
#'
#' Computes restricted log-likelihood using the specified backend.
#'
#' @param theta Parameter vector
#' @param tree A phylo object
#' @param Y Data frame with species in first column, traits in remaining columns
#' @param backend Which backend to use: "R", "armadillo", or "eigen"
#' @param Sigma_diag Logical, if TRUE Sigma is diagonal
#' @param B_diag Logical, if TRUE B is diagonal
#' @param Sigma_fixed Fixed Sigma matrix
#' @param B_fixed Fixed B matrix
#' @param precompute_list List of precomputed values
#' @param model Evolutionary model
#' @param evo_model_par_fixed Fixed evolutionary parameter
#'
#' @return Log-likelihood value
#' @export
restricted_loglik <- function(theta, tree, Y,
                               backend = c("R", "armadillo", "eigen"),
                               Sigma_diag = FALSE, B_diag = FALSE,
                               Sigma_fixed = NULL, B_fixed = NULL,
                               precompute_list = NULL,
                               model = "BM", evo_model_par_fixed = NULL) {
  
  backend <- match.arg(backend)
  
  if (backend == "R") {
    return(restricted_log_likelihood(
      theta = theta, tree = tree, Y = Y, ret_logL = TRUE,
      Sigma_diag = Sigma_diag, B_diag = B_diag,
      Sigma_fixed = Sigma_fixed, B_fixed = B_fixed,
      precompute_list = precompute_list,
      model = model, evo_model_par_fixed = evo_model_par_fixed
    ))
  }
  
  # For C++ backends, we need to preprocess data
  if (is.null(precompute_list)) {
    traits <- Y[, -1, drop = FALSE]
    m <- ncol(traits)
    species <- Y[, 1]
    nind <- nrow(Y)
    Nnode <- tree$Nnode + length(tree$tip.label)
    n_obs <- sum(!is.na(traits))
  } else {
    traits <- precompute_list$traits
    m <- precompute_list$m
    species <- precompute_list$species
    nind <- precompute_list$nind
    Nnode <- precompute_list$Nnode
    n_obs <- precompute_list$n_obs
  }
  
  # Convert theta to matrices
  pars <- theta2vals(
    theta = theta, m = m, Sigma_diag = Sigma_diag, B_diag = B_diag,
    Sigma_fixed = Sigma_fixed, B_fixed = B_fixed,
    model = model, evo_model_par_fixed = evo_model_par_fixed
  )
  Sigma <- pars$Sigma
  B <- pars$B
  evo_model_par <- pars$evo_model_par
  
  # Transform tree for non-BM models
  working_tree <- tree
  if (model != "BM" && is.null(evo_model_par_fixed) && !is.null(evo_model_par)) {
    working_tree <- dist_from_root(tree)
    working_tree <- transf_branch_lengths(working_tree, model, evo_model_par)
  }
  
  # Create species index (which tip node each observation belongs to)
  species_idx <- match(species, working_tree$tip.label)
  
  # Convert traits to matrix
  traits_mat <- as.matrix(traits)
  
  # Call appropriate C++ function
  if (backend == "armadillo") {
    return(restricted_loglik_arma(
      Sigma = Sigma,
      B = B,
      edge = working_tree$edge,
      edge_length = working_tree$edge.length,
      traits = traits_mat,
      species_idx = as.integer(species_idx),
      n_obs = as.integer(n_obs),
      Nnode = as.integer(Nnode)
    ))
  } else if (backend == "eigen") {
    return(restricted_loglik_eigen(
      Sigma = Sigma,
      B = B,
      edge = working_tree$edge,
      edge_length = working_tree$edge.length,
      traits = traits_mat,
      species_idx = as.integer(species_idx),
      n_obs = as.integer(n_obs),
      Nnode = as.integer(Nnode)
    ))
  }
}

#' Benchmark backends
#'
#' Compare performance of different backends on the same data.
#'
#' @param tree A phylo object
#' @param Y Data frame with traits
#' @param n_reps Number of repetitions for timing
#'
#' @return Data frame with timing results
#' @export
benchmark_backends <- function(tree, Y, n_reps = 10) {
  
  # Get initial parameters
  m <- ncol(Y) - 1
  n_params <- m * (m + 1)  # Sigma and B lower triangles
  theta <- rep(0, n_params)
  
  backends <- c("R", "armadillo", "eigen")
  results <- data.frame(
    backend = backends,
    mean_time_ms = NA_real_,
    logL = NA_real_
  )
  
  for (i in seq_along(backends)) {
    backend <- backends[i]
    
    # Warm up
    tryCatch({
      logL <- restricted_loglik(theta, tree, Y, backend = backend)
      results$logL[i] <- logL
      
      # Time it
      times <- numeric(n_reps)
      for (j in seq_len(n_reps)) {
        start <- Sys.time()
        restricted_loglik(theta, tree, Y, backend = backend)
        times[j] <- as.numeric(Sys.time() - start) * 1000
      }
      results$mean_time_ms[i] <- mean(times)
    }, error = function(e) {
      message("Backend '", backend, "' failed: ", e$message)
    })
  }
  
  results
}
