#' Restricted log-likelihood for phylogenetic mixed model
#'
#' Calculates the REML log-likelihood using a postorder tree traversal algorithm.
#' This is O(n) in the number of species, handling missing data and within-species
#' variation.
#'
#' @param theta Parameter vector
#' @param tree A phylo object
#' @param Y Data frame with species in first column, traits in remaining columns
#' @param ret_logL If TRUE, return only the log-likelihood value
#' @param Sigma_diag Logical, if TRUE Sigma is diagonal
#' @param B_diag Logical, if TRUE B is diagonal
#' @param Sigma_fixed Fixed Sigma matrix
#' @param B_fixed Fixed B matrix
#' @param precompute_list List of precomputed values for efficiency
#' @param model Evolutionary model
#' @param evo_model_par_fixed Fixed evolutionary parameter
#'
#' @return Log-likelihood value or list with additional information
#' @export
restricted_log_likelihood <- function(theta, tree, Y, ret_logL = TRUE,
                                      Sigma_diag = FALSE, B_diag = FALSE,
                                      Sigma_fixed = NULL, B_fixed = NULL,
                                      precompute_list = NULL,
                                      model = "BM", evo_model_par_fixed = NULL) {
  
  # Extract or compute precomputed values
  if (is.null(precompute_list)) {
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
  } else {
    observations <- precompute_list$observations
    nind <- precompute_list$nind
    Nnode <- precompute_list$Nnode
    species <- precompute_list$species
    traits <- precompute_list$traits
    m <- precompute_list$m
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
  
  # Transform tree for non-BM models during optimization
  working_tree <- tree
  if (model != "BM" && is.null(evo_model_par_fixed) && !is.null(evo_model_par)) {
    if (!is.null(precompute_list$original_tree)) {
      working_tree <- precompute_list$original_tree
      working_tree <- dist_from_root(working_tree)
      working_tree <- transf_branch_lengths(working_tree, model, evo_model_par)
    } else {
      working_tree <- dist_from_root(tree)
      working_tree <- transf_branch_lengths(working_tree, model, evo_model_par)
    }
  }
  
  # Initialize per-individual storage
  p_ind <- lapply(seq_len(nind), function(i) matrix(0, m, m))
  Vr_ind <- lapply(seq_len(nind), function(i) matrix(0, m, 1))
  Q_ind <- lapply(seq_len(nind), function(i) matrix(0, 1, 1))
  logW_ind <- numeric(nind)
  
  # Initialize per-node storage
  p <- lapply(seq_len(Nnode), function(i) matrix(0, m, m))
  Vr <- lapply(seq_len(Nnode), function(i) matrix(0, m, 1))
  Q <- lapply(seq_len(Nnode), function(i) matrix(0, 1, 1))
  logW <- numeric(Nnode)
  
  # Process observations for each species
  for (sp in names(observations)) {
    sp_node <- which(working_tree$tip.label == sp)
    obs_idx <- observations[[sp]]
    
    for (idx in obs_idx) {
      obs_traits <- which(!is.na(traits[idx, ]))
      if (length(obs_traits) == 0) next
      
      Yvec <- matrix(as.numeric(traits[idx, obs_traits]), ncol = 1)
      Ba <- B[obs_traits, obs_traits, drop = FALSE]
      # Return bad likelihood if singular
      Ba_inv <- tryCatch(
        solve(Ba),
        error = function(e) NULL
      )
      if (is.null(Ba_inv)) {
        if (ret_logL) return(-1e100)
        return(list(logL = -1e100))
      }

      p_ind[[idx]][obs_traits, obs_traits] <- Ba_inv
      Vr_ind[[idx]][obs_traits, 1] <- Ba_inv %*% Yvec
      Q_ind[[idx]] <- t(Yvec) %*% Ba_inv %*% Yvec
      logW_ind[idx] <- determinant(Ba, logarithm = TRUE)$modulus[1]
      
      p[[sp_node]] <- p[[sp_node]] + p_ind[[idx]]
      Vr[[sp_node]] <- Vr[[sp_node]] + Vr_ind[[idx]]
      Q[[sp_node]] <- Q[[sp_node]] + Q_ind[[idx]]
      logW[sp_node] <- logW[sp_node] + logW_ind[idx]
    }
  }
  
  # Postorder traversal (tips to root)
  for (i in seq_len(nrow(working_tree$edge))) {
    anc <- working_tree$edge[i, 1]
    des <- working_tree$edge[i, 2]
    edge_length <- working_tree$edge.length[i]
    T_edge <- edge_length * Sigma
    
    pA <- p[[des]]
    itpa <- diag(m) + T_edge %*% pA
    # Return bad likelihood if singular (tells optimizer to avoid this region)
    itpainv <- tryCatch(
      solve(itpa),
      error = function(e) {
        if (ret_logL) return(-1e100)
        return(NULL)
      }
    )
    if (is.null(itpainv) || (is.numeric(itpainv) && !is.matrix(itpainv) && length(itpainv) == 1)) {
      if (ret_logL) return(-1e100)
      return(list(logL = -1e100))
    }
    p[[des]] <- pA %*% itpainv
    logW[des] <- logW[des] + determinant(itpa, logarithm = TRUE)$modulus[1]
    Q[[des]] <- Q[[des]] - t(Vr[[des]]) %*% itpainv %*% T_edge %*% Vr[[des]]
    Vr[[des]] <- t(t(Vr[[des]]) %*% itpainv)
    
    p[[anc]] <- p[[anc]] + p[[des]]
    Vr[[anc]] <- Vr[[anc]] + Vr[[des]]
    Q[[anc]] <- Q[[anc]] + Q[[des]]
    logW[anc] <- logW[anc] + logW[des]
  }
  
  # Compute mean and log-likelihood at root
  # Note: anc is now the root node from the last edge processed
  root <- anc
  mu <- tryCatch(
    solve(p[[root]], Vr[[root]]),
    error = function(e) NULL
  )
  if (is.null(mu)) {
    if (ret_logL) return(-1e100)
    return(list(logL = -1e100))
  }

  logL <- -0.5 * (
    (n_obs - m) * log(2 * pi) +
      determinant(p[[root]], logarithm = TRUE)$modulus[1] +
      logW[root] +
      Q[[root]] -
      2 * t(mu) %*% Vr[[root]] +
      t(mu) %*% p[[root]] %*% mu
  )
  
  if (ret_logL) {
    return(as.numeric(logL))
  }
  
  list(
    logL = as.numeric(logL),
    mu = mu,
    logW = logW[root],
    Vr = Vr[[root]],
    Q = Q[[root]],
    p = p[[root]]
  )
}
