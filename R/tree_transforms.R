#' Tree Branch Length Transformations for Evolutionary Models
#'
#' Functions that transform phylogenetic tree branch lengths to implement
#' different models of trait evolution. Each model modifies the tree so
#' that BM on the transformed tree is equivalent to the specified model
#' on the original tree.
#'
#' @name tree-transforms
NULL


#' Transform Tree Branch Lengths for Evolutionary Model
#'
#' Dispatch function that applies the appropriate branch length transformation
#' based on the model name. The tree must have root distances precomputed
#' (via \code{\link{dist_from_root}}) for models that require them
#' (delta, EB, OU).
#'
#' @param tree A phylo object. Should have \code{root_dist} computed for
#'   delta, EB, and OU models.
#' @param model Character string: "BM", "lambda", "kappa", "delta", "EB", or "OU"
#' @param par Model parameter value
#'
#' @return Tree with transformed branch lengths
#'
#' @details
#' Parameter constraints:
#' \itemize{
#'   \item \strong{lambda}: par in [0, 1]. Scales internal branches; tips adjusted
#'     to preserve root-to-tip distance. lambda=1 is BM, lambda=0 is star.
#'   \item \strong{kappa}: par > 0. Raises each branch length to the power kappa.
#'     kappa=1 is BM. kappa<1 = punctuated, kappa>1 = gradualism.
#'   \item \strong{delta}: par > 0. Raises node heights to power delta, then
#'     recomputes branch lengths. delta=1 is BM. delta<1 = early change,
#'     delta>1 = late change.
#'   \item \strong{EB}: par <= 0. Early Burst / ACDC model. Multiplies branches by
#'     exp(rate * t_from_root). par=0 is BM, par<0 = declining rates.
#'   \item \strong{OU}: par >= 0. Ornstein-Uhlenbeck via Pagel's lambda
#'     approximation. Requires ultrametric tree. par=0 is BM.
#' }
#'
#' @export
transform_tree <- function(tree, model, par) {
  switch(model,
    "BM" = tree,
    "lambda" = transform_lambda(tree, par),
    "kappa" = transform_kappa(tree, par),
    "delta" = transform_delta(tree, par),
    "EB" = transform_eb(tree, par),
    "OU" = transform_ou(tree, par),
    stop("Unknown model: '", model, "'. Supported: BM, lambda, kappa, delta, EB, OU")
  )
}


#' Pagel's Lambda Transform
#'
#' Multiplies all internal branch lengths by lambda. For tip edges, adds
#' back the lost height to preserve root-to-tip distances. This effectively
#' controls the strength of phylogenetic signal.
#'
#' @param tree A phylo object with \code{root_dist} precomputed
#' @param lambda Numeric in [0, 1]
#'
#' @return Tree with transformed branch lengths
#'
#' @details
#' For each edge with descendant node d:
#' \itemize{
#'   \item Internal edges: new_length = lambda * old_length
#'   \item Tip edges: new_length = lambda * old_length + (1 - lambda) * root_dist[d]
#' }
#' At lambda=1, the tree is unchanged (pure BM).
#' At lambda=0, all tips connect directly to root (star tree).
#'
#' @keywords internal
transform_lambda <- function(tree, lambda) {
  if (lambda < 0 || lambda > 1) {
    stop("Lambda must be between 0 and 1, got ", lambda)
  }
  if (lambda == 1) return(tree)

  edge_length <- tree$edge.length
  des <- tree$edge[, 2]
  nspecies <- length(tree$tip.label)

  # Identify tip edges
  is_tip <- des <= nspecies

  # Multiply all branches by lambda
  edge_length <- edge_length * lambda

  # Add back root-to-tip distance for tip edges
  edge_length[is_tip] <- edge_length[is_tip] +
    (1 - lambda) * tree$root_dist[des[is_tip]]

  tree$edge.length <- edge_length
  tree
}


#' Kappa Transform
#'
#' Raises each branch length to the power kappa. This models punctuational
#' (kappa < 1) versus gradualist (kappa > 1) evolution.
#'
#' @param tree A phylo object
#' @param kappa Numeric > 0
#'
#' @return Tree with transformed branch lengths
#'
#' @details
#' For each edge: new_length = old_length^kappa
#'
#' At kappa=1, the tree is unchanged (BM).
#' At kappa=0, all branches have length 1 (pure punctuational model).
#'
#' @keywords internal
transform_kappa <- function(tree, kappa) {
  if (kappa <= 0) {
    stop("Kappa must be positive, got ", kappa)
  }
  if (kappa == 1) return(tree)

  tree$edge.length <- tree$edge.length^kappa
  tree
}


#' Delta Transform
#'
#' Raises node heights (root-to-node distances) to the power delta, then
#' recomputes branch lengths as differences in transformed heights. This models
#' whether evolution is concentrated early (delta < 1) or late (delta > 1)
#' in the tree.
#'
#' @param tree A phylo object with \code{root_dist} precomputed
#' @param delta Numeric > 0
#'
#' @return Tree with transformed branch lengths
#'
#' @details
#' For each node i: new_height[i] = old_height[i]^delta
#' For each edge (anc, des): new_length = new_height[des] - new_height[anc]
#'
#' The root has height 0, so root height is always 0 regardless of delta.
#' At delta=1, the tree is unchanged (BM).
#'
#' @keywords internal
transform_delta <- function(tree, delta) {
  if (delta <= 0) {
    stop("Delta must be positive, got ", delta)
  }
  if (delta == 1) return(tree)

  # Transform node heights
  new_heights <- tree$root_dist^delta

  # Recompute branch lengths from transformed heights
  anc <- tree$edge[, 1]
  des <- tree$edge[, 2]
  tree$edge.length <- new_heights[des] - new_heights[anc]

  # Guard against negative branch lengths from numerical issues
  tree$edge.length[tree$edge.length < 0] <- 0

  tree
}


#' Early Burst (EB / ACDC) Transform
#'
#' Multiplies each branch length by an exponential function of time from the
#' root, modeling declining (or increasing) rates of evolution. Also called
#' the ACDC (Accelerating-Decelerating) model.
#'
#' @param tree A phylo object with \code{root_dist} precomputed
#' @param rate Numeric <= 0 for Early Burst (rate > 0 would be Late Burst)
#'
#' @return Tree with transformed branch lengths
#'
#' @details
#' For each edge from ancestor (anc) to descendant (des):
#' new_length = (exp(rate * t_des) - exp(rate * t_anc)) / rate
#'
#' where t_anc and t_des are distances from root.
#'
#' At rate=0, the tree is unchanged (BM).
#' At rate<0, branches near the root are longer (early burst).
#' At rate>0, branches near the tips are longer (late burst).
#'
#' The implementation uses geiger's formula: integral of exp(rate * t) dt
#' over the branch interval [t_anc, t_des].
#'
#' @keywords internal
transform_eb <- function(tree, rate) {
  if (rate == 0) return(tree)

  anc <- tree$edge[, 1]
  des <- tree$edge[, 2]
  t_anc <- tree$root_dist[anc]
  t_des <- tree$root_dist[des]

  # Integral of exp(rate * t) from t_anc to t_des, divided by original length
  # = (exp(rate * t_des) - exp(rate * t_anc)) / rate
  tree$edge.length <- (exp(rate * t_des) - exp(rate * t_anc)) / rate

  # Guard against negative branch lengths
  tree$edge.length[tree$edge.length < 0] <- 0

  tree
}


#' Ornstein-Uhlenbeck Transform
#'
#' Transforms branch lengths using the OU-to-lambda approximation for
#' ultrametric trees. For ultrametric trees, OU with selection strength alpha
#' is equivalent to Pagel's lambda with lambda = exp(-2 * alpha * T),
#' where T is the tree height, combined with a variance rescaling.
#'
#' This is the "OUfixedRoot" parameterization used by geiger.
#'
#' @param tree A phylo object with \code{root_dist} precomputed
#' @param alpha Numeric >= 0 (selection strength)
#'
#' @return Tree with transformed branch lengths
#'
#' @details
#' For ultrametric trees, the OU covariance between species i and j is:
#'
#'   Cov(i,j) = sigma^2 / (2*alpha) * exp(-2*alpha*d_ij/2)
#'
#' This can be represented as BM on a tree with transformed branch lengths.
#' The transform uses the "OUfixedRoot" parameterization from geiger:
#'
#' For each edge from time t_anc to t_des (root-to-node distances):
#'   new_length = (1/(2*alpha)) * (exp(2*alpha*t_des) - exp(2*alpha*t_anc))
#'
#' The resulting tree, when used with BM, produces the OU covariance structure.
#'
#' At alpha=0, the tree is unchanged (BM limit).
#'
#' @keywords internal
transform_ou <- function(tree, alpha) {
  if (alpha < 0) {
    stop("OU alpha must be non-negative, got ", alpha)
  }
  if (alpha == 0) return(tree)

  # Warn for non-ultrametric trees -- the OUfixedRoot transform is only
  # an exact representation of OU on ultrametric trees
  if (requireNamespace("ape", quietly = TRUE)) {
    if (!ape::is.ultrametric(tree, tol = 0.01)) {
      warning("OU tree transform uses the OUfixedRoot approximation, which is ",
              "only exact for ultrametric trees. Results on non-ultrametric trees ",
              "may be inaccurate.", call. = FALSE)
    }
  }

  anc <- tree$edge[, 1]
  des <- tree$edge[, 2]
  t_anc <- tree$root_dist[anc]
  t_des <- tree$root_dist[des]

  # OUfixedRoot transform from geiger
  # new_length = (1/(2*alpha)) * (exp(2*alpha*t_des) - exp(2*alpha*t_anc))
  two_alpha <- 2 * alpha
  tree$edge.length <- (exp(two_alpha * t_des) - exp(two_alpha * t_anc)) / two_alpha

  # Guard against negative branch lengths
  tree$edge.length[tree$edge.length < 0] <- 0

  tree
}


#' Validate Evolutionary Model Parameter
#'
#' Checks that a model parameter is within valid bounds.
#'
#' @param model Character string: model name
#' @param par Numeric parameter value
#' @param par_name Character string: name of parameter (for error messages)
#'
#' @keywords internal
validate_evo_par <- function(model, par, par_name = "parameter") {
  switch(model,
    "lambda" = {
      if (par < 0 || par > 1)
        stop(par_name, ": lambda must be in [0, 1], got ", par)
    },
    "kappa" = {
      if (par <= 0)
        stop(par_name, ": kappa must be positive, got ", par)
    },
    "delta" = {
      if (par <= 0)
        stop(par_name, ": delta must be positive, got ", par)
    },
    "EB" = {
      if (par > 0)
        stop(par_name, ": EB rate must be <= 0, got ", par)
    },
    "OU" = {
      if (par < 0)
        stop(par_name, ": OU alpha must be >= 0, got ", par)
    }
  )
  invisible(TRUE)
}


#' Jacobian of Branch Length Transform w.r.t. Evolutionary Parameter
#'
#' Computes dt'_k / d(evo_param) for each edge in the tree.
#' Used in the chain rule for analytic gradients of the log-likelihood
#' with respect to evolutionary model parameters.
#'
#' @param tree A phylo object with \code{root_dist} precomputed
#'   (via \code{\link{dist_from_root}})
#' @param model Character string: "lambda", "kappa", "delta", "EB", or "OU"
#' @param par Current parameter value
#'
#' @return Numeric vector of length \code{nrow(tree$edge)}, where each element
#'   is the derivative of the transformed branch length with respect to the
#'   evolutionary parameter.
#'
#' @details
#' The derivatives for each model:
#' \describe{
#'   \item{lambda}{Internal: dt'/dlambda = t_orig. Tip: dt'/dlambda = t_orig - root_dist[tip]}
#'   \item{kappa}{dt'/dkappa = t_orig^kappa * log(t_orig) (per-branch, since t' = t^kappa)}
#'   \item{delta}{dt'/ddelta = d_des^delta * log(d_des) - d_anc^delta * log(d_anc)}
#'   \item{EB}{dt'/dr = [d_des*exp(r*d_des) - d_anc*exp(r*d_anc)]/r - [exp(r*d_des) - exp(r*d_anc)]/r^2}
#'   \item{OU}{dt'/dalpha = [2*d_des*exp(2a*d_des) - 2*d_anc*exp(2a*d_anc)]/(2a) - [exp(2a*d_des) - exp(2a*d_anc)]/(2a)^2}
#' }
#'
#' @keywords internal
transform_jacobian <- function(tree, model, par) {
  n_edges <- nrow(tree$edge)
  jac <- numeric(n_edges)

  switch(model,
    "lambda" = {
      des <- tree$edge[, 2]
      nspecies <- length(tree$tip.label)
      is_tip <- des <= nspecies

      # Internal edges: dt'/dlambda = t_orig (the original branch length)
      # But we need the ORIGINAL branch length before transform.
      # The original branch length is stored in tree$edge.length only if
      # the tree hasn't been transformed yet. Since we receive the original
      # tree (pre-transform), we use tree$edge.length directly.
      jac <- tree$edge.length

      # Tip edges: dt'/dlambda = t_orig - root_dist[tip]
      # where root_dist is the original root-to-tip distance
      jac[is_tip] <- tree$edge.length[is_tip] - tree$root_dist[des[is_tip]]
    },

    "kappa" = {
      # t' = t^kappa
      # dt'/dkappa = t^kappa * log(t)
      t_orig <- tree$edge.length
      # Handle zero-length branches: 0^kappa * log(0) = 0 (by convention)
      safe_t <- ifelse(t_orig > 0, t_orig, 1)  # avoid log(0)
      jac <- ifelse(t_orig > 0, t_orig^par * log(safe_t), 0)
    },

    "delta" = {
      # t' = d_des^delta - d_anc^delta
      # dt'/ddelta = d_des^delta * log(d_des) - d_anc^delta * log(d_anc)
      anc <- tree$edge[, 1]
      des <- tree$edge[, 2]
      d_anc <- tree$root_dist[anc]
      d_des <- tree$root_dist[des]

      # Handle d=0 (root has distance 0): 0^delta * log(0) = 0
      safe_anc <- ifelse(d_anc > 0, d_anc, 1)
      safe_des <- ifelse(d_des > 0, d_des, 1)
      jac <- ifelse(d_des > 0, d_des^par * log(safe_des), 0) -
             ifelse(d_anc > 0, d_anc^par * log(safe_anc), 0)
    },

    "EB" = {
      # t' = (exp(r*d_des) - exp(r*d_anc)) / r
      # dt'/dr = [d_des*exp(r*d_des) - d_anc*exp(r*d_anc)]/r
      #        - [exp(r*d_des) - exp(r*d_anc)]/r^2
      anc <- tree$edge[, 1]
      des <- tree$edge[, 2]
      d_anc <- tree$root_dist[anc]
      d_des <- tree$root_dist[des]
      r <- par

      if (abs(r) < 1e-10) {
        # Limit as r -> 0: dt'/dr = (d_des^2 - d_anc^2) / 2
        jac <- (d_des^2 - d_anc^2) / 2
      } else {
        exp_r_des <- exp(r * d_des)
        exp_r_anc <- exp(r * d_anc)
        jac <- (d_des * exp_r_des - d_anc * exp_r_anc) / r -
               (exp_r_des - exp_r_anc) / r^2
      }
    },

    "OU" = {
      # t' = (exp(2a*d_des) - exp(2a*d_anc)) / (2a)
      # Let u = 2a, so t' = (exp(u*d_des) - exp(u*d_anc)) / u
      # dt'/dalpha = (dt'/du) * (du/dalpha) = 2 * dt'/du
      # dt'/du = [d_des*exp(u*d_des) - d_anc*exp(u*d_anc)]/u
      #        - [exp(u*d_des) - exp(u*d_anc)]/u^2
      anc <- tree$edge[, 1]
      des <- tree$edge[, 2]
      d_anc <- tree$root_dist[anc]
      d_des <- tree$root_dist[des]
      alpha <- par

      if (abs(alpha) < 1e-10) {
        # Limit as alpha -> 0:
        # t' -> d_des - d_anc (original branch length)
        # dt'/dalpha -> 2*(d_des^2 - d_anc^2)/2 = d_des^2 - d_anc^2
        # (from L'Hopital or Taylor expansion)
        jac <- 2 * (d_des^2 - d_anc^2) / 2
      } else {
        u <- 2 * alpha
        exp_u_des <- exp(u * d_des)
        exp_u_anc <- exp(u * d_anc)
        # dt'/du
        dtdu <- (d_des * exp_u_des - d_anc * exp_u_anc) / u -
                (exp_u_des - exp_u_anc) / u^2
        # Chain rule: du/dalpha = 2
        jac <- 2 * dtdu
      }
    },

    stop("No Jacobian defined for model: ", model,
         ". BM has no evolutionary parameter.")
  )

  jac
}
