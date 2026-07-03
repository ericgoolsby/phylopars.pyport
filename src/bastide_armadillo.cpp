// [[Rcpp::depends(RcppArmadillo)]]
#define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// Helper: Moore-Penrose pseudoinverse with log pseudo-determinant
struct MPInv {
  arma::mat inv;
  double logd;
  int rank;
};

MPInv mp_inv(const arma::mat& x, double tol = 1e-8) {
  MPInv result;
  arma::vec s;
  arma::mat U, V;

  bool ok = arma::svd(U, s, V, x);
  if (!ok) {
    result.inv = arma::zeros<arma::mat>(x.n_cols, x.n_rows);
    result.logd = -1e100;
    result.rank = 0;
    return result;
  }

  double threshold = tol * s(0);
  arma::uvec positive = arma::find(s > threshold);

  if (positive.n_elem == s.n_elem) {
    // Full rank
    result.inv = V * arma::diagmat(1.0 / s) * U.t();
    result.logd = arma::sum(arma::log(s));
    result.rank = s.n_elem;
  } else if (positive.n_elem == 0) {
    // Zero matrix
    result.inv = arma::zeros<arma::mat>(x.n_cols, x.n_rows);
    result.logd = -1e100;
    result.rank = 0;
  } else {
    // Rank deficient
    arma::vec s_pos = s(positive);
    arma::mat U_pos = U.cols(positive);
    arma::mat V_pos = V.cols(positive);
    result.inv = V_pos * arma::diagmat(1.0 / s_pos) * U_pos.t();
    result.logd = arma::sum(arma::log(s_pos));
    result.rank = positive.n_elem;
  }

  return result;
}

// Helper: log pseudo-determinant
double plogdet(const arma::mat& x, double tol = 1e-8) {
  arma::vec s = arma::svd(x);
  double threshold = tol * s(0);
  arma::uvec positive = arma::find(s > threshold);

  if (positive.n_elem == 0) return -1e100;
  return arma::sum(arma::log(s(positive)));
}


//' Bastide Traversal (Forward Only) - Armadillo
//'
//' Computes ML log-likelihood using Bastide 2021 algorithm.
//' This is the forward pass only (tips to root).
//'
//' @param R Phylogenetic rate matrix (p x p)
//' @param S Residual covariance matrix (p x p)
//' @param edge Edge matrix from phylo object (n_edges x 2), 1-indexed
//' @param edge_length Edge lengths vector
//' @param traits Trait matrix (nind x p), NA for missing
//' @param species_idx Species index for each observation (1-indexed tip node)
//' @param n_tips Number of tip species
//' @param n_nodes Total number of nodes (tips + internal)
//'
//' @return ML log-likelihood value
//' @export
// [[Rcpp::export]]
double loglik_bastide_arma(const arma::mat& R,
                           const arma::mat& S,
                           const arma::imat& edge,
                           const arma::vec& edge_length,
                           const arma::mat& traits,
                           const arma::ivec& species_idx,
                           int n_tips,
                           int n_nodes) {

  int p = R.n_rows;
  int N = traits.n_rows;
  int n_edges = edge.n_rows;
  int m = n_nodes;
  int total = N + m;
  arma::mat Ip = arma::eye(p, p);

  // Augmented edge: tree edges + fake root edge
  // In R we add: c(m + 1, n + 1) which is (n_nodes + 1, n_tips + 1)
  // Here we work with 0-indexed

  // Edge lengths: 0 for obs, tree lengths, 0 for fake root
  arma::vec len(total + 1, arma::fill::zeros);
  for (int i = 0; i < n_edges; i++) {
    len(N + i) = edge_length(i);
  }

  // Initialize storage
  arma::mat m_mat(total + 1, p, arma::fill::zeros);
  arma::field<arma::mat> P(total + 1);
  arma::field<arma::mat> Pstar(total + 1);
  arma::field<arma::mat> delta(total + 1);
  arma::field<arma::mat> Sigma(total + 1);
  arma::field<arma::mat> c_mat(total + 1);
  arma::vec r(total + 1, arma::fill::zeros);
  arma::vec r1(total + 1, arma::fill::zeros);
  arma::vec r2(total + 1, arma::fill::zeros);

  for (int i = 0; i <= total; i++) {
    P(i) = arma::zeros<arma::mat>(p, p);
    Pstar(i) = arma::zeros<arma::mat>(p, p);
    delta(i) = arma::eye<arma::mat>(p, p);
    Sigma(i) = arma::zeros<arma::mat>(p, p);
    c_mat(i) = arma::zeros<arma::mat>(p, p);
  }

  // Map observations to species (0-indexed)
  arma::ivec node_obs = species_idx - 1;  // Convert to 0-indexed

  //-------------------------------------------------------------------------
  // FORWARD PASS (tips to root)
  //-------------------------------------------------------------------------

  for (int k = 0; k < total; k++) {
    int des, anc;

    if (k < N) {
      // Processing observation k
      des = k;
      anc = node_obs(k) + N;

      // c_mat = S for observations
      c_mat(des) = S;

      // Handle missing data
      arma::uvec obs_traits;
      std::vector<arma::uword> obs_vec;
      for (int j = 0; j < p; j++) {
        if (!std::isnan(traits(k, j))) {
          obs_vec.push_back(j);
          m_mat(des, j) = traits(k, j);
        } else {
          delta(des)(j, j) = 0.0;
        }
      }

      // Pstar = ginv(delta %*% S %*% delta)
      arma::mat dSd = delta(des) * S * delta(des);
      Pstar(des) = mp_inv(dSd).inv;

      // Likelihood contribution
      int rank_Pstar = 0;
      for (int j = 0; j < p; j++) {
        if (delta(des)(j, j) > 0.5) rank_Pstar++;
      }
      r(des) = -0.5 * rank_Pstar * std::log(2.0 * arma::datum::pi) + 0.5 * plogdet(Pstar(des));

    } else {
      // Processing internal branch
      int edge_idx = k - N;

      if (edge_idx < n_edges) {
        des = edge(edge_idx, 1) - 1 + N;  // 0-indexed node + N offset
        anc = edge(edge_idx, 0) - 1 + N;
      } else {
        // Fake root edge
        des = n_tips + N;  // Root node index (n+1 in 1-indexed = n_tips in 0-indexed)
        anc = m + N;       // Fake ancestor
      }

      // c_mat = len * R
      Sigma(des) = len(k) * R;
      c_mat(des) = Sigma(des);

      // Update m_mat: m = ginv(P) %*% m
      if (arma::accu(arma::abs(P(des))) > 0) {
        arma::mat P_inv = mp_inv(P(des)).inv;
        m_mat.row(des) = (P_inv * m_mat.row(des).t()).t();
      }

      if (k < (m + N - 1)) {
        // Not at root yet - compute Pstar
        if (len(k) > 0) {
          arma::mat Sigma_inv = mp_inv(Sigma(des)).inv;
          arma::mat sum_mat = P(des) + Sigma_inv;
          arma::mat invPSigma_P;
          bool ok = arma::solve(invPSigma_P, sum_mat, P(des));
          if (!ok) return -1e100;

          Pstar(des) = P(des) - P(des) * invPSigma_P;

          // Likelihood contribution
          double qform = arma::as_scalar(m_mat.row(des) * P(des) * m_mat.row(des).t());
          r(des) = r1(des) + 0.5 * qform - 0.5 * r2(des) + 0.5 * plogdet(Ip - invPSigma_P);
        }
      }
    }

    // Propagate to ancestor
    if (k < (m + N - 1)) {
      arma::mat ta_Pstar = Pstar(des);  // a = I for BM, so t(a) %*% Pstar %*% a = Pstar
      P(anc) = P(anc) + ta_Pstar;
      m_mat.row(anc) = m_mat.row(anc) + (ta_Pstar * m_mat.row(des).t()).t();
      r1(anc) = r1(anc) + r(des);
      r2(anc) = r2(anc) + arma::as_scalar(m_mat.row(des) * Pstar(des) * m_mat.row(des).t());
    } else {
      // At root: finalize likelihood
      double qform = arma::as_scalar(m_mat.row(des) * P(des) * m_mat.row(des).t());
      r(des) = r1(des) + 0.5 * qform - 0.5 * r2(des);
    }
  }

  // Root index
  int root_idx = n_tips + N;  // Root is n+1 in 1-indexed = n_tips in 0-indexed, + N offset

  // Final log-likelihood
  double logL = r(root_idx);

  if (!arma::is_finite(logL)) return -1e100;
  return logL;
}


//' Bastide Gradient Computation - Armadillo
//'
//' Computes likelihood and analytic gradients using Bastide 2021 algorithm.
//' Runs both forward pass (for likelihood) and backward pass (for gradients).
//'
//' @param R Phylogenetic rate matrix (p x p)
//' @param S Residual covariance matrix (p x p)
//' @param edge Edge matrix from phylo object (n_edges x 2), 1-indexed
//' @param edge_length Edge lengths vector
//' @param traits Trait matrix (nind x p), NA for missing
//' @param species_idx Species index for each observation (1-indexed tip node)
//' @param n_tips Number of tip species
//' @param n_nodes Total number of nodes (tips + internal)
//'
//' @return List with logL, mu, grad_R, grad_S
//' @export
// [[Rcpp::export]]
List bastide_gradient_arma(const arma::mat& R,
                           const arma::mat& S,
                           const arma::imat& edge,
                           const arma::vec& edge_length,
                           const arma::mat& traits,
                           const arma::ivec& species_idx,
                           int n_tips,
                           int n_nodes) {

  int p = R.n_rows;
  int N = traits.n_rows;
  int n_edges = edge.n_rows;
  int m = n_nodes;
  int total = N + m;
  arma::mat Ip = arma::eye(p, p);

  // Edge lengths
  arma::vec len(total + 1, arma::fill::zeros);
  for (int i = 0; i < n_edges; i++) {
    len(N + i) = edge_length(i);
  }

  // Initialize storage for forward pass
  arma::mat m_mat(total + 1, p, arma::fill::zeros);
  arma::field<arma::mat> P(total + 1);
  arma::field<arma::mat> Pstar(total + 1);
  arma::field<arma::mat> delta(total + 1);
  arma::field<arma::mat> Sigma(total + 1);
  arma::field<arma::mat> c_mat(total + 1);
  arma::vec r(total + 1, arma::fill::zeros);
  arma::vec r1(total + 1, arma::fill::zeros);
  arma::vec r2(total + 1, arma::fill::zeros);

  // Initialize storage for backward pass
  arma::mat qn(total + 1, p, arma::fill::zeros);
  arma::mat qnstar(total + 1, p, arma::fill::zeros);
  arma::mat M(total + 1, p, arma::fill::zeros);
  arma::mat Pstar_mkmk(total + 1, p, arma::fill::zeros);
  arma::field<arma::mat> Q(total + 1);
  arma::field<arma::mat> Qstar(total + 1);
  arma::field<arma::mat> V(total + 1);
  arma::field<arma::mat> Pstar_mk(total + 1);

  for (int i = 0; i <= total; i++) {
    P(i) = arma::zeros<arma::mat>(p, p);
    Pstar(i) = arma::zeros<arma::mat>(p, p);
    delta(i) = arma::eye<arma::mat>(p, p);
    Sigma(i) = arma::zeros<arma::mat>(p, p);
    c_mat(i) = arma::zeros<arma::mat>(p, p);
    Q(i) = arma::zeros<arma::mat>(p, p);
    Qstar(i) = arma::zeros<arma::mat>(p, p);
    V(i) = arma::zeros<arma::mat>(p, p);
    Pstar_mk(i) = arma::zeros<arma::mat>(p, p);
  }

  // Map observations to species (0-indexed)
  arma::ivec node_obs = species_idx - 1;

  // Build parent_nodes_of_observations (0-indexed)
  arma::ivec parent_nodes_of_obs = node_obs;

  // Build edge lookup for backward pass
  // edge_lookup(i) = index of edge with child i (0-indexed node)
  arma::ivec edge_child(n_edges);
  arma::ivec edge_parent(n_edges);
  for (int i = 0; i < n_edges; i++) {
    edge_child(i) = edge(i, 1) - 1;  // 0-indexed
    edge_parent(i) = edge(i, 0) - 1;
  }

  //-------------------------------------------------------------------------
  // FORWARD PASS
  //-------------------------------------------------------------------------

  for (int k = 0; k < total; k++) {
    int des, anc;

    if (k < N) {
      des = k;
      anc = node_obs(k) + N;
      c_mat(des) = S;

      for (int j = 0; j < p; j++) {
        if (!std::isnan(traits(k, j))) {
          m_mat(des, j) = traits(k, j);
        } else {
          delta(des)(j, j) = 0.0;
        }
      }

      arma::mat dSd = delta(des) * S * delta(des);
      Pstar(des) = mp_inv(dSd).inv;

      int rank_Pstar = 0;
      for (int j = 0; j < p; j++) {
        if (delta(des)(j, j) > 0.5) rank_Pstar++;
      }
      r(des) = -0.5 * rank_Pstar * std::log(2.0 * arma::datum::pi) + 0.5 * plogdet(Pstar(des));

    } else {
      int edge_idx = k - N;

      if (edge_idx < n_edges) {
        des = edge(edge_idx, 1) - 1 + N;
        anc = edge(edge_idx, 0) - 1 + N;
      } else {
        des = n_tips + N;
        anc = m + N;
      }

      Sigma(des) = len(k) * R;
      c_mat(des) = Sigma(des);

      if (arma::accu(arma::abs(P(des))) > 0) {
        arma::mat P_inv = mp_inv(P(des)).inv;
        m_mat.row(des) = (P_inv * m_mat.row(des).t()).t();
      }

      if (k < (m + N - 1) && len(k) > 0) {
        arma::mat Sigma_inv = mp_inv(Sigma(des)).inv;
        arma::mat sum_mat = P(des) + Sigma_inv;
        arma::mat invPSigma_P;
        bool ok = arma::solve(invPSigma_P, sum_mat, P(des));
        if (!ok) {
          return List::create(Named("logL") = -1e100);
        }

        Pstar(des) = P(des) - P(des) * invPSigma_P;
        double qform = arma::as_scalar(m_mat.row(des) * P(des) * m_mat.row(des).t());
        r(des) = r1(des) + 0.5 * qform - 0.5 * r2(des) + 0.5 * plogdet(Ip - invPSigma_P);
      }
    }

    if (k < (m + N - 1)) {
      arma::mat ta_Pstar = Pstar(des);
      P(anc) = P(anc) + ta_Pstar;
      m_mat.row(anc) = m_mat.row(anc) + (ta_Pstar * m_mat.row(des).t()).t();
      r1(anc) = r1(anc) + r(des);
      r2(anc) = r2(anc) + arma::as_scalar(m_mat.row(des) * Pstar(des) * m_mat.row(des).t());
    } else {
      double qform = arma::as_scalar(m_mat.row(des) * P(des) * m_mat.row(des).t());
      r(des) = r1(des) + 0.5 * qform - 0.5 * r2(des);
    }
  }

  int root_idx = n_tips + N;
  double logL = r(root_idx);
  arma::rowvec mu = m_mat.row(root_idx);

  if (!arma::is_finite(logL)) {
    return List::create(Named("logL") = -1e100);
  }

  //-------------------------------------------------------------------------
  // BACKWARD PASS
  //-------------------------------------------------------------------------

  for (int k = total - 1; k >= 0; k--) {
    int des, anc;

    if (k >= N) {
      int edge_idx = k - N;
      if (edge_idx < n_edges) {
        des = edge(edge_idx, 1) - 1 + N;
        anc = edge(edge_idx, 0) - 1 + N;
      } else {
        des = n_tips + N;
        anc = m + N;
      }
    } else {
      des = k;
      anc = node_obs(k) + N;
    }

    int anc_node = anc - N;  // 0-indexed node number of ancestor

    if (anc_node == n_tips) {
      // Direct descendant of root
      Q(des) = mp_inv(c_mat(des)).inv;
      qn.row(des) = mu;  // a=I, b=0 for BM

    } else if (k < total - 1) {
      // Find sibling nodes
      if (k >= N) {
        // Internal node: siblings are other tree edges from same parent
        int current_edge = k - N;
        for (int e = 0; e < n_edges; e++) {
          if (e != current_edge && edge_parent(e) == anc_node) {
            int sib = edge_child(e) + N;
            Pstar_mk(des) = Pstar_mk(des) + Pstar(sib);
            Pstar_mkmk.row(des) = Pstar_mkmk.row(des) + (Pstar(sib) * m_mat.row(sib).t()).t();
          }
        }
      } else {
        // Observation: siblings are other obs with same species
        int my_species = node_obs(k);
        for (int obs = 0; obs < N; obs++) {
          if (obs != k && node_obs(obs) == my_species) {
            Pstar_mk(des) = Pstar_mk(des) + Pstar(obs);
            Pstar_mkmk.row(des) = Pstar_mkmk.row(des) + (Pstar(obs) * m_mat.row(obs).t()).t();
          }
        }
      }

      // Compute Q and qn
      Qstar(des) = Pstar_mk(des) + Q(anc);
      arma::mat Qstar_inv = mp_inv(Qstar(des)).inv;
      qnstar.row(des) = (Qstar_inv * (Pstar_mkmk.row(des).t() + Q(anc) * qn.row(anc).t())).t();

      arma::mat temp = Qstar_inv + c_mat(des);
      Q(des) = mp_inv(temp).inv;
      qn.row(des) = qnstar.row(des);  // a=I, b=0 for BM
    }

    // Compute M and V
    if (k == total - 1) {
      M.row(des) = m_mat.row(des);
    } else if (k >= N) {
      // Internal node
      arma::mat PQ_sum = P(des) + Q(des);
      V(des) = mp_inv(PQ_sum).inv;
      M.row(des) = (V(des) * (P(des) * m_mat.row(des).t() + Q(des) * qn.row(des).t())).t();
    } else {
      // Observation
      int n_obs_traits = 0;
      for (int j = 0; j < p; j++) {
        if (delta(k)(j, j) > 0.5) n_obs_traits++;
      }

      if (n_obs_traits == p) {
        // Fully observed
        M.row(k) = traits.row(k);
      } else if (n_obs_traits > 0) {
        // Partially observed
        for (int j = 0; j < p; j++) {
          if (delta(k)(j, j) > 0.5) {
            M(k, j) = traits(k, j);
          } else {
            // Missing trait - use conditional
            M(k, j) = qn(k, j);
          }
        }
        // V for missing traits
        arma::uvec miss_idx;
        std::vector<arma::uword> miss_vec;
        for (int j = 0; j < p; j++) {
          if (delta(k)(j, j) < 0.5) miss_vec.push_back(j);
        }
        if (!miss_vec.empty()) {
          miss_idx = arma::conv_to<arma::uvec>::from(miss_vec);
          arma::mat Q_mm = Q(k).submat(miss_idx, miss_idx);
          V(k).submat(miss_idx, miss_idx) = mp_inv(Q_mm).inv;
        }
      } else {
        // Fully missing
        arma::mat PQ_sum = P(des) + Q(des);
        V(des) = mp_inv(PQ_sum).inv;
        M.row(des) = (V(des) * (P(des) * m_mat.row(des).t() + Q(des) * qn.row(des).t())).t();
      }
    }
  }

  //-------------------------------------------------------------------------
  // COMPUTE GRADIENTS
  //-------------------------------------------------------------------------

  arma::mat grad_R(p, p, arma::fill::zeros);
  arma::mat grad_S(p, p, arma::fill::zeros);

  // Gradient w.r.t. R (sum over internal edges)
  for (int k = N; k < total; k++) {
    int edge_idx = k - N;
    if (edge_idx >= n_edges) continue;  // Skip fake root edge

    int des = edge(edge_idx, 1) - 1 + N;
    double t_k = edge_length(edge_idx);

    if (t_k > 0) {
      arma::mat Q_des = Q(des);
      arma::rowvec M_des = M.row(des);
      arma::rowvec qn_des = qn.row(des);
      arma::mat V_des = V(des);

      arma::rowvec diff = M_des - qn_des;
      arma::mat outer_diff = diff.t() * diff;

      arma::mat Q_inv = mp_inv(Q_des).inv;
      arma::mat inner = Q_inv - outer_diff - V_des;
      grad_R = grad_R - 0.5 * t_k * Q_des * inner * Q_des;
    }
  }

  // Gradient w.r.t. S (sum over observations)
  for (int k = 0; k < N; k++) {
    arma::mat Q_k = Q(k);
    arma::rowvec M_k = M.row(k);
    arma::rowvec qn_k = qn.row(k);
    arma::mat V_k = V(k);

    arma::rowvec diff = M_k - qn_k;
    arma::mat outer_diff = diff.t() * diff;

    arma::mat Q_inv = mp_inv(Q_k).inv;
    arma::mat inner = Q_inv - outer_diff - V_k;
    grad_S = grad_S - 0.5 * Q_k * inner * Q_k;
  }

  return List::create(
    Named("logL") = logL,
    Named("mu") = mu,
    Named("grad_R") = grad_R,
    Named("grad_S") = grad_S
  );
}


//' EM Algorithm - Armadillo
//'
//' Full EM algorithm for phylogenetic mixed models.
//'
//' @param R_init Initial phylogenetic rate matrix
//' @param S_init Initial residual covariance matrix
//' @param mu_init Initial mean vector
//' @param edge Edge matrix (n_edges x 2), 1-indexed
//' @param edge_length Edge lengths
//' @param traits Trait matrix (nind x p)
//' @param species_idx Species index for each observation (1-indexed)
//' @param n_tips Number of tips
//' @param n_nodes Total nodes
//' @param max_iter Maximum iterations
//' @param tol Convergence tolerance
//' @param verbose Print progress
//'
//' @return List with R, S, mu, logL, converged, iterations
//' @export
// [[Rcpp::export]]
List em_fit_arma(const arma::mat& R_init,
                 const arma::mat& S_init,
                 const arma::vec& mu_init,
                 const arma::imat& edge,
                 const arma::vec& edge_length,
                 const arma::mat& traits,
                 const arma::ivec& species_idx,
                 int n_tips,
                 int n_nodes,
                 int max_iter = 100,
                 double tol = 1e-6,
                 bool verbose = false) {

  int p = R_init.n_rows;
  int N = traits.n_rows;
  int n_edges = edge.n_rows;
  int m = n_nodes;
  int total = N + m;
  arma::mat Ip = arma::eye(p, p);

  arma::mat R = R_init;
  arma::mat S = S_init;
  arma::vec mu = mu_init;

  arma::ivec node_obs = species_idx - 1;  // 0-indexed

  double logL_prev = -1e100;
  double logL = -1e100;
  bool converged = false;
  int iter = 0;

  for (iter = 0; iter < max_iter; iter++) {
    //-----------------------------------------------------------------------
    // E-STEP: Compute conditional expectations and covariances
    //-----------------------------------------------------------------------

    // Initialize storage
    arma::field<arma::mat> pA_orig(total);
    arma::field<arma::mat> pA_prop(total);
    arma::mat pY_orig(total, p, arma::fill::zeros);
    arma::mat pY_prop(total, p, arma::fill::zeros);
    arma::mat cond_exp(total, p, arma::fill::zeros);
    arma::mat exps(total, p, arma::fill::zeros);
    arma::field<arma::mat> vars(total);
    arma::field<arma::mat> covars(total);

    for (int i = 0; i < total; i++) {
      pA_orig(i) = arma::zeros<arma::mat>(p, p);
      pA_prop(i) = arma::zeros<arma::mat>(p, p);
      vars(i) = arma::zeros<arma::mat>(p, p);
      covars(i) = arma::zeros<arma::mat>(p, p);
    }

    // Forward pass: observations
    for (int i = 0; i < N; i++) {
      int species_node = node_obs(i);
      int anc_i = species_node + N;
      int des_i = i;

      std::vector<arma::uword> obs_vec;
      for (int j = 0; j < p; j++) {
        if (!std::isnan(traits(i, j))) {
          obs_vec.push_back(j);
        }
      }
      if (obs_vec.empty()) continue;

      arma::uvec obs_idx = arma::conv_to<arma::uvec>::from(obs_vec);

      if ((int)obs_vec.size() == p) {
        // Fully observed
        arma::mat S_inv;
        bool ok = arma::inv(S_inv, S);
        if (!ok) continue;

        pA_orig(des_i) = S_inv;
        pA_prop(des_i) = S_inv;
        cond_exp.row(des_i) = traits.row(i);
        pY_orig.row(des_i) = (S_inv * traits.row(i).t()).t();
        pY_prop.row(des_i) = pY_orig.row(des_i);
      } else {
        // Partially observed
        arma::mat S_sub = S.submat(obs_idx, obs_idx);
        arma::mat S_sub_inv;
        bool ok = arma::inv(S_sub_inv, S_sub);
        if (!ok) continue;

        pA_orig(des_i).submat(obs_idx, obs_idx) = S_sub_inv;
        pA_prop(des_i) = pA_orig(des_i);
        for (size_t k = 0; k < obs_vec.size(); k++) {
          cond_exp(des_i, obs_vec[k]) = traits(i, obs_vec[k]);
        }
        arma::vec y_obs(obs_vec.size());
        for (size_t k = 0; k < obs_vec.size(); k++) {
          y_obs(k) = traits(i, obs_vec[k]);
        }
        arma::vec pY_sub = S_sub_inv * y_obs;
        for (size_t k = 0; k < obs_vec.size(); k++) {
          pY_orig(des_i, obs_vec[k]) = pY_sub(k);
        }
        pY_prop.row(des_i) = pY_orig.row(des_i);
      }

      // Propagate to species node
      pA_orig(anc_i) = pA_orig(anc_i) + pA_prop(des_i);
      pA_prop(anc_i) = pA_orig(anc_i);
      pY_orig.row(anc_i) = pY_orig.row(anc_i) + pY_prop.row(des_i);
      pY_prop.row(anc_i) = pY_orig.row(anc_i);
    }

    // Forward pass: tree edges
    for (int i = 0; i < n_edges; i++) {
      int anc_i = edge(i, 0) - 1 + N;
      int des_i = edge(i, 1) - 1 + N;
      double t_e = edge_length(i);

      arma::mat Sigma_e = t_e * R;

      // Compute conditional expectation
      if (arma::accu(arma::abs(pA_orig(des_i))) > 0) {
        arma::mat pA_inv;
        bool ok = arma::inv(pA_inv, pA_orig(des_i));
        if (ok) {
          cond_exp.row(des_i) = (pA_inv * pY_orig.row(des_i).t()).t();
        }
      }

      if (arma::accu(arma::abs(pA_orig(des_i))) > 0) {
        arma::mat itpa = Ip + Sigma_e * pA_orig(des_i);
        arma::mat itpainv;
        bool ok = arma::inv(itpainv, itpa);
        if (ok) {
          pA_prop(des_i) = pA_orig(des_i) * itpainv;
          pY_prop.row(des_i) = (itpainv.t() * pY_orig.row(des_i).t()).t();
        }
      }

      // Propagate to parent
      pA_orig(anc_i) = pA_orig(anc_i) + pA_prop(des_i);
      pA_prop(anc_i) = pA_orig(anc_i);
      pY_orig.row(anc_i) = pY_orig.row(anc_i) + pY_prop.row(des_i);
      pY_prop.row(anc_i) = pY_orig.row(anc_i);
    }

    // Root
    int root_idx = n_tips + N;
    if (arma::accu(arma::abs(pA_orig(root_idx))) > 0) {
      arma::mat root_var;
      bool ok = arma::inv(root_var, pA_orig(root_idx));
      if (ok) {
        cond_exp.row(root_idx) = (root_var * pY_orig.row(root_idx).t()).t();
        vars(root_idx) = root_var;
      }
    }
    exps.row(root_idx) = mu.t();

    // Backward pass: tree edges
    for (int i = n_edges - 1; i >= 0; i--) {
      int anc_i = edge(i, 0) - 1 + N;
      int des_i = edge(i, 1) - 1 + N;
      double t_e = edge_length(i);

      arma::mat Sigma_e = t_e * R;

      if (arma::accu(arma::abs(pA_orig(des_i))) > 0) {
        arma::mat itpa = Ip + Sigma_e * pA_orig(des_i);
        arma::mat itpainv;
        bool ok = arma::inv(itpainv, itpa);
        if (ok) {
          covars(des_i) = (itpainv * vars(anc_i)).t();

          arma::mat new_p = pA_orig(des_i) * itpainv;
          exps.row(des_i) = (itpainv * exps.row(anc_i).t() +
                             Sigma_e * new_p * cond_exp.row(des_i).t()).t();

          vars(des_i) = itpainv * Sigma_e + itpainv * covars(des_i);
        }
      } else {
        exps.row(des_i) = exps.row(anc_i);
        covars(des_i) = vars(anc_i);
        vars(des_i) = vars(anc_i) + Sigma_e;
      }
    }

    // Backward pass: observations
    for (int i = 0; i < N; i++) {
      int anc_i = node_obs(i) + N;
      int des_i = i;

      std::vector<arma::uword> obs_vec;
      for (int j = 0; j < p; j++) {
        if (!std::isnan(traits(i, j))) {
          obs_vec.push_back(j);
        }
      }

      if ((int)obs_vec.size() == p) {
        // Fully observed
        exps.row(des_i) = traits.row(i);
        vars(des_i).zeros();
        covars(des_i).zeros();
      } else if (!obs_vec.empty()) {
        // Partially observed
        arma::uvec obs_idx = arma::conv_to<arma::uvec>::from(obs_vec);
        for (size_t k = 0; k < obs_vec.size(); k++) {
          exps(des_i, obs_vec[k]) = traits(i, obs_vec[k]);
        }
        // Missing traits inherit from parent
        for (int j = 0; j < p; j++) {
          if (std::isnan(traits(i, j))) {
            exps(des_i, j) = exps(anc_i, j);
          }
        }
      } else {
        // Fully missing
        exps.row(des_i) = exps.row(anc_i);
        vars(des_i) = S + vars(anc_i);
        covars(des_i) = vars(anc_i);
      }
    }

    //-----------------------------------------------------------------------
    // M-STEP: Update R, S, mu
    //-----------------------------------------------------------------------

    // Update mu
    std::vector<int> root_edges;
    for (int i = 0; i < n_edges; i++) {
      if (edge(i, 0) - 1 == n_tips) {
        root_edges.push_back(i);
      }
    }

    if (!root_edges.empty()) {
      double weight_sum = 0;
      arma::vec mu_new(p, arma::fill::zeros);
      for (size_t i = 0; i < root_edges.size(); i++) {
        int e = root_edges[i];
        double w = 1.0 / edge_length(e);
        int child_idx = edge(e, 1) - 1 + N;
        mu_new = mu_new + w * exps.row(child_idx).t();
        weight_sum += w;
      }
      mu = mu_new / weight_sum;
    }

    // Update R
    arma::mat sum_exp(p, p, arma::fill::zeros);
    arma::mat sum_var(p, p, arma::fill::zeros);

    for (int e = 0; e < n_edges; e++) {
      int parent_node = edge(e, 0) - 1;
      int child_node = edge(e, 1) - 1;
      double t_e = edge_length(e);

      int parent_idx = parent_node + N;
      int child_idx = child_node + N;

      arma::rowvec diff_exp = exps.row(child_idx) - exps.row(parent_idx);
      sum_exp = sum_exp + diff_exp.t() * diff_exp / t_e;

      arma::mat var_child = vars(child_idx);
      arma::mat var_parent = vars(parent_idx);
      arma::mat cov_cp = covars(child_idx);

      arma::mat var_diff = var_child + var_parent - cov_cp - cov_cp.t();
      sum_var = sum_var + var_diff / t_e;
    }

    R = (sum_exp + sum_var) / (m - 1);
    R = (R + R.t()) / 2.0;

    // Update S
    arma::mat sum_exp_ind(p, p, arma::fill::zeros);
    arma::mat sum_var_ind(p, p, arma::fill::zeros);

    for (int i = 0; i < N; i++) {
      int species_idx_i = node_obs(i) + N;

      arma::rowvec diff_exp = exps.row(species_idx_i) - exps.row(i);
      sum_exp_ind = sum_exp_ind + diff_exp.t() * diff_exp;

      arma::mat var_species = vars(species_idx_i);
      arma::mat var_obs = vars(i);
      arma::mat cov_so = covars(i);

      arma::mat var_diff = var_species + var_obs - cov_so - cov_so.t();
      sum_var_ind = sum_var_ind + var_diff;
    }

    S = (sum_exp_ind + sum_var_ind) / N;
    S = (S + S.t()) / 2.0;

    // Ensure positive definiteness
    arma::vec eigval_R;
    arma::mat eigvec_R;
    arma::eig_sym(eigval_R, eigvec_R, R);
    for (int i = 0; i < p; i++) {
      if (eigval_R(i) <= 1e-10) eigval_R(i) = 1e-10;
    }
    R = eigvec_R * arma::diagmat(eigval_R) * eigvec_R.t();

    arma::vec eigval_S;
    arma::mat eigvec_S;
    arma::eig_sym(eigval_S, eigvec_S, S);
    for (int i = 0; i < p; i++) {
      if (eigval_S(i) <= 1e-10) eigval_S(i) = 1e-10;
    }
    S = eigvec_S * arma::diagmat(eigval_S) * eigvec_S.t();

    // Compute log-likelihood
    logL = loglik_bastide_arma(R, S, edge, edge_length, traits, species_idx, n_tips, n_nodes);

    if (verbose) {
      Rcpp::Rcout << "EM iter " << (iter + 1) << ": logL = " << logL << std::endl;
    }

    if (std::abs(logL - logL_prev) < tol) {
      converged = true;
      break;
    }

    logL_prev = logL;
  }

  return List::create(
    Named("R") = R,
    Named("S") = S,
    Named("mu") = mu,
    Named("logL") = logL,
    Named("converged") = converged,
    Named("iterations") = iter + 1
  );
}
