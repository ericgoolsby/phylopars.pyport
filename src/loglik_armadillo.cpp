// [[Rcpp::depends(RcppArmadillo)]]
#define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

//' Restricted log-likelihood using RcppArmadillo
//'
//' Fast C++ implementation of the postorder tree traversal algorithm.
//'
//' @param Sigma Phylogenetic covariance matrix (m x m)
//' @param B Residual covariance matrix (m x m)
//' @param edge Edge matrix from phylo object (n_edges x 2), 1-indexed
//' @param edge_length Edge lengths vector
//' @param traits Trait matrix (nind x m), NA for missing
//' @param species_idx Species index for each observation (1-indexed tip node)
//' @param n_obs Total number of non-missing observations
//' @param Nnode Total number of nodes (tips + internal)
//'
//' @return Log-likelihood value, or -1e100 if singular
//' @export
// [[Rcpp::export]]
double restricted_loglik_arma(const arma::mat& Sigma,
                               const arma::mat& B,
                               const arma::imat& edge,
                               const arma::vec& edge_length,
                               const arma::mat& traits,
                               const arma::ivec& species_idx,
                               int n_obs,
                               int Nnode) {
  
  int m = Sigma.n_rows;
  int nind = traits.n_rows;
  int n_edges = edge.n_rows;
  
  // Initialize per-node storage using field (like list of matrices)
  arma::field<arma::mat> p(Nnode);
  arma::field<arma::vec> Vr(Nnode);
  arma::vec Q(Nnode, fill::zeros);
  arma::vec logW(Nnode, fill::zeros);
  
  for (int i = 0; i < Nnode; i++) {
    p(i) = arma::zeros<arma::mat>(m, m);
    Vr(i) = arma::zeros<arma::vec>(m);
  }
  
  // Process observations for each individual
  for (int idx = 0; idx < nind; idx++) {
    int sp_node = species_idx(idx) - 1;  // Convert to 0-indexed
    
    // Find which traits are observed for this individual
    std::vector<arma::uword> obs_vec;
    for (int j = 0; j < m; j++) {
      if (!std::isnan(traits(idx, j))) {
        obs_vec.push_back(j);
      }
    }
    if (obs_vec.empty()) continue;
    
    arma::uvec obs_traits = arma::conv_to<arma::uvec>::from(obs_vec);
    int n_obs_traits = obs_traits.n_elem;
    
    // Extract observed trait values and submatrix of B
    arma::vec Yvec(n_obs_traits);
    for (size_t k = 0; k < obs_traits.n_elem; k++) {
      Yvec(k) = traits(idx, obs_traits(k));
    }
    arma::mat Ba = B.submat(obs_traits, obs_traits);
    
    // Check condition number before inverting
    double rcond = arma::rcond(Ba);
    if (rcond < 1e-15 || !arma::is_finite(rcond)) {
      return -1e100;
    }
    
    arma::mat Ba_inv;
    bool success = arma::inv(Ba_inv, Ba);
    if (!success) return -1e100;
    
    // Update p, Vr, Q, logW for this observation
    for (size_t k = 0; k < obs_traits.n_elem; k++) {
      for (size_t l = 0; l < obs_traits.n_elem; l++) {
        p(sp_node)(obs_traits(k), obs_traits(l)) += Ba_inv(k, l);
      }
      Vr(sp_node)(obs_traits(k)) += arma::dot(Ba_inv.row(k), Yvec);
    }
    Q(sp_node) += arma::as_scalar(Yvec.t() * Ba_inv * Yvec);
    
    double log_det_val;
    double log_det_sign;
    arma::log_det(log_det_val, log_det_sign, Ba);
    if (!arma::is_finite(log_det_val)) return -1e100;
    logW(sp_node) += log_det_val;
  }
  
  // Postorder traversal (tips to root)
  int root = -1;
  for (int i = 0; i < n_edges; i++) {
    int anc = edge(i, 0) - 1;  // 0-indexed
    int des = edge(i, 1) - 1;
    double t = edge_length(i);
    arma::mat T_edge = t * Sigma;
    
    arma::mat pA = p(des);
    arma::mat itpa = arma::eye(m, m) + T_edge * pA;
    
    // Check condition number
    double rcond = arma::rcond(itpa);
    if (rcond < 1e-15 || !arma::is_finite(rcond)) {
      return -1e100;
    }
    
    arma::mat itpainv;
    bool inv_ok = arma::inv(itpainv, itpa);
    if (!inv_ok) return -1e100;
    
    p(des) = pA * itpainv;
    
    double log_det_val;
    double log_det_sign;
    arma::log_det(log_det_val, log_det_sign, itpa);
    if (!arma::is_finite(log_det_val)) return -1e100;
    logW(des) += log_det_val;
    
    Q(des) -= arma::as_scalar(Vr(des).t() * itpainv * T_edge * Vr(des));
    Vr(des) = itpainv.t() * Vr(des);
    
    p(anc) += p(des);
    Vr(anc) += Vr(des);
    Q(anc) += Q(des);
    logW(anc) += logW(des);
    
    root = anc;
  }
  
  // Check root p matrix condition
  double rcond_root = arma::rcond(p(root));
  if (rcond_root < 1e-15 || !arma::is_finite(rcond_root)) {
    return -1e100;
  }
  
  // Compute mean and log-likelihood at root
  arma::vec mu;
  bool solve_success = arma::solve(mu, p(root), Vr(root));
  if (!solve_success) return -1e100;
  
  double log_det_p;
  double log_det_sign;
  arma::log_det(log_det_p, log_det_sign, p(root));
  if (!arma::is_finite(log_det_p)) return -1e100;
  
  double logL = -0.5 * (
    (n_obs - m) * std::log(2.0 * arma::datum::pi) +
    log_det_p +
    logW(root) +
    Q(root) -
    2.0 * arma::dot(mu, Vr(root)) +
    arma::as_scalar(mu.t() * p(root) * mu)
  );
  
  if (!arma::is_finite(logL)) return -1e100;
  
  return logL;
}
