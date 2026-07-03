// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <cmath>
using namespace Rcpp;

typedef Eigen::MatrixXd MatrixXd;
typedef Eigen::VectorXd VectorXd;

//' Restricted log-likelihood using RcppEigen
//'
//' Fast C++ implementation using Eigen library.
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
double restricted_loglik_eigen(const Eigen::MatrixXd& Sigma,
                                const Eigen::MatrixXd& B,
                                const Eigen::MatrixXi& edge,
                                const Eigen::VectorXd& edge_length,
                                const Eigen::MatrixXd& traits,
                                const Eigen::VectorXi& species_idx,
                                int n_obs,
                                int Nnode) {
  
  int m = Sigma.rows();
  int nind = traits.rows();
  int n_edges = edge.rows();
  
  // Initialize per-node storage
  std::vector<MatrixXd> p(Nnode, MatrixXd::Zero(m, m));
  std::vector<VectorXd> Vr(Nnode, VectorXd::Zero(m));
  VectorXd Q = VectorXd::Zero(Nnode);
  VectorXd logW = VectorXd::Zero(Nnode);
  
  // Process observations for each individual
  for (int idx = 0; idx < nind; idx++) {
    int sp_node = species_idx(idx) - 1;  // Convert to 0-indexed
    
    // Find which traits are observed
    std::vector<int> obs_traits;
    for (int j = 0; j < m; j++) {
      if (!std::isnan(traits(idx, j))) {
        obs_traits.push_back(j);
      }
    }
    if (obs_traits.empty()) continue;
    
    int n_obs_traits = obs_traits.size();
    
    // Extract observed trait values and submatrix of B
    VectorXd Yvec(n_obs_traits);
    MatrixXd Ba(n_obs_traits, n_obs_traits);
    for (int k = 0; k < n_obs_traits; k++) {
      Yvec(k) = traits(idx, obs_traits[k]);
      for (int l = 0; l < n_obs_traits; l++) {
        Ba(k, l) = B(obs_traits[k], obs_traits[l]);
      }
    }
    
    // Use LLT (Cholesky) - fails cleanly if not positive definite
    Eigen::LLT<MatrixXd> llt_Ba(Ba);
    if (llt_Ba.info() != Eigen::Success) {
      return -1e100;
    }
    
    MatrixXd Ba_inv = llt_Ba.solve(MatrixXd::Identity(n_obs_traits, n_obs_traits));
    
    for (int k = 0; k < n_obs_traits; k++) {
      for (int l = 0; l < n_obs_traits; l++) {
        p[sp_node](obs_traits[k], obs_traits[l]) += Ba_inv(k, l);
      }
      Vr[sp_node](obs_traits[k]) += Ba_inv.row(k).dot(Yvec);
    }
    Q(sp_node) += Yvec.transpose() * Ba_inv * Yvec;
    
    // log det from Cholesky: 2 * sum(log(diag(L)))
    double log_det = 2.0 * llt_Ba.matrixL().toDenseMatrix().diagonal().array().log().sum();
    if (!std::isfinite(log_det)) return -1e100;
    logW(sp_node) += log_det;
  }
  
  // Postorder traversal (tips to root)
  int root = -1;
  for (int i = 0; i < n_edges; i++) {
    int anc = edge(i, 0) - 1;  // 0-indexed
    int des = edge(i, 1) - 1;
    double t = edge_length(i);
    MatrixXd T_edge = t * Sigma;
    
    MatrixXd pA = p[des];
    MatrixXd itpa = MatrixXd::Identity(m, m) + T_edge * pA;
    
    // Use PartialPivLU for general matrices
    Eigen::PartialPivLU<MatrixXd> lu_itpa(itpa);
    
    // Check for singularity via determinant
    double det = lu_itpa.determinant();
    if (std::abs(det) < 1e-100 || !std::isfinite(det)) {
      return -1e100;
    }
    
    MatrixXd itpainv = lu_itpa.inverse();
    
    p[des] = pA * itpainv;
    double log_det = std::log(std::abs(det));
    if (!std::isfinite(log_det)) return -1e100;
    logW(des) += log_det;
    
    Q(des) -= Vr[des].transpose() * itpainv * T_edge * Vr[des];
    Vr[des] = itpainv.transpose() * Vr[des];
    
    p[anc] += p[des];
    Vr[anc] += Vr[des];
    Q(anc) += Q(des);
    logW(anc) += logW(des);
    
    root = anc;
  }
  
  // Check root p matrix
  Eigen::PartialPivLU<MatrixXd> lu_p_root(p[root]);
  double det_root = lu_p_root.determinant();
  if (std::abs(det_root) < 1e-100 || !std::isfinite(det_root)) {
    return -1e100;
  }
  
  VectorXd mu = lu_p_root.solve(Vr[root]);
  
  double log_det_p = std::log(std::abs(det_root));
  if (!std::isfinite(log_det_p)) return -1e100;
  
  double logL = -0.5 * (
    (n_obs - m) * std::log(2.0 * M_PI) +
    log_det_p +
    logW(root) +
    Q(root) -
    2.0 * mu.dot(Vr[root]) +
    mu.transpose() * p[root] * mu
  );
  
  if (!std::isfinite(logL)) return -1e100;
  
  return logL;
}
