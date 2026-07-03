// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <cmath>
using namespace Rcpp;

// Helper: Moore-Penrose pseudoinverse with log pseudo-determinant
struct MPInvEigen {
  Eigen::MatrixXd inv;
  double logd;
  int rank;
};

MPInvEigen mp_inv_eigen(const Eigen::MatrixXd& x, double tol = 1e-8) {
  MPInvEigen result;

  Eigen::JacobiSVD<Eigen::MatrixXd> svd(x, Eigen::ComputeThinU | Eigen::ComputeThinV);
  Eigen::VectorXd s = svd.singularValues();

  double threshold = tol * s(0);
  int rank = 0;
  double logd = 0.0;

  for (int i = 0; i < s.size(); i++) {
    if (s(i) > threshold) {
      rank++;
      logd += std::log(s(i));
    }
  }

  if (rank == 0) {
    result.inv = Eigen::MatrixXd::Zero(x.cols(), x.rows());
    result.logd = -1e100;
    result.rank = 0;
    return result;
  }

  // Compute pseudoinverse
  Eigen::VectorXd s_inv = Eigen::VectorXd::Zero(s.size());
  for (int i = 0; i < s.size(); i++) {
    if (s(i) > threshold) {
      s_inv(i) = 1.0 / s(i);
    }
  }

  result.inv = svd.matrixV() * s_inv.asDiagonal() * svd.matrixU().transpose();
  result.logd = logd;
  result.rank = rank;

  return result;
}

// Helper: log pseudo-determinant
double plogdet_eigen(const Eigen::MatrixXd& x, double tol = 1e-8) {
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(x, Eigen::ComputeThinU | Eigen::ComputeThinV);
  Eigen::VectorXd s = svd.singularValues();

  double threshold = tol * s(0);
  double logd = 0.0;
  int count = 0;

  for (int i = 0; i < s.size(); i++) {
    if (s(i) > threshold) {
      logd += std::log(s(i));
      count++;
    }
  }

  if (count == 0) return -1e100;
  return logd;
}


//' Bastide Traversal (Forward Only) - Eigen
//'
//' Computes ML log-likelihood using Bastide 2021 algorithm (Eigen implementation).
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
double loglik_bastide_eigen(const Eigen::MatrixXd& R,
                            const Eigen::MatrixXd& S,
                            const Eigen::MatrixXi& edge,
                            const Eigen::VectorXd& edge_length,
                            const Eigen::MatrixXd& traits,
                            const Eigen::VectorXi& species_idx,
                            int n_tips,
                            int n_nodes) {

  int p = R.rows();
  int N = traits.rows();
  int n_edges = edge.rows();
  int m = n_nodes;
  int total = N + m;
  Eigen::MatrixXd Ip = Eigen::MatrixXd::Identity(p, p);

  // Edge lengths
  Eigen::VectorXd len = Eigen::VectorXd::Zero(total + 1);
  for (int i = 0; i < n_edges; i++) {
    len(N + i) = edge_length(i);
  }

  // Initialize storage
  Eigen::MatrixXd m_mat = Eigen::MatrixXd::Zero(total + 1, p);
  std::vector<Eigen::MatrixXd> P(total + 1, Eigen::MatrixXd::Zero(p, p));
  std::vector<Eigen::MatrixXd> Pstar(total + 1, Eigen::MatrixXd::Zero(p, p));
  std::vector<Eigen::MatrixXd> delta(total + 1, Eigen::MatrixXd::Identity(p, p));
  std::vector<Eigen::MatrixXd> Sigma(total + 1, Eigen::MatrixXd::Zero(p, p));
  std::vector<Eigen::MatrixXd> c_mat(total + 1, Eigen::MatrixXd::Zero(p, p));
  Eigen::VectorXd r = Eigen::VectorXd::Zero(total + 1);
  Eigen::VectorXd r1 = Eigen::VectorXd::Zero(total + 1);
  Eigen::VectorXd r2 = Eigen::VectorXd::Zero(total + 1);

  // Map observations to species (0-indexed)
  Eigen::VectorXi node_obs = species_idx.array() - 1;

  //-------------------------------------------------------------------------
  // FORWARD PASS
  //-------------------------------------------------------------------------

  for (int k = 0; k < total; k++) {
    int des, anc;

    if (k < N) {
      des = k;
      anc = node_obs(k) + N;
      c_mat[des] = S;

      for (int j = 0; j < p; j++) {
        if (!std::isnan(traits(k, j))) {
          m_mat(des, j) = traits(k, j);
        } else {
          delta[des](j, j) = 0.0;
        }
      }

      Eigen::MatrixXd dSd = delta[des] * S * delta[des];
      Pstar[des] = mp_inv_eigen(dSd).inv;

      int rank_Pstar = 0;
      for (int j = 0; j < p; j++) {
        if (delta[des](j, j) > 0.5) rank_Pstar++;
      }
      r(des) = -0.5 * rank_Pstar * std::log(2.0 * M_PI) + 0.5 * plogdet_eigen(Pstar[des]);

    } else {
      int edge_idx = k - N;

      if (edge_idx < n_edges) {
        des = edge(edge_idx, 1) - 1 + N;
        anc = edge(edge_idx, 0) - 1 + N;
      } else {
        des = n_tips + N;
        anc = m + N;
      }

      Sigma[des] = len(k) * R;
      c_mat[des] = Sigma[des];

      if (P[des].norm() > 0) {
        Eigen::MatrixXd P_inv = mp_inv_eigen(P[des]).inv;
        m_mat.row(des) = (P_inv * m_mat.row(des).transpose()).transpose();
      }

      if (k < (m + N - 1) && len(k) > 0) {
        Eigen::MatrixXd Sigma_inv = mp_inv_eigen(Sigma[des]).inv;
        Eigen::MatrixXd sum_mat = P[des] + Sigma_inv;
        Eigen::MatrixXd invPSigma_P = sum_mat.ldlt().solve(P[des]);

        Pstar[des] = P[des] - P[des] * invPSigma_P;
        double qform = (m_mat.row(des) * P[des] * m_mat.row(des).transpose())(0, 0);
        r(des) = r1(des) + 0.5 * qform - 0.5 * r2(des) + 0.5 * plogdet_eigen(Ip - invPSigma_P);
      }
    }

    if (k < (m + N - 1)) {
      Eigen::MatrixXd ta_Pstar = Pstar[des];
      P[anc] = P[anc] + ta_Pstar;
      m_mat.row(anc) = m_mat.row(anc) + (ta_Pstar * m_mat.row(des).transpose()).transpose();
      r1(anc) = r1(anc) + r(des);
      r2(anc) = r2(anc) + (m_mat.row(des) * Pstar[des] * m_mat.row(des).transpose())(0, 0);
    } else {
      double qform = (m_mat.row(des) * P[des] * m_mat.row(des).transpose())(0, 0);
      r(des) = r1(des) + 0.5 * qform - 0.5 * r2(des);
    }
  }

  int root_idx = n_tips + N;
  double logL = r(root_idx);

  if (!std::isfinite(logL)) return -1e100;
  return logL;
}


//' Bastide Gradient Computation - Eigen
//'
//' Computes likelihood and analytic gradients (Eigen implementation).
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
List bastide_gradient_eigen(const Eigen::MatrixXd& R,
                            const Eigen::MatrixXd& S,
                            const Eigen::MatrixXi& edge,
                            const Eigen::VectorXd& edge_length,
                            const Eigen::MatrixXd& traits,
                            const Eigen::VectorXi& species_idx,
                            int n_tips,
                            int n_nodes) {

  int p = R.rows();
  int N = traits.rows();
  int n_edges = edge.rows();
  int m = n_nodes;
  int total = N + m;
  Eigen::MatrixXd Ip = Eigen::MatrixXd::Identity(p, p);

  // Edge lengths
  Eigen::VectorXd len = Eigen::VectorXd::Zero(total + 1);
  for (int i = 0; i < n_edges; i++) {
    len(N + i) = edge_length(i);
  }

  // Initialize storage
  Eigen::MatrixXd m_mat = Eigen::MatrixXd::Zero(total + 1, p);
  std::vector<Eigen::MatrixXd> P(total + 1, Eigen::MatrixXd::Zero(p, p));
  std::vector<Eigen::MatrixXd> Pstar(total + 1, Eigen::MatrixXd::Zero(p, p));
  std::vector<Eigen::MatrixXd> delta(total + 1, Eigen::MatrixXd::Identity(p, p));
  std::vector<Eigen::MatrixXd> Sigma(total + 1, Eigen::MatrixXd::Zero(p, p));
  std::vector<Eigen::MatrixXd> c_mat(total + 1, Eigen::MatrixXd::Zero(p, p));
  Eigen::VectorXd r = Eigen::VectorXd::Zero(total + 1);
  Eigen::VectorXd r1 = Eigen::VectorXd::Zero(total + 1);
  Eigen::VectorXd r2 = Eigen::VectorXd::Zero(total + 1);

  // Backward pass storage
  Eigen::MatrixXd qn = Eigen::MatrixXd::Zero(total + 1, p);
  Eigen::MatrixXd qnstar = Eigen::MatrixXd::Zero(total + 1, p);
  Eigen::MatrixXd M = Eigen::MatrixXd::Zero(total + 1, p);
  Eigen::MatrixXd Pstar_mkmk = Eigen::MatrixXd::Zero(total + 1, p);
  std::vector<Eigen::MatrixXd> Q(total + 1, Eigen::MatrixXd::Zero(p, p));
  std::vector<Eigen::MatrixXd> Qstar(total + 1, Eigen::MatrixXd::Zero(p, p));
  std::vector<Eigen::MatrixXd> V(total + 1, Eigen::MatrixXd::Zero(p, p));
  std::vector<Eigen::MatrixXd> Pstar_mk(total + 1, Eigen::MatrixXd::Zero(p, p));

  Eigen::VectorXi node_obs = species_idx.array() - 1;

  //-------------------------------------------------------------------------
  // FORWARD PASS
  //-------------------------------------------------------------------------

  for (int k = 0; k < total; k++) {
    int des, anc;

    if (k < N) {
      des = k;
      anc = node_obs(k) + N;
      c_mat[des] = S;

      for (int j = 0; j < p; j++) {
        if (!std::isnan(traits(k, j))) {
          m_mat(des, j) = traits(k, j);
        } else {
          delta[des](j, j) = 0.0;
        }
      }

      Eigen::MatrixXd dSd = delta[des] * S * delta[des];
      Pstar[des] = mp_inv_eigen(dSd).inv;

      int rank_Pstar = 0;
      for (int j = 0; j < p; j++) {
        if (delta[des](j, j) > 0.5) rank_Pstar++;
      }
      r(des) = -0.5 * rank_Pstar * std::log(2.0 * M_PI) + 0.5 * plogdet_eigen(Pstar[des]);

    } else {
      int edge_idx = k - N;

      if (edge_idx < n_edges) {
        des = edge(edge_idx, 1) - 1 + N;
        anc = edge(edge_idx, 0) - 1 + N;
      } else {
        des = n_tips + N;
        anc = m + N;
      }

      Sigma[des] = len(k) * R;
      c_mat[des] = Sigma[des];

      if (P[des].norm() > 0) {
        Eigen::MatrixXd P_inv = mp_inv_eigen(P[des]).inv;
        m_mat.row(des) = (P_inv * m_mat.row(des).transpose()).transpose();
      }

      if (k < (m + N - 1) && len(k) > 0) {
        Eigen::MatrixXd Sigma_inv = mp_inv_eigen(Sigma[des]).inv;
        Eigen::MatrixXd sum_mat = P[des] + Sigma_inv;
        Eigen::MatrixXd invPSigma_P = sum_mat.ldlt().solve(P[des]);

        Pstar[des] = P[des] - P[des] * invPSigma_P;
        double qform = (m_mat.row(des) * P[des] * m_mat.row(des).transpose())(0, 0);
        r(des) = r1(des) + 0.5 * qform - 0.5 * r2(des) + 0.5 * plogdet_eigen(Ip - invPSigma_P);
      }
    }

    if (k < (m + N - 1)) {
      Eigen::MatrixXd ta_Pstar = Pstar[des];
      P[anc] = P[anc] + ta_Pstar;
      m_mat.row(anc) = m_mat.row(anc) + (ta_Pstar * m_mat.row(des).transpose()).transpose();
      r1(anc) = r1(anc) + r(des);
      r2(anc) = r2(anc) + (m_mat.row(des) * Pstar[des] * m_mat.row(des).transpose())(0, 0);
    } else {
      double qform = (m_mat.row(des) * P[des] * m_mat.row(des).transpose())(0, 0);
      r(des) = r1(des) + 0.5 * qform - 0.5 * r2(des);
    }
  }

  int root_idx = n_tips + N;
  double logL = r(root_idx);
  Eigen::RowVectorXd mu = m_mat.row(root_idx);

  if (!std::isfinite(logL)) {
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

    int anc_node = anc - N;

    if (anc_node == n_tips) {
      Q[des] = mp_inv_eigen(c_mat[des]).inv;
      qn.row(des) = mu;

    } else if (k < total - 1) {
      // Find siblings
      if (k >= N) {
        int current_edge = k - N;
        for (int e = 0; e < n_edges; e++) {
          if (e != current_edge && (edge(e, 0) - 1) == anc_node) {
            int sib = edge(e, 1) - 1 + N;
            Pstar_mk[des] = Pstar_mk[des] + Pstar[sib];
            Pstar_mkmk.row(des) = Pstar_mkmk.row(des) + (Pstar[sib] * m_mat.row(sib).transpose()).transpose();
          }
        }
      } else {
        int my_species = node_obs(k);
        for (int obs = 0; obs < N; obs++) {
          if (obs != k && node_obs(obs) == my_species) {
            Pstar_mk[des] = Pstar_mk[des] + Pstar[obs];
            Pstar_mkmk.row(des) = Pstar_mkmk.row(des) + (Pstar[obs] * m_mat.row(obs).transpose()).transpose();
          }
        }
      }

      Qstar[des] = Pstar_mk[des] + Q[anc];
      Eigen::MatrixXd Qstar_inv = mp_inv_eigen(Qstar[des]).inv;
      Eigen::VectorXd temp = Pstar_mkmk.row(des).transpose() + Q[anc] * qn.row(anc).transpose();
      qnstar.row(des) = (Qstar_inv * temp).transpose();

      Eigen::MatrixXd temp2 = Qstar_inv + c_mat[des];
      Q[des] = mp_inv_eigen(temp2).inv;
      qn.row(des) = qnstar.row(des);
    }

    // Compute M and V
    if (k == total - 1) {
      M.row(des) = m_mat.row(des);
    } else if (k >= N) {
      Eigen::MatrixXd PQ_sum = P[des] + Q[des];
      V[des] = mp_inv_eigen(PQ_sum).inv;
      Eigen::VectorXd temp = P[des] * m_mat.row(des).transpose() + Q[des] * qn.row(des).transpose();
      M.row(des) = (V[des] * temp).transpose();
    } else {
      int n_obs_traits = 0;
      for (int j = 0; j < p; j++) {
        if (delta[k](j, j) > 0.5) n_obs_traits++;
      }

      if (n_obs_traits == p) {
        M.row(k) = traits.row(k);
        // V[k] stays zero for fully observed
      } else if (n_obs_traits > 0) {
        // Partially observed - set M and compute V for missing traits
        std::vector<int> miss_vec;
        for (int j = 0; j < p; j++) {
          if (delta[k](j, j) > 0.5) {
            M(k, j) = traits(k, j);
          } else {
            M(k, j) = qn(k, j);
            miss_vec.push_back(j);
          }
        }
        // Compute V for missing traits
        if (!miss_vec.empty()) {
          int n_miss = miss_vec.size();
          Eigen::MatrixXd Q_mm(n_miss, n_miss);
          for (int i = 0; i < n_miss; i++) {
            for (int j = 0; j < n_miss; j++) {
              Q_mm(i, j) = Q[k](miss_vec[i], miss_vec[j]);
            }
          }
          Eigen::MatrixXd Q_mm_inv = mp_inv_eigen(Q_mm).inv;
          for (int i = 0; i < n_miss; i++) {
            for (int j = 0; j < n_miss; j++) {
              V[k](miss_vec[i], miss_vec[j]) = Q_mm_inv(i, j);
            }
          }
        }
      } else {
        Eigen::MatrixXd PQ_sum = P[des] + Q[des];
        V[des] = mp_inv_eigen(PQ_sum).inv;
        Eigen::VectorXd temp = P[des] * m_mat.row(des).transpose() + Q[des] * qn.row(des).transpose();
        M.row(des) = (V[des] * temp).transpose();
      }
    }
  }

  //-------------------------------------------------------------------------
  // COMPUTE GRADIENTS
  //-------------------------------------------------------------------------

  Eigen::MatrixXd grad_R = Eigen::MatrixXd::Zero(p, p);
  Eigen::MatrixXd grad_S = Eigen::MatrixXd::Zero(p, p);

  // Gradient w.r.t. R
  for (int k = N; k < total; k++) {
    int edge_idx = k - N;
    if (edge_idx >= n_edges) continue;

    int des = edge(edge_idx, 1) - 1 + N;
    double t_k = edge_length(edge_idx);

    if (t_k > 0) {
      Eigen::MatrixXd Q_des = Q[des];
      Eigen::RowVectorXd M_des = M.row(des);
      Eigen::RowVectorXd qn_des = qn.row(des);
      Eigen::MatrixXd V_des = V[des];

      Eigen::RowVectorXd diff = M_des - qn_des;
      Eigen::MatrixXd outer_diff = diff.transpose() * diff;

      Eigen::MatrixXd Q_inv = mp_inv_eigen(Q_des).inv;
      Eigen::MatrixXd inner = Q_inv - outer_diff - V_des;
      grad_R = grad_R - 0.5 * t_k * Q_des * inner * Q_des;
    }
  }

  // Gradient w.r.t. S
  for (int k = 0; k < N; k++) {
    Eigen::MatrixXd Q_k = Q[k];
    Eigen::RowVectorXd M_k = M.row(k);
    Eigen::RowVectorXd qn_k = qn.row(k);
    Eigen::MatrixXd V_k = V[k];

    Eigen::RowVectorXd diff = M_k - qn_k;
    Eigen::MatrixXd outer_diff = diff.transpose() * diff;

    Eigen::MatrixXd Q_inv = mp_inv_eigen(Q_k).inv;
    Eigen::MatrixXd inner = Q_inv - outer_diff - V_k;
    grad_S = grad_S - 0.5 * Q_k * inner * Q_k;
  }

  return List::create(
    Named("logL") = logL,
    Named("mu") = mu,
    Named("grad_R") = grad_R,
    Named("grad_S") = grad_S
  );
}


//' EM Algorithm - Eigen
//'
//' Full EM algorithm (Eigen implementation).
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
List em_fit_eigen(const Eigen::MatrixXd& R_init,
                  const Eigen::MatrixXd& S_init,
                  const Eigen::VectorXd& mu_init,
                  const Eigen::MatrixXi& edge,
                  const Eigen::VectorXd& edge_length,
                  const Eigen::MatrixXd& traits,
                  const Eigen::VectorXi& species_idx,
                  int n_tips,
                  int n_nodes,
                  int max_iter = 100,
                  double tol = 1e-6,
                  bool verbose = false) {

  int p = R_init.rows();
  int N = traits.rows();
  int n_edges = edge.rows();
  int m = n_nodes;
  int total = N + m;
  Eigen::MatrixXd Ip = Eigen::MatrixXd::Identity(p, p);

  Eigen::MatrixXd R = R_init;
  Eigen::MatrixXd S = S_init;
  Eigen::VectorXd mu = mu_init;

  Eigen::VectorXi node_obs = species_idx.array() - 1;

  double logL_prev = -1e100;
  double logL = -1e100;
  bool converged = false;
  int iter = 0;

  for (iter = 0; iter < max_iter; iter++) {
    // Initialize storage
    std::vector<Eigen::MatrixXd> pA_orig(total, Eigen::MatrixXd::Zero(p, p));
    std::vector<Eigen::MatrixXd> pA_prop(total, Eigen::MatrixXd::Zero(p, p));
    Eigen::MatrixXd pY_orig = Eigen::MatrixXd::Zero(total, p);
    Eigen::MatrixXd pY_prop = Eigen::MatrixXd::Zero(total, p);
    Eigen::MatrixXd cond_exp = Eigen::MatrixXd::Zero(total, p);
    Eigen::MatrixXd exps = Eigen::MatrixXd::Zero(total, p);
    std::vector<Eigen::MatrixXd> vars(total, Eigen::MatrixXd::Zero(p, p));
    std::vector<Eigen::MatrixXd> covars(total, Eigen::MatrixXd::Zero(p, p));

    //-----------------------------------------------------------------------
    // E-STEP: Forward pass - observations
    //-----------------------------------------------------------------------
    for (int i = 0; i < N; i++) {
      int anc_i = node_obs(i) + N;
      int des_i = i;

      std::vector<int> obs_vec;
      for (int j = 0; j < p; j++) {
        if (!std::isnan(traits(i, j))) {
          obs_vec.push_back(j);
        }
      }
      if (obs_vec.empty()) continue;

      if ((int)obs_vec.size() == p) {
        Eigen::MatrixXd S_inv = S.ldlt().solve(Ip);
        pA_orig[des_i] = S_inv;
        pA_prop[des_i] = S_inv;
        cond_exp.row(des_i) = traits.row(i);
        pY_orig.row(des_i) = (S_inv * traits.row(i).transpose()).transpose();
        pY_prop.row(des_i) = pY_orig.row(des_i);
      } else {
        // Extract submatrix for observed traits
        int n_obs = obs_vec.size();
        Eigen::MatrixXd S_sub(n_obs, n_obs);
        Eigen::VectorXd y_obs(n_obs);
        for (int k = 0; k < n_obs; k++) {
          y_obs(k) = traits(i, obs_vec[k]);
          for (int l = 0; l < n_obs; l++) {
            S_sub(k, l) = S(obs_vec[k], obs_vec[l]);
          }
        }
        Eigen::MatrixXd S_sub_inv = S_sub.ldlt().solve(Eigen::MatrixXd::Identity(n_obs, n_obs));
        Eigen::VectorXd pY_sub = S_sub_inv * y_obs;

        for (int k = 0; k < n_obs; k++) {
          for (int l = 0; l < n_obs; l++) {
            pA_orig[des_i](obs_vec[k], obs_vec[l]) = S_sub_inv(k, l);
          }
          cond_exp(des_i, obs_vec[k]) = traits(i, obs_vec[k]);
          pY_orig(des_i, obs_vec[k]) = pY_sub(k);
        }
        pA_prop[des_i] = pA_orig[des_i];
        pY_prop.row(des_i) = pY_orig.row(des_i);
      }

      pA_orig[anc_i] = pA_orig[anc_i] + pA_prop[des_i];
      pA_prop[anc_i] = pA_orig[anc_i];
      pY_orig.row(anc_i) = pY_orig.row(anc_i) + pY_prop.row(des_i);
      pY_prop.row(anc_i) = pY_orig.row(anc_i);
    }

    // Forward pass - tree edges
    for (int i = 0; i < n_edges; i++) {
      int anc_i = edge(i, 0) - 1 + N;
      int des_i = edge(i, 1) - 1 + N;
      double t_e = edge_length(i);

      Eigen::MatrixXd Sigma_e = t_e * R;

      if (pA_orig[des_i].norm() > 0) {
        Eigen::MatrixXd pA_inv = pA_orig[des_i].ldlt().solve(Ip);
        cond_exp.row(des_i) = (pA_inv * pY_orig.row(des_i).transpose()).transpose();

        // itpa = I + Sigma_e * P is NOT symmetric (product of two SPD matrices
        // is generally not symmetric), so use partialPivLu instead of ldlt
        Eigen::MatrixXd itpa = Ip + Sigma_e * pA_orig[des_i];
        Eigen::MatrixXd itpainv = itpa.partialPivLu().solve(Ip);

        pA_prop[des_i] = pA_orig[des_i] * itpainv;
        pY_prop.row(des_i) = (itpainv.transpose() * pY_orig.row(des_i).transpose()).transpose();
      }

      pA_orig[anc_i] = pA_orig[anc_i] + pA_prop[des_i];
      pA_prop[anc_i] = pA_orig[anc_i];
      pY_orig.row(anc_i) = pY_orig.row(anc_i) + pY_prop.row(des_i);
      pY_prop.row(anc_i) = pY_orig.row(anc_i);
    }

    // Root
    int root_idx = n_tips + N;
    if (pA_orig[root_idx].norm() > 0) {
      Eigen::MatrixXd root_var = pA_orig[root_idx].ldlt().solve(Ip);
      cond_exp.row(root_idx) = (root_var * pY_orig.row(root_idx).transpose()).transpose();
      vars[root_idx] = root_var;
    }
    exps.row(root_idx) = mu.transpose();

    // Backward pass - tree edges
    for (int i = n_edges - 1; i >= 0; i--) {
      int anc_i = edge(i, 0) - 1 + N;
      int des_i = edge(i, 1) - 1 + N;
      double t_e = edge_length(i);

      Eigen::MatrixXd Sigma_e = t_e * R;

      if (pA_orig[des_i].norm() > 0) {
        // itpa = I + Sigma_e * P is NOT symmetric, use partialPivLu
        Eigen::MatrixXd itpa = Ip + Sigma_e * pA_orig[des_i];
        Eigen::MatrixXd itpainv = itpa.partialPivLu().solve(Ip);

        covars[des_i] = (itpainv * vars[anc_i]).transpose();

        Eigen::MatrixXd new_p = pA_orig[des_i] * itpainv;
        exps.row(des_i) = (itpainv * exps.row(anc_i).transpose() +
                          Sigma_e * new_p * cond_exp.row(des_i).transpose()).transpose();

        vars[des_i] = itpainv * Sigma_e + itpainv * covars[des_i];
      } else {
        exps.row(des_i) = exps.row(anc_i);
        covars[des_i] = vars[anc_i];
        vars[des_i] = vars[anc_i] + Sigma_e;
      }
    }

    // Backward pass - observations
    for (int i = 0; i < N; i++) {
      int anc_i = node_obs(i) + N;
      int des_i = i;

      std::vector<int> obs_vec;
      for (int j = 0; j < p; j++) {
        if (!std::isnan(traits(i, j))) {
          obs_vec.push_back(j);
        }
      }

      if ((int)obs_vec.size() == p) {
        exps.row(des_i) = traits.row(i);
        vars[des_i].setZero();
        covars[des_i].setZero();
      } else if (!obs_vec.empty()) {
        for (int k = 0; k < (int)obs_vec.size(); k++) {
          exps(des_i, obs_vec[k]) = traits(i, obs_vec[k]);
        }
        for (int j = 0; j < p; j++) {
          if (std::isnan(traits(i, j))) {
            exps(des_i, j) = exps(anc_i, j);
          }
        }
      } else {
        exps.row(des_i) = exps.row(anc_i);
        vars[des_i] = S + vars[anc_i];
        covars[des_i] = vars[anc_i];
      }
    }

    //-----------------------------------------------------------------------
    // M-STEP
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
      Eigen::VectorXd mu_new = Eigen::VectorXd::Zero(p);
      for (size_t i = 0; i < root_edges.size(); i++) {
        int e = root_edges[i];
        double w = 1.0 / edge_length(e);
        int child_idx = edge(e, 1) - 1 + N;
        mu_new = mu_new + w * exps.row(child_idx).transpose();
        weight_sum += w;
      }
      mu = mu_new / weight_sum;
    }

    // Update R
    Eigen::MatrixXd sum_exp = Eigen::MatrixXd::Zero(p, p);
    Eigen::MatrixXd sum_var = Eigen::MatrixXd::Zero(p, p);

    for (int e = 0; e < n_edges; e++) {
      int parent_node = edge(e, 0) - 1;
      int child_node = edge(e, 1) - 1;
      double t_e = edge_length(e);

      int parent_idx = parent_node + N;
      int child_idx = child_node + N;

      Eigen::RowVectorXd diff_exp = exps.row(child_idx) - exps.row(parent_idx);
      sum_exp = sum_exp + diff_exp.transpose() * diff_exp / t_e;

      Eigen::MatrixXd var_child = vars[child_idx];
      Eigen::MatrixXd var_parent = vars[parent_idx];
      Eigen::MatrixXd cov_cp = covars[child_idx];

      Eigen::MatrixXd var_diff = var_child + var_parent - cov_cp - cov_cp.transpose();
      sum_var = sum_var + var_diff / t_e;
    }

    R = (sum_exp + sum_var) / (m - 1);
    R = (R + R.transpose()) / 2.0;

    // Update S
    Eigen::MatrixXd sum_exp_ind = Eigen::MatrixXd::Zero(p, p);
    Eigen::MatrixXd sum_var_ind = Eigen::MatrixXd::Zero(p, p);

    for (int i = 0; i < N; i++) {
      int species_idx_i = node_obs(i) + N;

      Eigen::RowVectorXd diff_exp = exps.row(species_idx_i) - exps.row(i);
      sum_exp_ind = sum_exp_ind + diff_exp.transpose() * diff_exp;

      Eigen::MatrixXd var_species = vars[species_idx_i];
      Eigen::MatrixXd var_obs = vars[i];
      Eigen::MatrixXd cov_so = covars[i];

      Eigen::MatrixXd var_diff = var_species + var_obs - cov_so - cov_so.transpose();
      sum_var_ind = sum_var_ind + var_diff;
    }

    S = (sum_exp_ind + sum_var_ind) / N;
    S = (S + S.transpose()) / 2.0;

    // Ensure positive definiteness
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig_R(R);
    Eigen::VectorXd eigval_R = eig_R.eigenvalues();
    for (int i = 0; i < p; i++) {
      if (eigval_R(i) <= 1e-10) eigval_R(i) = 1e-10;
    }
    R = eig_R.eigenvectors() * eigval_R.asDiagonal() * eig_R.eigenvectors().transpose();

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig_S(S);
    Eigen::VectorXd eigval_S = eig_S.eigenvalues();
    for (int i = 0; i < p; i++) {
      if (eigval_S(i) <= 1e-10) eigval_S(i) = 1e-10;
    }
    S = eig_S.eigenvectors() * eigval_S.asDiagonal() * eig_S.eigenvectors().transpose();

    // Compute log-likelihood
    logL = loglik_bastide_eigen(R, S, edge, edge_length, traits, species_idx, n_tips, n_nodes);

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
