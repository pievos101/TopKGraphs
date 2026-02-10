#include <Rcpp.h>
#include <unordered_set>
#include <random>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerMatrix walk_with_jaccard_fast(
    List adj_list,
    int start_node,
    int walk_depth = 20,
    int n_iter = 50,
    double eps = 1e-3
) {
  int n = adj_list.size();
  double alpha = 0.0;
  double beta = 0.0;

  // Precompute degrees and neighbor sets
  std::vector<int> deg(n);
  std::vector<std::unordered_set<int>> neighbors_set(n);

  for (int i = 0; i < n; i++) {
    IntegerVector nei = adj_list[i];
    deg[i] = nei.size();
    neighbors_set[i].insert(nei.begin(), nei.end());
  }

  std::unordered_set<int> start_neighbors(neighbors_set[start_node - 1].begin(),
                                          neighbors_set[start_node - 1].end());

  // Random number generator
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis(0.0, 1.0);

  IntegerMatrix WALKS(walk_depth + 1, n_iter);

  for (int iter = 0; iter < n_iter; iter++) {
    std::vector<int> walk(walk_depth + 1);
    walk[0] = start_node;
    int current = start_node;

    for (int i = 1; i <= walk_depth; i++) {

      if (dis(gen) < alpha) {  // restart
        current = start_node;
        walk[i] = current;
        continue;
      }

      auto &cur_neigh = neighbors_set[current - 1];
      int m = cur_neigh.size();
      if (m == 0) {  // no neighbors, restart
        current = start_node;
        walk[i] = current;
        continue;
      }

      std::vector<int> neigh_vec(cur_neigh.begin(), cur_neigh.end());
      std::vector<double> weights(m);

      // Compute weights
      for (int j = 0; j < m; j++) {
        int v = neigh_vec[j];
        auto &v_neigh = neighbors_set[v - 1];

        // intersection size
        int inter_len = 0;
        for (int u : v_neigh)
          if (start_neighbors.count(u)) inter_len++;

        int union_len = v_neigh.size() + start_neighbors.size() - inter_len;
        double jaccard = (union_len == 0) ? 0.0 : (double) inter_len / union_len;

        int d = std::max(deg[v - 1], 1);
        weights[j] = (1.0 / std::pow(d, beta)) * (jaccard + eps);
      }

      // Normalize weights
      double sum_w = 0.0;
      for (double w : weights) sum_w += w;
      if (sum_w <= 0.0 || !R_finite(sum_w)) sum_w = 1.0;
      for (double &w : weights) w /= sum_w;

      // Weighted sample
      std::discrete_distribution<> dd(weights.begin(), weights.end());
      current = neigh_vec[dd(gen)];

      walk[i] = current;
    }

    for (int row = 0; row <= walk_depth; row++)
      WALKS(row, iter) = walk[row];
  }

  return WALKS;
}



// [[Rcpp::export]]
IntegerMatrix walk_with_jaccard_degree_safe_cpp(
    List adj_list,
    int start_node,
    int walk_depth = 20,  
    int n_iter = 50,       // NEW: number of walks
    double eps = 1e-3
) {

  double alpha = 0.0;
  double beta = 0.0;
  int n = adj_list.size();

  // degrees
  IntegerVector deg(n);
  for (int i = 0; i < n; i++) {
    IntegerVector nei = adj_list[i];
    deg[i] = nei.size();
  }

  // matrix to store all walks: rows = walk_depth + 1, cols = n_iter
  IntegerMatrix WALKS(walk_depth + 1, n_iter);

  // Precompute start neighbors
  IntegerVector start_neighbors = adj_list[start_node - 1];

  for (int iter = 0; iter < n_iter; iter++) {
    IntegerVector walk(walk_depth + 1);
    walk[0] = start_node;
    int current = start_node;

    for (int i = 1; i <= walk_depth; i++) {

      // restart
      if (R::runif(0.0, 1.0) < alpha) {
        current = start_node;
        walk[i] = current;
        continue;
      }

      IntegerVector neighbors = adj_list[current - 1];

      // if no neighbors, restart
      if (neighbors.size() == 0) {
        current = start_node;
        walk[i] = current;
        continue;
      }

      int m = neighbors.size();
      NumericVector weights(m);

      // compute weights
      for (int j = 0; j < m; j++) {
        int v = neighbors[j];
        IntegerVector nv = adj_list[v - 1];

        // Jaccard similarity
        int inter_len = 0;
        for (int a : nv) {
          for (int b : start_neighbors) {
            if (a == b) {
              inter_len++;
              break;
            }
          }
        }

        int union_len = nv.size() + start_neighbors.size() - inter_len;
        double jaccard = (union_len == 0) ? 0.0 : (double) inter_len / union_len;

        int d = std::max(deg[v - 1], 1);
        weights[j] = (1.0 / std::pow(d, beta)) * (jaccard + eps);
      }

      // normalize weights
      double sum_w = sum(weights);
      if (!R_finite(sum_w) || sum_w <= 0) {
        weights = NumericVector(m, 1.0 / m);
      } else {
        weights = weights / sum_w;
      }

      // sample
      if (m == 1) {
        current = neighbors[0];
      } else {
        current = Rcpp::sample(neighbors, 1, false, weights)[0];
      }

      walk[i] = current;
    }

    // store this walk in the matrix
    for (int row = 0; row <= walk_depth; row++) {
      WALKS(row, iter) = walk[row];
    }
  }

  return WALKS;
}









// [[Rcpp::export]]
IntegerVector walk_with_jaccard_degree_safe_cpp_old(
    List adj_list,
    int start_node,
    int walk_depth = 20,  
    double eps = 1e-3
) {

  double alpha = 0.0;
  double beta = 0.0;
  int n = adj_list.size();

  // degrees
  IntegerVector deg(n);
  for (int i = 0; i < n; i++) {
    IntegerVector nei = adj_list[i];
    deg[i] = nei.size();
  }

  IntegerVector walk(walk_depth + 1);
  walk[0] = start_node;
  int current = start_node;

  IntegerVector start_neighbors = adj_list[start_node - 1];

  for (int i = 1; i <= walk_depth; i++) {

    // restart
    if (R::runif(0.0, 1.0) < alpha) {
      current = start_node;
      walk[i] = current;
      continue;
    }

    IntegerVector neighbors = adj_list[current - 1];

    // if no neighbors, restart
    if (neighbors.size() == 0) {
      current = start_node;
      walk[i] = current;
      continue;
    }

    int m = neighbors.size();
    NumericVector weights(m);

    // compute weights
    for (int j = 0; j < m; j++) {
      int v = neighbors[j];
      IntegerVector nv = adj_list[v - 1];

      // Jaccard similarity
      int inter_len = 0;
      for (int a : nv) {
        for (int b : start_neighbors) {
          if (a == b) {
            inter_len++;
            break;
          }
        }
      }

      int union_len = nv.size() + start_neighbors.size() - inter_len;
      double jaccard = (union_len == 0) ? 0.0 : (double) inter_len / union_len;

      int d = std::max(deg[v - 1], 1);
      weights[j] = (1.0 / std::pow(d, beta)) * (jaccard + eps);
    }

    // normalize weights
    double sum_w = sum(weights);
    if (!R_finite(sum_w) || sum_w <= 0) {
      weights = NumericVector(m, 1.0 / m);
    } else {
      weights = weights / sum_w;
    }

    // sample
    if (m == 1) {
      current = neighbors[0];
    } else {
      current = Rcpp::sample(neighbors, 1, false, weights)[0];
    }

    walk[i] = current;
  }

  return walk;
}