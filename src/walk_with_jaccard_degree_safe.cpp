#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector walk_with_jaccard_degree_safe_cpp(
    List adj_list,
    int start_node,
    int walk_depth = 20,
    double alpha = 0.3,
    double beta = 2.0,
    double eps = 1e-3
) {

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