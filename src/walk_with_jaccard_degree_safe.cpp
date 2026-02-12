#include <Rcpp.h>
#include <unordered_set>
#include <random>
using namespace Rcpp;

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
NumericVector stationary_jaccard_walk_cpp(
    List adj_list,
    int start_node,
    int max_iter = 100,
    double tol = 1e-8,
    double eps = 1e-3
) {
    int n = adj_list.size();

    // precompute start neighbors
    IntegerVector start_neighbors = adj_list[start_node - 1];

    // initial probability vector (start node)
    NumericVector x(n, 0.0);
    x[start_node - 1] = 1.0;

    NumericVector x_next(n);

    for (int iter = 0; iter < max_iter; iter++) {
        x_next = NumericVector(n, 0.0);

        for (int u = 0; u < n; u++) {
            if (x[u] < 1e-12) continue; // skip negligible probabilities

            IntegerVector neighbors = adj_list[u];
            int m = neighbors.size();

            if (m == 0) {
                // no neighbors → stay at current node
                x_next[u] += x[u];
                continue;
            }

            // compute Jaccard + eps weights
            NumericVector weights(m);
            for (int j = 0; j < m; j++) {
                int v = neighbors[j] - 1;
                IntegerVector nv = adj_list[v];

                int inter_len = 0;
                for (int a : nv) {
                    for (int b : start_neighbors) {
                        if (a == b) { inter_len++; break; }
                    }
                }
                int union_len = nv.size() + start_neighbors.size() - inter_len;
                double jaccard = (union_len == 0) ? 0.0 : (double)inter_len / union_len;

                weights[j] = jaccard + eps;
            }

            // normalize weights
            double sum_w = sum(weights);
            if (!R_finite(sum_w) || sum_w <= 0) {
                weights = NumericVector(m, 1.0 / m);
            } else {
                weights = weights / sum_w;
            }

            // propagate probabilities
            for (int j = 0; j < m; j++) {
                int v = neighbors[j] - 1;
                x_next[v] += x[u] * weights[j];
            }
        }

        // check convergence
        double diff = sum(abs(x_next - x));
        x = x_next;
        if (diff < tol) break;
    }

    return x;  // long-term visitation probabilities
}