#include <Rcpp.h>
#include <vector>
#include <unordered_set>
#include <algorithm>
#include <cmath>

using namespace Rcpp;

// [[Rcpp::export]]
IntegerMatrix walk_with_jaccard_degree_safe_cpp(
    List adj_list,
    int start_node,
    int walk_depth = 20,
    int n_iter = 50,
    double eps = 1e-3
) {

    const int n = adj_list.size();

    if (n == 0)
        stop("adj_list must not be empty.");

    if (start_node < 1 || start_node > n)
        stop("Invalid start_node.");

    if (walk_depth < 0 || n_iter < 0)
        stop("walk_depth and n_iter must be nonnegative.");

    if (!std::isfinite(eps) || eps <= 0.0)
        stop("eps must be finite and positive.");

    const int source = start_node - 1;

    // =====================================================
    // 1. Convert R adjacency lists to C++ vectors
    // =====================================================

    std::vector<std::vector<int>> adj(n);

    std::vector<int> degree(n);

    for (int i = 0; i < n; ++i) {

        IntegerVector neighbors = adj_list[i];

        degree[i] = neighbors.size();

        adj[i].reserve(neighbors.size());

        for (int v : neighbors) {

            if (v < 1 || v > n)
                stop("Invalid node ID in adjacency list.");

            adj[i].push_back(v - 1);
        }
    }

    // =====================================================
    // 2. Cache source-neighbor membership
    // =====================================================

    std::unordered_set<int> source_neighbors(
        adj[source].begin(),
        adj[source].end()
    );

    const int source_degree = degree[source];

    // =====================================================
    // 3. Jaccard cache
    // =====================================================

    std::vector<double> jaccard_cache(n, 0.0);

    std::vector<unsigned char> jaccard_ready(n, 0);

    auto get_jaccard = [&](int v) -> double {

        if (jaccard_ready[v])
            return jaccard_cache[v];

        int intersection = 0;

        for (int a : adj[v]) {

            if (source_neighbors.find(a) !=
                source_neighbors.end()) {

                ++intersection;
            }
        }

        const int union_size =
            degree[v] + source_degree - intersection;

        const double jaccard =
            (union_size == 0)
            ? 0.0
            : static_cast<double>(intersection) /
              static_cast<double>(union_size);

        jaccard_cache[v] = jaccard;

        jaccard_ready[v] = 1;

        return jaccard;
    };

    // =====================================================
    // 4. Transition probability cache
    // =====================================================

    std::vector<std::vector<double>> transition_cache(n);

    std::vector<unsigned char> transition_ready(n, 0);

    auto get_transition = [&](int u)
        -> const std::vector<double>& {

        if (transition_ready[u])
            return transition_cache[u];

        const int m = degree[u];

        std::vector<double>& cumulative =
            transition_cache[u];

        cumulative.resize(m);

        double total_weight = 0.0;

        // Calculate Jaccard weights once
        for (int j = 0; j < m; ++j) {

            const int v = adj[u][j];

            const double weight =
                get_jaccard(v) + eps;

            cumulative[j] = weight;

            total_weight += weight;
        }

        // Convert weights into cumulative probabilities
        if (!std::isfinite(total_weight) ||
            total_weight <= 0.0) {

            for (int j = 0; j < m; ++j) {

                cumulative[j] =
                    static_cast<double>(j + 1) / m;
            }

        } else {

            double running_sum = 0.0;

            for (int j = 0; j < m; ++j) {

                running_sum +=
                    cumulative[j] / total_weight;

                cumulative[j] = running_sum;
            }
        }

        if (m > 0)
            cumulative[m - 1] = 1.0;

        transition_ready[u] = 1;

        return transition_cache[u];
    };

    // =====================================================
    // 5. Generate random walks
    // =====================================================

    IntegerMatrix WALKS(walk_depth + 1, n_iter);

    for (int iter = 0; iter < n_iter; ++iter) {

        int current = source;

        WALKS(0, iter) = start_node;

        for (int t = 1; t <= walk_depth; ++t) {

            const int m = degree[current];

            // Preserve original isolated-node behavior
            if (m == 0) {

                current = source;

                WALKS(t, iter) = current + 1;

                continue;
            }

            // Only one neighbor: deterministic transition
            if (m == 1) {

                current = adj[current][0];

                WALKS(t, iter) = current + 1;

                continue;
            }

            // Get cached cumulative probabilities
            const std::vector<double>& cumulative =
                get_transition(current);

            // Sample from the cached distribution
            const double r = R::runif(0.0, 1.0);

            auto it = std::lower_bound(
                cumulative.begin(),
                cumulative.end(),
                r
            );

            int index = static_cast<int>(
                it - cumulative.begin()
            );

            if (index >= m)
                index = m - 1;

            // Move to selected neighbor
            current = adj[current][index];

            WALKS(t, iter) = current + 1;
        }
    }

    return WALKS;
}


#include <Rcpp.h>
#include <unordered_set>
#include <random>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerMatrix walk_with_jaccard_degree_safe_cpp_old(
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