#ifndef GFLOW_KNN_H_
#define GFLOW_KNN_H_

#include <vector>

// Leak-proof result: flat row-major buffers inside std::vector.
struct knn_result_t {
  int n = 0;                 // number of points
  int k = 0;                 // neighbors per point
  std::vector<int> indices;        // size n*k, neighbors of i at [i*k ... i*k+k-1]
  std::vector<double> distances;   // size n*k, same layout
};

// C++ helper: compute kNN over X as vector-of-rows
knn_result_t kNN(const std::vector<std::vector<double>>& X, int k);

#endif // GFLOW_KNN_H_
