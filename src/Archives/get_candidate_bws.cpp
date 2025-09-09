#include "get_candidate_bws.h"

/**
* @brief Computes candidate bandwidth values from composite paths using a quantile-based strategy
*
* Computes min_bw (smallest distance that ensures at least one path has min_path_size vertices),
* max_bw (largest distance in any path), and generates n_bws candidate bandwidth values
* using quantile-based subdivision of the [min_bw, max_bw] interval.
*
* @param composite_paths Vector of composite paths
* @param min_path_size Minimum required path size
* @param n_bws Number of bandwidth values to generate
*
* @return std::tuple<double, double, std::vector<double>> Returns (min_bw, max_bw, candidate_bws)
*/
std::tuple<double, double, std::vector<double>> uniform_grid_graph_t::get_candidate_bws(
   const std::vector<compose_path_t>& composite_paths,
   size_t min_path_size,
   size_t n_bws) const {

   // Collect all unique distances from composite paths
   std::set<double> all_distances;
   for (const auto& path : composite_paths) {
       all_distances.insert(path.dist_to_ref_vertex.begin(),
                          path.dist_to_ref_vertex.end());
   }

   // Convert to vector for easier manipulation
   std::vector<double> sorted_distances(all_distances.begin(), all_distances.end());

   // Find max_bw - it's the largest distance in any path
   double max_bw = sorted_distances.back();

   // Find min_bw - smallest distance that ensures at least one path has min_path_size vertices
   double min_bw = max_bw;  // Start with max_bw and work down
   for (const auto& path : composite_paths) {
       // Count how many vertices are within each distance
       for (size_t i = 0; i < path.vertices.size(); ++i) {
           double dist = path.dist_to_ref_vertex[i];

           // Count vertices within this distance
           size_t vertices_within_dist = std::count_if(
               path.dist_to_ref_vertex.begin(),
               path.dist_to_ref_vertex.end(),
               [dist](double d) { return d <= dist; }
           );

           // If we have enough vertices and this distance is smaller than current min_bw
           if (vertices_within_dist >= min_path_size && dist < min_bw) {
               min_bw = dist;
           }
       }
   }

   // Filter distances to be within [min_bw, max_bw]
   auto it_min = std::lower_bound(sorted_distances.begin(), sorted_distances.end(), min_bw);
   std::vector<double> filtered_distances(it_min, sorted_distances.end());

   // Generate candidate bandwidths using quantile strategy
   std::vector<double> candidate_bws;
   candidate_bws.reserve(n_bws);

   // Always include min_bw and max_bw
   candidate_bws.push_back(min_bw);

   // Calculate indices for internal quantiles
   if (n_bws > 2) {
       size_t n_distances = filtered_distances.size();
       for (size_t i = 1; i < n_bws - 1; ++i) {
           // Calculate index for this quantile
           size_t index = (i * (n_distances - 1)) / (n_bws - 1);
           candidate_bws.push_back(filtered_distances[index]);
       }
   }

   // Add max_bw if we haven't already
   if (candidate_bws.back() != max_bw) {
       candidate_bws.push_back(max_bw);
   }

   return {min_bw, max_bw, candidate_bws};
}
