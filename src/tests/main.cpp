// src/tests/main.cpp
#include "test_helpers.hpp"
#include "riem_dcx.hpp"
#include <R.h>
#include <Rinternals.h>
#include <numeric>
#include <algorithm>

// Test fixture with friend access
class riem_dcx_test_fixture {
public:
    static void setup_neighbors(riem_dcx_t& dcx,
                                const std::vector<std::unordered_set<index_t>>& neighbors) {
        dcx.neighbor_sets = neighbors;
    }

    static void call_initialize_reference_measure(
        riem_dcx_t& dcx,
        const std::vector<std::vector<index_t>>& knn_neighbors,
        const std::vector<std::vector<double>>& knn_distances,
        bool use_counting,
        double normalization) {
        dcx.initialize_reference_measure(knn_neighbors, knn_distances,
                                        use_counting, normalization);
    }

    static void call_compute_edge_mass_matrix(riem_dcx_t& dcx) {
        dcx.compute_edge_mass_matrix();
    }

    static const std::vector<double>& get_reference_measure(const riem_dcx_t& dcx) {
        return dcx.reference_measure;
    }
};

void init_r() {
    static bool initialized = false;
    if (!initialized) {
        initialized = true;
    }
}

TEST(reference_measure_counting) {
    init_r();
    riem_dcx_t dcx;
    dcx.init_dims(1, {10, 0});

    // Setup kNN data
    std::vector<std::vector<index_t>> knn_neighbors(10);
    std::vector<std::vector<double>> knn_distances(10);

    for (size_t i = 0; i < 10; ++i) {
        knn_neighbors[i] = {static_cast<index_t>(i)};
        knn_distances[i] = {0.0};
    }

    riem_dcx_test_fixture::call_initialize_reference_measure(
        dcx, knn_neighbors, knn_distances, true, 1.0
    );

    const auto& ref_measure = riem_dcx_test_fixture::get_reference_measure(dcx);

    double sum = std::accumulate(ref_measure.begin(), ref_measure.end(), 0.0);
    ASSERT_NEAR(sum, 10.0, 1e-10);

    bool all_positive = std::all_of(ref_measure.begin(), ref_measure.end(),
                                    [](double v) { return v > 0; });
    ASSERT_TRUE(all_positive);
}

TEST(mass_matrix_symmetry) {
    init_r();
    riem_dcx_t dcx;
    dcx.init_dims(2, {3, 3, 1});  // 3 vertices, 3 edges, 1 triangle

    // Build 3-vertex triangle
    dcx.S[0].simplex_verts = {{0}, {1}, {2}};
    dcx.S[1].simplex_verts = {{0,1}, {0,2}, {1,2}};
    dcx.S[2].simplex_verts = {{0,1,2}};

    // Setup densities
    dcx.rho.rho[0] = vec_t::Ones(3);
    dcx.rho.rho[1] = vec_t::Ones(3);

    // Setup neighbors
    std::vector<std::unordered_set<index_t>> neighbors = {
        {0,1,2}, {0,1,2}, {0,1,2}
    };
    riem_dcx_test_fixture::setup_neighbors(dcx, neighbors);

    riem_dcx_test_fixture::call_compute_edge_mass_matrix(dcx);

    // Check symmetry
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            ASSERT_NEAR(dcx.g.M[1].coeff(i,j),
                       dcx.g.M[1].coeff(j,i), 1e-12);
        }
    }
}

int main() {
    return run_all_tests();
}
