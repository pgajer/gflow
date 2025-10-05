#include "riem_dcx.hpp"
#include <R.h>

// ============================================================
// INTERNAL HELPERS (anonymous namespace)
// ============================================================

namespace {

double apply_filter_function(
    double lambda,
    double eta,
    filter_type_t filter_type
) {
    // ... implementation ...
}

} // anonymous namespace

// ============================================================
// INITIALIZATION HELPERS (private member functions)
// ============================================================

void riem_dcx_t::initialize_reference_measure(
    const std::vector<std::vector<index_t>>& knn_neighbors,
    const std::vector<std::vector<double>>& knn_distances,
    index_t k,
    bool use_counting_measure
) {
    // ... implementation ...
}

void riem_dcx_t::compute_initial_densities(
    const std::vector<std::unordered_set<index_t>>& neighbor_sets
) {
    // ... implementation ...
}

// ... other initialization helpers ...

// ============================================================
// ITERATION HELPERS (public/private member functions)
// ============================================================

vec_t riem_dcx_t::apply_damped_heat_diffusion(
    const vec_t& rho_current,
    double t,
    double beta
) {
    // ... implementation ...
}

void riem_dcx_t::update_edge_densities_from_vertices() {
    // ... implementation ...
}

void riem_dcx_t::apply_response_coherence_modulation(
    const vec_t& y_hat,
    double gamma
) {
    // ... implementation ...
}

gcv_result_t riem_dcx_t::smooth_response_via_spectral_filter(
    const vec_t& y,
    int n_eigenpairs,
    filter_type_t filter_type
) {
    // ... implementation ...
}

convergence_status_t riem_dcx_t::check_convergence(
    const vec_t& y_hat_prev,
    const vec_t& y_hat_curr,
    const std::vector<vec_t>& rho_prev,
    const std::vector<vec_t>& rho_curr,
    double epsilon_y,
    double epsilon_rho,
    int iteration,
    int max_iterations
) {
    // ... implementation ...
}

// ============================================================
// MAIN REGRESSION METHOD
// ============================================================

void riem_dcx_t::fit_knn_riem_graph_regression(
    const spmat_t& X,
    const vec_t& y,
    index_t k,
    bool use_counting_measure,
    double density_normalization,
    double t_diffusion,
    double beta_damping,
    double gamma_modulation,
    int n_eigenpairs,
    filter_type_t filter_type,
    double epsilon_y,
    double epsilon_rho,
    int max_iterations
) {
    // ================================================================
    // PART I: INITIALIZATION
    // ================================================================

    // Phase 1: Build 1-skeleton geometry
    // ... (similar to current build_nerve_from_knn but focused on 1-skeleton)

    // Phase 2: Initialize reference measure
    // ... calls this->initialize_reference_measure()

    // Phase 3: Compute initial densities
    // ... calls this->compute_initial_densities()

    // Phase 4: Build initial metric
    // ... calls this->initialize_metric_from_density()

    // Phase 5: Assemble initial Laplacian
    // ... calls this->assemble_operators()

    // Phase 6: Initial response smoothing
    auto gcv_result = this->smooth_response_via_spectral_filter(
        y, n_eigenpairs, filter_type
    );
    vec_t y_hat_curr = gcv_result.y_hat;
    sig.y_hat_hist.clear();
    sig.y_hat_hist.push_back(y_hat_curr);

    // ================================================================
    // PART II: ITERATIVE REFINEMENT
    // ================================================================

    vec_t y_hat_prev;
    std::vector<vec_t> rho_prev;

    for (int iter = 1; iter <= max_iterations; ++iter) {
        y_hat_prev = y_hat_curr;
        rho_prev = rho.rho;

        // Step 1: Density diffusion
        rho.rho[0] = this->apply_damped_heat_diffusion(
            rho.rho[0], t_diffusion, beta_damping
        );

        // Step 2: Edge density update
        this->update_edge_densities_from_vertices();

        // Step 3: Response-coherence modulation
        this->apply_response_coherence_modulation(y_hat_curr, gamma_modulation);

        // Step 4: Metric update
        this->update_metric_from_density();

        // Step 5: Laplacian reassembly
        this->assemble_operators();

        // Step 6: Response smoothing
        gcv_result = this->smooth_response_via_spectral_filter(
            y, n_eigenpairs, filter_type
        );
        y_hat_curr = gcv_result.y_hat;
        sig.y_hat_hist.push_back(y_hat_curr);

        // Step 7: Convergence check
        auto status = this->check_convergence(
            y_hat_prev, y_hat_curr,
            rho_prev, rho.rho,
            epsilon_y, epsilon_rho,
            iter, max_iterations
        );

        Rprintf("Iteration %d: response_change=%.6f, density_change=%.6f\n",
               iter, status.response_change, status.max_density_change);

        if (status.converged) {
            Rprintf("%s\n", status.message.c_str());
            break;
        }
    }

    // ================================================================
    // PART III: FINALIZATION
    // ================================================================

    // Store original response
    sig.y = y;

    // Final state already in place
}
