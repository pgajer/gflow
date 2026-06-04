#!/usr/bin/env Rscript

repo_root <- normalizePath(getwd(), mustWork = TRUE)
if (!file.exists(file.path(repo_root, "DESCRIPTION"))) {
  stop("Run from the gflow repository root.")
}

out_dir <- file.path(repo_root, "split_audit", "sa3_compiled_code")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

target_methods <- data.frame(
  method_cluster = c(
    "kernel_local_polynomial",
    "malps",
    "lpl_tf",
    "slpl_tf",
    "ssrhe"
  ),
  r_files = c(
    "R/kernel_local_polynomial_cv.R",
    "R/malps.R",
    "R/lpl_tf.R",
    "R/slpl_tf.R",
    "R/ssrhe_hessian_energy.R"
  ),
  public_entrypoints = c(
    "kernel.local.polynomial.cv;predict.kernel.local.polynomial.cv;print.kernel.local.polynomial.cv",
    "fit.malps;refit.malps;predict.malps;bootstrap.malps;malps.gcv;malps.smoother.matrix",
    "lpl.tf.operator;fit.lpl.tf;refit.lpl.tf;predict.lpl_tf;print.lpl_tf_operator",
    "slpl.tf.operator;fit.slpl.tf;refit.slpl.tf;predict.slpl_tf;print.slpl_tf_operator",
    "ssrhe.hessian.operator;fit.ssrhe.hessian.regression;fit.ssrhe.hessian.regression.cv;fit.ssrhe.hessian.regression.gcv;fit.ssrhe.hessian.l1.regression"
  ),
  stringsAsFactors = FALSE
)
write.csv(target_methods, file.path(out_dir, "sa3_target_method_cluster.csv"), row.names = FALSE)

direct_native <- data.frame(
  method_cluster = c(
    "kernel_local_polynomial",
    "kernel_local_polynomial",
    "kernel_local_polynomial",
    "lpl_tf",
    "ssrhe"
  ),
  r_file = c(
    "R/kernel_local_polynomial_cv.R",
    "R/kernel_local_polynomial_cv.R",
    "R/kernel_local_polynomial_cv.R",
    "R/lpl_tf.R",
    "R/ssrhe_hessian_energy.R"
  ),
  r_symbol = c(
    "rcpp_kernel_local_polynomial_cv_coordinates",
    "rcpp_kernel_local_polynomial_predict_coordinates",
    "rcpp_local_pca_chart",
    "rcpp_local_pca_chart",
    "S_ssrhe_hessian_operator"
  ),
  native_symbol = c(
    "_gflow_rcpp_kernel_local_polynomial_cv_coordinates",
    "_gflow_rcpp_kernel_local_polynomial_predict_coordinates",
    "_gflow_rcpp_local_pca_chart",
    "_gflow_rcpp_local_pca_chart",
    "S_ssrhe_hessian_operator"
  ),
  call_type = c("Rcpp wrapper", "Rcpp wrapper", "Rcpp wrapper", "Rcpp wrapper", "direct .Call"),
  stringsAsFactors = FALSE
)
write.csv(direct_native, file.path(out_dir, "sa3_target_direct_native_calls.csv"), row.names = FALSE)

compiled_files <- data.frame(
  compiled_file = c(
    "src/kernel_local_polynomial_cv_rcpp.cpp",
    "src/local_pca_charts.cpp",
    "src/local_pca_charts.hpp",
    "src/local_pca_charts_rcpp.cpp",
    "src/ssrhe_hessian_energy.cpp",
    "src/ssrhe_hessian_energy_r.h",
    "src/RcppExports.cpp",
    "src/init.c"
  ),
  used_by_cluster = c(
    "kernel_local_polynomial",
    "kernel_local_polynomial;lpl_tf;ssrhe",
    "kernel_local_polynomial;lpl_tf;ssrhe",
    "kernel_local_polynomial;lpl_tf",
    "ssrhe",
    "ssrhe",
    "kernel_local_polynomial;lpl_tf;shared_generated",
    "kernel_local_polynomial;lpl_tf;ssrhe;shared_generated"
  ),
  main_external_or_internal_dependencies = c(
    "Rcpp;R_ext/Lapack;ANN",
    "Rcpp;Eigen",
    "Eigen",
    "Rcpp;Eigen;local_pca_charts.hpp",
    "Rcpp;Eigen;local_pca_charts.hpp",
    "Rinternals",
    "Rcpp;generated",
    "R;Rinternals;R_ext/Rdynload;many package headers"
  ),
  current_sa3_recommendation = c(
    "Move with geosmooth cluster unless kernel smoother remains in gflow temporarily.",
    "Keep in gflow as shared geometry service for first split; expose stable R-level/API contract to geosmooth.",
    "Keep in gflow with local_pca_charts.cpp for first split.",
    "Keep in gflow for first split; export/import stable wrapper, or duplicate only if geosmooth must be independent.",
    "Move with geosmooth cluster if SSRHE moves; it is response-smoother code using shared local PCA.",
    "Move with SSRHE implementation.",
    "Regenerate separately in each package after split; do not move by hand.",
    "Regenerate/hand-maintain separately in each package after split; do not move by hand."
  ),
  stringsAsFactors = FALSE
)
write.csv(compiled_files, file.path(out_dir, "sa3_target_compiled_files.csv"), row.names = FALSE)

native_registration <- data.frame(
  native_symbol = c(
    "_gflow_rcpp_kernel_local_polynomial_cv_coordinates",
    "_gflow_rcpp_kernel_local_polynomial_predict_coordinates",
    "_gflow_rcpp_local_pca_chart",
    "S_ssrhe_hessian_operator"
  ),
  registration_file = "src/init.c",
  generated_or_manual = c("Rcpp generated plus manual init registration",
                          "Rcpp generated plus manual init registration",
                          "Rcpp generated plus manual init registration",
                          "manual CallMethods registration"),
  arity = c(6L, 6L, 9L, 13L),
  future_registration_action = c(
    "Register in geosmooth if kernel C++ backend moves.",
    "Register in geosmooth if kernel C++ backend moves.",
    "Do not move initially if local PCA remains in gflow; provide an exported/importable wrapper.",
    "Register in geosmooth if SSRHE moves."
  ),
  stringsAsFactors = FALSE
)
write.csv(native_registration, file.path(out_dir, "sa3_target_native_registration.csv"), row.names = FALSE)

later_geosmooth_native <- data.frame(
  r_file = c(
    "R/harmonic_extension.R",
    "R/harmonic_extension.R",
    "R/harmonic_smoother.R",
    "R/harmonic_smoother.R",
    "R/mean_shift_smoother.R",
    "R/metric_graph_lowpass.R",
    "R/metric_graph_lowpass.R",
    "R/riem_dcx_regression.R",
    "R/ulogit.R"
  ),
  native_symbol = c(
    "S_compute_harmonic_extension",
    "S_select_max_density_trajectory",
    "S_perform_harmonic_smoothing",
    "S_harmonic_smoother",
    "S_mean_shift_data_smoother*",
    "S_metric_graph_lowpass_operator",
    "S_metric_graph_lowpass_spectrum",
    "S_fit_rdgraph_regression;S_compute_hop_extremp_radii_batch",
    "S_ulogit;S_eigen_ulogit"
  ),
  sa3_status = c(
    "Outside first cluster move; audit later.",
    "Outside first cluster move; audit later.",
    "Outside first cluster move; audit later.",
    "Outside first cluster move; audit later.",
    "Outside first cluster move; consider gflowx or geosmooth later.",
    "Outside first cluster move; graph spectral smoother needs separate decision.",
    "Outside first cluster move; graph spectral smoother needs separate decision.",
    "Outside first cluster move; RDGraph is larger than first cluster.",
    "Outside first cluster move; likely gflowx or legacy geosmooth."
  ),
  stringsAsFactors = FALSE
)
write.csv(later_geosmooth_native, file.path(out_dir, "sa3_later_geosmooth_native_symbols.csv"), row.names = FALSE)

cat("SA3 compiled-code audit tables written to ", out_dir, "\n", sep = "")
for (f in list.files(out_dir, full.names = FALSE)) {
  cat(" - ", f, "\n", sep = "")
}
