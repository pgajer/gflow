fit.knn.riem.graph.regression <- function(
                                          X,                             # Feature matrix
                                          y,                             # Response vector
                                          k = NULL,                      # Auto-select if NULL
                                          method = c("counting", "adaptive"),
                                          t = NULL,                      # Auto-select based on spectral gap
                                          beta = NULL,                   # Auto-select as 0.1/t
                                          gamma = 1.0,                   # Default Cauchy kernel
                                          n.eigenpairs = 200,
                                          filter = c("heat_kernel", "tikhonov"),
                                          epsilon.y = 1e-4,
                                          epsilon.rho = 1e-4,
                                          max.iter = 50
                                          ) {
    ## ... parameter validation and auto-selection ...

    result <- .Call(
        "S_fit_knn_riem_graph_regression",
        X,
        as.double(y),
        as.integer(k),
        params,
        PACKAGE = "gflow"
    )

    class(result) <- c("knn_riem_fit", "riem_fit")
    result
}
