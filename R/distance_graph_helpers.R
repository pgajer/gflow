.dgh_as_distance_matrix <- function(D, name = "D") {
    D <- as.matrix(D)
    if (!is.numeric(D)) {
        stop(name, " must be a numeric matrix.")
    }
    if (nrow(D) != ncol(D)) {
        stop(name, " must be square.")
    }
    if (nrow(D) < 2L) {
        stop(name, " must have at least two rows.")
    }
    if (any(!is.finite(D))) {
        stop(name, " cannot contain NA/NaN/Inf.")
    }
    if (any(D < 0)) {
        stop(name, " cannot contain negative distances.")
    }
    if (max(abs(D - t(D))) > sqrt(.Machine$double.eps)) {
        stop(name, " must be symmetric.")
    }
    if (any(abs(diag(D)) > sqrt(.Machine$double.eps))) {
        stop("diag(", name, ") must be zero.")
    }
    storage.mode(D) <- "double"
    D
}

.dgh_validate_k <- function(k, n) {
    if (!is.numeric(k) || length(k) != 1L || !is.finite(k) ||
        k != floor(k) || k < 1L || k >= n) {
        stop("k must be a positive integer smaller than nrow(D).")
    }
    as.integer(k)
}

.dgh_knn_sets_from_dist <- function(D, k) {
    D <- .dgh_as_distance_matrix(D)
    n <- nrow(D)
    k <- .dgh_validate_k(k, n)

    lapply(seq_len(n), function(i) {
        ord <- order(D[i, ], seq_len(n), decreasing = FALSE)
        ord <- ord[ord != i]
        ord[seq_len(k)]
    })
}

.dgh_edge_matrix_to_lists <- function(edge.matrix, weights, n.vertices) {
    adj.list <- vector("list", n.vertices)
    weight.list <- vector("list", n.vertices)
    for (i in seq_len(n.vertices)) {
        adj.list[[i]] <- integer(0)
        weight.list[[i]] <- numeric(0)
    }
    if (!is.null(edge.matrix) && nrow(edge.matrix) > 0L) {
        for (r in seq_len(nrow(edge.matrix))) {
            i <- edge.matrix[r, 1L]
            j <- edge.matrix[r, 2L]
            w <- weights[r]
            adj.list[[i]] <- c(adj.list[[i]], j)
            adj.list[[j]] <- c(adj.list[[j]], i)
            weight.list[[i]] <- c(weight.list[[i]], w)
            weight.list[[j]] <- c(weight.list[[j]], w)
        }
        for (i in seq_len(n.vertices)) {
            if (length(adj.list[[i]]) > 1L) {
                ord <- order(adj.list[[i]])
                adj.list[[i]] <- adj.list[[i]][ord]
                weight.list[[i]] <- weight.list[[i]][ord]
            }
        }
    }
    list(adj.list = adj.list, weight.list = weight.list)
}

.dgh_graph_result <- function(edge.matrix,
                              weights,
                              n.vertices,
                              type,
                              k,
                              extra = list()) {
    if (is.null(edge.matrix) || length(edge.matrix) == 0L) {
        edge.matrix <- matrix(integer(0), ncol = 2L)
    }
    edge.matrix <- matrix(as.integer(edge.matrix), ncol = 2L)
    weights <- as.numeric(weights)
    lists <- .dgh_edge_matrix_to_lists(edge.matrix, weights, n.vertices)
    out <- c(
        list(
            type = type,
            k = as.integer(k),
            n.vertices = as.integer(n.vertices),
            edge.matrix = edge.matrix,
            edge.weight = weights,
            adj.list = lists$adj.list,
            weight.list = lists$weight.list,
            n.edges = nrow(edge.matrix)
        ),
        extra
    )
    class(out) <- c(paste0(type, "_dist_graph"), "list")
    out
}

.dgh_mknn_graph_from_dist <- function(D, k) {
    D <- .dgh_as_distance_matrix(D)
    n <- nrow(D)
    k <- .dgh_validate_k(k, n)
    nn <- .dgh_knn_sets_from_dist(D, k)

    edges <- list()
    m <- 0L
    for (i in seq_len(n - 1L)) {
        for (j in seq.int(i + 1L, n)) {
            if (j %in% nn[[i]] && i %in% nn[[j]]) {
                m <- m + 1L
                edges[[m]] <- c(i, j)
            }
        }
    }

    edge.matrix <- if (m == 0L) matrix(integer(0), ncol = 2L) else do.call(rbind, edges)
    weights <- if (nrow(edge.matrix) == 0L) numeric(0) else D[edge.matrix]
    .dgh_graph_result(edge.matrix, weights, n, "mknn", k,
                      extra = list(knn.sets = nn, weight.type = "distance"))
}

.dgh_iknn_graph_from_dist <- function(D, k) {
    D <- .dgh_as_distance_matrix(D)
    n <- nrow(D)
    k <- .dgh_validate_k(k, n)
    nn <- .dgh_knn_sets_from_dist(D, k)

    edges <- list()
    direct.distance <- numeric()
    shared.detour <- numeric()
    intersection.size <- integer()
    m <- 0L
    for (i in seq_len(n - 1L)) {
        for (j in seq.int(i + 1L, n)) {
            common <- intersect(nn[[i]], nn[[j]])
            if (length(common) > 0L) {
                m <- m + 1L
                edges[[m]] <- c(i, j)
                direct.distance[m] <- D[i, j]
                shared.detour[m] <- min(D[i, common] + D[j, common])
                intersection.size[m] <- length(common)
            }
        }
    }

    edge.matrix <- if (m == 0L) matrix(integer(0), ncol = 2L) else do.call(rbind, edges)
    .dgh_graph_result(
        edge.matrix,
        shared.detour,
        n,
        "iknn",
        k,
        extra = list(
            knn.sets = nn,
            weight.type = "shared_neighbor_detour",
            direct.distance = direct.distance,
            intersection.size = intersection.size
        )
    )
}

.dgh_phate_support_graph_from_dist <- function(D,
                                               k,
                                               decay = 40,
                                               threshold = 1e-4,
                                               bandwidth.scale = 1,
                                               symmetrize = c("mean", "max")) {
    D <- .dgh_as_distance_matrix(D)
    n <- nrow(D)
    k <- .dgh_validate_k(k, n)
    symmetrize <- match.arg(symmetrize)

    if (!is.numeric(decay) || length(decay) != 1L ||
        !is.finite(decay) || decay <= 0) {
        stop("decay must be a single positive finite number.")
    }
    if (!is.numeric(threshold) || length(threshold) != 1L ||
        !is.finite(threshold) || threshold < 0 || threshold >= 1) {
        stop("threshold must be a single finite number in [0, 1).")
    }
    if (!is.numeric(bandwidth.scale) || length(bandwidth.scale) != 1L ||
        !is.finite(bandwidth.scale) || bandwidth.scale <= 0) {
        stop("bandwidth.scale must be a single positive finite number.")
    }

    nn <- .dgh_knn_sets_from_dist(D, k)
    sigma <- numeric(n)
    for (i in seq_len(n)) {
        row <- D[i, ]
        row[i] <- Inf
        sigma[i] <- sort(row, partial = k)[k] * bandwidth.scale
    }
    sigma <- pmax(sigma, .Machine$double.eps)

    A <- matrix(0, nrow = n, ncol = n)
    diag(A) <- 1
    for (i in seq_len(n)) {
        for (j in seq_len(n)) {
            if (i == j) {
                next
            }
            val <- exp(-((D[i, j] / sigma[i])^decay))
            if (is.finite(val) && val >= threshold) {
                A[i, j] <- val
            }
        }
    }

    K <- if (symmetrize == "max") {
        pmax(A, t(A))
    } else {
        0.5 * (A + t(A))
    }
    row.sum <- rowSums(K)
    P <- K / row.sum
    P[!is.finite(P)] <- 0

    idx <- which(upper.tri(K) & K > 0, arr.ind = TRUE)
    edge.matrix <- if (nrow(idx) == 0L) matrix(integer(0), ncol = 2L) else idx
    weights <- if (nrow(edge.matrix) == 0L) numeric(0) else K[edge.matrix]

    .dgh_graph_result(
        edge.matrix,
        weights,
        n,
        "phate_support",
        k,
        extra = list(
            knn.sets = nn,
            weight.type = "symmetrized_affinity",
            decay = as.numeric(decay),
            threshold = as.numeric(threshold),
            bandwidth.scale = as.numeric(bandwidth.scale),
            symmetrize = symmetrize,
            sigma = sigma,
            A = A,
            K = K,
            P = P,
            n.retained.directed = rowSums(A > 0)
        )
    )
}
