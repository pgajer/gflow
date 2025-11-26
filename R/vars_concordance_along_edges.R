## Suite of tools testing concordance of how two variable change along graph edges

#' Permutation test for edge-wise concordance (CORRECTED)
#'
#' @param delta.y Vector of y differences across edges
#' @param delta.z Vector of z differences across edges
#' @param n.perm Number of permutations (default 10000)
#' @param n.cores Number of cores for parallel processing (default: all available)
#' @return List with test statistic, p-value, and permutation distribution
test.edge.concordance <- function(delta.y, delta.z, n.perm = 10000,
                                 n.cores = parallel::detectCores() - 1) {

    require(parallel)

    ## Observed test statistics
    obs.concordant <- mean((delta.y > 0 & delta.z > 0) |
                          (delta.y < 0 & delta.z < 0))

    obs.phi <- cor(delta.y > 0, delta.z > 0)

    ## Parallel permutation
    cl <- makeCluster(n.cores)
    clusterExport(cl, c("delta.y", "delta.z"), envir = environment())

    perm.results <- parLapply(cl, 1:n.perm, function(i) {
        z.perm <- sample(delta.z)
        concordant <- mean((delta.y > 0 & z.perm > 0) |
                          (delta.y < 0 & z.perm < 0))
        phi <- cor(delta.y > 0, z.perm > 0)
        c(concordant = concordant, phi = phi)
    })

    stopCluster(cl)

    ## Extract results
    perm.matrix <- do.call(rbind, perm.results)
    perm.concordant <- perm.matrix[, "concordant"]
    perm.phi <- perm.matrix[, "phi"]

    ## CORRECTED p-values:
    ## Test whether observed is extreme relative to permutation distribution
    ## Two-sided test
    if (obs.concordant >= median(perm.concordant)) {
        ## Observed is on high end
        p.concordant <- 2 * mean(perm.concordant >= obs.concordant)
    } else {
        ## Observed is on low end
        p.concordant <- 2 * mean(perm.concordant <= obs.concordant)
    }
    p.concordant <- min(p.concordant, 1.0)  # Cap at 1

    ## Similar for phi coefficient
    if (obs.phi >= median(perm.phi)) {
        p.phi <- 2 * mean(perm.phi >= obs.phi)
    } else {
        p.phi <- 2 * mean(perm.phi <= obs.phi)
    }
    p.phi <- min(p.phi, 1.0)

    ## One-sided p-values (for reference)
    p.concordant.greater <- mean(perm.concordant >= obs.concordant)
    p.concordant.less <- mean(perm.concordant <= obs.concordant)

    ## Effect size: standardized difference from null expectation
    null.mean <- mean(perm.concordant)
    null.sd <- sd(perm.concordant)
    z.score <- (obs.concordant - null.mean) / null.sd

    list(
        obs.concordant = obs.concordant,
        obs.phi = obs.phi,
        null.mean = null.mean,
        null.sd = null.sd,
        z.score = z.score,
        p.value.two.sided = p.concordant,
        p.value.greater = p.concordant.greater,
        p.value.less = p.concordant.less,
        p.value.phi = p.phi,
        perm.concordant = perm.concordant,
        perm.phi = perm.phi,
        n.edges = length(delta.y),
        n.perm = n.perm
    )
}

## Enhanced diagnostic function
print.concordance.test <- function(result) {
    cat("\n")
    cat("=" , rep("=", 70), "\n", sep = "")
    cat("  EDGE-WISE CONCORDANCE TEST\n")
    cat("=" , rep("=", 70), "\n", sep = "")
    cat(sprintf("Number of edges: %d\n", result$n.edges))
    cat(sprintf("Number of permutations: %d\n", result$n.perm))
    cat("\n")
    cat("Observed concordance: ", sprintf("%.6f\n", result$obs.concordant))
    cat("Null mean (permuted): ", sprintf("%.6f\n", result$null.mean))
    cat("Null SD (permuted):   ", sprintf("%.6f\n", result$null.sd))
    cat("Z-score:              ", sprintf("%.4f\n", result$z.score))
    cat("\n")
    cat("P-values:\n")
    cat(sprintf("  Two-sided:         %.6f %s\n",
                result$p.value.two.sided,
                ifelse(result$p.value.two.sided < 0.001, "***",
                      ifelse(result$p.value.two.sided < 0.01, "**",
                            ifelse(result$p.value.two.sided < 0.05, "*", "")))))
    cat(sprintf("  One-sided (>):     %.6f\n", result$p.value.greater))
    cat(sprintf("  One-sided (<):     %.6f\n", result$p.value.less))
    cat(sprintf("  Phi coefficient:   %.6f (p = %.6f)\n",
                result$obs.phi, result$p.value.phi))
    cat("\n")

    ## Interpretation
    if (result$obs.concordant > result$null.mean) {
        cat("Interpretation: CONCORDANCE (same direction changes)\n")
        if (result$p.value.two.sided < 0.05) {
            cat("  Significantly more concordant than expected by chance\n")
        }
    } else if (result$obs.concordant < result$null.mean) {
        cat("Interpretation: DISCORDANCE (opposite direction changes)\n")
        if (result$p.value.two.sided < 0.05) {
            cat("  Significantly more discordant than expected by chance\n")
        }
    } else {
        cat("Interpretation: No evidence of association\n")
    }
    cat("=" , rep("=", 70), "\n", sep = "")
    cat("\n")

    invisible(result)
}

## Visualization function
plot.concordance.test <- function(result, main = "Permutation Test Results") {

    par(mfrow = c(1, 2))

    ## Plot 1: Histogram of permutation distribution
    hist(result$perm.concordant, breaks = 50,
         xlim = range(c(result$perm.concordant, result$obs.concordant, 0.5)),
         main = "Permutation Distribution",
         xlab = "Concordance Proportion",
         col = "lightgray", border = "white")
    abline(v = result$obs.concordant, col = "red", lwd = 3, lty = 1)
    abline(v = result$null.mean, col = "blue", lwd = 2, lty = 2)
    abline(v = 0.5, col = "darkgreen", lwd = 1, lty = 3)

    legend("topright",
           legend = c("Observed", "Null mean", "0.5 (independence)"),
           col = c("red", "blue", "darkgreen"),
           lwd = c(3, 2, 1),
           lty = c(1, 2, 3))

    ## Add text annotations
    y.pos <- par("usr")[4] * 0.9
    text(result$obs.concordant, y.pos,
         sprintf("Obs: %.4f\np=%.4f",
                result$obs.concordant,
                result$p.value.two.sided),
         pos = 4, col = "red")

    ## Plot 2: QQ plot to check if observed is in tail
    qqnorm(result$perm.concordant, main = "Q-Q Plot (Permutation Dist)")
    qqline(result$perm.concordant, col = "blue")
    points(qnorm(rank(result$obs.concordant) / (length(result$perm.concordant) + 1)),
           result$obs.concordant, col = "red", pch = 19, cex = 2)

    par(mfrow = c(1, 1))
}


#' Chi-square test for 2x2 concordance table
test.concordance.chisq <- function(delta.y, delta.z) {
    tbl <- table(delta.y > 0, delta.z > 0)
    chi.result <- chisq.test(tbl)

    ## Effect size: Cramér's V
    n <- sum(tbl)
    cramer.v <- sqrt(chi.result$statistic / n)

    ## Odds ratio (how much more likely concordance vs discordance)
    or <- (tbl[1,1] * tbl[2,2]) / (tbl[1,2] * tbl[2,1])

    list(
        table = tbl,
        chi.sq = chi.result$statistic,
        p.value = chi.result$p.value,
        cramer.v = as.numeric(cramer.v),
        odds.ratio = or,
        concordance = (tbl[1,1] + tbl[2,2]) / n
    )
}

#' Bayesian test: Is concordance proportion > 0.5?
#'
#' Uses Beta-Binomial model with weakly informative prior
test.concordance.bayes <- function(delta.y, delta.z,
                                   prior.a = 1, prior.b = 1,
                                   n.samples = 10000) {

    ## Count concordant and discordant edges
    concordant <- sum((delta.y > 0 & delta.z > 0) |
                     (delta.y < 0 & delta.z < 0))
    n <- length(delta.y)

    ## Posterior: Beta(a + concordant, b + n - concordant)
    post.a <- prior.a + concordant
    post.b <- prior.b + (n - concordant)

    ## Sample from posterior
    theta.samples <- rbeta(n.samples, post.a, post.b)

    ## Posterior probability that concordance > 0.5
    p.greater.half <- mean(theta.samples > 0.5)

    ## 95% credible interval
    ci <- quantile(theta.samples, c(0.025, 0.975))

    ## Bayes factor for H1: theta > 0.5 vs H0: theta = 0.5
    ## Using Savage-Dickey density ratio approximation
    prior.density.at.half <- dbeta(0.5, prior.a, prior.b)
    post.density.at.half <- dbeta(0.5, post.a, post.b)
    bf10 <- post.density.at.half / prior.density.at.half

    list(
        n.concordant = concordant,
        n.total = n,
        obs.proportion = concordant / n,
        posterior.mean = post.a / (post.a + post.b),
        posterior.samples = theta.samples,
        credible.interval = ci,
        prob.greater.half = p.greater.half,
        bayes.factor = bf10
    )
}

#' Extract concordance and discordance subgraphs
#'
#' @param adj.list Graph adjacency list
#' @param weight.list Edge weights (lengths)
#' @param vertex.delta.y List of y-differences for each vertex
#' @param vertex.delta.z List of z-differences for each vertex
#' @param min.component.size Minimum size for "signal" components (default 5)
#' @return List containing concordance and discordance subgraphs with components
extract.concordance.subgraphs <- function(adj.list, weight.list,
                                         vertex.delta.y, vertex.delta.z,
                                         min.component.size = 5) {

    n.vertices <- length(adj.list)

    ## Initialize edge lists for concordance and discordance subgraphs
    concordant.edges <- list()
    discordant.edges <- list()

    concordant.adj.list <- vector("list", n.vertices)
    concordant.weight.list <- vector("list", n.vertices)
    discordant.adj.list <- vector("list", n.vertices)
    discordant.weight.list <- vector("list", n.vertices)

    ## Mark edges as concordant or discordant
    for (v in 1:n.vertices) {
        neighbors <- adj.list[[v]]
        if (length(neighbors) == 0) next

        delta.y <- vertex.delta.y[[v]]
        delta.z <- vertex.delta.z[[v]]

        ## Concordant: both same sign (both positive or both negative)
        concordant.mask <- (delta.y > 0 & delta.z > 0) | (delta.y < 0 & delta.z < 0)

        ## Discordant: opposite signs
        discordant.mask <- (delta.y > 0 & delta.z < 0) | (delta.y < 0 & delta.z > 0)

        ## Build subgraph adjacency lists
        if (any(concordant.mask)) {
            concordant.adj.list[[v]] <- neighbors[concordant.mask]
            concordant.weight.list[[v]] <- weight.list[[v]][concordant.mask]

            ## Store edges (avoid duplicates by only storing v < u)
            for (i in which(concordant.mask)) {
                u <- neighbors[i]
                if (v < u) {
                    concordant.edges[[length(concordant.edges) + 1]] <-
                        list(v = v, u = u,
                             weight = weight.list[[v]][i],
                             delta.y = delta.y[i],
                             delta.z = delta.z[i])
                }
            }
        }

        if (any(discordant.mask)) {
            discordant.adj.list[[v]] <- neighbors[discordant.mask]
            discordant.weight.list[[v]] <- weight.list[[v]][discordant.mask]

            for (i in which(discordant.mask)) {
                u <- neighbors[i]
                if (v < u) {
                    discordant.edges[[length(discordant.edges) + 1]] <-
                        list(v = v, u = u,
                             weight = weight.list[[v]][i],
                             delta.y = delta.y[i],
                             delta.z = delta.z[i])
                }
            }
        }
    }

    ## Find connected components
    concordant.components <- find.connected.components(concordant.adj.list)
    discordant.components <- find.connected.components(discordant.adj.list)

    ## Classify components by size
    concordant.large <- which(
        sapply(concordant.components$component.sizes, function(s) s >= min.component.size)
    )
    discordant.large <- which(
        sapply(discordant.components$component.sizes, function(s) s >= min.component.size)
    )

    list(
        concordant = list(
            adj.list = concordant.adj.list,
            weight.list = concordant.weight.list,
            edges = concordant.edges,
            components = concordant.components,
            large.components = concordant.large,
            n.edges = length(concordant.edges)
        ),
        discordant = list(
            adj.list = discordant.adj.list,
            weight.list = discordant.weight.list,
            edges = discordant.edges,
            components = discordant.components,
            large.components = discordant.large,
            n.edges = length(discordant.edges)
        ),
        min.component.size = min.component.size
    )
}

#' Find connected components in a graph
#'
#' @param adj.list Graph adjacency list
#' @return List with component assignments and sizes
find.connected.components <- function(adj.list) {

    n.vertices <- length(adj.list)
    component.id <- rep(NA, n.vertices)
    current.component <- 0

    ## BFS for each unvisited vertex
    for (start in 1:n.vertices) {
        if (!is.na(component.id[start])) next
        if (length(adj.list[[start]]) == 0) next  # Isolated vertex

        current.component <- current.component + 1
        queue <- start
        component.id[start] <- current.component

        while (length(queue) > 0) {
            v <- queue[1]
            queue <- queue[-1]

            neighbors <- adj.list[[v]]
            for (u in neighbors) {
                if (is.na(component.id[u])) {
                    component.id[u] <- current.component
                    queue <- c(queue, u)
                }
            }
        }
    }

    ## Component sizes and member lists
    components.list <- split(which(!is.na(component.id)),
                            component.id[!is.na(component.id)])
    component.sizes <- sapply(components.list, length)

    list(
        component.id = component.id,
        n.components = current.component,
        components = components.list,
        component.sizes = component.sizes
    )
}

#' Permutation test for concordance component structure
#'
#' Tests whether observed large components are more extensive than expected
test.concordance.components <- function(vertex.delta.y, vertex.delta.z,
                                       adj.list, weight.list,
                                       observed.subgraphs,
                                       n.perm = 1000) {

    ## Observed statistics
    obs.n.large.concordant <- length(observed.subgraphs$concordant$large.components)
    obs.max.concordant.size <- max(observed.subgraphs$concordant$components$component.sizes,
                                   na.rm = TRUE)
    obs.total.concordant.vertices <- sum(
        observed.subgraphs$concordant$components$component.sizes[
            observed.subgraphs$concordant$large.components
        ]
    )

    ## Permutation distribution
    perm.n.large <- integer(n.perm)
    perm.max.size <- integer(n.perm)
    perm.total.vertices <- integer(n.perm)

    for (i in 1:n.perm) {
        ## Permute z while keeping y fixed (preserves graph structure)
        vertex.delta.z.perm <- lapply(vertex.delta.z, sample)

        perm.subgraphs <- extract.concordance.subgraphs(
            adj.list, weight.list,
            vertex.delta.y, vertex.delta.z.perm,
            min.component.size = observed.subgraphs$min.component.size
        )

        perm.n.large[i] <- length(perm.subgraphs$concordant$large.components)
        perm.max.size[i] <- max(perm.subgraphs$concordant$components$component.sizes,
                               na.rm = TRUE)
        perm.total.vertices[i] <- sum(
            perm.subgraphs$concordant$components$component.sizes[
                perm.subgraphs$concordant$large.components
            ]
        )
    }

    ## P-values
    p.n.large <- mean(perm.n.large >= obs.n.large.concordant)
    p.max.size <- mean(perm.max.size >= obs.max.concordant.size)
    p.total.vertices <- mean(perm.total.vertices >= obs.total.concordant.vertices)

    list(
        observed = list(
            n.large.components = obs.n.large.concordant,
            max.component.size = obs.max.concordant.size,
            total.vertices.in.large = obs.total.concordant.vertices
        ),
        p.values = list(
            n.large = p.n.large,
            max.size = p.max.size,
            total.vertices = p.total.vertices
        ),
        permutation = list(
            n.large = perm.n.large,
            max.size = perm.max.size,
            total.vertices = perm.total.vertices
        )
    )
}

#' Find paths through concordance subgraph
#'
#' Identifies trajectories where y and z consistently move together
find.concordant.paths <- function(concordant.subgraph,
                                 start.vertices, end.vertices,
                                 max.path.length = 10) {

    paths <- list()

    for (start in start.vertices) {
        for (end in end.vertices) {
            ## Find all paths from start to end in concordant subgraph
            ## (using BFS or DFS with path tracking)
            path <- find.path.in.subgraph(
                concordant.subgraph$adj.list,
                start, end, max.path.length
            )

            if (!is.null(path)) {
                paths[[length(paths) + 1]] <- path
            }
        }
    }

    paths
}

#' Characterize concordance components
#'
#' Identifies shared features of vertices in large components
characterize.components <- function(subgraphs,
                                   X,  # Feature matrix
                                   y,  # Outcome
                                   clinical.data) {

    large.components <- subgraphs$concordant$components$components[
        subgraphs$concordant$large.components
    ]

    results <- lapply(1:length(large.components), function(comp.idx) {
        vertices <- large.components[[comp.idx]]

        ## Summarize features in this component
        list(
            component.id = comp.idx,
            size = length(vertices),

            ## Microbiome features
            mean.abundances = colMeans(X[vertices, ]),

            ## Outcome distribution
            y.mean = mean(y[vertices]),
            y.range = range(y[vertices]),

            ## Clinical features
            clinical.summary = summarize.clinical(clinical.data[vertices, ]),

            ## Topological features
            density = compute.subgraph.density(
                subgraphs$concordant$adj.list[vertices]
            ),
            diameter = compute.diameter(
                subgraphs$concordant$adj.list[vertices]
            )
        )
    })

    results
}

#' Relate concordance components to gradient flow basins
#'
#' Tests whether concordance structure aligns with basin boundaries
relate.concordance.to.basins <- function(subgraphs, basin.assignments,
                                        basin.extrema) {

    concordant.components <- subgraphs$concordant$components$components

    ## For each large component, determine basin membership
    component.basin.overlap <- lapply(concordant.components, function(vertices) {
        basin.counts <- table(basin.assignments[vertices])

        list(
            vertices = vertices,
            primary.basin = as.numeric(names(which.max(basin.counts))),
            basin.purity = max(basin.counts) / length(vertices),
            n.basins.spanned = length(basin.counts),
            crosses.basin.boundary = length(basin.counts) > 1
        )
    })

    ## Components that cross basin boundaries are particularly interesting:
    ## They represent regions where y and z change concordantly
    ## DESPITE gradient flow pointing in different directions

    boundary.crossing.components <- which(
        sapply(component.basin.overlap, function(x) x$crosses.basin.boundary)
    )

    list(
        component.basin.overlap = component.basin.overlap,
        boundary.crossing = boundary.crossing.components
    )
}

#' Build multi-feature concordance network
#'
#' For each pair of features, create concordance subgraph
#' Then find meta-structure across all feature pairs
build.concordance.network <- function(adj.list, weight.list,
                                     y, Z,  # Multiple features
                                     vertex.delta.y, vertex.delta.Z) {

    n.features <- ncol(Z)

    ## Compute concordance subgraph for each feature
    feature.subgraphs <- lapply(1:n.features, function(j) {
        extract.concordance.subgraphs(
            adj.list, weight.list,
            vertex.delta.y, vertex.delta.Z[[j]]
        )
    })

    ## Build feature concordance matrix: which features have overlapping
    ## concordant components?
    feature.concordance.matrix <- matrix(0, n.features, n.features)

    for (i in 1:(n.features-1)) {
        for (j in (i+1):n.features) {
            ## Jaccard similarity of large component vertices
            vertices.i <- unlist(feature.subgraphs[[i]]$concordant$components$components[
                feature.subgraphs[[i]]$concordant$large.components
            ])
            vertices.j <- unlist(feature.subgraphs[[j]]$concordant$components$components[
                feature.subgraphs[[j]]$concordant$large.components
            ])

            overlap <- length(intersect(vertices.i, vertices.j))
            union <- length(union(vertices.i, vertices.j))

            feature.concordance.matrix[i, j] <- overlap / union
            feature.concordance.matrix[j, i] <- feature.concordance.matrix[i, j]
        }
    }

    ## Identify feature modules: groups of features with shared
    ## concordance structure
    feature.modules <- cluster::pam(
        as.dist(1 - feature.concordance.matrix),
        k = determine.optimal.k(feature.concordance.matrix)
    )

    list(
        feature.subgraphs = feature.subgraphs,
        feature.concordance.matrix = feature.concordance.matrix,
        feature.modules = feature.modules
    )
}

#' Visualize concordance subgraph with components highlighted
plot.concordance.subgraph <- function(subgraphs,
                                     graph.layout,  # 2D or 3D coordinates
                                     y.values,
                                     highlight.large = TRUE) {

    require(igraph)

    ## Build igraph object
    edges.df <- do.call(rbind, lapply(subgraphs$concordant$edges, function(e) {
        data.frame(from = e$v, to = e$u)
    }))

    g <- graph_from_data_frame(edges.df, directed = FALSE)

    ## Color vertices by component
    V(g)$component <- subgraphs$concordant$components$component.id[as.numeric(V(g)$name)]

    if (highlight.large) {
        ## Large components get distinct colors
        large.comp.ids <- names(subgraphs$concordant$components$components)[
            subgraphs$concordant$large.components
        ]

        V(g)$color <- ifelse(
            V(g)$component %in% large.comp.ids,
            rainbow(length(large.comp.ids))[match(V(g)$component, large.comp.ids)],
            "gray80"
        )
    }

    ## Size vertices by outcome value
    V(g)$size <- scales::rescale(y.values[as.numeric(V(g)$name)], to = c(2, 10))

    plot(g,
         layout = graph.layout[as.numeric(V(g)$name), ],
         main = "Concordance Subgraph (Large Components Highlighted)")

    ## Add legend
    legend("topright",
           legend = c("Large components", "Small components"),
           pch = 21,
           pt.bg = c(rainbow(length(large.comp.ids))[1], "gray80"),
           pt.cex = 2)
}


#' Find dense modules in concordance subgraph using community detection
#'
#' @param concordant.subgraph The concordant subgraph from extract.concordance.subgraphs()
#' @param method Community detection method ("louvain", "leiden", "infomap", "walktrap")
#' @param min.module.size Minimum module size to consider
#' @return List with module assignments and quality metrics
find.concordance.modules <- function(concordant.subgraph,
                                    method = "louvain",
                                    min.module.size = 10,
                                    resolution = 1.0) {

    require(igraph)

    ## Build igraph object from concordant edges
    if (length(concordant.subgraph$edges) == 0) {
        return(list(n.modules = 0, modules = list()))
    }

    edges.df <- do.call(rbind, lapply(concordant.subgraph$edges, function(e) {
        data.frame(from = e$v, to = e$u, weight = 1/e$weight)  # Weight by inverse distance
    }))

    g <- graph_from_data_frame(edges.df, directed = FALSE)

    ## Apply community detection algorithm
    communities <- switch(method,
        "louvain" = cluster_louvain(g, weights = E(g)$weight,
                                     resolution = resolution),
        "leiden" = {
            ## Requires leidenalg Python package via reticulate
            ## Or use igraph's implementation if available
            cluster_leiden(g, weights = E(g)$weight,
                          resolution_parameter = resolution)
        },
        "infomap" = cluster_infomap(g, e.weights = E(g)$weight),
        "walktrap" = cluster_walktrap(g, weights = E(g)$weight, steps = 4),
        "fast_greedy" = cluster_fast_greedy(g, weights = E(g)$weight),
        stop("Unknown method: ", method)
    )

    ## Extract module assignments for all vertices
    n.vertices <- max(as.numeric(V(g)$name))
    module.id <- rep(NA, n.vertices)
    module.id[as.numeric(V(g)$name)] <- membership(communities)

    ## Get modules as vertex lists
    modules.list <- split(as.numeric(V(g)$name), membership(communities))
    module.sizes <- sapply(modules.list, length)

    ## Filter by size
    large.modules <- which(module.sizes >= min.module.size)

    ## Compute module quality metrics
    module.metrics <- lapply(large.modules, function(mod.idx) {
        vertices <- modules.list[[mod.idx]]

        compute.module.metrics(
            vertices,
            concordant.subgraph$adj.list,
            concordant.subgraph$weight.list
        )
    })

    list(
        module.id = module.id,
        modules = modules.list[large.modules],
        module.sizes = module.sizes[large.modules],
        module.metrics = module.metrics,
        modularity = modularity(communities),
        n.modules = length(large.modules),
        method = method
    )
}

#' Compute quality metrics for a module
compute.module.metrics <- function(vertices, adj.list, weight.list) {

    ## Internal edges: both endpoints in module
    n.internal <- 0
    total.internal.weight <- 0

    ## Boundary edges: one endpoint in module, one outside
    n.boundary <- 0
    total.boundary.weight <- 0

    for (v in vertices) {
        neighbors <- adj.list[[v]]
        weights <- weight.list[[v]]

        for (i in seq_along(neighbors)) {
            u <- neighbors[i]
            w <- weights[i]

            if (u %in% vertices) {
                ## Internal edge (count each edge once)
                if (v < u) {
                    n.internal <- n.internal + 1
                    total.internal.weight <- total.internal.weight + w
                }
            } else {
                ## Boundary edge
                n.boundary <- n.boundary + 1
                total.boundary.weight <- total.boundary.weight + w
            }
        }
    }

    ## Possible internal edges
    n.possible.internal <- length(vertices) * (length(vertices) - 1) / 2

    ## Density: actual internal edges / possible internal edges
    density <- n.internal / n.possible.internal

    ## Conductance: boundary edges / (internal + boundary edges)
    ## Lower conductance = better separated module
    conductance <- n.boundary / (n.internal + n.boundary)

    ## Cut ratio: boundary edges / (size * (n - size))
    n.total <- length(adj.list)
    cut.ratio <- n.boundary / (length(vertices) * (n.total - length(vertices)))

    list(
        size = length(vertices),
        n.internal.edges = n.internal,
        n.boundary.edges = n.boundary,
        density = density,
        conductance = conductance,
        cut.ratio = cut.ratio,
        internal.weight = total.internal.weight,
        boundary.weight = total.boundary.weight
    )
}

#' K-core decomposition of concordance subgraph
#'
#' Finds nested sequence of dense subgraphs
#' Higher k-cores = denser, more robust concordance signal
find.concordance.kcores <- function(concordant.subgraph,
                                   min.k = 5) {

    require(igraph)

    ## Build igraph
    edges.df <- do.call(rbind, lapply(concordant.subgraph$edges, function(e) {
        data.frame(from = e$v, to = e$u)
    }))

    g <- graph_from_data_frame(edges.df, directed = FALSE)

    ## Compute k-core decomposition
    coreness <- coreness(g)

    ## Map back to full vertex set
    n.vertices <- max(as.numeric(V(g)$name))
    vertex.coreness <- rep(0, n.vertices)
    vertex.coreness[as.numeric(V(g)$name)] <- coreness

    ## Extract k-cores for k >= min.k
    kcores <- list()
    max.k <- max(coreness)

    for (k in min.k:max.k) {
        vertices.in.kcore <- as.numeric(V(g)$name)[coreness >= k]

        if (length(vertices.in.kcore) > 0) {
            kcores[[as.character(k)]] <- list(
                k = k,
                vertices = vertices.in.kcore,
                size = length(vertices.in.kcore)
            )
        }
    }

    list(
        vertex.coreness = vertex.coreness,
        kcores = kcores,
        max.k = max.k
    )
}

#' Find dense modules using local density peaks
#'
#' Identifies regions of locally high concordant edge density
find.density.based.modules <- function(concordant.subgraph,
                                      min.module.size = 10,
                                      density.threshold.quantile = 0.9) {

    n.vertices <- length(concordant.subgraph$adj.list)

    ## Compute local density: number of concordant neighbors
    local.density <- sapply(1:n.vertices, function(v) {
        length(concordant.subgraph$adj.list[[v]])
    })

    ## Also compute "second-order" density: density of neighbors
    second.order.density <- sapply(1:n.vertices, function(v) {
        neighbors <- concordant.subgraph$adj.list[[v]]
        if (length(neighbors) == 0) return(0)
        mean(local.density[neighbors])
    })

    ## Combined density score
    density.score <- local.density * second.order.density

    ## Find density peaks (vertices with higher density than all neighbors)
    density.peaks <- sapply(1:n.vertices, function(v) {
        neighbors <- concordant.subgraph$adj.list[[v]]
        if (length(neighbors) == 0) return(FALSE)

        ## v is a peak if its density is highest among neighbors
        density.score[v] > max(density.score[neighbors])
    })

    peak.vertices <- which(density.peaks &
                          density.score > quantile(density.score,
                                                  density.threshold.quantile))

    ## Assign each vertex to nearest peak
    module.assignment <- assign.to.nearest.peak(
        concordant.subgraph$adj.list,
        peak.vertices,
        density.score
    )

    ## Extract modules
    modules.list <- split(which(!is.na(module.assignment)),
                         module.assignment[!is.na(module.assignment)])

    ## Filter by size
    module.sizes <- sapply(modules.list, length)
    large.modules <- modules.list[module.sizes >= min.module.size]

    list(
        modules = large.modules,
        peak.vertices = peak.vertices,
        density.score = density.score,
        module.assignment = module.assignment
    )
}

## Helper: assign vertices to nearest density peak
assign.to.nearest.peak <- function(adj.list, peaks, density.score) {

    n.vertices <- length(adj.list)
    assignment <- rep(NA, n.vertices)

    ## BFS from each peak
    for (peak.idx in seq_along(peaks)) {
        peak <- peaks[peak.idx]
        assignment[peak] <- peak.idx

        queue <- peak
        visited <- logical(n.vertices)
        visited[peak] <- TRUE

        while (length(queue) > 0) {
            v <- queue[1]
            queue <- queue[-1]

            neighbors <- adj.list[[v]]
            for (u in neighbors) {
                if (!visited[u]) {
                    visited[u] <- TRUE

                    ## Assign to this peak if density is decreasing
                    if (density.score[u] <= density.score[v]) {
                        assignment[u] <- peak.idx
                        queue <- c(queue, u)
                    }
                }
            }
        }
    }

    assignment
}

#' Test module concordance significance
#'
#' Permutation test: are modules more densely concordant than random?
test.module.significance <- function(modules,
                                    vertex.delta.y, vertex.delta.z,
                                    adj.list, weight.list,
                                    n.perm = 1000) {

    ## For each module, compute observed concordance metrics
    observed.metrics <- lapply(modules, function(vertices) {

        ## Extract edges within module
        module.edges.y <- list()
        module.edges.z <- list()

        for (v in vertices) {
            neighbors <- adj.list[[v]]
            mask <- neighbors %in% vertices

            module.edges.y <- c(module.edges.y, list(vertex.delta.y[[v]][mask]))
            module.edges.z <- c(module.edges.z, list(vertex.delta.z[[v]][mask]))
        }

        all.y <- unlist(module.edges.y)
        all.z <- unlist(module.edges.z)

        list(
            concordance = mean((all.y > 0 & all.z > 0) |
                             (all.y < 0 & all.z < 0)),
            n.edges = length(all.y),
            phi = cor(all.y > 0, all.z > 0)
        )
    })

    ## Permutation test for each module
    module.pvalues <- lapply(seq_along(modules), function(mod.idx) {

        obs <- observed.metrics[[mod.idx]]
        vertices <- modules[[mod.idx]]

        ## Permute z within the entire graph, recompute concordance in module
        perm.concordance <- replicate(n.perm, {

            ## Permute z
            vertex.delta.z.perm <- lapply(vertex.delta.z, sample)

            ## Recompute module concordance
            module.edges.y <- list()
            module.edges.z.perm <- list()

            for (v in vertices) {
                neighbors <- adj.list[[v]]
                mask <- neighbors %in% vertices

                module.edges.y <- c(module.edges.y, list(vertex.delta.y[[v]][mask]))
                module.edges.z.perm <- c(module.edges.z.perm,
                                        list(vertex.delta.z.perm[[v]][mask]))
            }

            all.y <- unlist(module.edges.y)
            all.z.perm <- unlist(module.edges.z.perm)

            mean((all.y > 0 & all.z.perm > 0) |
                 (all.y < 0 & all.z.perm < 0))
        })

        ## Two-sided p-value
        if (obs$concordance >= median(perm.concordance)) {
            p.val <- 2 * mean(perm.concordance >= obs$concordance)
        } else {
            p.val <- 2 * mean(perm.concordance <= obs$concordance)
        }

        list(
            observed = obs,
            p.value = min(p.val, 1.0),
            perm.dist = perm.concordance
        )
    })

    module.pvalues
}
