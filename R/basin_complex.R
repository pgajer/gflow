# Canonical basin-complex schema and Phase B constructor.

.basin.complex.methods <- c(
    "trajectory_flow",
    "superlevel_merge_tree",
    "geodesic_reachability",
    "rtcb",
    "overlap_cell_complex"
)

.basin.complex.directions <- c("max", "min", "both")

.stop.basin.complex <- function(message,
                                argument = NULL,
                                class = "gflow_basin_input_error",
                                details = NULL) {
    condition <- structure(
        list(
            message = message,
            call = NULL,
            argument = argument,
            details = details
        ),
        class = c(class, "error", "condition")
    )
    stop(condition)
}

.basin.assert.named.list <- function(x, argument) {
    if (!is.list(x) || is.object(x)) {
        .stop.basin.complex(
            sprintf("'%s' must be an ordinary named list.", argument),
            argument
        )
    }
    if (length(x) > 0L) {
        nms <- names(x)
        if (is.null(nms) || anyNA(nms) || any(!nzchar(nms)) ||
            anyDuplicated(nms)) {
            .stop.basin.complex(
                sprintf("'%s' must have unique, nonempty names.", argument),
                argument
            )
        }
    }
    x
}

.basin.merge.params <- function(defaults, supplied, argument) {
    supplied <- .basin.assert.named.list(supplied, argument)
    unknown <- setdiff(names(supplied), names(defaults))
    if (length(unknown) > 0L) {
        .stop.basin.complex(
            sprintf(
                "Unknown %s field%s: %s.",
                argument,
                if (length(unknown) == 1L) "" else "s",
                paste(unknown, collapse = ", ")
            ),
            argument,
            details = list(unknown = unknown)
        )
    }
    out <- defaults
    for (name in names(supplied)) {
        out[name] <- supplied[name]
    }
    out
}

.basin.assert.logical <- function(x, argument) {
    if (!is.logical(x) || length(x) != 1L || is.na(x)) {
        .stop.basin.complex(
            sprintf("'%s' must be TRUE or FALSE.", argument),
            argument
        )
    }
    x
}

.basin.assert.number <- function(x,
                                 argument,
                                 lower = -Inf,
                                 upper = Inf,
                                 lower.open = FALSE,
                                 upper.open = FALSE,
                                 finite = TRUE) {
    ok <- is.numeric(x) && length(x) == 1L && !is.na(x)
    if (ok && finite) {
        ok <- is.finite(x)
    }
    if (ok) {
        ok <- if (lower.open) x > lower else x >= lower
    }
    if (ok) {
        ok <- if (upper.open) x < upper else x <= upper
    }
    if (!ok) {
        interval <- paste0(
            if (lower.open) "(" else "[",
            format(lower),
            ", ",
            format(upper),
            if (upper.open) ")" else "]"
        )
        .stop.basin.complex(
            sprintf("'%s' must be a numeric scalar in %s.", argument, interval),
            argument
        )
    }
    as.numeric(x)
}

.basin.assert.integer <- function(x,
                                  argument,
                                  lower = -Inf,
                                  upper = Inf) {
    ok <- is.numeric(x) && length(x) == 1L && is.finite(x) &&
        !is.na(x) && x == floor(x) && x >= lower && x <= upper
    if (!ok) {
        .stop.basin.complex(
            sprintf(
                "'%s' must be a whole-number scalar between %s and %s.",
                argument,
                format(lower),
                format(upper)
            ),
            argument
        )
    }
    as.integer(x)
}

.basin.assert.choice <- function(x, choices, argument) {
    if (length(x) > 1L && identical(as.character(x), choices)) {
        x <- choices[[1L]]
    }
    if (!is.character(x) || length(x) != 1L || is.na(x) ||
        !x %in% choices) {
        .stop.basin.complex(
            sprintf(
                "'%s' must be one of: %s.",
                argument,
                paste(choices, collapse = ", ")
            ),
            argument
        )
    }
    x
}

.basin.assert.vector <- function(x,
                                 argument,
                                 n,
                                 nonnegative = FALSE,
                                 positive.total = FALSE) {
    if (!is.numeric(x) || is.object(x) || length(x) != n ||
        any(!is.finite(x))) {
        .stop.basin.complex(
            sprintf(
                "'%s' must be a finite numeric vector of length %d.",
                argument,
                n
            ),
            argument
        )
    }
    if (nonnegative && any(x < 0)) {
        .stop.basin.complex(
            sprintf("'%s' must contain only nonnegative values.", argument),
            argument
        )
    }
    if (positive.total && sum(x) <= 0) {
        .stop.basin.complex(
            sprintf("'%s' must have a positive total.", argument),
            argument
        )
    }
    as.numeric(x)
}

.basin.graph.defaults <- function() {
    list(
        edge.length.symmetry.tolerance = sqrt(.Machine$double.eps)
    )
}

.basin.validate.graph <- function(adj.list,
                                  edge.length.list,
                                  graph.params) {
    if (!is.list(adj.list) || is.object(adj.list)) {
        .stop.basin.complex("'adj.list' must be an ordinary list.", "adj.list")
    }
    if (!is.list(edge.length.list) || is.object(edge.length.list)) {
        .stop.basin.complex(
            "'edge.length.list' must be an ordinary list.",
            "edge.length.list"
        )
    }
    n <- length(adj.list)
    if (n == 0L) {
        .stop.basin.complex("'adj.list' must describe at least one vertex.", "adj.list")
    }
    if (length(edge.length.list) != n) {
        .stop.basin.complex(
            "'adj.list' and 'edge.length.list' must have the same length.",
            "edge.length.list"
        )
    }

    graph.params <- .basin.merge.params(
        .basin.graph.defaults(),
        graph.params,
        "graph.params"
    )
    tolerance <- .basin.assert.number(
        graph.params$edge.length.symmetry.tolerance,
        "graph.params$edge.length.symmetry.tolerance",
        lower = 0
    )
    graph.params$edge.length.symmetry.tolerance <- tolerance

    normalized.adj <- vector("list", n)
    normalized.length <- vector("list", n)
    names(normalized.adj) <- names(adj.list)
    names(normalized.length) <- names(edge.length.list)

    for (vertex in seq_len(n)) {
        neighbors <- adj.list[[vertex]]
        lengths <- edge.length.list[[vertex]]

        if (!is.numeric(neighbors) || is.object(neighbors) ||
            any(!is.finite(neighbors)) ||
            any(neighbors != floor(neighbors))) {
            .stop.basin.complex(
                sprintf(
                    "'adj.list[[%d]]' must contain finite whole-number vertex ids.",
                    vertex
                ),
                "adj.list",
                details = list(vertex = vertex)
            )
        }
        neighbors <- as.integer(neighbors)
        if (any(neighbors < 1L | neighbors > n)) {
            .stop.basin.complex(
                sprintf("'adj.list[[%d]]' contains an out-of-range vertex id.", vertex),
                "adj.list",
                details = list(vertex = vertex)
            )
        }
        if (any(neighbors == vertex)) {
            .stop.basin.complex(
                sprintf("'adj.list[[%d]]' contains a self-loop.", vertex),
                "adj.list",
                details = list(vertex = vertex)
            )
        }
        if (anyDuplicated(neighbors)) {
            .stop.basin.complex(
                sprintf("'adj.list[[%d]]' contains duplicate neighbors.", vertex),
                "adj.list",
                details = list(vertex = vertex)
            )
        }
        if (!is.numeric(lengths) || is.object(lengths) ||
            length(lengths) != length(neighbors) ||
            any(!is.finite(lengths)) || any(lengths < 0)) {
            .stop.basin.complex(
                sprintf(
                    "'edge.length.list[[%d]]' must contain one finite nonnegative length per neighbor.",
                    vertex
                ),
                "edge.length.list",
                details = list(vertex = vertex)
            )
        }
        normalized.adj[[vertex]] <- neighbors
        normalized.length[[vertex]] <- as.numeric(lengths)
    }

    for (vertex in seq_len(n)) {
        neighbors <- normalized.adj[[vertex]]
        for (edge.index in seq_along(neighbors)) {
            neighbor <- neighbors[[edge.index]]
            reverse.index <- which(normalized.adj[[neighbor]] == vertex)
            if (length(reverse.index) != 1L) {
                .stop.basin.complex(
                    sprintf(
                        "Edge %d -> %d does not have exactly one reverse edge.",
                        vertex,
                        neighbor
                    ),
                    "adj.list",
                    details = list(vertex = vertex, neighbor = neighbor)
                )
            }
            difference <- abs(
                normalized.length[[vertex]][[edge.index]] -
                    normalized.length[[neighbor]][[reverse.index]]
            )
            if (difference > tolerance) {
                .stop.basin.complex(
                    sprintf(
                        "Reverse edge lengths for %d <-> %d differ by more than the symmetry tolerance.",
                        vertex,
                        neighbor
                    ),
                    "edge.length.list",
                    details = list(
                        vertex = vertex,
                        neighbor = neighbor,
                        difference = difference,
                        tolerance = tolerance
                    )
                )
            }
        }
    }

    component <- integer(n)
    n.components <- 0L
    for (start in seq_len(n)) {
        if (component[[start]] != 0L) {
            next
        }
        n.components <- n.components + 1L
        queue <- start
        component[[start]] <- n.components
        cursor <- 1L
        while (cursor <= length(queue)) {
            vertex <- queue[[cursor]]
            cursor <- cursor + 1L
            unseen <- normalized.adj[[vertex]][
                component[normalized.adj[[vertex]]] == 0L
            ]
            if (length(unseen) > 0L) {
                unseen <- unique(unseen)
                component[unseen] <- n.components
                queue <- c(queue, unseen)
            }
        }
    }

    list(
        adj.list = normalized.adj,
        edge.length.list = normalized.length,
        graph.params = graph.params,
        validation = list(
            valid = TRUE,
            n.vertices = n,
            n.edges = as.integer(sum(lengths(normalized.adj)) / 2L),
            n.components = n.components,
            component = component,
            isolated.vertices = which(lengths(normalized.adj) == 0L),
            has.zero.length.edges = any(
                unlist(normalized.length, use.names = FALSE) == 0
            ),
            undirected = TRUE,
            edge.length.symmetric = TRUE
        )
    )
}

.basin.method.defaults <- function(method, n.vertices) {
    switch(
        method,
        trajectory_flow = list(
            modulation = "CLOSEST",
            edge.length.quantile.thld = 0.9,
            long.edge.fallback = "allow_and_flag",
            store.trajectories = TRUE,
            max.trajectory.length = 10000L,
            symmetric.seeding = TRUE,
            tie.breaking = FALSE,
            tie.seed = NULL,
            primary.assignment.policy = "backend_primary"
        ),
        superlevel_merge_tree = list(
            plateau.tolerance = 0,
            elder.rule = "birth_then_representative_vertex",
            primary.assignment.policy = "elder_at_merge",
            root.death.policy = "component_opposite_extreme"
        ),
        geodesic_reachability = list(
            edge.length.quantile.thld = 0.9,
            store.trajectories = FALSE,
            primary.assignment.policy = "none"
        ),
        rtcb = list(
            edge.length.quantile.thld = 0.9,
            store.trajectories = FALSE,
            n.min = max(20L, as.integer(ceiling(sqrt(n.vertices)))),
            m.min = 0.7,
            q.min = NULL,
            run.max = 2L,
            tau0 = 0.1,
            kappa = 1.5,
            k.max = 24L,
            h.max = 1000L,
            d.path.max = Inf,
            eta.step = 0,
            epsilon.d = 1e-12,
            eps.num = 1e-12,
            sink.prune = FALSE,
            sink.prune.max.iter = 5L,
            max.paths.per.sink = 64L,
            primary.assignment.policy = "backend_primary"
        ),
        overlap_cell_complex = list(
            basin.merge.overlap.thld = 0.1,
            min.asc.desc.cell.size.thld = 1L,
            min.asc.asc.cell.size.thld = 1L,
            min.desc.desc.cell.size.thld = 1L,
            cell.graph.params = list(),
            primary.assignment.policy = "none"
        )
    )
}

.basin.resolve.method.params <- function(method,
                                         method.params,
                                         n.vertices,
                                         vertex.density) {
    defaults <- .basin.method.defaults(method, n.vertices)
    if (method == "rtcb" && "n.min" %in% names(method.params) &&
        is.null(method.params$n.min)) {
        method.params$n.min <- defaults$n.min
    }
    params <- .basin.merge.params(defaults, method.params, "method.params")

    if (method == "trajectory_flow") {
        params$modulation <- .basin.assert.choice(
            params$modulation,
            c("CLOSEST", "NONE", "DENSITY", "EDGELEN", "DENSITY_EDGELEN"),
            "method.params$modulation"
        )
        params$edge.length.quantile.thld <- .basin.assert.number(
            params$edge.length.quantile.thld,
            "method.params$edge.length.quantile.thld",
            lower = 0,
            upper = 1,
            lower.open = TRUE
        )
        params$long.edge.fallback <- .basin.assert.choice(
            params$long.edge.fallback,
            c("allow_and_flag", "allow", "forbid"),
            "method.params$long.edge.fallback"
        )
        params$store.trajectories <- .basin.assert.logical(
            params$store.trajectories,
            "method.params$store.trajectories"
        )
        params$max.trajectory.length <- .basin.assert.integer(
            params$max.trajectory.length,
            "method.params$max.trajectory.length",
            lower = 1
        )
        params$symmetric.seeding <- .basin.assert.logical(
            params$symmetric.seeding,
            "method.params$symmetric.seeding"
        )
        params$tie.breaking <- .basin.assert.logical(
            params$tie.breaking,
            "method.params$tie.breaking"
        )
        if (params$tie.breaking) {
            if (is.null(params$tie.seed)) {
                .stop.basin.complex(
                    "'method.params$tie.seed' is required when tie breaking is enabled.",
                    "method.params$tie.seed"
                )
            }
            params$tie.seed <- .basin.assert.integer(
                params$tie.seed,
                "method.params$tie.seed",
                lower = -.Machine$integer.max,
                upper = .Machine$integer.max
            )
        } else if (!is.null(params$tie.seed)) {
            .stop.basin.complex(
                "'method.params$tie.seed' must be NULL when tie breaking is disabled.",
                "method.params$tie.seed"
            )
        }
        params$primary.assignment.policy <- .basin.assert.choice(
            params$primary.assignment.policy,
            c("backend_primary", "largest_membership_then_extremum", "none"),
            "method.params$primary.assignment.policy"
        )
        if (params$modulation %in% c("DENSITY", "DENSITY_EDGELEN") &&
            is.null(vertex.density)) {
            .stop.basin.complex(
                sprintf(
                    "'vertex.density' is required for modulation '%s'.",
                    params$modulation
                ),
                "vertex.density"
            )
        }
    } else if (method == "superlevel_merge_tree") {
        params$plateau.tolerance <- .basin.assert.number(
            params$plateau.tolerance,
            "method.params$plateau.tolerance",
            lower = 0
        )
        if (params$plateau.tolerance != 0) {
            .stop.basin.complex(
                "Only exact plateau contraction is supported; 'plateau.tolerance' must be zero.",
                "method.params$plateau.tolerance"
            )
        }
        params$elder.rule <- .basin.assert.choice(
            params$elder.rule,
            "birth_then_representative_vertex",
            "method.params$elder.rule"
        )
        params$primary.assignment.policy <- .basin.assert.choice(
            params$primary.assignment.policy,
            "elder_at_merge",
            "method.params$primary.assignment.policy"
        )
        params$root.death.policy <- .basin.assert.choice(
            params$root.death.policy,
            "component_opposite_extreme",
            "method.params$root.death.policy"
        )
    } else if (method == "geodesic_reachability") {
        params$edge.length.quantile.thld <- .basin.assert.number(
            params$edge.length.quantile.thld,
            "method.params$edge.length.quantile.thld",
            lower = 0,
            upper = 1,
            lower.open = TRUE
        )
        params$store.trajectories <- .basin.assert.logical(
            params$store.trajectories,
            "method.params$store.trajectories"
        )
        params$primary.assignment.policy <- .basin.assert.choice(
            params$primary.assignment.policy,
            c("none", "largest_membership_then_extremum"),
            "method.params$primary.assignment.policy"
        )
    } else if (method == "rtcb") {
        params$edge.length.quantile.thld <- .basin.assert.number(
            params$edge.length.quantile.thld,
            "method.params$edge.length.quantile.thld",
            lower = 0,
            upper = 1,
            lower.open = TRUE
        )
        params$store.trajectories <- .basin.assert.logical(
            params$store.trajectories,
            "method.params$store.trajectories"
        )
        params$n.min <- .basin.assert.integer(
            params$n.min,
            "method.params$n.min",
            lower = 1
        )
        params$m.min <- .basin.assert.number(
            params$m.min,
            "method.params$m.min",
            lower = -1,
            upper = 1,
            lower.open = TRUE,
            upper.open = TRUE
        )
        if (!is.null(params$q.min)) {
            params$q.min <- .basin.assert.number(
                params$q.min,
                "method.params$q.min",
                lower = 0,
                upper = 1,
                lower.open = TRUE
            )
        }
        params$run.max <- .basin.assert.integer(
            params$run.max,
            "method.params$run.max",
            lower = 0
        )
        params$tau0 <- .basin.assert.number(
            params$tau0,
            "method.params$tau0",
            lower = 0
        )
        params$kappa <- .basin.assert.number(
            params$kappa,
            "method.params$kappa",
            lower = 0,
            lower.open = TRUE
        )
        for (name in c(
            "k.max", "h.max", "sink.prune.max.iter", "max.paths.per.sink"
        )) {
            params[[name]] <- .basin.assert.integer(
                params[[name]],
                paste0("method.params$", name),
                lower = 1
            )
        }
        params$d.path.max <- .basin.assert.number(
            params$d.path.max,
            "method.params$d.path.max",
            lower = 0,
            lower.open = TRUE,
            finite = FALSE
        )
        params$eta.step <- .basin.assert.number(
            params$eta.step,
            "method.params$eta.step",
            lower = 0
        )
        for (name in c("epsilon.d", "eps.num")) {
            params[[name]] <- .basin.assert.number(
                params[[name]],
                paste0("method.params$", name),
                lower = 0,
                lower.open = TRUE
            )
        }
        params$sink.prune <- .basin.assert.logical(
            params$sink.prune,
            "method.params$sink.prune"
        )
        params$primary.assignment.policy <- .basin.assert.choice(
            params$primary.assignment.policy,
            c("backend_primary", "none"),
            "method.params$primary.assignment.policy"
        )
    } else {
        params$basin.merge.overlap.thld <- .basin.assert.number(
            params$basin.merge.overlap.thld,
            "method.params$basin.merge.overlap.thld",
            lower = 0,
            upper = 1
        )
        for (name in c(
            "min.asc.desc.cell.size.thld",
            "min.asc.asc.cell.size.thld",
            "min.desc.desc.cell.size.thld"
        )) {
            params[[name]] <- .basin.assert.integer(
                params[[name]],
                paste0("method.params$", name),
                lower = 1
            )
        }
        params$cell.graph.params <- .basin.assert.named.list(
            params$cell.graph.params,
            "method.params$cell.graph.params"
        )
        params$primary.assignment.policy <- .basin.assert.choice(
            params$primary.assignment.policy,
            "none",
            "method.params$primary.assignment.policy"
        )
    }
    params
}

.basin.simplify.defaults <- function() {
    list(
        relative.value = list(
            enabled = FALSE,
            min.relative.value.max = 1.1,
            max.relative.value.min = 0.9
        ),
        maxima.clustering = list(
            enabled = FALSE,
            overlap.threshold = 0.15
        ),
        minima.clustering = list(
            enabled = FALSE,
            overlap.threshold = 0.15
        ),
        geometric.filter = list(
            enabled = FALSE,
            mean.neighbor.distance.percentile.max = 0.9,
            mean.hop.distance.percentile.max = 0.9,
            degree.percentile.max = 0.9,
            hop.k = 2L
        ),
        support.filter = list(
            enabled = FALSE,
            min.basin.size = 0L,
            min.trajectory.count = 0L,
            min.basin.mass = 0
        ),
        expansion = list(
            enabled = FALSE,
            policy = "nearest_retained_basin"
        )
    )
}

.basin.resolve.simplify.params <- function(method, simplify.params) {
    defaults <- .basin.simplify.defaults()
    supplied <- .basin.merge.params(defaults, simplify.params, "simplify.params")
    resolved <- defaults
    for (stage in names(defaults)) {
        resolved[[stage]] <- .basin.merge.params(
            defaults[[stage]],
            supplied[[stage]],
            paste0("simplify.params$", stage)
        )
        resolved[[stage]]$enabled <- .basin.assert.logical(
            resolved[[stage]]$enabled,
            paste0("simplify.params$", stage, "$enabled")
        )
    }

    for (name in c("min.relative.value.max", "max.relative.value.min")) {
        resolved$relative.value[[name]] <- .basin.assert.number(
            resolved$relative.value[[name]],
            paste0("simplify.params$relative.value$", name),
            lower = 0
        )
    }
    for (stage in c("maxima.clustering", "minima.clustering")) {
        resolved[[stage]]$overlap.threshold <- .basin.assert.number(
            resolved[[stage]]$overlap.threshold,
            paste0("simplify.params$", stage, "$overlap.threshold"),
            lower = 0,
            upper = 1
        )
    }
    for (name in c(
        "mean.neighbor.distance.percentile.max",
        "mean.hop.distance.percentile.max",
        "degree.percentile.max"
    )) {
        resolved$geometric.filter[[name]] <- .basin.assert.number(
            resolved$geometric.filter[[name]],
            paste0("simplify.params$geometric.filter$", name),
            lower = 0,
            upper = 1
        )
    }
    resolved$geometric.filter$hop.k <- .basin.assert.integer(
        resolved$geometric.filter$hop.k,
        "simplify.params$geometric.filter$hop.k",
        lower = 1
    )
    for (name in c("min.basin.size", "min.trajectory.count")) {
        resolved$support.filter[[name]] <- .basin.assert.integer(
            resolved$support.filter[[name]],
            paste0("simplify.params$support.filter$", name),
            lower = 0
        )
    }
    resolved$support.filter$min.basin.mass <- .basin.assert.number(
        resolved$support.filter$min.basin.mass,
        "simplify.params$support.filter$min.basin.mass",
        lower = 0,
        upper = 1
    )
    resolved$expansion$policy <- .basin.assert.choice(
        resolved$expansion$policy,
        "nearest_retained_basin",
        "simplify.params$expansion$policy"
    )

    applicable <- list(
        relative.value = c("trajectory_flow", "geodesic_reachability", "rtcb"),
        maxima.clustering = c("trajectory_flow", "geodesic_reachability", "rtcb"),
        minima.clustering = c("trajectory_flow", "geodesic_reachability", "rtcb"),
        geometric.filter = c("trajectory_flow", "geodesic_reachability", "rtcb"),
        support.filter = c(
            "trajectory_flow",
            "superlevel_merge_tree",
            "geodesic_reachability",
            "rtcb"
        ),
        expansion = c("geodesic_reachability", "rtcb")
    )
    for (stage in names(applicable)) {
        if (!method %in% applicable[[stage]] &&
            !identical(resolved[[stage]], defaults[[stage]])) {
            .stop.basin.complex(
                sprintf(
                    "Refinement stage '%s' is not applicable to method '%s'; disabled nondefault settings are also rejected.",
                    stage,
                    method
                ),
                paste0("simplify.params$", stage)
            )
        }
    }
    resolved
}

.empty.extrema.table <- function() {
    structure(
        list(
            extremum.id = character(),
            type = character(),
            representative.vertex = integer(),
            extremum.value = numeric(),
            plateau.id = character(),
            plateau.vertices = I(list()),
            n.plateau.vertices = integer(),
            is.retained = logical(),
            retention.status = character()
        ),
        class = "data.frame",
        row.names = integer()
    )
}

.empty.basin.table <- function() {
    structure(
        list(
            basin.id = character(),
            extremum.id = character(),
            type = character(),
            method = character(),
            retained = logical(),
            extremum.vertex = integer(),
            extremum.value = numeric(),
            birth.level = numeric(),
            death.level = numeric(),
            persistence = numeric(),
            parent.basin.id = character(),
            raw.support.vertices = I(list()),
            raw.support.size = integer(),
            raw.support.mass = numeric(),
            retained.support.vertices = I(list()),
            retained.support.size = integer(),
            retained.support.mass = numeric(),
            primary.support.vertices = I(list()),
            primary.support.size = integer(),
            primary.support.mass = numeric(),
            raw.allocated.mass = numeric(),
            assignment.status = character(),
            retention.status = character()
        ),
        class = "data.frame",
        row.names = integer()
    )
}

.empty.membership.table <- function() {
    structure(
        list(
            vertex = integer(),
            direction = character(),
            basin.id = character(),
            membership.weight = numeric(),
            membership.status = character(),
            source.stage = character(),
            is.primary = logical()
        ),
        class = "data.frame",
        row.names = integer()
    )
}

.empty.assignment.table <- function() {
    structure(
        list(
            vertex = integer(),
            direction = character(),
            basin.id = character(),
            assignment.weight = numeric(),
            assignment.status = character(),
            assignment.policy = character(),
            root.vertex = integer(),
            next.vertex = integer()
        ),
        class = "data.frame",
        row.names = integer()
    )
}

.empty.merge.table <- function() {
    structure(
        list(
            event.id = character(),
            direction = character(),
            losing.basin.id = character(),
            surviving.basin.id = character(),
            merge.plateau.id = character(),
            merge.vertices = I(list()),
            merge.level = numeric(),
            birth.level = numeric(),
            death.level = numeric(),
            persistence = numeric(),
            event.status = character()
        ),
        class = "data.frame",
        row.names = integer()
    )
}

.empty.diagnostics.table <- function() {
    structure(
        list(
            stage = character(),
            method = character(),
            condition.class = character(),
            message = character(),
            affected.vertices = I(list())
        ),
        class = "data.frame",
        row.names = integer()
    )
}

.empty.refinement.stages <- function() {
    structure(
        list(
            stage = character(),
            status = character(),
            parameters = I(list()),
            input.basin.count.max = integer(),
            input.basin.count.min = integer(),
            output.basin.count.max = integer(),
            output.basin.count.min = integer(),
            input.membership.count = integer(),
            output.membership.count = integer(),
            retained.basin.ids = I(list()),
            removed.basin.ids = I(list()),
            warnings = I(list()),
            diagnostics = I(list()),
            summary.snapshot = I(list())
        ),
        class = "data.frame",
        row.names = integer()
    )
}

.basin.diagnostic <- function(stage,
                              method,
                              condition.class,
                              message,
                              affected.vertices = integer()) {
    structure(
        list(
            stage = as.character(stage),
            method = as.character(method),
            condition.class = as.character(condition.class),
            message = as.character(message),
            affected.vertices = I(list(as.integer(affected.vertices)))
        ),
        class = "data.frame",
        row.names = 1L
    )
}

.basin.table.schemas <- function() {
    list(
        extrema = vapply(.empty.extrema.table(), typeof, character(1)),
        basin.table = vapply(.empty.basin.table(), typeof, character(1)),
        membership = vapply(.empty.membership.table(), typeof, character(1)),
        assignment = vapply(.empty.assignment.table(), typeof, character(1)),
        merge.table = vapply(.empty.merge.table(), typeof, character(1)),
        diagnostics = vapply(.empty.diagnostics.table(), typeof, character(1)),
        refinement.stages = vapply(.empty.refinement.stages(), typeof, character(1))
    )
}

.validate.basin.table.schema <- function(x, schema, name) {
    if (!is.data.frame(x) || !identical(names(x), names(schema)) ||
        !identical(
            unname(vapply(x, typeof, character(1))),
            unname(schema)
        )) {
        .stop.basin.complex(
            sprintf("'%s' does not match the canonical table schema.", name),
            name,
            class = "gflow_basin_schema_error"
        )
    }
    invisible(TRUE)
}

.validate.successful.basin.tables <- function(object) {
    requested <- .basin.requested.directions(object$direction)
    expected <- expand.grid(
        vertex = seq_len(object$n.vertices),
        direction = requested,
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
    )
    expected.key <- paste(expected$vertex, expected$direction, sep = "\r")
    assignment.key <- paste(
        object$assignment$vertex,
        object$assignment$direction,
        sep = "\r"
    )
    if (anyDuplicated(assignment.key) ||
        !setequal(assignment.key, expected.key)) {
        .stop.basin.complex(
            "Successful assignment rows must uniquely cover every requested vertex and direction.",
            "object$assignment",
            class = "gflow_basin_schema_error"
        )
    }
    if (any(!object$membership$direction %in% requested) ||
        any(object$membership$vertex < 1L |
            object$membership$vertex > object$n.vertices) ||
        any(!is.finite(object$membership$membership.weight)) ||
        any(object$membership$membership.weight < 0)) {
        .stop.basin.complex(
            "Membership rows contain an invalid direction, vertex, or weight.",
            "object$membership",
            class = "gflow_basin_schema_error"
        )
    }
    basin.key <- paste(
        object$basin.table$basin.id,
        object$basin.table$type,
        sep = "\r"
    )
    membership.basin.key <- paste(
        object$membership$basin.id,
        object$membership$direction,
        sep = "\r"
    )
    if (any(!membership.basin.key %in% basin.key)) {
        .stop.basin.complex(
            "Membership rows reference an unknown basin or direction.",
            "object$membership$basin.id",
            class = "gflow_basin_schema_error"
        )
    }
    if (nrow(object$membership) > 0L) {
        membership.group <- paste(
            object$membership$vertex,
            object$membership$direction,
            sep = "\r"
        )
        weight.sums <- tapply(
            object$membership$membership.weight,
            membership.group,
            sum
        )
        if (any(abs(weight.sums - 1) > 1e-12)) {
            .stop.basin.complex(
                "Membership weights must sum to one for every covered vertex and direction.",
                "object$membership$membership.weight",
                class = "gflow_basin_schema_error"
            )
        }
    }

    valid.status <- c(
        assigned = 1,
        unassigned = 0,
        filtered = 0,
        not_applicable = NA_real_,
        failed = NA_real_
    )
    if (any(!object$assignment$assignment.status %in% names(valid.status))) {
        .stop.basin.complex(
            "Assignment rows contain an unknown status.",
            "object$assignment$assignment.status",
            class = "gflow_basin_schema_error"
        )
    }
    expected.weight <- unname(
        valid.status[object$assignment$assignment.status]
    )
    weight.matches <- (is.na(expected.weight) &
        is.na(object$assignment$assignment.weight)) |
        (!is.na(expected.weight) &
            object$assignment$assignment.weight == expected.weight)
    if (any(!weight.matches)) {
        .stop.basin.complex(
            "Assignment weights do not match assignment statuses.",
            "object$assignment$assignment.weight",
            class = "gflow_basin_schema_error"
        )
    }
    assigned <- object$assignment$assignment.status == "assigned"
    assigned.key <- paste(
        object$assignment$basin.id[assigned],
        object$assignment$direction[assigned],
        sep = "\r"
    )
    if (any(is.na(object$assignment$basin.id[assigned])) ||
        any(!assigned.key %in% basin.key)) {
        .stop.basin.complex(
            "Assigned rows must reference a basin in the same direction.",
            "object$assignment$basin.id",
            class = "gflow_basin_schema_error"
        )
    }
    invisible(TRUE)
}

.basin.package.version <- function() {
    tryCatch(
        as.character(utils::packageVersion("gflow")),
        error = function(e) NA_character_
    )
}

.basin.or <- function(x, y) {
    if (is.null(x)) y else x
}

.new.basin.complex <- function(method,
                               direction,
                               graph.input,
                               field,
                               parameters,
                               provenance,
                               status = c("ok", "partial", "failed"),
                               extrema = .empty.extrema.table(),
                               basin.table = .empty.basin.table(),
                               membership = .empty.membership.table(),
                               assignment = .empty.assignment.table(),
                               merge.table = NULL,
                               trajectory.forest = NULL,
                               cell.complex = NULL,
                               diagnostics = .empty.diagnostics.table(),
                               warnings = character(),
                               raw = list()) {
    status <- .basin.assert.choice(
        status,
        c("ok", "partial", "failed"),
        "status"
    )
    object <- list(
        method = method,
        direction = direction,
        n.vertices = length(graph.input$adj.list),
        graph.input = graph.input,
        field = field,
        parameters = parameters,
        provenance = provenance,
        extrema = extrema,
        basin.table = basin.table,
        membership = membership,
        assignment = assignment,
        merge.table = merge.table,
        trajectory.forest = trajectory.forest,
        cell.complex = cell.complex,
        diagnostics = diagnostics,
        warnings = as.character(warnings),
        status = status,
        raw = raw
    )
    class(object) <- c(
        paste0("basin_complex_", method),
        "basin_complex",
        "gflow_basin_complex",
        "list"
    )
    .validate.basin.complex(object)
    object
}

.validate.basin.complex <- function(object) {
    if (!is.list(object) || !inherits(object, "basin_complex")) {
        .stop.basin.complex(
            "'object' must inherit from 'basin_complex'.",
            "object",
            class = "gflow_basin_schema_error"
        )
    }
    required <- c(
        "method", "direction", "n.vertices", "graph.input", "field",
        "parameters", "provenance", "extrema", "basin.table", "membership",
        "assignment", "merge.table", "trajectory.forest", "cell.complex",
        "diagnostics", "warnings", "status", "raw"
    )
    if (!all(required %in% names(object))) {
        .stop.basin.complex(
            "The basin complex is missing required top-level fields.",
            "object",
            class = "gflow_basin_schema_error",
            details = list(missing = setdiff(required, names(object)))
        )
    }
    method <- .basin.assert.choice(
        object$method,
        .basin.complex.methods,
        "object$method"
    )
    .basin.assert.choice(
        object$direction,
        .basin.complex.directions,
        "object$direction"
    )
    .basin.assert.choice(
        object$status,
        c("ok", "partial", "failed"),
        "object$status"
    )
    expected.class <- c(
        paste0("basin_complex_", method),
        "basin_complex",
        "gflow_basin_complex",
        "list"
    )
    if (!identical(class(object), expected.class)) {
        .stop.basin.complex(
            "Canonical basin-complex classes are missing or out of order.",
            "object",
            class = "gflow_basin_schema_error"
        )
    }
    if (!is.numeric(object$n.vertices) || length(object$n.vertices) != 1L ||
        object$n.vertices < 1L ||
        object$n.vertices != length(object$graph.input$adj.list)) {
        .stop.basin.complex(
            "'n.vertices' is inconsistent with 'graph.input'.",
            "object$n.vertices",
            class = "gflow_basin_schema_error"
        )
    }
    if (length(object$field$input.values) != object$n.vertices ||
        length(object$field$construction.values) != object$n.vertices) {
        .stop.basin.complex(
            "Stored field vectors are inconsistent with 'n.vertices'.",
            "object$field",
            class = "gflow_basin_schema_error"
        )
    }
    if (!is.null(object$field$vertex.mass.normalized)) {
        mass <- object$field$vertex.mass.normalized
        if (length(mass) != object$n.vertices || any(!is.finite(mass)) ||
            any(mass < 0) ||
            !isTRUE(all.equal(sum(mass), 1, tolerance = 1e-12))) {
            .stop.basin.complex(
                "Stored normalized vertex mass is invalid.",
                "object$field$vertex.mass.normalized",
                class = "gflow_basin_schema_error"
            )
        }
    }

    schemas <- .basin.table.schemas()
    for (name in c(
        "extrema", "basin.table", "membership", "assignment", "diagnostics"
    )) {
        .validate.basin.table.schema(object[[name]], schemas[[name]], name)
    }
    .validate.basin.table.schema(
        object$provenance$refinement.stages,
        schemas$refinement.stages,
        "provenance$refinement.stages"
    )
    if (!is.null(object$merge.table)) {
        .validate.basin.table.schema(
            object$merge.table,
            schemas$merge.table,
            "merge.table"
        )
    }
    if (object$status == "failed" && nrow(object$diagnostics) == 0L) {
        .stop.basin.complex(
            "A failed basin complex must contain at least one diagnostic.",
            "object$diagnostics",
            class = "gflow_basin_schema_error"
        )
    }
    if (object$status == "ok") {
        requested <- if (object$direction == "both") {
            c("max", "min")
        } else {
            object$direction
        }
        expected.rows <- object$n.vertices * length(requested)
        if (nrow(object$assignment) != expected.rows) {
            .stop.basin.complex(
                "A successful basin complex must have one assignment row per requested vertex and direction.",
                "object$assignment",
                class = "gflow_basin_schema_error"
            )
        }
        .validate.successful.basin.tables(object)
    }
    invisible(TRUE)
}

.new.failed.basin.complex <- function(method,
                                      direction,
                                      graph.input,
                                      field,
                                      parameters,
                                      condition.class,
                                      message,
                                      affected.vertices = integer()) {
    merge.table <- if (method == "superlevel_merge_tree") {
        .empty.merge.table()
    } else {
        NULL
    }
    trajectory.forest <- if (method == "trajectory_flow") {
        list(
            modulation = parameters$method.params$modulation,
            long.edge.fallback = parameters$method.params$long.edge.fallback,
            next.vertex = integer(),
            root.vertex = integer(),
            trajectories = list(),
            trajectory.status = "failed",
            tie.breaking = parameters$method.params$tie.breaking,
            tie.seed = parameters$method.params$tie.seed
        )
    } else {
        NULL
    }
    cell.complex <- if (method == "overlap_cell_complex") {
        list(
            cells = list(),
            basin.intersections = list(),
            complex.graph = NULL,
            cell.summary = data.frame(),
            cluster.assignments = list(),
            cluster.mappings = list(),
            simplified.field = numeric()
        )
    } else {
        NULL
    }
    provenance <- list(
        package.version = .basin.package.version(),
        method.backend = NA_character_,
        construction = list(
            status = "not_implemented",
            completed = FALSE,
            membership.weight.policy = NA_character_
        ),
        refinement.stages = .empty.refinement.stages(),
        random.seed = .basin.or(parameters$method.params$tie.seed, NULL)
    )
    .new.basin.complex(
        method = method,
        direction = direction,
        graph.input = graph.input,
        field = field,
        parameters = parameters,
        provenance = provenance,
        status = "failed",
        merge.table = merge.table,
        trajectory.forest = trajectory.forest,
        cell.complex = cell.complex,
        diagnostics = .basin.diagnostic(
            stage = "construction",
            method = method,
            condition.class = condition.class,
            message = message,
            affected.vertices = affected.vertices
        )
    )
}

#' Create a Canonical Basin Complex
#'
#' Creates the canonical schema used by all basin-construction methods in
#' \pkg{gflow}. Trajectory-flow and geodesic-reachability methods adapt the
#' corresponding legacy backends; methods whose adapters are not yet available
#' return a structured recoverable-failure object.
#'
#' @param adj.list Undirected graph adjacency list using 1-based vertex ids.
#' @param edge.length.list Numeric edge-length vectors parallel to
#'   `adj.list`.
#' @param field Finite scalar field over graph vertices.
#' @param method Basin-construction family.
#' @param direction One of `"max"`, `"min"`, or `"both"`.
#' @param vertex.mass Optional finite nonnegative empirical or quadrature mass.
#'   The total must be positive. No uniform fallback is applied when omitted.
#' @param vertex.density Optional finite nonnegative vertex density used by
#'   density-modulated trajectory methods.
#' @param graph.params Named graph-validation parameters.
#' @param method.params Named parameters for the selected construction method.
#' @param simplify.params Named post-construction refinement parameters.
#' @param verbose Logical scalar controlling backend progress output.
#'
#' @return A `basin_complex`. Implemented adapters return `status = "ok"`.
#'   A valid input for a method whose adapter is unavailable returns
#'   `status = "failed"` with a structured
#'   `gflow_basin_backend_not_implemented` diagnostic.
#'
#' @export
create.basin.complex <- function(
    adj.list,
    edge.length.list,
    field,
    method = c(
        "trajectory_flow",
        "superlevel_merge_tree",
        "geodesic_reachability",
        "rtcb",
        "overlap_cell_complex"
    ),
    direction = "both",
    vertex.mass = NULL,
    vertex.density = NULL,
    graph.params = list(
        edge.length.symmetry.tolerance = sqrt(.Machine$double.eps)
    ),
    method.params = list(),
    simplify.params = list(),
    verbose = FALSE
) {
    method <- .basin.assert.choice(
        method,
        .basin.complex.methods,
        "method"
    )
    direction <- .basin.assert.choice(
        direction,
        .basin.complex.directions,
        "direction"
    )
    verbose <- .basin.assert.logical(verbose, "verbose")
    if (method == "overlap_cell_complex" && direction != "both") {
        .stop.basin.complex(
            "'overlap_cell_complex' requires direction = 'both'.",
            "direction"
        )
    }

    graph <- .basin.validate.graph(
        adj.list,
        edge.length.list,
        graph.params
    )
    n <- length(graph$adj.list)
    field <- .basin.assert.vector(field, "field", n)
    if (!is.null(vertex.mass)) {
        vertex.mass <- .basin.assert.vector(
            vertex.mass,
            "vertex.mass",
            n,
            nonnegative = TRUE,
            positive.total = TRUE
        )
        mass.total <- sum(vertex.mass)
        normalized.mass <- vertex.mass / mass.total
    } else {
        mass.total <- NA_real_
        normalized.mass <- NULL
    }
    if (!is.null(vertex.density)) {
        vertex.density <- .basin.assert.vector(
            vertex.density,
            "vertex.density",
            n,
            nonnegative = TRUE
        )
    }

    resolved.method.params <- .basin.resolve.method.params(
        method,
        method.params,
        n,
        vertex.density
    )
    resolved.simplify.params <- .basin.resolve.simplify.params(
        method,
        simplify.params
    )
    tie.breaking <- identical(resolved.method.params$tie.breaking, TRUE)
    tie.seed <- .basin.or(resolved.method.params$tie.seed, NULL)
    field.record <- list(
        input.values = field,
        construction.values = field,
        vertex.mass.input = vertex.mass,
        vertex.mass.normalized = normalized.mass,
        vertex.mass.input.total = mass.total,
        vertex.density = vertex.density,
        tie.policy = list(
            enabled = tie.breaking,
            seed = tie.seed,
            applied = FALSE,
            perturbation.scale = 0
        )
    )
    graph.input <- list(
        adj.list = graph$adj.list,
        edge.length.list = graph$edge.length.list,
        graph.params = graph$graph.params,
        validation = graph$validation
    )
    parameters <- list(
        method.params = resolved.method.params,
        simplify.params = resolved.simplify.params,
        verbose = verbose
    )
    if (method == "geodesic_reachability") {
        return(.create.geodesic.reachability.complex(
            direction = direction,
            graph.input = graph.input,
            field = field.record,
            parameters = parameters
        ))
    }
    if (method == "trajectory_flow") {
        return(.create.trajectory.flow.complex(
            direction = direction,
            graph.input = graph.input,
            field = field.record,
            parameters = parameters,
            verbose = verbose
        ))
    }
    .new.failed.basin.complex(
        method = method,
        direction = direction,
        graph.input = graph.input,
        field = field.record,
        parameters = parameters,
        condition.class = "gflow_basin_backend_not_implemented",
        message = sprintf(
            "The canonical '%s' backend is not implemented in Phase B.",
            method
        )
    )
}

.basin.accessor.object <- function(object) {
    if (!inherits(object, "basin_complex")) {
        .stop.basin.complex(
            "'object' must inherit from 'basin_complex'.",
            "object"
        )
    }
    .validate.basin.complex(object)
    object
}

#' Extract the Canonical Basin Table
#'
#' @param object A canonical `basin_complex`.
#'
#' @return A canonical basin table.
#' @export
get.basin.table <- function(object) {
    .basin.accessor.object(object)$basin.table
}

#' Extract Canonical Basin Membership
#'
#' @inheritParams get.basin.table
#' @return A canonical raw-membership table.
#' @export
get.basin.membership <- function(object) {
    .basin.accessor.object(object)$membership
}

#' Extract Canonical Primary Basin Assignments
#'
#' @inheritParams get.basin.table
#' @return A canonical primary-assignment table.
#' @export
get.basin.assignment <- function(object) {
    .basin.accessor.object(object)$assignment
}

.basin.unsupported.accessor <- function(name, object, required) {
    required <- .basin.assert.logical(required, "required")
    if (required) {
        .stop.basin.complex(
            sprintf(
                "Accessor '%s' is not available for basin method '%s'.",
                name,
                object$method
            ),
            "object"
        )
    }
    NULL
}

#' Extract a Canonical Basin Merge Tree
#'
#' @inheritParams get.basin.table
#' @param required If `TRUE`, error when the selected method has no merge tree.
#' @return A canonical merge table or `NULL`.
#' @export
get.basin.merge.tree <- function(object, required = FALSE) {
    object <- .basin.accessor.object(object)
    if (is.null(object$merge.table)) {
        return(.basin.unsupported.accessor(
            "get.basin.merge.tree",
            object,
            required
        ))
    }
    object$merge.table
}

#' Extract a Canonical Basin Trajectory Forest
#'
#' @inheritParams get.basin.merge.tree
#' @return A trajectory-forest record or `NULL`.
#' @export
get.basin.trajectory.forest <- function(object, required = FALSE) {
    object <- .basin.accessor.object(object)
    if (is.null(object$trajectory.forest)) {
        return(.basin.unsupported.accessor(
            "get.basin.trajectory.forest",
            object,
            required
        ))
    }
    object$trajectory.forest
}

#' Extract Canonical Basin Cells
#'
#' @inheritParams get.basin.merge.tree
#' @return A canonical cell-complex record or `NULL`.
#' @export
get.basin.cells <- function(object, required = FALSE) {
    object <- .basin.accessor.object(object)
    if (is.null(object$cell.complex)) {
        return(.basin.unsupported.accessor(
            "get.basin.cells",
            object,
            required
        ))
    }
    object$cell.complex
}

#' Convert an Object to a Canonical Basin Complex
#'
#' @param object Object to convert.
#' @param ... Additional arguments passed to methods.
#'
#' @return A canonical `basin_complex`.
#' @export
as.basin.complex <- function(object, ...) {
    UseMethod("as.basin.complex")
}

#' @method as.basin.complex basin_complex
#' @export
as.basin.complex.basin_complex <- function(object, ...) {
    .validate.basin.complex(object)
    object
}

#' @method as.basin.complex default
#' @export
as.basin.complex.default <- function(object, ...) {
    .stop.basin.complex(
        sprintf(
            "No basin-complex adapter is implemented for class: %s.",
            paste(class(object), collapse = "/")
        ),
        "object",
        class = "gflow_basin_conversion_error"
    )
}

#' @export
print.basin_complex <- function(x, ...) {
    cat("<basin_complex>\n")
    cat("  method:    ", x$method, "\n", sep = "")
    cat("  direction: ", x$direction, "\n", sep = "")
    cat("  vertices:  ", x$n.vertices, "\n", sep = "")
    cat("  status:    ", x$status, "\n", sep = "")
    cat("  basins:    ", nrow(x$basin.table), "\n", sep = "")
    if (x$status != "ok" && nrow(x$diagnostics) > 0L) {
        cat("  diagnostic:", x$diagnostics$message[[1L]], "\n")
    }
    invisible(x)
}

#' @export
summary.basin_complex <- function(object, ...) {
    structure(
        list(
            method = object$method,
            direction = object$direction,
            status = object$status,
            n.vertices = object$n.vertices,
            n.components = object$graph.input$validation$n.components,
            n.basins = nrow(object$basin.table),
            n.memberships = nrow(object$membership),
            n.assignments = nrow(object$assignment),
            n.diagnostics = nrow(object$diagnostics),
            has.vertex.mass = !is.null(object$field$vertex.mass.normalized)
        ),
        class = c("summary.basin_complex", "list")
    )
}

#' @export
print.summary.basin_complex <- function(x, ...) {
    cat("Canonical Basin Complex Summary\n")
    cat("  Method: ", x$method, "\n", sep = "")
    cat("  Direction: ", x$direction, "\n", sep = "")
    cat("  Status: ", x$status, "\n", sep = "")
    cat(
        "  Vertices/components: ",
        x$n.vertices,
        "/",
        x$n.components,
        "\n",
        sep = ""
    )
    cat(
        "  Basins/memberships/assignments: ",
        x$n.basins,
        "/",
        x$n.memberships,
        "/",
        x$n.assignments,
        "\n",
        sep = ""
    )
    invisible(x)
}

#' @export
plot.basin_complex <- function(x,
                               xlab = "Vertex",
                               ylab = "Field",
                               main = "Canonical Basin Complex Field",
                               ...) {
    graphics::plot(
        seq_len(x$n.vertices),
        x$field$input.values,
        xlab = xlab,
        ylab = ylab,
        main = main,
        ...
    )
    invisible(x)
}

#' @export
as.data.frame.basin_complex <- function(x,
                                        row.names = NULL,
                                        optional = FALSE,
                                        ...) {
    as.data.frame(
        x$basin.table,
        row.names = row.names,
        optional = optional,
        ...
    )
}

#' @export
basin.basin_complex <- function(object,
                                id,
                                stage = c("retained", "raw", "primary"),
                                type = c("full", "vertices"),
                                ...) {
    stage <- match.arg(stage)
    type <- match.arg(type)
    if (!is.character(id) || length(id) != 1L || is.na(id)) {
        .stop.basin.complex("'id' must be one basin id.", "id")
    }
    index <- which(object$basin.table$basin.id == id)
    if (length(index) != 1L) {
        .stop.basin.complex(
            sprintf("Basin id '%s' was not found exactly once.", id),
            "id"
        )
    }
    if (type == "full") {
        return(object$basin.table[index, , drop = FALSE])
    }
    object$basin.table[[paste0(stage, ".support.vertices")]][[index]]
}

#' @export
basins.basin_complex <- function(object,
                                 stage = c("retained", "raw", "primary"),
                                 direction = NULL,
                                 ...) {
    stage <- match.arg(stage)
    table <- object$basin.table
    if (!is.null(direction)) {
        direction <- .basin.assert.choice(
            direction,
            c("max", "min"),
            "direction"
        )
        table <- table[table$type == direction, , drop = FALSE]
    }
    vertices <- unclass(table[[paste0(stage, ".support.vertices")]])
    if (length(vertices) > 0L) {
        names(vertices) <- table$basin.id
    }
    vertices
}
