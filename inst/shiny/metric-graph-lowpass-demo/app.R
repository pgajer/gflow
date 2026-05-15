if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("The shiny package is required to run this demo.", call. = FALSE)
}
if (!requireNamespace("gflow", quietly = TRUE)) {
    stop("The gflow package must be installed or loaded to run this demo.",
         call. = FALSE)
}

mixture.weights <- function(weights, n.components) {
    weights <- as.numeric(weights[seq_len(n.components)])
    weights[!is.finite(weights) | weights < 0] <- 0
    if (sum(weights) <= 0) {
        rep(1 / n.components, n.components)
    } else {
        weights / sum(weights)
    }
}

gaussian.mixture <- function(x, means, sds, weights, scale.to.max = TRUE) {
    n.components <- length(weights)
    y <- numeric(length(x))
    for (j in seq_len(n.components)) {
        y <- y + weights[j] * stats::dnorm(x, mean = means[j], sd = sds[j])
    }
    if (isTRUE(scale.to.max) && max(y) > 0) {
        y <- y / max(y)
    }
    y
}

sample.unit.interval <- function(n, type, mean, sd) {
    if (identical(type, "uniform")) {
        return(sort(stats::runif(n)))
    }
    lower.p <- stats::pnorm(0, mean = mean, sd = sd)
    upper.p <- stats::pnorm(1, mean = mean, sd = sd)
    sort(stats::qnorm(stats::runif(n, lower.p, upper.p), mean = mean, sd = sd))
}

rlaplace <- function(n, scale) {
    u <- stats::runif(n) - 0.5
    -scale * sign(u) * log1p(-2 * abs(u))
}

make.synthetic.data <- function(input, resample.count) {
    n.components <- as.integer(input$n_components)
    means <- c(input$mean1, input$mean2, input$mean3)[seq_len(n.components)]
    sds <- pmax(c(input$sd1, input$sd2, input$sd3)[seq_len(n.components)], 1e-6)
    weights <- mixture.weights(c(input$weight1, input$weight2, input$weight3),
                               n.components)

    seed <- as.integer(input$seed) + as.integer(resample.count)
    set.seed(seed)

    x <- sample.unit.interval(
        n = as.integer(input$sample_size),
        type = input$sampling_distribution,
        mean = input$sampling_mean,
        sd = max(input$sampling_sd, 1e-6)
    )
    y.truth.sample <- gaussian.mixture(
        x = x,
        means = means,
        sds = sds,
        weights = weights,
        scale.to.max = input$scale_truth
    )
    noise.scale <- max(input$noise_scale, 0)
    noise <- switch(
        input$noise_type,
        gaussian = stats::rnorm(length(x), sd = noise.scale),
        laplace = rlaplace(length(x), scale = noise.scale)
    )
    y <- y.truth.sample + noise

    x.grid <- seq(0, 1, length.out = 401)
    y.truth.grid <- gaussian.mixture(
        x = x.grid,
        means = means,
        sds = sds,
        weights = weights,
        scale.to.max = input$scale_truth
    )

    list(
        x = x,
        y = y,
        y.truth.sample = y.truth.sample,
        x.grid = x.grid,
        y.truth.grid = y.truth.grid,
        means = means,
        sds = sds,
        weights = weights,
        seed = seed
    )
}

build.knn.graph.1d <- function(x, k, edge.length.floor = 1e-8) {
    n <- length(x)
    k <- min(max(as.integer(k), 1L), n - 1L)
    ord <- order(x)
    rank.of <- integer(n)
    rank.of[ord] <- seq_len(n)

    edges <- matrix(integer(0), ncol = 2)
    for (i in seq_len(n)) {
        rank.i <- rank.of[i]
        candidate.ranks <- seq.int(max(1L, rank.i - k), min(n, rank.i + k))
        candidates <- ord[candidate.ranks]
        candidates <- candidates[candidates != i]
        if (length(candidates) > k) {
            candidates <- candidates[order(abs(x[candidates] - x[i]))[seq_len(k)]]
        }
        if (length(candidates) > 0) {
            edges <- rbind(edges, cbind(rep(i, length(candidates)), candidates))
        }
    }

    edges <- t(apply(edges, 1L, sort))
    edges <- unique(edges)
    adj <- vector("list", n)
    weights <- vector("list", n)
    for (i in seq_len(n)) {
        adj[[i]] <- integer(0)
        weights[[i]] <- numeric(0)
    }
    for (e in seq_len(nrow(edges))) {
        i <- edges[e, 1L]
        j <- edges[e, 2L]
        len <- max(abs(x[i] - x[j]), edge.length.floor)
        adj[[i]] <- c(adj[[i]], j)
        weights[[i]] <- c(weights[[i]], len)
        adj[[j]] <- c(adj[[j]], i)
        weights[[j]] <- c(weights[[j]], len)
    }

    for (i in seq_len(n)) {
        if (length(adj[[i]]) > 0) {
            o <- order(adj[[i]])
            adj[[i]] <- as.integer(adj[[i]][o])
            weights[[i]] <- as.numeric(weights[[i]][o])
        }
    }
    list(adj.list = adj, weight.list = weights, edge.count = nrow(edges))
}

eta.grid.from.input <- function(input) {
    if (identical(input$eta_mode, "auto")) {
        return(NULL)
    }
    min.eta <- max(input$eta_min, .Machine$double.eps)
    max.eta <- max(input$eta_max, min.eta * 1.0001)
    exp(seq(log(min.eta), log(max.eta), length.out = as.integer(input$eta_count)))
}

fit.metric.lowpass.demo <- function(data, input) {
    graph <- build.knn.graph.1d(
        data$x,
        k = input$knn_k,
        edge.length.floor = input$edge_length_floor
    )

    gflow::fit.metric.graph.lowpass(
        adj.list = graph$adj.list,
        weight.list = graph$weight.list,
        y = data$y,
        conductance.rule = input$conductance_rule,
        conductance.epsilon = input$conductance_epsilon,
        conductance.alpha = input$conductance_alpha,
        conductance.sigma = if (isTRUE(input$use_auto_sigma)) NULL else input$conductance_sigma,
        conductance.sigma.rule = input$sigma_rule,
        conductance.sigma.quantile = input$sigma_quantile,
        conductance.local.k = input$local_k,
        n.eigenpairs = min(as.integer(input$n_eigenpairs), length(data$x) - 1L),
        filter.type = input$filter_type,
        eta.grid = eta.grid.from.input(input),
        n.candidates = as.integer(input$n_candidates),
        eigen.solver = input$eigen_solver,
        dense.eigen.threshold = as.integer(input$dense_eigen_threshold),
        dense.fallback.threshold = as.integer(input$dense_fallback_threshold),
        dense.fallback = input$dense_fallback,
        verbose = FALSE
    )
}

range.finite <- function(...) {
    values <- unlist(list(...), use.names = FALSE)
    values <- values[is.finite(values)]
    if (length(values) == 0) {
        c(0, 1)
    } else {
        range(values)
    }
}

data.signature <- function(input, resample.count) {
    paste(
        input$n_components,
        input$mean1, input$mean2, input$mean3,
        input$sd1, input$sd2, input$sd3,
        input$weight1, input$weight2, input$weight3,
        input$scale_truth,
        input$sample_size,
        input$sampling_distribution,
        input$sampling_mean,
        input$sampling_sd,
        input$seed,
        input$noise_type,
        input$noise_scale,
        resample.count,
        sep = "|"
    )
}

model.signature <- function(input) {
    paste(
        input$knn_k,
        input$edge_length_floor,
        input$conductance_rule,
        input$conductance_epsilon,
        input$conductance_alpha,
        input$use_auto_sigma,
        input$conductance_sigma,
        input$sigma_rule,
        input$sigma_quantile,
        input$local_k,
        input$n_eigenpairs,
        input$filter_type,
        input$eta_mode,
        input$eta_min,
        input$eta_max,
        input$eta_count,
        input$n_candidates,
        input$eigen_solver,
        input$dense_eigen_threshold,
        input$dense_fallback,
        input$dense_fallback_threshold,
        sep = "|"
    )
}

ui <- shiny::fluidPage(
    shiny::tags$head(
        shiny::tags$style(shiny::HTML("
            body { background: #f7f7f4; color: #1f2933; }
            .control-card {
                background: #ffffff;
                border: 1px solid #d9ded7;
                border-radius: 6px;
                padding: 14px 16px;
                margin-bottom: 14px;
                box-shadow: 0 1px 2px rgba(15, 23, 42, 0.04);
            }
            .control-card h3 {
                margin-top: 0;
                margin-bottom: 12px;
                font-size: 18px;
                font-weight: 650;
            }
            .btn-primary {
                background-color: #1f6f8b;
                border-color: #1f6f8b;
            }
            .status-line {
                font-weight: 600;
                margin: 8px 0 12px;
            }
            .plot-card {
                background: #ffffff;
                border: 1px solid #d9ded7;
                border-radius: 6px;
                padding: 12px;
            }
            pre {
                background: #101820;
                color: #f4f7f5;
                border-radius: 6px;
            }
        "))
    ),
    shiny::titlePanel("Metric Graph Low-Pass Regression: 1D Gaussian Mixture"),
    shiny::fluidRow(
        shiny::column(
            width = 4,
            shiny::div(
                class = "control-card",
                shiny::h3("Gaussian Mixture"),
                shiny::selectInput("n_components", "Components", choices = 1:3, selected = 2),
                shiny::checkboxInput("scale_truth", "Scale truth curve to maximum 1", TRUE),
                shiny::sliderInput("mean1", "Component 1 mean", 0, 1, 0.28, step = 0.01),
                shiny::sliderInput("sd1", "Component 1 standard deviation", 0.005, 0.25, 0.055, step = 0.005),
                shiny::sliderInput("weight1", "Component 1 weight", 0, 2, 1, step = 0.05),
                shiny::conditionalPanel(
                    "Number(input.n_components) >= 2",
                    shiny::sliderInput("mean2", "Component 2 mean", 0, 1, 0.68, step = 0.01),
                    shiny::sliderInput("sd2", "Component 2 standard deviation", 0.005, 0.25, 0.08, step = 0.005),
                    shiny::sliderInput("weight2", "Component 2 weight", 0, 2, 0.8, step = 0.05)
                ),
                shiny::conditionalPanel(
                    "Number(input.n_components) >= 3",
                    shiny::sliderInput("mean3", "Component 3 mean", 0, 1, 0.48, step = 0.01),
                    shiny::sliderInput("sd3", "Component 3 standard deviation", 0.005, 0.25, 0.04, step = 0.005),
                    shiny::sliderInput("weight3", "Component 3 weight", 0, 2, 0.5, step = 0.05)
                )
            ),
            shiny::div(
                class = "control-card",
                shiny::h3("Sampling and Noise"),
                shiny::sliderInput("sample_size", "Sample size", 20, 1200, 250, step = 10),
                shiny::selectInput(
                    "sampling_distribution",
                    "Sampling distribution on [0, 1]",
                    choices = c("uniform", "gaussian"),
                    selected = "uniform"
                ),
                shiny::conditionalPanel(
                    "input.sampling_distribution == 'gaussian'",
                    shiny::sliderInput("sampling_mean", "Sampling mean", 0, 1, 0.5, step = 0.01),
                    shiny::sliderInput("sampling_sd", "Sampling standard deviation", 0.02, 0.5, 0.18, step = 0.01)
                ),
                shiny::selectInput("noise_type", "Noise type", choices = c("gaussian", "laplace")),
                shiny::sliderInput("noise_scale", "Noise sd / Laplace scale", 0, 1, 0.08, step = 0.01),
                shiny::numericInput("seed", "Random seed", value = 2026, min = 1, step = 1),
                shiny::actionButton("resample", "Resample", class = "btn-default")
            )
        ),
        shiny::column(
            width = 4,
            shiny::div(
                class = "control-card",
                shiny::h3("Graph Construction"),
                shiny::sliderInput("knn_k", "Symmetric 1D kNN graph: k", 1, 30, 5, step = 1),
                shiny::numericInput("edge_length_floor", "Edge length floor", value = 1e-8, min = 1e-12, step = 1e-8),
                shiny::helpText("Graph edges use metric lengths abs(x_i - x_j). These lengths are passed to fit.metric.graph.lowpass(), not precomputed conductances.")
            ),
            shiny::div(
                class = "control-card",
                shiny::h3("fit.metric.graph.lowpass() Parameters"),
                shiny::selectInput(
                    "conductance_rule",
                    "Conductance rule",
                    choices = c(
                        "inverse.length.power",
                        "exp.length",
                        "exp.length.squared",
                        "self.tuned.gaussian"
                    ),
                    selected = "inverse.length.power"
                ),
                shiny::numericInput("conductance_epsilon", "Conductance epsilon", value = 1e-8, min = 1e-12, step = 1e-8),
                shiny::conditionalPanel(
                    "input.conductance_rule == 'inverse.length.power'",
                    shiny::sliderInput("conductance_alpha", "Inverse length power alpha", 0.1, 4, 1, step = 0.1)
                ),
                shiny::conditionalPanel(
                    "input.conductance_rule == 'exp.length' || input.conductance_rule == 'exp.length.squared'",
                    shiny::checkboxInput("use_auto_sigma", "Use automatic sigma", TRUE),
                    shiny::conditionalPanel(
                        "!input.use_auto_sigma",
                        shiny::numericInput("conductance_sigma", "Conductance sigma", value = 0.08, min = 1e-6, step = 0.01)
                    ),
                    shiny::selectInput("sigma_rule", "Automatic sigma rule", choices = c("edge.quantile", "median")),
                    shiny::sliderInput("sigma_quantile", "Sigma edge quantile", 0.05, 0.95, 0.75, step = 0.05)
                ),
                shiny::conditionalPanel(
                    "input.conductance_rule == 'self.tuned.gaussian'",
                    shiny::sliderInput("local_k", "Self-tuned local k", 1, 30, 5, step = 1)
                ),
                shiny::selectInput(
                    "filter_type",
                    "Filter type",
                    choices = c("heat_kernel", "tikhonov", "cubic_spline",
                                "gaussian", "exponential", "butterworth"),
                    selected = "heat_kernel"
                ),
                shiny::sliderInput("n_eigenpairs", "Eigenpairs", 5, 250, 50, step = 5),
                shiny::selectInput("eta_mode", "Eta grid", choices = c("auto", "manual"), selected = "auto"),
                shiny::conditionalPanel(
                    "input.eta_mode == 'manual'",
                    shiny::numericInput("eta_min", "Eta min", value = 1e-4, min = 1e-12, step = 1e-4),
                    shiny::numericInput("eta_max", "Eta max", value = 1e2, min = 1e-12, step = 1),
                    shiny::sliderInput("eta_count", "Eta count", 5, 120, 40, step = 5)
                ),
                shiny::conditionalPanel(
                    "input.eta_mode == 'auto'",
                    shiny::sliderInput("n_candidates", "Auto eta candidates", 5, 120, 40, step = 5)
                ),
                shiny::helpText("Phase 1 uses the unnormalized weighted Laplacian L = D - C."),
                shiny::actionButton("fit", "Fit metric low-pass", class = "btn-primary")
            )
        ),
        shiny::column(
            width = 4,
            shiny::div(
                class = "control-card",
                shiny::h3("Solver"),
                shiny::selectInput("eigen_solver", "Eigen solver", choices = c("auto", "sparse", "dense")),
                shiny::numericInput("dense_eigen_threshold", "Dense eigen threshold", value = 200, min = 1, step = 10),
                shiny::selectInput("dense_fallback", "Dense fallback", choices = c("auto", "never", "always")),
                shiny::numericInput("dense_fallback_threshold", "Dense fallback threshold", value = 5000, min = 1, step = 100)
            ),
            shiny::div(
                class = "control-card",
                shiny::h3("Status"),
                shiny::uiOutput("fit_status"),
                shiny::verbatimTextOutput("diagnostics")
            )
        )
    ),
    shiny::fluidRow(
        shiny::column(
            width = 12,
            shiny::div(
                class = "plot-card",
                shiny::plotOutput("mixture_plot", height = "620px")
            )
        )
    )
)

server <- function(input, output, session) {
    resample.count <- shiny::reactiveVal(0L)
    fit.state <- shiny::reactiveVal(NULL)
    fit.error <- shiny::reactiveVal(NULL)

    shiny::observeEvent(input$resample, {
        resample.count(resample.count() + 1L)
    })

    synthetic.data <- shiny::reactive({
        make.synthetic.data(input, resample.count())
    })

    current.data.signature <- shiny::reactive({
        data.signature(input, resample.count())
    })

    current.model.signature <- shiny::reactive({
        model.signature(input)
    })

    shiny::observeEvent(input$fit, {
        fit.error(NULL)
        data <- synthetic.data()
        shiny::withProgress(message = "Fitting metric graph low-pass model", value = 0.25, {
            result <- tryCatch(
                {
                    shiny::incProgress(0.5)
                    fit.metric.lowpass.demo(data, input)
                },
                error = function(e) e
            )
            shiny::incProgress(0.25)
            if (inherits(result, "error")) {
                fit.state(NULL)
                fit.error(conditionMessage(result))
            } else {
                fit.state(list(
                    fit = result,
                    data.signature = current.data.signature(),
                    model.signature = current.model.signature()
                ))
            }
        })
    })

    fit.is.current <- shiny::reactive({
        state <- fit.state()
        !is.null(state) &&
            identical(state$data.signature, current.data.signature()) &&
            identical(state$model.signature, current.model.signature())
    })

    output$fit_status <- shiny::renderUI({
        if (!is.null(fit.error())) {
            return(shiny::div(class = "status-line", style = "color: #b42318;",
                              paste("Fit failed:", fit.error())))
        }
        state <- fit.state()
        if (is.null(state)) {
            return(shiny::div(class = "status-line", "No fit has been run yet."))
        }
        if (!fit.is.current()) {
            return(shiny::div(class = "status-line", style = "color: #9a6700;",
                              "Current fit is stale. Press Fit metric low-pass."))
        }
        shiny::div(class = "status-line", style = "color: #1a7f37;",
                   "Fit is current.")
    })

    output$mixture_plot <- shiny::renderPlot({
        data <- synthetic.data()
        state <- fit.state()
        use.fit <- fit.is.current()

        fitted.values <- if (use.fit) state$fit$fitted.values else NULL
        ylim <- range.finite(data$y.truth.grid, data$y, fitted.values)
        pad <- diff(ylim) * 0.08
        if (!is.finite(pad) || pad == 0) {
            pad <- 0.1
        }
        ylim <- ylim + c(-pad, pad)

        graphics::plot(
            data$x.grid,
            data$y.truth.grid,
            type = "l",
            lwd = 3,
            col = "#c62828",
            xlim = c(0, 1),
            ylim = ylim,
            xlab = "x",
            ylab = "response",
            main = "Metric graph low-pass fit on 1D Gaussian mixture data"
        )
        graphics::points(data$x, data$y, pch = 16, cex = 0.65, col = "#111111")
        if (isTRUE(use.fit)) {
            o <- order(data$x)
            graphics::lines(data$x[o], fitted.values[o], col = "#1565c0", lwd = 3)
        }
        graphics::legend(
            "topright",
            legend = c("noise-free mixture", "observed data",
                       if (isTRUE(use.fit)) "metric graph low-pass fit" else NULL),
            col = c("#c62828", "#111111", if (isTRUE(use.fit)) "#1565c0" else NULL),
            lwd = c(3, NA, if (isTRUE(use.fit)) 3 else NULL),
            pch = c(NA, 16, if (isTRUE(use.fit)) NA else NULL),
            bty = "n"
        )
    })

    output$diagnostics <- shiny::renderPrint({
        data <- synthetic.data()
        graph <- build.knn.graph.1d(
            data$x,
            k = input$knn_k,
            edge.length.floor = input$edge_length_floor
        )
        edge.lengths <- unlist(graph$weight.list, use.names = FALSE)

        cat("Data\n")
        cat("  n:", length(data$x), "\n")
        cat("  seed:", data$seed, "\n")
        cat("  mixture weights:", paste(signif(data$weights, 3), collapse = ", "), "\n")
        cat("  graph edges:", graph$edge.count, "\n")
        cat("  edge length summary:\n")
        print(summary(edge.lengths))

        if (!is.null(fit.error())) {
            cat("\nLast fit error:\n")
            cat(fit.error(), "\n")
            return(invisible())
        }

        state <- fit.state()
        if (is.null(state)) {
            cat("\nNo fit diagnostics yet.\n")
            return(invisible())
        }
        if (!fit.is.current()) {
            cat("\nFit diagnostics are stale for the current controls.\n")
            return(invisible())
        }

        fit <- state$fit
        conductances <- fit$operator$edge.table$conductance
        lengths <- fit$operator$edge.table$length

        cat("\nFit\n")
        cat("  eta:", signif(fit$gcv$eta.optimal, 6), "\n")
        cat("  GCV:", signif(fit$gcv$gcv.optimal, 6), "\n")
        cat("  effective df:", signif(fit$gcv$effective.df, 6), "\n")
        cat("  eigen backend:", fit$spectral$backend, "\n")
        cat("  conductance rule:", fit$conductance$rule, "\n")
        cat("  fitted range:", paste(signif(range(fit$fitted.values), 4), collapse = " to "), "\n")
        cat("  operator edge length summary:\n")
        print(summary(lengths))
        cat("  conductance summary:\n")
        print(summary(conductances))
    })
}

shiny::shinyApp(ui = ui, server = server)
