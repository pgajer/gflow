repo.dir <- "/Users/pgajer/current_projects/gflow"
report.dir <- file.path(repo.dir, "dev/metric-graph-lowpass-benchmarks")
out.dir <- file.path(report.dir, "reports")
dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)

if (requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(repo.dir, quiet = TRUE)
} else {
  library(gflow)
}

source(file.path(report.dir, "multi_scenario_1d_mixture_helpers.R"))

path_graph_to_conductance_weights <- function(graph,
                                              conductance.rule = "exp.length",
                                              conductance.sigma.rule = "edge.quantile",
                                              conductance.sigma.quantile = 0.75,
                                              conductance.floor = 1e-4,
                                              normalize = TRUE) {
  op <- metric.graph.lowpass.operator(
    adj.list = graph$adj_list,
    weight.list = graph$weight_list,
    conductance.rule = conductance.rule,
    conductance.sigma.rule = conductance.sigma.rule,
    conductance.sigma.quantile = conductance.sigma.quantile,
    laplacian.type = "unnormalized",
    return.sparse = FALSE,
    verbose = FALSE
  )

  conductance <- pmax(op$edge.table$conductance, conductance.floor)
  if (isTRUE(normalize)) {
    conductance <- conductance / stats::median(conductance)
  }

  n <- length(graph$adj_list)
  out <- vector("list", n)
  for (i in seq_len(n)) out[[i]] <- numeric()
  for (e in seq_len(nrow(op$edge.table))) {
    u <- op$edge.table$u[e]
    v <- op$edge.table$v[e]
    cval <- conductance[e]
    out[[u]] <- c(out[[u]], cval)
    out[[v]] <- c(out[[v]], cval)
  }
  lapply(out, as.double)
}

score_subset_prediction <- function(dataset.row, method, fit.out, x, y, truth,
                                    lambda = NA_real_, weight.rule = NA_character_,
                                    conductance.source = NA_character_,
                                    trend.order = NA_integer_) {
  if (inherits(fit.out$result, "error")) {
    return(data.frame(
      dataset.row,
      method = method,
      status = "error",
      trend.order = trend.order,
      rmse.truth = NA_real_,
      rmse.observed = NA_real_,
      mae.truth = NA_real_,
      runtime.sec = fit.out$elapsed,
      lambda = lambda,
      weight.rule = weight.rule,
      conductance.source = conductance.source,
      message = conditionMessage(fit.out$result),
      stringsAsFactors = FALSE
    ))
  }

  yhat <- as.numeric(fit.out$result$fitted.values)
  if (is.na(lambda) && !is.null(fit.out$result$lambda)) {
    lambda <- fit.out$result$lambda
  }
  data.frame(
    dataset.row,
    method = method,
    status = "ok",
    trend.order = trend.order,
    rmse.truth = sqrt(mean((yhat - truth)^2)),
    rmse.observed = sqrt(mean((yhat - y)^2)),
    mae.truth = mean(abs(yhat - truth)),
    runtime.sec = fit.out$elapsed,
    lambda = lambda,
    weight.rule = weight.rule,
    conductance.source = conductance.source,
    message = "",
    stringsAsFactors = FALSE
  )
}

prediction_frame <- function(dataset.row, x, y, truth, method, fit.out) {
  if (inherits(fit.out$result, "error")) return(NULL)
  data.frame(
    dataset.id = dataset.row$dataset.id,
    shape.id = dataset.row$shape.id,
    shape.name = dataset.row$shape.name,
    variant.id = dataset.row$variant.id,
    variant.name = dataset.row$variant.name,
    replicate = dataset.row$replicate,
    x = x,
    y = y,
    truth = truth,
    method = method,
    fitted = as.numeric(fit.out$result$fitted.values),
    stringsAsFactors = FALSE
  )
}

fit_graph_trend_subset <- function(graph, y, order, weight.rule, weight.list,
                                   lambda.grid, foldid, maxsteps = 300L,
                                   approx = FALSE) {
  time_expression(
    fit.graph.trend.filtering(
      adj.list = graph$adj_list,
      weight.list = weight.list,
      y = y,
      order = order,
      lambda.grid = lambda.grid,
      lambda.selection = "cv",
      weight.rule = weight.rule,
      foldid = foldid,
      maxsteps = maxsteps,
      approx = approx,
      verbose = FALSE
    )
  )
}

fit_quadratic_trendfilter_subset <- function(x, y) {
  fit_genlasso_trendfilter_cv(x, y)
}

run_graph_trend_filtering_subset <- function(n = 250L,
                                             shape.ids = c("S01", "S03", "S05", "S09",
                                                           "S10", "S11", "S14", "S16"),
                                             variant.ids = c("V1", "V2", "V3"),
                                             replicates = 1L,
                                             lambda.grid.order0 = 10^seq(0.7, -2, length.out = 10L),
                                             lambda.grid.higher = 10^seq(0.7, -1, length.out = 6L),
                                             progress = TRUE) {
  shapes <- shape_catalog()
  variants <- variant_catalog()
  scenario.grid <- make_scenario_grid(
    replicates = replicates,
    shape.ids = shape.ids,
    variant.ids = variant.ids
  )
  scenario.grid <- scenario.grid[order(scenario.grid$shape.id, scenario.grid$variant.id,
                                       scenario.grid$replicate), ]
  rownames(scenario.grid) <- NULL

  graph.fit.grid <- data.frame(
    order = c(0L, 0L, 0L, 1L, 2L, 2L),
    weight.rule = c("unit", "conductance", "sqrt.conductance",
                    "unit", "unit", "sqrt.conductance"),
    conductance.source = c("none",
                           "exp.length sigma=edge.q0.75 floor=1e-4 median-normalized",
                           "exp.length sigma=edge.q0.75 floor=1e-4 median-normalized",
                           "none",
                           "none",
                           "exp.length sigma=edge.q0.75 floor=1e-4 median-normalized"),
    stringsAsFactors = FALSE
  )

  results <- vector("list", nrow(scenario.grid) * (nrow(graph.fit.grid) + 1L))
  predictions <- vector("list", nrow(scenario.grid) * (nrow(graph.fit.grid) + 1L))
  idx <- 0L
  pidx <- 0L

  for (s in seq_len(nrow(scenario.grid))) {
    if (isTRUE(progress)) {
      message(sprintf("[%d/%d] %s", s, nrow(scenario.grid), scenario.grid$dataset.id[s]))
    }
    shape <- shapes[match(scenario.grid$shape.id[s], shapes$shape.id), ]
    variant <- variants[match(scenario.grid$variant.id[s], variants$variant.id), ]
    dat <- make_dataset(shape, variant, scenario.grid$replicate[s], n = n)
    graph <- make_path_graph(dat$x)
    foldid <- rep(seq_len(5L), length.out = n)

    quad <- fit_quadratic_trendfilter_subset(dat$x, dat$y)
    idx <- idx + 1L
    results[[idx]] <- score_subset_prediction(
      scenario.grid[s, ], "genlasso.trendfilter.ord2.cv", quad,
      dat$x, dat$y, dat$truth,
      lambda = if (!inherits(quad$result, "error")) quad$result$eta else NA_real_,
      trend.order = 2L
    )
    pidx <- pidx + 1L
    predictions[[pidx]] <- prediction_frame(
      scenario.grid[s, ], dat$x, dat$y, dat$truth,
      "genlasso.trendfilter.ord2.cv", quad
    )

    conductance.weights <- path_graph_to_conductance_weights(graph)

    for (gg in seq_len(nrow(graph.fit.grid))) {
      spec <- graph.fit.grid[gg, ]
      weight.list <- if (identical(spec$conductance.source, "none")) {
        graph$weight_list
      } else {
        conductance.weights
      }
      lambda.grid <- if (spec$order == 0L) lambda.grid.order0 else lambda.grid.higher
      fit.out <- fit_graph_trend_subset(
        graph = graph,
        y = dat$y,
        order = spec$order,
        weight.rule = spec$weight.rule,
        weight.list = weight.list,
        lambda.grid = lambda.grid,
        foldid = foldid,
        maxsteps = if (spec$order == 0L) 300L else 300L,
        approx = spec$order > 0L
      )
      method <- paste0("fit.graph.trend.filtering.order", spec$order, ".", spec$weight.rule)
      idx <- idx + 1L
      results[[idx]] <- score_subset_prediction(
        scenario.grid[s, ], method, fit.out,
        dat$x, dat$y, dat$truth,
        weight.rule = spec$weight.rule,
        conductance.source = spec$conductance.source,
        trend.order = spec$order
      )
      pidx <- pidx + 1L
      predictions[[pidx]] <- prediction_frame(
        scenario.grid[s, ], dat$x, dat$y, dat$truth,
        method, fit.out
      )
    }
  }

  list(
    results = do.call(rbind, results[seq_len(idx)]),
    predictions = do.call(rbind, predictions[seq_len(pidx)]),
    scenario.grid = scenario.grid,
    lambda.grid.order0 = lambda.grid.order0,
    lambda.grid.higher = lambda.grid.higher,
    graph.fit.grid = graph.fit.grid
  )
}

out <- run_graph_trend_filtering_subset(progress = TRUE)

csv.file <- file.path(out.dir, "graph_trend_filtering_subset_results.csv")
pred.file <- file.path(out.dir, "graph_trend_filtering_subset_predictions.csv")
rds.file <- file.path(out.dir, "graph_trend_filtering_subset_results.rds")
write.csv(out$results, csv.file, row.names = FALSE)
write.csv(out$predictions, pred.file, row.names = FALSE)
saveRDS(out, rds.file)

message(sprintf("Wrote %s", csv.file))
message(sprintf("Wrote %s", pred.file))
message(sprintf("Wrote %s", rds.file))
message(sprintf("Rows: %d; predictions: %d; datasets: %d; errors: %d",
                nrow(out$results), nrow(out$predictions),
                length(unique(out$results$dataset.id)),
                sum(out$results$status != "ok")))
