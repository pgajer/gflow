#!/usr/bin/env Rscript

## Focused PHATE-graph embedding benchmark.
##
## This development report tests whether the PHATE graph support, embedded by
## weighted GRIP or by PHATE diffusion-potential MDS variants, recovers the
## original 3D geometry of the paraboloid and indefinite quadratic graph
## examples.

if (requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(".", quiet = TRUE)
} else {
  library(gflow)
}

required <- c("grip", "smacof", "ggplot2")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0L) {
  stop("Missing required package(s): ", paste(missing, collapse = ", "))
}

out.dir <- file.path("dev", "phate-graph-mds-comparison", "report")
fig.dir <- file.path(out.dir, "figures")
dir.create(fig.dir, recursive = TRUE, showWarnings = FALSE)

csv.path <- file.path(out.dir, "phate_graph_mds_geometry_metrics.csv")
html.path <- file.path(out.dir, "phate_graph_mds_geometry_report.html")

datasets <- c("paraboloid", "saddle")
n.values <- c(50L, 100L, 200L, 300L)
k.values <- 3L:10L
decay.values <- c(10L, 20L, 40L, 80L)
embedding.methods <- c("weighted_grip", "classic", "metric", "nonmetric")

html.escape <- function(x) {
  x <- as.character(x)
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub("\"", "&quot;", x, fixed = TRUE)
  x
}

make.dataset <- function(dataset, n) {
  seed <- switch(dataset,
                 paraboloid = 21000L,
                 saddle = 31000L,
                 41000L) + as.integer(n)
  set.seed(seed)
  theta <- stats::runif(n, 0, 2 * pi)
  radius <- sqrt(stats::runif(n, 0, 1))
  x <- radius * cos(theta)
  y <- radius * sin(theta)
  z <- if (identical(dataset, "paraboloid")) x^2 + y^2 else x^2 - y^2
  list(X = cbind(x, y, z), seed = seed)
}

dataset.label <- function(dataset) {
  switch(dataset,
         paraboloid = "2D paraboloid graph: z = x^2 + y^2",
         saddle = "2D indefinite quadratic graph: z = x^2 - y^2",
         dataset)
}

edges.from.kernel <- function(K) {
  K <- as.matrix(K)
  idx <- which(upper.tri(K) & K > 0, arr.ind = TRUE)
  if (nrow(idx) == 0L) {
    return(matrix(integer(0), ncol = 2L))
  }
  idx <- idx[order(idx[, 1], idx[, 2]), , drop = FALSE]
  cbind(idx[, 1] - 1L, idx[, 2] - 1L)
}

edge.lengths <- function(X, edges) {
  if (is.null(edges) || nrow(edges) == 0L) {
    return(numeric(0))
  }
  sqrt(rowSums((X[edges[, 1] + 1L, , drop = FALSE] -
                  X[edges[, 2] + 1L, , drop = FALSE])^2))
}

component.stats <- function(n, edges) {
  parent <- seq_len(n)
  size <- rep.int(1L, n)
  find <- function(x) {
    while (parent[x] != x) {
      parent[x] <<- parent[parent[x]]
      x <- parent[x]
    }
    x
  }
  unite <- function(a, b) {
    ra <- find(a)
    rb <- find(b)
    if (ra == rb) return()
    if (size[ra] < size[rb]) {
      tmp <- ra
      ra <- rb
      rb <- tmp
    }
    parent[rb] <<- ra
    size[ra] <<- size[ra] + size[rb]
  }
  if (!is.null(edges) && nrow(edges) > 0L) {
    for (i in seq_len(nrow(edges))) {
      unite(edges[i, 1] + 1L, edges[i, 2] + 1L)
    }
  }
  roots <- vapply(seq_len(n), find, integer(1))
  tab <- table(roots)
  list(n_components = length(tab),
       largest_component_fraction = as.numeric(max(tab)) / n)
}

pairwise <- function(X) {
  as.numeric(stats::dist(X))
}

embedding.metrics <- function(ref, obs) {
  ref <- as.matrix(ref)
  obs <- as.matrix(obs)
  d.ref <- pairwise(ref)
  d.obs <- pairwise(obs)
  dist.denom <- sum(d.obs^2)
  dist.scale <- if (is.finite(dist.denom) && dist.denom > 0) {
    sum(d.ref * d.obs) / dist.denom
  } else {
    NA_real_
  }
  rel.dist.error <- if (is.finite(dist.scale)) {
    sqrt(sum((d.ref - dist.scale * d.obs)^2) / sum(d.ref^2))
  } else {
    NA_real_
  }

  A <- scale(ref, center = TRUE, scale = FALSE)
  B <- scale(obs, center = TRUE, scale = FALSE)
  sv <- svd(t(B) %*% A)
  rotation <- sv$u %*% t(sv$v)
  B.rot <- B %*% rotation
  proc.denom <- sum(B.rot^2)
  proc.scale <- if (is.finite(proc.denom) && proc.denom > 0) {
    sum(A * B.rot) / proc.denom
  } else {
    NA_real_
  }
  procrustes.error <- if (is.finite(proc.scale)) {
    sqrt(sum((A - proc.scale * B.rot)^2) / sum(A^2))
  } else {
    NA_real_
  }

  data.frame(
    rel_dist_error = rel.dist.error,
    procrustes_error = procrustes.error,
    pearson = suppressWarnings(stats::cor(d.ref, d.obs, method = "pearson")),
    spearman = suppressWarnings(stats::cor(d.ref, d.obs, method = "spearman"))
  )
}

weighted.grip.embedding <- function(X, edges, seed) {
  grip::grip.layout.weighted(
    edges = edges + 1L,
    edge_weights = edge.lengths(X, edges),
    n = nrow(X),
    dim = 3L,
    rounds = 24L,
    final_rounds = 32L,
    num_init = min(8L, nrow(X)),
    num_nbrs = min(16L, nrow(X) - 1L),
    coarse_repulsion_factor = 0.8,
    coarse_repulsion_sample = min(12L, nrow(X)),
    coarse_repulsion_exact_below = 48L,
    seed = seed,
    disconnected = "components"
  )
}

run.with.warnings <- function(expr) {
  warnings <- character()
  value <- withCallingHandlers(
    expr,
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(value = value, warnings = unique(warnings))
}

safe.embedding <- function(expr) {
  tryCatch(
    {
      out <- run.with.warnings(expr)
      list(ok = TRUE, value = out$value, warnings = out$warnings, error = NA_character_)
    },
    error = function(e) {
      list(ok = FALSE, value = NULL, warnings = character(), error = conditionMessage(e))
    }
  )
}

append.rows <- function(rows, path) {
  df <- do.call(rbind, rows)
  utils::write.table(df, path, sep = ",", row.names = FALSE,
                     col.names = !file.exists(path), append = file.exists(path),
                     qmethod = "double")
}

if (file.exists(csv.path)) {
  file.remove(csv.path)
}

message("Writing incremental metrics to ", normalizePath(csv.path, mustWork = FALSE))
case.counter <- 0L

for (dataset in datasets) {
  for (n in n.values) {
    ds <- make.dataset(dataset, n)
    X <- ds$X
    for (k in k.values) {
      for (decay in decay.values) {
        case.counter <- case.counter + 1L
        case.id <- paste(dataset, paste0("n", n), paste0("k", k),
                         paste0("decay", decay), sep = "_")
        message(sprintf("[%03d] %s", case.counter, case.id))

        core.out <- safe.embedding(phate.core(
          X = X,
          k = k,
          alpha.decay = decay,
          kernel.mode = "phate",
          vne.method = "phate",
          potential.mode = "phate",
          gamma = 1,
          compute.D.pot = TRUE,
          verbose = FALSE
        ))

        if (!isTRUE(core.out$ok)) {
          rows <- lapply(embedding.methods, function(method) {
            data.frame(
              dataset = dataset, n = n, k = k, decay = decay,
              case_id = case.id, method = method, status = "core_error",
              rel_dist_error = NA_real_, procrustes_error = NA_real_,
              pearson = NA_real_, spearman = NA_real_,
              t = NA_integer_, edge_count = NA_integer_,
              n_components = NA_integer_,
              largest_component_fraction = NA_real_,
              warnings = NA_character_, error = core.out$error,
              stringsAsFactors = FALSE
            )
          })
          append.rows(rows, csv.path)
          next
        }

        core <- core.out$value
        edges <- edges.from.kernel(core$K)
        comps <- component.stats(n, edges)

        rows <- list()

        grip.out <- safe.embedding(weighted.grip.embedding(
          X = X,
          edges = edges,
          seed = ds$seed + 7000L + 100L * k + decay
        ))
        if (isTRUE(grip.out$ok)) {
          m <- embedding.metrics(X, grip.out$value)
          rows[[length(rows) + 1L]] <- data.frame(
            dataset = dataset, n = n, k = k, decay = decay,
            case_id = case.id, method = "weighted_grip", status = "ok",
            m,
            t = core$t, edge_count = nrow(edges),
            n_components = comps$n_components,
            largest_component_fraction = comps$largest_component_fraction,
            warnings = paste(c(core.out$warnings, grip.out$warnings), collapse = " | "),
            error = NA_character_,
            stringsAsFactors = FALSE
          )
        } else {
          rows[[length(rows) + 1L]] <- data.frame(
            dataset = dataset, n = n, k = k, decay = decay,
            case_id = case.id, method = "weighted_grip", status = "embed_error",
            rel_dist_error = NA_real_, procrustes_error = NA_real_,
            pearson = NA_real_, spearman = NA_real_,
            t = core$t, edge_count = nrow(edges),
            n_components = comps$n_components,
            largest_component_fraction = comps$largest_component_fraction,
            warnings = paste(core.out$warnings, collapse = " | "),
            error = grip.out$error,
            stringsAsFactors = FALSE
          )
        }

        for (method in c("classic", "metric", "nonmetric")) {
          emb.out <- safe.embedding(phate.embed(
            core = core,
            ndim = 3L,
            method = method,
            maxit = 200L,
            tol = 1e-3,
            seed = ds$seed + 9000L + 100L * k + decay +
              match(method, c("classic", "metric", "nonmetric")),
            verbose = FALSE
          ))

          if (isTRUE(emb.out$ok)) {
            emb <- emb.out$value
            m <- embedding.metrics(X, emb$embedding)
            rows[[length(rows) + 1L]] <- data.frame(
              dataset = dataset, n = n, k = k, decay = decay,
              case_id = case.id, method = method, status = "ok",
              m,
              t = core$t, edge_count = nrow(edges),
              n_components = comps$n_components,
              largest_component_fraction = comps$largest_component_fraction,
              warnings = paste(c(core.out$warnings, emb.out$warnings), collapse = " | "),
              error = NA_character_,
              stringsAsFactors = FALSE
            )
          } else {
            rows[[length(rows) + 1L]] <- data.frame(
              dataset = dataset, n = n, k = k, decay = decay,
              case_id = case.id, method = method, status = "embed_error",
              rel_dist_error = NA_real_, procrustes_error = NA_real_,
              pearson = NA_real_, spearman = NA_real_,
              t = core$t, edge_count = nrow(edges),
              n_components = comps$n_components,
              largest_component_fraction = comps$largest_component_fraction,
              warnings = paste(core.out$warnings, collapse = " | "),
              error = emb.out$error,
              stringsAsFactors = FALSE
            )
          }
        }

        append.rows(rows, csv.path)
      }
    }
  }
}

metrics <- utils::read.csv(csv.path, stringsAsFactors = FALSE)
ok <- metrics[metrics$status == "ok", , drop = FALSE]

method.order <- c("weighted_grip", "classic", "metric", "nonmetric")
ok$method <- factor(ok$method, levels = method.order)

case.best.rel <- ok[order(ok$case_id, ok$rel_dist_error), ]
case.best.rel <- case.best.rel[!duplicated(case.best.rel$case_id), ]
case.best.proc <- ok[order(ok$case_id, ok$procrustes_error), ]
case.best.proc <- case.best.proc[!duplicated(case.best.proc$case_id), ]

summ <- aggregate(
  cbind(rel_dist_error, procrustes_error, pearson, spearman) ~ dataset + method,
  ok,
  function(x) c(mean = mean(x), median = stats::median(x), max = max(x))
)
summ <- do.call(data.frame, summ)
names(summ) <- gsub("\\.", "_", names(summ))

win.rel <- as.data.frame(table(case.best.rel$dataset, case.best.rel$method),
                         stringsAsFactors = FALSE)
names(win.rel) <- c("dataset", "method", "rel_dist_wins")
win.proc <- as.data.frame(table(case.best.proc$dataset, case.best.proc$method),
                          stringsAsFactors = FALSE)
names(win.proc) <- c("dataset", "method", "procrustes_wins")
win <- merge(win.rel, win.proc, by = c("dataset", "method"), all = TRUE)
win$rel_dist_wins[is.na(win$rel_dist_wins)] <- 0L
win$procrustes_wins[is.na(win$procrustes_wins)] <- 0L

win.by.n <- as.data.frame(table(case.best.rel$dataset, case.best.rel$n,
                                case.best.rel$method),
                          stringsAsFactors = FALSE)
names(win.by.n) <- c("dataset", "n", "method", "rel_dist_wins")
win.by.n <- win.by.n[win.by.n$rel_dist_wins > 0L, , drop = FALSE]
win.by.n <- win.by.n[order(win.by.n$dataset, as.integer(as.character(win.by.n$n)),
                           -win.by.n$rel_dist_wins), ]

win.by.decay <- as.data.frame(table(case.best.rel$dataset, case.best.rel$decay,
                                    case.best.rel$method),
                              stringsAsFactors = FALSE)
names(win.by.decay) <- c("dataset", "decay", "method", "rel_dist_wins")
win.by.decay <- win.by.decay[win.by.decay$rel_dist_wins > 0L, , drop = FALSE]
win.by.decay <- win.by.decay[order(win.by.decay$dataset,
                                   as.integer(as.character(win.by.decay$decay)),
                                   -win.by.decay$rel_dist_wins), ]

by.decay <- aggregate(
  cbind(rel_dist_error, procrustes_error) ~ dataset + decay + method,
  ok,
  mean
)

utils::write.csv(summ, file.path(out.dir, "phate_graph_mds_geometry_summary.csv"),
                 row.names = FALSE)
utils::write.csv(win, file.path(out.dir, "phate_graph_mds_geometry_wins.csv"),
                 row.names = FALSE)
utils::write.csv(win.by.n, file.path(out.dir, "phate_graph_mds_geometry_wins_by_n.csv"),
                 row.names = FALSE)
utils::write.csv(win.by.decay, file.path(out.dir, "phate_graph_mds_geometry_wins_by_decay.csv"),
                 row.names = FALSE)
utils::write.csv(by.decay, file.path(out.dir, "phate_graph_mds_geometry_by_decay.csv"),
                 row.names = FALSE)

plot.method.box <- ggplot2::ggplot(ok, ggplot2::aes(x = method, y = rel_dist_error, fill = method)) +
  ggplot2::geom_boxplot(outlier.alpha = 0.35) +
  ggplot2::facet_wrap(~ dataset, scales = "free_y") +
  ggplot2::labs(x = NULL, y = "relative distance error",
                title = "Geometry recovery error by embedding method") +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(legend.position = "none",
                 axis.text.x = ggplot2::element_text(angle = 25, hjust = 1))
box.path <- file.path(fig.dir, "rel_dist_error_boxplot.png")
ggplot2::ggsave(box.path, plot.method.box, width = 8.5, height = 4.8, dpi = 160)

plot.decay <- ggplot2::ggplot(by.decay, ggplot2::aes(x = decay, y = rel_dist_error,
                                                     color = method, group = method)) +
  ggplot2::geom_line(linewidth = 0.8) +
  ggplot2::geom_point(size = 1.8) +
  ggplot2::facet_wrap(~ dataset, scales = "free_y") +
  ggplot2::scale_x_continuous(breaks = decay.values) +
  ggplot2::labs(x = "PHATE decay", y = "mean relative distance error",
                title = "Mean geometry error across n and k by decay") +
  ggplot2::theme_minimal(base_size = 12)
decay.path <- file.path(fig.dir, "rel_dist_error_by_decay.png")
ggplot2::ggsave(decay.path, plot.decay, width = 8.5, height = 4.8, dpi = 160)

format.num <- function(x, digits = 4L) {
  ifelse(is.na(x), "", formatC(x, digits = digits, format = "f"))
}

html.table <- function(df, digits = 4L, max.rows = Inf) {
  if (nrow(df) > max.rows) {
    df <- utils::head(df, max.rows)
  }
  df.out <- df
  for (nm in names(df.out)) {
    if (is.numeric(df.out[[nm]])) {
      df.out[[nm]] <- format.num(df.out[[nm]], digits = digits)
    }
  }
  header <- paste0("<tr>", paste(sprintf("<th>%s</th>", html.escape(names(df.out))),
                                 collapse = ""), "</tr>")
  rows <- apply(df.out, 1L, function(row) {
    paste0("<tr>", paste(sprintf("<td>%s</td>", html.escape(row)), collapse = ""), "</tr>")
  })
  paste0("<table><thead>", header, "</thead><tbody>",
         paste(rows, collapse = "\n"), "</tbody></table>")
}

top.win <- win[order(win$dataset, -win$rel_dist_wins, -win$procrustes_wins), ]
summary.display <- summ[, c(
  "dataset", "method",
  "rel_dist_error_mean", "rel_dist_error_median", "rel_dist_error_max",
  "procrustes_error_mean", "procrustes_error_median", "procrustes_error_max",
  "pearson_mean", "spearman_mean"
)]
summary.display <- summary.display[order(summary.display$dataset,
                                         summary.display$rel_dist_error_mean), ]

worst <- ok[order(ok$rel_dist_error, decreasing = TRUE), c(
  "dataset", "n", "k", "decay", "method", "rel_dist_error",
  "procrustes_error", "pearson", "spearman", "edge_count", "t"
)]

exceptions <- case.best.rel[case.best.rel$method != "weighted_grip", c(
  "dataset", "n", "k", "decay", "method", "rel_dist_error",
  "procrustes_error", "pearson", "spearman", "edge_count", "t"
)]

best.by.case <- merge(
  case.best.rel[, c("case_id", "dataset", "n", "k", "decay",
                   "method", "rel_dist_error", "procrustes_error")],
  case.best.proc[, c("case_id", "method", "procrustes_error")],
  by = "case_id",
  suffixes = c("_best_rel", "_best_proc")
)
best.by.case <- best.by.case[order(best.by.case$dataset, best.by.case$n,
                                   best.by.case$k, best.by.case$decay), ]

total.cases <- length(unique(ok$case_id))
grip.rel.wins <- sum(case.best.rel$method == "weighted_grip")
grip.proc.wins <- sum(case.best.proc$method == "weighted_grip")

css <- paste(
  "body{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif;margin:34px;line-height:1.48;color:#222}",
  "h1,h2{color:#16324f}h3{color:#243b53;margin-top:24px}",
  ".note{background:#f7f9fb;border-left:4px solid #3B6EA8;padding:12px 14px;margin:16px 0}",
  "table{border-collapse:collapse;width:100%;font-size:13px;margin:14px 0}",
  "th,td{border:1px solid #ddd;padding:5px 7px;text-align:left}th{background:#f2f5f8}",
  "code{background:#f4f4f4;padding:1px 3px}.fig{max-width:980px;width:100%;border:1px solid #ddd}",
  ".math-block{overflow-x:auto;margin:12px 0}",
  sep = "\n"
)

html <- paste0(
  "<!doctype html><html><head><meta charset=\"utf-8\">",
  "<title>PHATE Graph MDS Geometry Benchmark</title>",
  "<script>window.MathJax={tex:{inlineMath:[[\"$\",\"$\"],[\"\\\\(\",\"\\\\)\"]],displayMath:[[\"\\\\[\",\"\\\\]\"],[\"$$\",\"$$\"]]}};</script>",
  "<script defer src=\"https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-chtml.js\"></script>",
  "<style>", css, "</style></head><body>",
  "<h1>PHATE Graph MDS Geometry Benchmark</h1>",
  "<p><strong>Generated:</strong> ", html.escape(format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")), "</p>",
  "<div class=\"note\">This report tests the observation that, on the 2D paraboloid and indefinite quadratic graph datasets, the PHATE graph support embedded by weighted GRIP appears closest to the original sampled surface geometry. The comparison uses the same PHATE graph in every row and varies only the embedding method.</div>",
  "<h2>Question</h2>",
  "<p>For each dataset, sample size \\(n\\), PHATE bandwidth parameter \\(k\\), and decay \\(a\\), PHATE builds an adaptive affinity graph. The tested methods are:</p>",
  "<ul>",
  "<li><code>weighted_grip</code>: weighted GRIP layout of the PHATE graph support using original Euclidean edge lengths \\(\\ell_{ij}=\\lVert x_i-x_j\\rVert_2\\).</li>",
  "<li><code>classic</code>: PHATE diffusion-potential distances embedded by classical MDS, <code>stats::cmdscale(add = TRUE)</code>.</li>",
  "<li><code>metric</code>: PHATE diffusion-potential distances embedded by metric SMACOF, <code>smacof::smacofSym(type = \"ratio\")</code>, initialized by classical MDS.</li>",
  "<li><code>nonmetric</code>: PHATE diffusion-potential distances embedded by ordinal/nonmetric SMACOF, <code>smacof::smacofSym(type = \"ordinal\")</code>, initialized by classical MDS.</li>",
  "</ul>",
  "<h2>Metrics</h2>",
  "<p>The reference geometry is the original observed 3D surface coordinates \\(X\\). For an embedding \\(Y\\), relative distance error compares all pairwise distances after optimal scalar rescaling:</p>",
  "<div class=\"math-block\">\\[\\operatorname{rel\\_dist\\_error}(X,Y)=\\left(\\frac{\\sum_{i&lt;j}(d^X_{ij}-s_Dd^Y_{ij})^2}{\\sum_{i&lt;j}(d^X_{ij})^2}\\right)^{1/2},\\qquad s_D=\\frac{\\sum_{i&lt;j}d^X_{ij}d^Y_{ij}}{\\sum_{i&lt;j}(d^Y_{ij})^2}.\\]</div>",
  "<p>Procrustes error compares coordinates after centering, optimal rotation/reflection, and optimal scalar rescaling:</p>",
  "<div class=\"math-block\">\\[\\operatorname{procrustes\\_error}(X,Y)=\\frac{\\lVert X_c-s_PY_cR^*\\rVert_F}{\\lVert X_c\\rVert_F}.\\]</div>",
  "<p>Lower values are better. Errors below about 0.05 are visually very small, 0.05-0.10 are small, and values above 0.10 usually correspond to visible geometric distortion for these smooth 3D surface examples.</p>",
  "<h2>Headline Result</h2>",
  "<p>The benchmark ran ", total.cases, " PHATE graph cases: two datasets, four sample sizes, eight values of \\(k\\), and four decay values. Weighted GRIP won by relative distance error in <strong>", grip.rel.wins, "</strong> cases and by Procrustes error in <strong>", grip.proc.wins, "</strong> cases.</p>",
  "<p>The observation is therefore strongly supported, but not literally universal. The exceptions are concentrated in sparse low-\\(k\\) graph settings, especially paraboloid cases with \\(k=3\\) or \\(k=4\\), and a small number of saddle cases at \\(n=50\\).</p>",
  "<h2>Winner Counts</h2>",
  html.table(top.win, digits = 0L),
  "<h2>Winner Counts By Sample Size</h2>",
  "<p>These counts use relative distance error as the criterion. Each row counts the best method among the four tested embeddings for fixed dataset, \\(n\\), \\(k\\), and decay.</p>",
  html.table(win.by.n, digits = 0L),
  "<h2>Winner Counts By PHATE Decay</h2>",
  "<p>These counts check the suspicion that the weighted GRIP advantage holds across decay values. It mostly does, with a small number of MDS wins at every paraboloid decay and very few saddle exceptions.</p>",
  html.table(win.by.decay, digits = 0L),
  "<h2>Summary By Dataset And Method</h2>",
  html.table(summary.display, digits = 4L),
  "<h2>Relative Distance Error Distributions</h2>",
  "<p><img class=\"fig\" src=\"figures/", basename(box.path), "\" alt=\"Relative distance error boxplot\"></p>",
  "<h2>Decay Sensitivity</h2>",
  "<p><img class=\"fig\" src=\"figures/", basename(decay.path), "\" alt=\"Relative distance error by decay\"></p>",
  "<h2>Worst Individual Method Fits</h2>",
  html.table(utils::head(worst, 24L), digits = 4L),
  "<h2>Cases Where Weighted GRIP Did Not Win</h2>",
  "<p>This table lists every case where one of the PHATE diffusion-potential MDS methods had lower relative distance error than weighted GRIP on the same PHATE graph.</p>",
  html.table(exceptions, digits = 4L),
  "<h2>Per-Case Best Methods</h2>",
  "<p>The full per-case best-method table is written to the CSV files. The first rows are shown here for quick inspection.</p>",
  html.table(utils::head(best.by.case, 40L), digits = 4L),
  "<h2>Generated Files</h2>",
  "<ul>",
  "<li><code>", html.escape(normalizePath(csv.path, mustWork = FALSE)), "</code></li>",
  "<li><code>", html.escape(normalizePath(file.path(out.dir, "phate_graph_mds_geometry_summary.csv"), mustWork = FALSE)), "</code></li>",
  "<li><code>", html.escape(normalizePath(file.path(out.dir, "phate_graph_mds_geometry_wins.csv"), mustWork = FALSE)), "</code></li>",
  "<li><code>", html.escape(normalizePath(file.path(out.dir, "phate_graph_mds_geometry_wins_by_n.csv"), mustWork = FALSE)), "</code></li>",
  "<li><code>", html.escape(normalizePath(file.path(out.dir, "phate_graph_mds_geometry_wins_by_decay.csv"), mustWork = FALSE)), "</code></li>",
  "<li><code>", html.escape(normalizePath(file.path(out.dir, "phate_graph_mds_geometry_by_decay.csv"), mustWork = FALSE)), "</code></li>",
  "</ul>",
  "</body></html>"
)

writeLines(html, html.path)
message("Wrote ", normalizePath(html.path, mustWork = FALSE))
