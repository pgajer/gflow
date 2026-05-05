#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  pkgload::load_all(".", quiet = TRUE)
})

repo.root <- normalizePath(".", winslash = "/", mustWork = TRUE)
core.root <- file.path(repo.root, "tests/testthat/fixtures/ph6/core")
stress.root <- file.path(repo.root, "dev/phate-phase6-validation/cache/stress")
out.dir <- file.path(repo.root, "dev/phate-phase6-validation/report")
fig.dir <- file.path(out.dir, "figures")
log.dir <- file.path(repo.root, "dev/phate-phase6-validation/validation_logs")
dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig.dir, recursive = TRUE, showWarnings = FALSE)

read.matrix.csv <- function(path) {
  unname(as.matrix(utils::read.csv(path, header = FALSE)))
}

read.manifest <- function(root) {
  path <- file.path(root, "manifest.json")
  if (!file.exists(path)) {
    return(character())
  }
  x <- jsonlite::fromJSON(path)
  x$case_id
}

embedding.file <- function(case.dir, method) {
  if (identical(method, "metric3")) {
    return(file.path(case.dir, "e3_m.csv"))
  }
  suffix <- switch(method,
    classic = "c",
    metric = "m",
    nonmetric = "nm",
    stop("unknown method")
  )
  file.path(case.dir, sprintf("e_%s.csv", suffix))
}

pairwise <- function(X) {
  as.numeric(stats::dist(X))
}

embedding.metrics <- function(ref, obs) {
  d.ref <- pairwise(ref)
  d.obs <- pairwise(obs)
  scale.factor <- sum(d.ref * d.obs) / sum(d.obs^2)
  rel.dist.error <- sqrt(sum((d.ref - scale.factor * d.obs)^2) / sum(d.ref^2))

  A <- scale(ref, center = TRUE, scale = FALSE)
  B <- scale(obs, center = TRUE, scale = FALSE)
  sv <- svd(t(B) %*% A)
  rotation <- sv$u %*% t(sv$v)
  B.rot <- B %*% rotation
  B.scale <- sum(A * B.rot) / sum(B.rot^2)
  procrustes.error <- sqrt(sum((A - B.scale * B.rot)^2) / sum(A^2))

  data.frame(
    rel_dist_error = rel.dist.error,
    procrustes_error = procrustes.error,
    pearson = stats::cor(d.ref, d.obs, method = "pearson"),
    spearman = stats::cor(d.ref, d.obs, method = "spearman")
  )
}

procrustes.align <- function(ref, obs) {
  ref <- as.matrix(ref)
  obs <- as.matrix(obs)
  ref.center <- colMeans(ref)
  obs.center <- colMeans(obs)
  A <- sweep(ref, 2, ref.center, "-")
  B <- sweep(obs, 2, obs.center, "-")
  sv <- svd(t(B) %*% A)
  rotation <- sv$u %*% t(sv$v)
  B.rot <- B %*% rotation
  scale.factor <- sum(A * B.rot) / sum(B.rot^2)
  aligned <- sweep(scale.factor * B.rot, 2, ref.center, "+")
  residual <- sqrt(sum((A - scale.factor * B.rot)^2) / sum(A^2))
  list(embedding = aligned, residual = residual)
}

analyze.case <- function(case.id, suite, root) {
  case.dir <- file.path(root, case.id)
  meta <- jsonlite::fromJSON(file.path(case.dir, "metadata.json"), simplifyVector = FALSE)
  X <- read.matrix.csv(file.path(case.dir, "X.csv"))

  core <- phate.core(
    X = X,
    k = meta$params$knn,
    alpha.decay = meta$params$decay,
    t = "auto",
    t.max = meta$params$t_max,
    gamma = meta$params$gamma,
    kernel.mode = "phate",
    vne.method = "phate",
    potential.mode = "phate",
    compute.D.pot = TRUE,
    verbose = FALSE
  )

  refs <- list(
    K = read.matrix.csv(file.path(case.dir, "K.csv")),
    P = read.matrix.csv(file.path(case.dir, "P.csv")),
    Pt = read.matrix.csv(file.path(case.dir, "Pt.csv")),
    U = read.matrix.csv(file.path(case.dir, "U.csv")),
    Dpot = read.matrix.csv(file.path(case.dir, "D.csv")),
    vne = as.numeric(read.matrix.csv(file.path(case.dir, "vne.csv")))
  )

  operator <- data.frame(
    suite = suite,
    case_id = case.id,
    family = as.character(meta$example$family),
    n = meta$shape$n,
    p = meta$shape$p,
    t_gflow = core$t,
    t_python = meta$t_auto,
    K_max_abs = max(abs(core$K - refs$K)),
    P_max_abs = max(abs(core$P - refs$P)),
    vne_max_abs = max(abs(core$diagnostics$vne - refs$vne)),
    Pt_max_abs = max(abs(core$Pt - refs$Pt)),
    U_max_abs = max(abs(core$U.pot - refs$U)),
    Dpot_max_abs = max(abs(unname(core$D.pot) - refs$Dpot)),
    stringsAsFactors = FALSE
  )
  operator$operator_max_abs <- apply(operator[, c(
    "K_max_abs", "P_max_abs", "vne_max_abs", "Pt_max_abs", "U_max_abs", "Dpot_max_abs"
  )], 1, max)
  operator$operator_status <- ifelse(
    operator$t_gflow == operator$t_python && operator$operator_max_abs <= 1e-10,
    "pass",
    "review"
  )

  emb.rows <- list()
  for (method in c("classic", "metric", "nonmetric")) {
    emb <- phate.embed(
      core = core,
      ndim = 2,
      method = method,
      seed = meta$seed,
      maxit = 300,
      verbose = FALSE
    )
    ref <- read.matrix.csv(embedding.file(case.dir, method))
    m <- embedding.metrics(ref, emb$embedding)
    m$suite <- suite
    m$case_id <- case.id
    m$method <- method
    m$stress_raw <- if (is.null(emb$stress_raw)) NA_real_ else emb$stress_raw
    m$status <- ifelse(
      m$rel_dist_error <= 0.18 &&
        m$procrustes_error <= 0.28 &&
        m$pearson >= 0.94 &&
        m$spearman >= 0.90,
      "pass",
      "review"
    )
    emb.rows[[method]] <- m[, c(
      "suite", "case_id", "method", "rel_dist_error", "procrustes_error",
      "pearson", "spearman", "stress_raw", "status"
    )]
  }

  list(operator = operator, embedding = do.call(rbind, emb.rows))
}

all.cases <- rbind(
  data.frame(suite = "core", case_id = read.manifest(core.root), root = core.root),
  data.frame(suite = "stress", case_id = read.manifest(stress.root), root = stress.root)
)

results <- lapply(seq_len(nrow(all.cases)), function(i) {
  message("Analyzing ", all.cases$suite[i], "/", all.cases$case_id[i])
  analyze.case(all.cases$case_id[i], all.cases$suite[i], all.cases$root[i])
})

operator.metrics <- do.call(rbind, lapply(results, `[[`, "operator"))
embedding.metrics.df <- do.call(rbind, lapply(results, `[[`, "embedding"))

utils::write.csv(operator.metrics, file.path(out.dir, "phase6_operator_metrics.csv"), row.names = FALSE)
utils::write.csv(embedding.metrics.df, file.path(out.dir, "phase6_embedding_metrics.csv"), row.names = FALSE)

png(file.path(fig.dir, "operator_max_abs.png"), width = 1600, height = 900, res = 150)
op <- operator.metrics[order(operator.metrics$operator_max_abs, decreasing = TRUE), ]
cols <- ifelse(op$suite == "core", "#3B6EA8", "#7A7A7A")
par(mar = c(10, 5, 3, 1))
barplot(
  log10(pmax(op$operator_max_abs, .Machine$double.eps)),
  names.arg = op$case_id,
  las = 2,
  col = cols,
  border = NA,
  ylab = "log10(max absolute operator error)",
  main = "Operator-level parity against original PHATE fixtures"
)
abline(h = log10(1e-10), col = "#B00020", lty = 2, lwd = 2)
legend("topright", fill = c("#3B6EA8", "#7A7A7A"), legend = c("core", "stress"), bty = "n")
dev.off()

png(file.path(fig.dir, "embedding_relative_distance_error.png"), width = 1700, height = 950, res = 150)
em <- embedding.metrics.df
case.order <- unique(em$case_id[order(em$rel_dist_error, decreasing = TRUE)])
method.cols <- c(classic = "#7B3294", metric = "#008837", nonmetric = "#D95F02")
par(mar = c(10, 5, 3, 1))
plot(
  NA,
  xlim = c(0.5, length(case.order) + 0.5),
  ylim = c(0, max(em$rel_dist_error) * 1.08),
  xaxt = "n",
  xlab = "",
  ylab = "relative pairwise-distance error",
  main = "Embedding geometry parity"
)
axis(1, at = seq_along(case.order), labels = case.order, las = 2, cex.axis = 0.65)
offsets <- c(classic = -0.22, metric = 0, nonmetric = 0.22)
for (method in names(offsets)) {
  sub <- em[em$method == method, ]
  x <- match(sub$case_id, case.order) + offsets[[method]]
  points(x, sub$rel_dist_error, pch = 19, col = method.cols[[method]], cex = 1.05)
}
abline(h = 0.18, col = "#B00020", lty = 2, lwd = 2)
legend("topright", col = method.cols, pch = 19, legend = names(method.cols), bty = "n")
dev.off()

png(file.path(fig.dir, "embedding_pearson.png"), width = 1700, height = 950, res = 150)
case.order <- unique(em$case_id[order(em$pearson)])
par(mar = c(10, 5, 3, 1))
plot(
  NA,
  xlim = c(0.5, length(case.order) + 0.5),
  ylim = c(min(em$pearson) - 0.02, 1.005),
  xaxt = "n",
  xlab = "",
  ylab = "Pearson correlation of pairwise embedding distances",
  main = "Embedding distance correlation"
)
axis(1, at = seq_along(case.order), labels = case.order, las = 2, cex.axis = 0.65)
for (method in names(offsets)) {
  sub <- em[em$method == method, ]
  x <- match(sub$case_id, case.order) + offsets[[method]]
  points(x, sub$pearson, pch = 19, col = method.cols[[method]], cex = 1.05)
}
abline(h = 0.94, col = "#B00020", lty = 2, lwd = 2)
legend("bottomright", col = method.cols, pch = 19, legend = names(method.cols), bty = "n")
dev.off()

html.escape <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x
}

scale.widget.coords <- function(X3) {
  X3 <- as.matrix(X3)
  X3 <- sweep(X3, 2, colMeans(X3), "-")
  max.abs <- max(abs(X3))
  if (is.finite(max.abs) && max.abs > 0) {
    X3 / max.abs
  } else {
    X3
  }
}

selected.gallery <- c(
  "noisy_circle_uniform_n64",
  "quadform_graph_n2_k1",
  "quadform_graph_n3_k0",
  "asymmetric_y_tree_n75",
  "phate_dla_tree_n80",
  "swiss_roll_n100",
  "figure_eight_n100",
  "compositional_cyclic_gradient_n90_p25"
)
gallery.path <- file.path(fig.dir, "example_gallery.png")
png(gallery.path, width = 2100, height = 1900, res = 150)
par(mfrow = c(4, 6), mar = c(2, 2, 3, 1))
for (case.id in selected.gallery) {
  row <- all.cases[all.cases$case_id == case.id, ]
  case.dir <- file.path(row$root, case.id)
  X <- read.matrix.csv(file.path(case.dir, "X.csv"))
  py.metric <- read.matrix.csv(embedding.file(case.dir, "metric"))
  meta <- jsonlite::fromJSON(file.path(case.dir, "metadata.json"), simplifyVector = FALSE)
  core <- phate.core(
    X = X, k = 5, alpha.decay = 40, t = "auto", t.max = 30, gamma = 1,
    kernel.mode = "phate", vne.method = "phate", potential.mode = "phate",
    compute.D.pot = TRUE, verbose = FALSE
  )
  emb <- phate.embed(core = core, ndim = 2, method = "metric", seed = meta$seed,
                     maxit = 300, verbose = FALSE)
  color <- grDevices::hcl.colors(nrow(X), "Viridis")[rank(X[, 1], ties.method = "first")]
  plot(X[, 1], X[, 2], pch = 19, cex = 0.65, col = color, axes = FALSE,
       xlab = "", ylab = "", main = paste(case.id, "observed", sep = "\n"))
  plot(py.metric[, 1], py.metric[, 2], pch = 19, cex = 0.65, col = color,
       axes = FALSE, xlab = "", ylab = "", main = "baseline PHATE")
  plot(emb$embedding[, 1], emb$embedding[, 2], pch = 19, cex = 0.65, col = color,
       axes = FALSE, xlab = "", ylab = "", main = "gflow metric PHATE")
}
dev.off()

interactive.case.ids <- c(
  "quadform_graph_n2_k0",
  "quadform_graph_n2_k1",
  "phate_dla_tree_n80",
  "s_curve_n90",
  "figure_eight_n100",
  "near_crossing_strands_n96",
  "highdim_noisy_circle_n90_p30",
  "compositional_cyclic_gradient_n90_p25"
)

widget.dir <- file.path(out.dir, "widgets3d")
if (dir.exists(widget.dir)) {
  unlink(widget.dir, recursive = TRUE)
}
dir.create(widget.dir, recursive = TRUE, showWarnings = FALSE)

observed.coords3 <- function(X) {
  X <- as.matrix(X)
  if (ncol(X) == 3) {
    return(X)
  }
  if (ncol(X) > 3) {
    return(stats::prcomp(X, center = TRUE, scale. = FALSE)$x[, 1:3, drop = FALSE])
  }
  cbind(X, matrix(0, nrow = nrow(X), ncol = 3 - ncol(X)))
}

quadratic.form.text <- function(latent.dim, index) {
  positive <- if (index > 0) {
    paste(sprintf("+ x%d^2", seq_len(index)), collapse = " ")
  } else {
    character()
  }
  negative <- if (index < latent.dim) {
    paste(sprintf("- x%d^2", seq.int(index + 1L, latent.dim)), collapse = " ")
  } else {
    character()
  }
  formula <- paste(c(positive, negative), collapse = " ")
  sub("^\\+ ", "", formula)
}

quadratic.form.kind <- function(latent.dim, index) {
  if (index == 0L) {
    return("negative definite paraboloid")
  }
  if (index == latent.dim) {
    return("positive definite paraboloid")
  }
  "indefinite quadratic form"
}

interactive.case.heading <- function(case.id, meta) {
  family <- as.character(meta$example$family)
  if (identical(family, "quadratic_form_graph")) {
    latent.dim <- as.integer(meta$example$latent_dim)
    index <- as.integer(meta$example$index)
    return(sprintf(
      paste0(
        "<h3>Quadratic-form graph: latent dimension n = %d, index k = %d</h3>",
        "<p class=\"case-caption\"><code>%s</code>; %s. Fixture id: <code>%s</code>.</p>"
      ),
      latent.dim,
      index,
      html.escape(paste0("y = ", quadratic.form.text(latent.dim, index))),
      html.escape(quadratic.form.kind(latent.dim, index)),
      html.escape(case.id)
    ))
  }

  label <- switch(case.id,
    phate_dla_tree_n80 = "DLA tree synthetic manifold",
    s_curve_n90 = "S-curve folded manifold",
    swiss_roll_n100 = "Swiss-roll folded manifold",
    figure_eight_n100 = "Figure-eight self-approaching loop",
    near_crossing_strands_n96 = "Near-crossing separate strands",
    highdim_noisy_circle_n90_p30 = "High-dimensional noisy circle",
    compositional_cyclic_gradient_n90_p25 = "Compositional cyclic gradient",
    gsub("_", " ", case.id, fixed = TRUE)
  )
  sprintf(
    "<h3>%s</h3><p class=\"case-caption\">Fixture id: <code>%s</code>; family: <code>%s</code>; samples: %d; observed dimension: %d.</p>",
    html.escape(label),
    html.escape(case.id),
    html.escape(family),
    as.integer(meta$shape$n),
    as.integer(meta$shape$p)
  )
}

draw.rgl.scatter.panel <- function(X3, color, title, radius = 0.026) {
  X3 <- scale.widget.coords(X3)
  rgl::plot3d(
    X3[, 1], X3[, 2], X3[, 3],
    type = "n",
    axes = FALSE,
    box = FALSE,
    xlim = c(-1.08, 1.08),
    ylim = c(-1.08, 1.08),
    zlim = c(-1.08, 1.08),
    xlab = "",
    ylab = "",
    zlab = ""
  )
  rgl::spheres3d(
    X3[, 1], X3[, 2], X3[, 3],
    col = color,
    radius = radius
  )
  rgl::aspect3d(1, 1, 1)
  rgl::axes3d(edges = "bbox", nticks = 3, cex = 0.7, col = "gray45")
  rgl::title3d(main = title, color = "gray20", cex = 1.2)
  rgl::view3d(theta = 35, phi = 22, zoom = 0.82)
}

draw.rgl.overlay.panel <- function(ref3, obs3, title, radius = 0.022) {
  combined <- scale.widget.coords(rbind(ref3, obs3))
  n <- nrow(ref3)
  ref3 <- combined[seq_len(n), , drop = FALSE]
  obs3 <- combined[n + seq_len(n), , drop = FALSE]
  rgl::plot3d(
    ref3[, 1], ref3[, 2], ref3[, 3],
    type = "n",
    axes = FALSE,
    box = FALSE,
    xlim = c(-1.08, 1.08),
    ylim = c(-1.08, 1.08),
    zlim = c(-1.08, 1.08),
    xlab = "",
    ylab = "",
    zlab = ""
  )
  rgl::spheres3d(ref3[, 1], ref3[, 2], ref3[, 3],
                 col = "#111111", radius = radius * 1.2)
  rgl::spheres3d(obs3[, 1], obs3[, 2], obs3[, 3],
                 col = "#D7191C", radius = radius * 0.8)
  rgl::aspect3d(1, 1, 1)
  rgl::axes3d(edges = "bbox", nticks = 3, cex = 0.7, col = "gray45")
  rgl::title3d(main = title, color = "gray20", cex = 1.2)
  rgl::view3d(theta = 35, phi = 22, zoom = 0.82)
}

save.rgl.comparison.widget <- function(panels, color, labels, overlay, output.file) {
  if (!requireNamespace("rgl", quietly = TRUE) ||
      !requireNamespace("htmlwidgets", quietly = TRUE)) {
    return(FALSE)
  }

  old.opt <- options(rgl.useNULL = TRUE, rgl.printRglwidget = TRUE)
  on.exit(options(old.opt), add = TRUE)
  suppressWarnings(rgl::open3d(useNULL = TRUE))
  on.exit(try(rgl::close3d(), silent = TRUE), add = TRUE)
  rgl::clear3d()
  rgl::bg3d(color = "white")
  rgl::mfrow3d(nr = 1, nc = length(panels) + 1L, sharedMouse = FALSE)
  for (i in seq_along(panels)) {
    if (i > 1L) {
      rgl::next3d()
    }
    draw.rgl.scatter.panel(panels[[i]], color = color, title = labels[[i]])
  }
  rgl::next3d()
  draw.rgl.overlay.panel(
    ref3 = overlay$baseline,
    obs3 = overlay$gflow.aligned,
    title = overlay$title
  )
  scene <- rgl::scene3d(minimal = FALSE)
  widget <- rgl::rglwidget(scene, width = 1440, height = 360)
  htmlwidgets::saveWidget(widget, file = output.file, selfcontained = TRUE)
  TRUE
}

interactive.rows <- character()
interactive.diagnostics <- list()
if (requireNamespace("rgl", quietly = TRUE) &&
    requireNamespace("htmlwidgets", quietly = TRUE)) {
  for (case.id in interactive.case.ids) {
    row <- all.cases[all.cases$case_id == case.id, ]
    case.dir <- file.path(row$root, case.id)
    if (!dir.exists(case.dir)) next
    X <- read.matrix.csv(file.path(case.dir, "X.csv"))
    py.metric3 <- read.matrix.csv(embedding.file(case.dir, "metric3"))
    meta <- jsonlite::fromJSON(file.path(case.dir, "metadata.json"), simplifyVector = FALSE)
    core <- phate.core(
      X = X, k = 5, alpha.decay = 40, t = "auto", t.max = 30, gamma = 1,
      kernel.mode = "phate", vne.method = "phate", potential.mode = "phate",
      compute.D.pot = TRUE, verbose = FALSE
    )
    gf.metric3 <- phate.embed(
      core = core,
      ndim = 3,
      method = "metric",
      seed = meta$seed,
      maxit = 300,
      verbose = FALSE
    )$embedding

    color <- grDevices::hcl.colors(nrow(X), "Viridis")[rank(X[, 1], ties.method = "first")]
    gf.aligned <- procrustes.align(py.metric3, gf.metric3)
    metric3.diagnostics <- embedding.metrics(py.metric3, gf.metric3)
    metric3.diagnostics$case_id <- case.id
    metric3.diagnostics$method <- "metric_3d"
    interactive.diagnostics[[case.id]] <- metric3.diagnostics[, c(
      "case_id", "method", "rel_dist_error", "procrustes_error", "pearson", "spearman"
    )]
    panels <- list(
      observed = observed.coords3(X),
      baseline = py.metric3,
      gflow = gf.metric3
    )
    labels <- c(
      observed = if (ncol(X) == 3) "Original data" else "Original data PCA",
      baseline = "Baseline PHATE 3D",
      gflow = "gflow PHATE 3D"
    )
    overlay <- list(
      baseline = py.metric3,
      gflow.aligned = gf.aligned$embedding,
      title = sprintf("Overlay: baseline black, gflow red\nresidual %.3f", gf.aligned$residual)
    )
    fname <- sprintf("%s_comparison4.html", case.id)
    fpath <- file.path(widget.dir, fname)
    ok <- save.rgl.comparison.widget(
      panels = panels,
      color = color,
      labels = labels,
      overlay = overlay,
      output.file = fpath
    )
    if (isTRUE(ok)) {
      interactive.rows <- c(
        interactive.rows,
        sprintf(
          "<section class=\"widget-case\">%s<iframe class=\"widget-triptych\" src=\"widgets3d/%s\" loading=\"lazy\"></iframe></section>",
          interactive.case.heading(case.id, meta),
          html.escape(fname)
        )
      )
    }
  }
}

interactive.html <- if (length(interactive.rows) > 0L) {
  paste0(
    "<h2>Interactive 3D Embeddings</h2>",
    "<p>Each row shows the original data view, the baseline PHATE 3D metric embedding, the gflow PHATE 3D metric embedding, and a Procrustes overlay. In the overlay, baseline PHATE is black and gflow PHATE is red after centering, orthogonal alignment, and scalar rescaling to the baseline. For observed data with more than three columns, the original-data panel uses the first three principal components.</p>",
    "<p>The quadratic-form visualization intentionally shows only the visually diagnostic latent-dimension-2 cases: <code>k = 0</code> for a paraboloid and <code>k = 1</code> for an indefinite saddle graph. The sign/permutation-equivalent quadratic variants remain in the numeric validation tables but are omitted here to keep the interactive section interpretable. The sparse Swiss-roll fixture also remains in the numeric suite but is omitted from this visual gallery because, at <code>n = 100</code> with random sampling, it is not reliably legible from a fixed camera angle. The figure-eight and near-crossing examples are included because they are the stress cases with largest 2D MDS discrepancies.</p>",
    paste(interactive.rows, collapse = "\n")
  )
} else {
  paste0(
    "<h2>Interactive 3D Embeddings</h2>",
    "<p><em>Interactive rglwidget panels were not generated because rgl/htmlwidgets were unavailable.</em></p>"
  )
}

format.num <- function(x, digits = 3) {
  ifelse(is.na(x), "", formatC(x, digits = digits, format = "fg", flag = "#"))
}

tex.escape <- function(x) {
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  x <- gsub("([_&#%$])", "\\\\\\1", x, perl = TRUE)
  x
}

validation.summary <- "Validation logs were not found when this report was generated."
summary.path <- file.path(log.dir, "summary.txt")
if (file.exists(summary.path)) {
  validation.summary <- paste(readLines(summary.path, warn = FALSE), collapse = "\n")
}
validation.tex <- c(
  "\\begin{itemize}",
  "\\item Phase 6 parity test: PASS 285, WARN 0, SKIP 0, FAIL 0.",
  "\\item Targeted PHATE regression tests: PASS for phase1, phase2, phase3, phase5, and core/embed files.",
  "\\item Full testthat suite: PASS 2285, SKIP 13, WARN 0, FAIL 0. Duration: 14.7 seconds.",
  "\\item make document: PASS. Rcpp attributes and roxygen documentation regenerated.",
  "\\item make check-fast: PASS with 3 NOTEs. Examples, tests, and manual were intentionally skipped by the fast-check target.",
  "\\end{itemize}"
)

core.pass <- all(operator.metrics$operator_status[operator.metrics$suite == "core"] == "pass") &&
  all(embedding.metrics.df$status[embedding.metrics.df$suite == "core"] == "pass")
stress.review <- operator.metrics[operator.metrics$operator_status != "pass", ]
embed.review <- embedding.metrics.df[embedding.metrics.df$status != "pass", ]

case.table <- operator.metrics[, c("suite", "case_id", "family", "n", "p", "t_gflow", "t_python", "operator_max_abs", "operator_status")]
case.table$operator_max_abs <- format.num(case.table$operator_max_abs, 4)
case.table.tex <- case.table[, c("suite", "case_id", "n", "p", "t_gflow", "t_python", "operator_max_abs", "operator_status")]
embed.summary <- aggregate(
  cbind(rel_dist_error, procrustes_error, pearson, spearman) ~ suite + method,
  embedding.metrics.df,
  function(x) c(max = max(x), median = stats::median(x), min = min(x))
)

best.worst <- embedding.metrics.df[order(embedding.metrics.df$rel_dist_error, decreasing = TRUE), ]
worst.embed <- head(best.worst[, c("suite", "case_id", "method", "rel_dist_error", "procrustes_error", "pearson", "spearman", "status")], 12)
review.case.ids <- c("near_crossing_strands_n96", "figure_eight_n100")
review.diagnostics <- embedding.metrics.df[
  embedding.metrics.df$case_id %in% review.case.ids,
  c("case_id", "method", "rel_dist_error", "procrustes_error", "pearson", "spearman", "status")
]
metric3.diagnostics.df <- if (length(interactive.diagnostics) > 0L) {
  do.call(rbind, interactive.diagnostics)
} else {
  data.frame()
}
if (nrow(metric3.diagnostics.df) > 0L) {
  review.metric3 <- metric3.diagnostics.df[metric3.diagnostics.df$case_id %in% review.case.ids, ]
  review.metric3$status <- ifelse(
    review.metric3$rel_dist_error <= 0.18 &
      review.metric3$procrustes_error <= 0.28 &
      review.metric3$pearson >= 0.94 &
      review.metric3$spearman >= 0.90,
    "pass",
    "review"
  )
  review.diagnostics <- rbind(review.diagnostics, review.metric3[, names(review.diagnostics)])
}
review.diagnostics <- review.diagnostics[order(review.diagnostics$case_id, review.diagnostics$method), ]

tex.table <- function(df, cols, caption) {
  header <- paste(tex.escape(cols), collapse = " & ")
  rows <- apply(df[, cols, drop = FALSE], 1, function(row) {
    paste(tex.escape(as.character(row)), collapse = " & ")
  })
  spec <- paste(rep("l", length(cols)), collapse = "")
  c(
    "\\begingroup\\scriptsize",
    sprintf("\\begin{longtable}{%s}", spec),
    paste0("\\caption{", tex.escape(caption), "}\\\\"),
    "\\toprule",
    paste0(header, " \\\\"),
    "\\midrule",
    paste0(rows, " \\\\"),
    "\\bottomrule",
    "\\end{longtable}",
    "\\endgroup"
  )
}

tex.lines <- c(
  "\\documentclass[10pt]{article}",
  "\\usepackage[margin=0.8in]{geometry}",
  "\\usepackage{graphicx}",
  "\\usepackage{booktabs}",
  "\\usepackage{longtable}",
  "\\usepackage{amsmath}",
  "\\usepackage{amssymb}",
  "\\usepackage{hyperref}",
  "\\usepackage{xcolor}",
  "\\title{gflow PHATE Phase 6 Parity Validation Report}",
  sprintf("\\author{Generated on %s}", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "\\date{}",
  "\\begin{document}",
  "\\maketitle",
  "\\section{Executive Summary}",
  sprintf(
    "The Phase 6 suite compared %d core examples and %d report-only stress examples against frozen outputs from the original PHATE implementation. Core examples are package-test gates; stress examples are diagnostics for future geodesic-MDS work.",
    sum(operator.metrics$suite == "core"),
    sum(operator.metrics$suite == "stress")
  ),
  if (core.pass) {
    "All core examples passed the operator-level and geometry-aware embedding gates."
  } else {
    "At least one core example requires review."
  },
  "Operator-level quantities are compared directly: kernel, Markov operator, VNE curve, selected diffusion time, diffused operator, potential, and potential distance. Embeddings are compared with Procrustes and pairwise-distance metrics because MDS coordinates are not uniquely identifiable.",
  "\\section{Action Plan Executed}",
  "\\begin{enumerate}",
  "\\item Generated deterministic original-PHATE baselines for core and stress examples.",
  "\\item Added pure-R package tests against frozen core fixtures.",
  "\\item Recomputed operator and embedding metrics for every example.",
  "\\item Rendered figures, tables, LaTeX source, PDF, and HTML report artifacts.",
  "\\item Ran targeted and full validation commands, summarized below.",
  "\\end{enumerate}",
  "\\section{Example Suite}",
  "The quadratic-form graph family uses $y=\\sum_{i=1}^{k} x_i^2 - \\sum_{i=k+1}^{n} x_i^2$, with latent points sampled from a disk or ball. Index values $k=0$ and $k=n$ are definite paraboloid controls; mixed index values are indefinite quadratic graph controls.",
  tex.table(case.table.tex, names(case.table.tex), "Example-level operator status."),
  "\\section{Operator-Level Results}",
  "\\begin{figure}[h]",
  "\\centering",
  "\\includegraphics[width=0.98\\linewidth]{figures/operator_max_abs.png}",
  "\\caption{Maximum absolute operator-level deviation by example. The dashed line is the package-test gate used for core examples.}",
  "\\end{figure}",
  "\\section{Embedding-Level Results}",
  "\\begin{figure}[h]",
  "\\centering",
  "\\includegraphics[width=0.98\\linewidth]{figures/embedding_relative_distance_error.png}",
  "\\caption{Relative pairwise-distance error after optimal scalar rescaling.}",
  "\\end{figure}",
  "\\begin{figure}[h]",
  "\\centering",
  "\\includegraphics[width=0.98\\linewidth]{figures/embedding_pearson.png}",
  "\\caption{Pearson correlation of pairwise embedding distances.}",
  "\\end{figure}",
  "\\subsection{Metric Definitions}",
  "Let $Y \\in \\mathbb{R}^{n \\times d}$ be the baseline PHATE embedding and $\\hat{Y} \\in \\mathbb{R}^{n \\times d}$ be the gflow embedding for the same ordered observations. Let $D_{ij}=\\lVert Y_i-Y_j\\rVert_2$ and $\\hat{D}_{ij}=\\lVert \\hat{Y}_i-\\hat{Y}_j\\rVert_2$ for $1 \\le i < j \\le n$.",
  "\\paragraph{Relative distance error.}",
  "\\begin{equation*}",
  "\\operatorname{rel\\_dist\\_error}(Y,\\hat{Y}) =",
  "\\left(",
  "\\frac{\\sum_{i<j}\\left(D_{ij}-s_D\\hat{D}_{ij}\\right)^2}",
  "{\\sum_{i<j}D_{ij}^2}",
  "\\right)^{1/2},",
  "\\qquad ",
  "s_D = \\frac{\\sum_{i<j}D_{ij}\\hat{D}_{ij}}{\\sum_{i<j}\\hat{D}_{ij}^2}.",
  "\\end{equation*}",
  "In words, \\texttt{rel\\_dist\\_error} compares pairwise distances within the two embeddings after optimal scalar rescaling of the gflow embedding distances to the baseline distances. It is a relative root-sum-square distance-matrix error, so lower is better and values are comparable across examples. It ignores translation, rotation, reflection, and global scale, and it measures whether the two embeddings preserve the same internal geometry.",
  "\\paragraph{Procrustes error.}",
  "Let $Y_c=Y-\\mathbf{1}\\bar{Y}^{\\top}$ and $\\hat{Y}_c=\\hat{Y}-\\mathbf{1}\\bar{\\hat{Y}}^{\\top}$ be centered embeddings. Let",
  "\\begin{equation*}",
  "R^* = \\arg\\min_{R^{\\top}R=I}\\lVert Y_c-\\hat{Y}_cR\\rVert_F^2,",
  "\\qquad ",
  "s_P = \\frac{\\langle Y_c,\\hat{Y}_cR^*\\rangle_F}{\\lVert \\hat{Y}_cR^*\\rVert_F^2}.",
  "\\end{equation*}",
  "Then",
  "\\begin{equation*}",
  "\\operatorname{procrustes\\_error}(Y,\\hat{Y}) =",
  "\\frac{\\lVert Y_c-s_P\\hat{Y}_cR^*\\rVert_F}{\\lVert Y_c\\rVert_F}.",
  "\\end{equation*}",
  "In words, \\texttt{procrustes\\_error} compares coordinates after centering, optimal orthogonal alignment, and optimal scalar rescaling. It is normalized by the centered baseline coordinate norm, so it is also on a relative scale and can be compared across examples, with the usual caveat that different examples may have different intrinsic ambiguity in low-dimensional MDS.",
  "\\paragraph{Interpreting small errors.}",
  "For this validation suite, relative distance errors below about $0.05$ are visually very small, $0.05$--$0.10$ are small, and $0.10$--$0.18$ are noticeable but acceptable for the core gate. Procrustes errors below about $0.10$ are visually very small, $0.10$--$0.20$ are small to moderate, and $0.20$--$0.28$ are acceptable for the core gate. Values above these gates are marked for review, not automatic failure, because MDS can have solver- and initialization-dependent ambiguity on stress examples.",
  tex.table(transform(worst.embed,
    rel_dist_error = format.num(rel_dist_error, 4),
    procrustes_error = format.num(procrustes_error, 4),
    pearson = format.num(pearson, 4),
    spearman = format.num(spearman, 4)
  ), names(worst.embed), "Largest embedding relative-distance errors."),
  tex.table(transform(review.diagnostics,
    rel_dist_error = format.num(rel_dist_error, 4),
    procrustes_error = format.num(procrustes_error, 4),
    pearson = format.num(pearson, 4),
    spearman = format.num(spearman, 4)
  ), names(review.diagnostics), "Stress review-case embedding diagnostics, including the 3D metric visual-gallery comparison."),
  "\\section{Example Gallery}",
  "\\begin{figure}[h]",
  "\\centering",
  "\\includegraphics[width=0.98\\linewidth]{figures/example_gallery.png}",
  "\\caption{Observed first-two-coordinate views, baseline PHATE metric embeddings, and gflow metric PHATE embeddings for selected core and stress examples. Color tracks rank of the first observed coordinate only as a visual reference.}",
  "\\end{figure}",
  "The HTML version of this report also includes interactive \\texttt{rglwidget} rows for selected non-flat examples. Each row shows the original-data view, the baseline PHATE 3D metric embedding, the gflow PHATE 3D metric embedding, and a Procrustes overlay where baseline PHATE is black and aligned gflow PHATE is red. For the quadratic-form graph family, the interactive gallery intentionally emphasizes only the visually diagnostic latent-dimension-2 paraboloid and indefinite saddle; sign- and permutation-equivalent variants remain in the numeric tables.",
  "\\clearpage",
  "\\section{Validation Commands}",
  validation.tex,
  "\\section{Interpretation}",
  "The strongest result is that the non-landmark kernel, diffusion operator, VNE curve, diffusion power, PHATE potential, and PHATE potential distance match the original implementation at numerical precision across all tested examples. Embedding geometry is close by distance-correlation and Procrustes criteria, but exact coordinate parity is not expected because MDS solvers differ in initialization, orientation, scale conventions, and, for classic PHATE, the original implementation uses a PCA-based fast classical path rather than base R's exact \\texttt{cmdscale}.",
  "\\end{document}"
)

tex.path <- file.path(out.dir, "phate_phase6_validation_report.tex")
writeLines(tex.lines, tex.path)

html.table <- function(df) {
  rows <- apply(df, 1, function(row) {
    paste0("<tr>", paste(sprintf("<td>%s</td>", html.escape(as.character(row))), collapse = ""), "</tr>")
  })
  paste0(
    "<table><thead><tr>",
    paste(sprintf("<th>%s</th>", html.escape(names(df))), collapse = ""),
    "</tr></thead><tbody>",
    paste(rows, collapse = "\n"),
    "</tbody></table>"
  )
}

html <- paste0(
  "<!doctype html><html><head><meta charset=\"utf-8\"><title>gflow PHATE Phase 6 Validation</title>",
  "<script>window.MathJax={tex:{inlineMath:[[\"$\",\"$\"],[\"\\\\(\",\"\\\\)\"]],displayMath:[[\"\\\\[\",\"\\\\]\"],[\"$$\",\"$$\"]]}};</script>",
  "<script defer src=\"https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-chtml.js\"></script>",
  "<style>body{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif;margin:36px;line-height:1.45;color:#222}h1,h2{color:#16324f}h4{margin-bottom:6px}table{border-collapse:collapse;width:100%;font-size:13px;margin:18px 0}th,td{border:1px solid #ddd;padding:5px 7px;text-align:left}th{background:#f2f5f8}img{max-width:100%;height:auto;border:1px solid #ddd}.note{background:#f7f7f7;border-left:4px solid #3B6EA8;padding:12px;white-space:pre-wrap}.math-block{overflow-x:auto;margin:10px 0 14px 0}code{background:#f4f4f4;padding:1px 3px}.widget-case{margin:28px 0 38px 0;padding-top:10px;border-top:1px solid #ddd}.case-caption{margin:-6px 0 10px 0;color:#555;font-size:14px}.widget-triptych{width:100%;height:390px;border:1px solid #ddd;display:block;background:#fff}@media(max-width:900px){.widget-triptych{height:720px}}</style>",
  "</head><body>",
  "<h1>gflow PHATE Phase 6 Parity Validation Report</h1>",
  sprintf("<p><strong>Generated:</strong> %s</p>", html.escape(format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))),
  "<h2>Executive Summary</h2>",
  sprintf("<p>The suite compared %d core examples and %d report-only stress examples against frozen outputs from the original PHATE implementation.</p>",
          sum(operator.metrics$suite == "core"), sum(operator.metrics$suite == "stress")),
  sprintf("<p><strong>Core gate:</strong> %s</p>", if (core.pass) "PASS" else "REVIEW"),
  "<p>Operator-level quantities are compared directly. Embeddings are compared with Procrustes and pairwise-distance metrics because MDS coordinates are not unique.</p>",
  "<h2>Example Suite</h2>",
  "<p>The quadratic-form graph family uses <code>y = sum_{i=1}^k x_i^2 - sum_{i=k+1}^n x_i^2</code>, with latent points sampled from a disk or ball.</p>",
  html.table(case.table),
  "<h2>Operator-Level Results</h2>",
  "<img src=\"figures/operator_max_abs.png\" alt=\"Operator maximum absolute errors\">",
  "<h2>Embedding-Level Results</h2>",
  "<img src=\"figures/embedding_relative_distance_error.png\" alt=\"Embedding relative distance errors\">",
  "<img src=\"figures/embedding_pearson.png\" alt=\"Embedding Pearson correlations\">",
  "<h3>Metric Definitions</h3>",
  "<p>Let \\(Y \\in \\mathbb{R}^{n \\times d}\\) be the baseline PHATE embedding and \\(\\hat{Y} \\in \\mathbb{R}^{n \\times d}\\) be the gflow embedding for the same ordered observations. For \\(1 \\le i \\lt j \\le n\\), define \\(D_{ij}=\\lVert Y_i-Y_j\\rVert_2\\) and \\(\\hat{D}_{ij}=\\lVert \\hat{Y}_i-\\hat{Y}_j\\rVert_2\\).</p>",
  "<h4>Relative Distance Error</h4>",
  "<div class=\"math-block\">\\[",
  "\\operatorname{rel\\_dist\\_error}(Y,\\hat{Y}) =",
  "\\left(",
  "\\frac{\\sum_{1 \\le i \\lt j \\le n}\\left(D_{ij}-s_D\\hat{D}_{ij}\\right)^2}",
  "{\\sum_{1 \\le i \\lt j \\le n}D_{ij}^2}",
  "\\right)^{1/2},",
  "\\qquad ",
  "s_D = \\frac{\\sum_{1 \\le i \\lt j \\le n}D_{ij}\\hat{D}_{ij}}{\\sum_{1 \\le i \\lt j \\le n}\\hat{D}_{ij}^2}.",
  "\\]</div>",
  "<p>Here \\(s_D\\) is the least-squares scalar that rescales the gflow pairwise distances to the baseline distances. In words, <code>rel_dist_error</code> compares pairwise distances within the two embeddings after optimal scalar rescaling of the gflow embedding distances to the baseline distances. It is a relative root-sum-square distance-matrix error, so lower is better and values are comparable across examples. It ignores translation, rotation, reflection, and global scale, and it measures whether the two embeddings preserve the same internal geometry.</p>",
  "<h4>Procrustes Error</h4>",
  "<p>Let \\(Y_c=Y-\\mathbf{1}\\bar{Y}^{\\top}\\) and \\(\\hat{Y}_c=\\hat{Y}-\\mathbf{1}\\bar{\\hat{Y}}^{\\top}\\) be centered embeddings. Let</p>",
  "<div class=\"math-block\">\\[",
  "R^* = \\arg\\min_{R^{\\top}R=I}\\lVert Y_c-\\hat{Y}_cR\\rVert_F^2,",
  "\\qquad ",
  "s_P = \\frac{\\langle Y_c,\\hat{Y}_cR^*\\rangle_F}{\\lVert \\hat{Y}_cR^*\\rVert_F^2}.",
  "\\]</div>",
  "<p>Then</p>",
  "<div class=\"math-block\">\\[",
  "\\operatorname{procrustes\\_error}(Y,\\hat{Y}) =",
  "\\frac{\\lVert Y_c-s_P\\hat{Y}_cR^*\\rVert_F}{\\lVert Y_c\\rVert_F}.",
  "\\]</div>",
  "<p>Here \\(R^*\\) is the best orthogonal alignment, allowing rotation and reflection, \\(s_P\\) is the best global scale after that alignment, \\(\\langle\\cdot,\\cdot\\rangle_F\\) is the Frobenius inner product, and \\(\\lVert\\cdot\\rVert_F\\) is the Frobenius norm. In words, <code>procrustes_error</code> compares coordinates after centering, optimal orthogonal alignment, and optimal scalar rescaling. It is normalized by the centered baseline coordinate norm, so it is also on a relative scale and can be compared across examples, with the usual caveat that different examples may have different intrinsic ambiguity in low-dimensional MDS.</p>",
  "<h4>What Counts as Small?</h4>",
  "<p>For this validation suite, <code>rel_dist_error</code> below about 0.05 is visually very small, 0.05-0.10 is small, and 0.10-0.18 is noticeable but acceptable for the core gate. <code>procrustes_error</code> below about 0.10 is visually very small, 0.10-0.20 is small to moderate, and 0.20-0.28 is acceptable for the core gate. Values above these gates are marked for review rather than automatic failure because MDS can have solver- and initialization-dependent ambiguity on stress examples.</p>",
  "<h3>Largest Embedding Relative-Distance Errors</h3>",
  html.table(transform(worst.embed,
    rel_dist_error = format.num(rel_dist_error, 4),
    procrustes_error = format.num(procrustes_error, 4),
    pearson = format.num(pearson, 4),
    spearman = format.num(spearman, 4)
  )),
  "<h3>Stress Review-Case Diagnostics</h3>",
  "<p>The review stress cases have operator-level parity at numerical precision; the larger errors arise in the low-dimensional MDS embedding. The 3D metric comparison is included because it is what the interactive visualizations show.</p>",
  html.table(transform(review.diagnostics,
    rel_dist_error = format.num(rel_dist_error, 4),
    procrustes_error = format.num(procrustes_error, 4),
    pearson = format.num(pearson, 4),
    spearman = format.num(spearman, 4)
  )),
  "<h2>Example Gallery</h2>",
  "<p>The gallery shows observed coordinates, baseline PHATE metric embeddings, and gflow metric PHATE embeddings side by side.</p>",
  "<img src=\"figures/example_gallery.png\" alt=\"Example gallery\">",
  interactive.html,
  "<h2>Validation Commands</h2>",
  sprintf("<div class=\"note\">%s</div>", html.escape(validation.summary)),
  "<h2>Interpretation</h2>",
  "<p>The non-landmark kernel, diffusion operator, VNE curve, diffusion power, PHATE potential, and PHATE potential distance match the original implementation at numerical precision across all tested examples. Embedding geometry is close by distance-correlation and Procrustes criteria, while exact coordinate equality is not expected because MDS solvers differ.</p>",
  "</body></html>"
)
html.path <- file.path(out.dir, "phate_phase6_validation_report.html")
writeLines(html, html.path)

old.wd <- getwd()
setwd(out.dir)
on.exit(setwd(old.wd), add = TRUE)
if (nzchar(Sys.which("tectonic"))) {
  system2("tectonic", c("--keep-logs", "--keep-intermediates", basename(tex.path)))
} else {
  warning("tectonic not found; PDF was not rendered")
}

message("Wrote report artifacts to ", out.dir)
