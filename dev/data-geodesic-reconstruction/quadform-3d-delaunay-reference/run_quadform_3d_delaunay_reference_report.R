#!/usr/bin/env Rscript

parse.args <- function(args) {
    out <- list()
    for (arg in args) {
        if (!grepl("^--", arg)) next
        kv <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1L]]
        key <- kv[[1L]]
        val <- if (length(kv) > 1L) paste(kv[-1L], collapse = "=") else TRUE
        out[[key]] <- val
    }
    out
}

`%||%` <- function(x, y) if (is.null(x)) y else x

html.escape <- function(x) {
    x <- gsub("&", "&amp;", x, fixed = TRUE)
    x <- gsub("<", "&lt;", x, fixed = TRUE)
    x <- gsub(">", "&gt;", x, fixed = TRUE)
    x
}

html.table <- function(df, digits = 4) {
    df <- as.data.frame(df)
    for (nm in names(df)) {
        if (is.numeric(df[[nm]])) {
            df[[nm]] <- signif(df[[nm]], digits)
        }
    }
    header <- paste(sprintf("<th>%s</th>", html.escape(names(df))), collapse = "")
    rows <- apply(df, 1L, function(row) {
        paste0("<tr>", paste(sprintf("<td>%s</td>", html.escape(row)), collapse = ""), "</tr>")
    })
    paste0("<table><thead><tr>", header, "</tr></thead><tbody>",
           paste(rows, collapse = "\n"), "</tbody></table>")
}

sample.parameter.3d <- function(n, radius, shape, seed) {
    set.seed(seed)
    if (identical(shape, "ball")) {
        dirs <- matrix(rnorm(n * 3L), ncol = 3L)
        dirs <- dirs / sqrt(rowSums(dirs^2))
        return(dirs * (radius * runif(n)^(1 / 3)))
    }
    matrix(runif(n * 3L, -radius, radius), ncol = 3L)
}

rel.error.summary <- function(D, D.ref) {
    gflow::summarize.isometry.deviation(D, D.ref, scale = TRUE)
}

args <- parse.args(commandArgs(trailingOnly = TRUE))
out.dir <- args[["output.dir"]]
if (is.null(out.dir)) {
    out.dir <- file.path(
        "dev", "data-geodesic-reconstruction", "quadform-3d-delaunay-reference"
    )
}
n.ref <- args[["n.ref"]]
if (is.null(n.ref)) {
    n.ref <- c(300L, 600L, 1200L)
} else {
    n.ref <- as.integer(strsplit(n.ref, ",", fixed = TRUE)[[1L]])
}
sample.n <- as.integer(args[["sample.n"]] %||% 35L)
candidate.multiplier <- as.numeric(args[["candidate.multiplier"]] %||% 4)
boundary.fraction <- as.numeric(args[["boundary.fraction"]] %||% 0.25)
edge.length.factor <- as.numeric(args[["edge.length.factor"]] %||% 3)

if (!requireNamespace("pkgload", quietly = TRUE)) {
    stop("The 'pkgload' package is required to run this report script.")
}
pkgload::load_all(".", quiet = TRUE)
if (!requireNamespace("geometry", quietly = TRUE)) {
    stop("The optional 'geometry' package is required. Install it to run this report.")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("The optional 'ggplot2' package is required to render figures.")
}

report.dir <- file.path(out.dir, "report")
fig.dir <- file.path(report.dir, "figures")
dir.create(fig.dir, recursive = TRUE, showWarnings = FALSE)
csv.dir <- file.path(out.dir, "results")
dir.create(csv.dir, recursive = TRUE, showWarnings = FALSE)

cases <- data.frame(
    case_id = c("ball_index3_paraboloid", "ball_index1_mixed", "cube_index2_mixed"),
    label = c(
        "Ball domain, index k = 3 (elliptic paraboloid in parameter dimension 3)",
        "Ball domain, index k = 1 (mixed-sign quadratic hypersurface)",
        "Cube domain, index k = 2 (mixed-sign quadratic hypersurface)"
    ),
    domain_shape = c("ball", "ball", "cube"),
    index_k = c(3L, 1L, 2L),
    seed = c(101L, 102L, 103L),
    stringsAsFactors = FALSE
)

results <- list()
payloads <- list()
for (case.row in seq_len(nrow(cases))) {
    case <- cases[case.row, ]
    X <- sample.parameter.3d(sample.n, radius = 1, shape = case$domain_shape,
                             seed = case$seed)
    for (nr in n.ref) {
        message(sprintf("case=%s n.ref=%d", case$case_id, nr))
        elapsed <- system.time({
            ref <- gflow::quadform.delaunay.geodesic.distances(
                X,
                index.k = case$index_k,
                domain.radius = 1,
                domain.shape = case$domain_shape,
                n.ref = nr,
                seed = case$seed + nr,
                candidate.multiplier = candidate.multiplier,
                boundary.fraction = boundary.fraction,
                edge.length.factor = edge.length.factor
            )
        })
        key <- paste(case$case_id, nr, sep = "::")
        payloads[[key]] <- ref
        results[[length(results) + 1L]] <- data.frame(
            case_id = case$case_id,
            label = case$label,
            domain_shape = case$domain_shape,
            index_k = case$index_k,
            n_ref_target = nr,
            n_reference_vertices = ref$n_reference_vertices,
            n_delaunay_edges = ref$n_delaunay_edges,
            n_edges = ref$n_edges,
            epsilon = ref$epsilon,
            filter_factor_used = ref$filter_factor_used,
            runtime_sec = unname(elapsed[["elapsed"]]),
            stringsAsFactors = FALSE
        )
    }
}
run.summary <- do.call(rbind, results)

error.rows <- list()
for (case_id in unique(run.summary$case_id)) {
    case.rows <- run.summary[run.summary$case_id == case_id, ]
    finest <- max(case.rows$n_ref_target)
    D.finest <- payloads[[paste(case_id, finest, sep = "::")]]$distances
    for (nr in case.rows$n_ref_target) {
        ref <- payloads[[paste(case_id, nr, sep = "::")]]
        err <- rel.error.summary(ref$distances, D.finest)
        error.rows[[length(error.rows) + 1L]] <- cbind(
            case_id = case_id,
            n_ref_target = nr,
            reference_n_ref_target = finest,
            err,
            stringsAsFactors = FALSE
        )
    }
}
error.summary <- do.call(rbind, error.rows)
summary.out <- merge(run.summary, error.summary,
                     by = c("case_id", "n_ref_target"), all.x = TRUE)
write.csv(run.summary, file.path(csv.dir, "delaunay_reference_run_summary.csv"),
          row.names = FALSE)
write.csv(error.summary, file.path(csv.dir, "delaunay_reference_error_summary.csv"),
          row.names = FALSE)

plot.df <- error.summary[error.summary$rel_rms_error > 0, , drop = FALSE]
p1 <- ggplot2::ggplot(plot.df, ggplot2::aes(
    x = n_ref_target, y = rel_rms_error, color = case_id, group = case_id
)) +
    ggplot2::geom_line() +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_x_continuous(trans = "log10") +
    ggplot2::scale_y_continuous(trans = "log10") +
    ggplot2::labs(
        x = "Target reference vertices",
        y = "Relative RMS error vs finest run",
        color = "Case",
        title = "Delaunay reference convergence"
    ) +
    ggplot2::theme_minimal(base_size = 12)
fig1 <- file.path(fig.dir, "convergence_rel_rms.png")
ggplot2::ggsave(fig1, p1, width = 8, height = 5, dpi = 150)

p2 <- ggplot2::ggplot(run.summary, ggplot2::aes(
    x = n_reference_vertices, y = runtime_sec, color = case_id, group = case_id
)) +
    ggplot2::geom_line() +
    ggplot2::geom_point(size = 2) +
    ggplot2::labs(
        x = "Actual reference vertices",
        y = "Runtime (seconds)",
        color = "Case",
        title = "Reference-oracle runtime"
    ) +
    ggplot2::theme_minimal(base_size = 12)
fig2 <- file.path(fig.dir, "runtime_by_vertices.png")
ggplot2::ggsave(fig2, p2, width = 8, height = 5, dpi = 150)

html <- paste0(
"<!doctype html><html><head><meta charset='utf-8'>",
"<title>3D Quadform Delaunay Reference Geodesics</title>",
"<script defer src='https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js'></script>",
"<style>",
"body{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif;max-width:1120px;margin:32px auto;line-height:1.5;color:#202124;padding:0 24px}",
"h1,h2{line-height:1.2} table{border-collapse:collapse;width:100%;font-size:13px;margin:16px 0 28px}",
"th,td{border:1px solid #d6d6d6;padding:6px 8px;text-align:right} th{text-align:left;background:#f5f5f5}",
"td:first-child,th:first-child{text-align:left}.note{background:#f7fafc;border-left:4px solid #3b82f6;padding:12px 16px;margin:16px 0}",
".fig{max-width:900px;width:100%;border:1px solid #ddd}",
"</style></head><body>",
"<h1>3D Quadratic Hypersurface Reference Geodesics</h1>",
"<p><strong>Build timestamp:</strong> ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "</p>",
"<div class='note'>This report validates the first experimental reference-geodesic oracle for parameter dimension 3. ",
"The method constructs an approximate epsilon-net in the 3D parameter domain, includes the sample points exactly, ",
"adds boundary candidates, builds the Delaunay one-skeleton, weights its edges by exact quadratic-hypersurface segment lengths, ",
"and computes shortest-path distances between sample vertices.</div>",
"<h2>Why This Exists</h2>",
"<p>The 2D quadform benchmark uses a regular parameter grid. In parameter dimension 3, a tensor grid grows too quickly and ",
"is a poor geometric primitive for an intrinsic reference oracle. The Delaunay one-skeleton over a quasi-uniform epsilon-net ",
"is a more geometric discretization: it adapts to the domain as a point cloud while retaining a triangulation-derived local graph.</p>",
"<p>For a quadratic graph hypersurface ",
"<span>q(x)=\\(\\sum_{i=1}^{k} c_i x_i^2 - \\sum_{i=k+1}^{3} c_i x_i^2\\)</span>, each Delaunay edge \\([u,v]\\) is weighted by the exact length ",
"of the embedded straight parameter segment \\(x(t)=u+t(v-u)\\): ",
"\\(\\int_0^1 \\sqrt{\\|v-u\\|_2^2 + (a+bt)^2}\\,dt\\). ",
"Shortest-path distances on this weighted reference graph approximate surface geodesic distances.</p>",
"<h2>Convergence Check</h2>",
"<p>Each curve compares a run at target reference size \\(N_{ref}\\) with the finest run for the same case after optimal scalar calibration. ",
"Lower relative RMS error indicates that the reference discretization is stabilizing as the epsilon-net is refined.</p>",
sprintf("<img class='fig' src='figures/%s' alt='Convergence relative RMS error'>", basename(fig1)),
"<h2>Runtime</h2>",
"<p>The timing includes epsilon-net construction, Delaunay tessellation, edge filtering/weighting, and shortest paths between the sample vertices.</p>",
sprintf("<img class='fig' src='figures/%s' alt='Runtime by reference vertices'>", basename(fig2)),
"<h2>Run Summary</h2>",
html.table(run.summary),
"<h2>Error Summary</h2>",
html.table(error.summary),
"<h2>Interpretation</h2>",
"<p>This is a readiness report, not yet the final high-dimensional benchmark. The quantities to watch are monotone stabilization of relative error, ",
"finite connected reference graphs after long-edge filtering, and runtime scaling. The current implementation intentionally uses Qhull through ",
"the optional R package <code>geometry</code>; the exact edge-length kernel is already implemented in C++ and can be reused if the Delaunay backend ",
"is later moved fully into package C++.</p>",
"</body></html>"
)
html.file <- file.path(report.dir, "quadform_3d_delaunay_reference_report.html")
writeLines(html, html.file)
message("Wrote ", normalizePath(html.file))
