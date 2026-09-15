check.log <- readLines("gflow.Rcheck/00check.log", warn = FALSE)
if (any(grepl("^Status:.*(ERROR|WARNING)", check.log))) {
    stop("The built-package check has errors or warnings; inspect 00check.log.")
}
library(gflow)
cat(R.version.string, "\ndgraphs:", as.character(packageVersion("dgraphs")), "\n")
info <- .Call("S_gflow_openmp_diag", PACKAGE = "gflow")
print(info)
if (Sys.getenv("GFLOW_DISABLE_OPENMP") == "1") stopifnot(!info$openmp_compiled)
if (Sys.getenv("GFLOW_BUILD_PROFILE") == "dev") stopifnot(info$openmp_compiled)
# End-to-end graph construction via the dependency's public interface.
if ("dgraph" %in% getNamespaceExports("dgraphs")) {
    graph <- dgraphs::dgraph(list(2L,c(1L,3L),2L),list(1,c(1,1),1))
    adj <- dgraphs::graph.adjacency(graph); lengths <- dgraphs::graph.lengths(graph)
} else {
    graph <- dgraphs::create.path.graph(list(2L,c(1L,3L),2L),list(1,c(1,1),1),h=1)
    adj <- graph$adj.list; lengths <- graph$edge.length.list
}
b <- create.basin.complex(adj,lengths,c(0,2,1),method="trajectory_flow",direction="both")
stopifnot(b$status == "ok", nrow(get.basin.table(b)) > 0,
          nrow(get.basin.membership(b)) > 0,
          all(is.finite(lcor(adj,lengths,c(0,2,1),c(1,0,2),hop.radius=2))))
cat("Installed dependency-to-basin workflow passed.\n")
