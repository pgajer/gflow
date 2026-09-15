# Public accessors bridge test fixtures across supported dependency representations.
.test.graph.adj <- function(g, stage = "final") {
    if (inherits(g,"dgraph")) return(getExportedValue("dgraphs","graph.adjacency")(g,stage))
    key <- if(stage == "final") "adj_list" else paste0(gsub(".","_",stage,fixed=TRUE),"_adj_list")
    g[[key]]
}
.test.graph.lengths <- function(g, stage = "final") {
    if (inherits(g,"dgraph")) return(getExportedValue("dgraphs","graph.lengths")(g,stage))
    key <- if(stage == "final") "weight_list" else paste0(gsub(".","_",stage,fixed=TRUE),"_weight_list")
    g[[key]]
}
.test.shortest.paths <- function(a,w,vertices) {
    n <- length(a); d <- matrix(Inf,n,n); diag(d) <- 0
    for(i in seq_len(n)) d[i,a[[i]]] <- w[[i]]
    for(k in seq_len(n)) d <- pmin(d,outer(d[,k],d[k,],`+`))
    d[vertices,vertices,drop=FALSE]
}
.test.circle.graph <- function(n, ...) {
    a <- lapply(seq_len(n),function(i) as.integer(c((i-2L)%%n+1L,i%%n+1L)))
    list(adj.list=a,weight.list=lapply(a,function(v) rep(2*sin(pi/n),length(v))))
}
.test.summary.graph <- function(g) {
    if (.dgraphs.object.api() && !inherits(g,"dgraph"))
        getExportedValue("dgraphs","dgraph")(g$adj_list,g$weight_list)
    else g
}
.test.graph.summary <- function(name) function(...) {
    args <- list(...)
    if(name == "compute.stability.metrics") {
        gs <- args[[1]]
        gs$geom_pruned_graphs <- lapply(gs$geom_pruned_graphs,.test.summary.graph)
        attr(gs,"k.values") <- seq.int(attr(gs,"kmin"),attr(gs,"kmax"))
        args[[1]] <- gs
    } else if(name == "compute.graph.summary.stability") {
        args[[1]] <- lapply(args[[1]],.test.summary.graph)
    } else {
        args[[1]] <- .test.summary.graph(args[[1]])
        if(name == "graph.summary.divergence") args[[2]] <- .test.summary.graph(args[[2]])
    }
    do.call(getExportedValue("dgraphs",name),args)
}
# Corruption tests deliberately alter the active payload, not obsolete fields.
.test.corrupt.graph <- function(g,adj=g$adj_list,lengths=g$weight_list) {
    if(inherits(g,"dgraph")) {
        g$stages$final$adj.list <- adj
        g$stages$final$length.list <- lengths
    } else {
        g$adj_list <- adj;g$weight_list <- lengths
    }
    g
}
