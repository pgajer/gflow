# Internal compatibility helper used by the protected graphMScx plot method.
#
# Generic graph manipulation belongs to dgraphs. This small conversion remains
# private because plot.graphMScx is part of the protected basin/complex closure.
convert.edge.label.list.to.edge.label.vector <- function(edge.label.list,
                                                         rm.duplicates = TRUE) {
    edge.label <- c()
    for (i in seq_along(edge.label.list)) {
        for (label in edge.label.list[[i]]) {
            edge.label <- c(edge.label, label)
        }
    }

    if (rm.duplicates) {
        edge.label <- unique(edge.label)
    }

    edge.label
}
