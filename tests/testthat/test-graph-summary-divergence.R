test_that("compute.graph.summary.pmf returns expected degree PMF", {
  g <- list(
    adj_list = list(c(2L), c(1L, 3L), c(2L))
  )

  pmf <- compute.graph.summary.pmf(g, summary = "degree_distribution", return.details = FALSE)

  expect_equal(unname(pmf[c("1", "2")]), c(2/3, 1/3), tolerance = 1e-8)
})

test_that("component_size_distribution is vertex-weighted", {
  g <- list(
    adj_list = list(c(2L), c(1L), integer(0))
  )

  pmf <- compute.graph.summary.pmf(g, summary = "component_size_distribution", return.details = FALSE)

  expect_equal(unname(pmf[c("1", "2")]), c(1/3, 2/3), tolerance = 1e-8)
})

test_that("neighborhood_label_distribution summarizes undirected label pairs", {
  g <- list(
    adj_list = list(c(2L), c(1L, 3L), c(2L))
  )
  labels <- c("A", "A", "B")

  pmf <- compute.graph.summary.pmf(
    g,
    summary = "neighborhood_label_distribution",
    labels = labels,
    return.details = FALSE
  )

  expect_equal(unname(pmf[c("A|A", "A|B")]), c(0.5, 0.5), tolerance = 1e-8)
})

test_that("edge_weight_distribution uses shared bins and returns zero divergence for identical graphs", {
  g <- list(
    adj_list = list(c(2L), c(1L, 3L), c(2L)),
    weight_list = list(c(0.2), c(0.2, 0.8), c(0.8))
  )

  div <- graph.summary.divergence(
    g,
    g,
    summary = "edge_weight_distribution",
    summary.args = list(n.bins = 4L),
    return.details = TRUE
  )

  expect_equal(div$value, 0, tolerance = 1e-12)
  expect_equal(sum(div$pmf1), 1, tolerance = 1e-12)
  expect_equal(sum(div$pmf2), 1, tolerance = 1e-12)
})

test_that("compute.graph.summary.stability preserves k values", {
  g1 <- list(adj_list = list(c(2L), c(1L, 3L), c(2L)))
  g2 <- list(adj_list = list(c(2L, 3L), c(1L), c(1L)))
  g3 <- list(adj_list = list(c(2L), c(1L), integer(0)))

  stab <- compute.graph.summary.stability(
    graphs = list(g1, g2, g3),
    summary = "degree_distribution",
    k.values = 3:5,
    return.details = TRUE
  )

  expect_equal(stab$k.values, 3:5)
  expect_equal(stab$k.transition.values, 3:4)
  expect_length(stab$values, 2)
})

test_that("compute.stability.metrics keeps legacy js.div field while using summary stability", {
  g1 <- list(adj_list = list(c(2L), c(1L, 3L), c(2L)))
  g2 <- list(adj_list = list(c(2L, 3L), c(1L), c(1L)))
  g3 <- list(adj_list = list(c(2L), c(1L), integer(0)))

  graphs <- structure(
    list(
      geom_pruned_graphs = list(g1, g2, g3),
      k_statistics = NULL
    ),
    class = "iknn_graphs",
    kmin = 3L,
    kmax = 5L
  )

  stab <- compute.stability.metrics(graphs, graph.type = "geom")
  deg.stab <- compute.graph.summary.stability(
    graphs = list(g1, g2, g3),
    summary = "degree_distribution",
    k.values = 3:5,
    return.details = TRUE
  )

  expect_true("js.div" %in% names(stab))
  expect_true("summary.stability" %in% names(stab))
  expect_equal(stab$js.div, deg.stab$values, tolerance = 1e-12)
  expect_equal(stab$summary.stability$degree_distribution$values, deg.stab$values, tolerance = 1e-12)
})

test_that("compute.stability.metrics supports alternate summary families through the higher-level API", {
  g1 <- list(adj_list = list(c(2L), c(1L, 3L), c(2L)))
  g2 <- list(adj_list = list(c(2L), c(1L), integer(0)))
  g3 <- list(adj_list = list(c(2L), c(1L, 3L), c(2L)))
  labels <- c("A", "A", "B")

  graphs <- structure(
    list(
      geom_pruned_graphs = list(g1, g2, g3),
      k_statistics = NULL
    ),
    class = "iknn_graphs",
    kmin = 3L,
    kmax = 5L
  )

  stab <- compute.stability.metrics(
    graphs,
    graph.type = "geom",
    summary = "neighborhood_label_distribution",
    labels = labels
  )

  label.stab <- compute.graph.summary.stability(
    graphs = list(g1, g2, g3),
    summary = "neighborhood_label_distribution",
    labels = labels,
    k.values = 3:5,
    return.details = TRUE
  )

  expect_equal(stab$summary, "neighborhood_label_distribution")
  expect_equal(stab$divergence, "js")
  expect_equal(stab$js.div, label.stab$values, tolerance = 1e-12)
  expect_equal(
    stab$summary.stability$neighborhood_label_distribution$values,
    label.stab$values,
    tolerance = 1e-12
  )
})
