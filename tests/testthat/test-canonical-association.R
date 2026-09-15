association.fixture <- function(y = c(0, 3, 1, 2, 0), ...) {
    a <- list(2L, c(1L, 3L), c(2L, 4L), c(3L, 5L), 4L)
    create.basin.complex(a, lapply(a, function(v) rep(1, length(v))), y,
                        method = "trajectory_flow", direction = "both", ...)
}

test_that("canonical association matches independently enumerated archived supports", {
    y <- c(0, 3, 1, 2, 0)
    b <- association.fixture(y)
    pack <- function(e, v) list(vertex = e, value = y[e], hop_idx = 0L,
                                basin_df = cbind(v, rep(0L, length(v))))
    archived <- structure(list(n_vertices = 5L,
        lmax_basins = list(pack(2L, 1:3), pack(4L, 3:5)),
        lmin_basins = list(pack(1L, 1:2), pack(3L, 2:4), pack(5L, 4:5))),
        class = "basins_of_attraction")
    old <- gfcor(y, y, archived, archived, vertex.mass = rep(.2, 5))
    now <- gfcor(y, y, b, b, support.stage = "raw")
    meta <- now$canonical
    now$canonical <- NULL
    expect_equal(now, old)
    expect_equal(sum(meta$vertex.mass), 1)
    expect_identical(meta$y$vertex.id, as.character(1:5))
    expect_length(meta$y$basin.ids$max, 2)
    expect_equal(gfcor(y, y, b, b, support.stage = "retained")$global, old$global)
    m <- gfassoc.membership(b, support.stage = "raw")
    expect_equal(m$max_membership[[3]], c(.5, .5))
    p <- gfassoc.polarity(y, m)
    expect_equal(now$global$A_pol, mean(p$polarity^2))
    expect_lt(now$global$A_pol, 1)
    reversed <- gfcor(y, -y, b, association.fixture(-y), support.stage = "raw")
    expect_equal(reversed$global$A_pol, -now$global$A_pol, tolerance = 1e-9)
    overlap <- gfassoc.overlap(m, m, vertex.mass = 1:5)
    expect_equal(overlap$total_mass, 1)
    expect_equal(sum(overlap$O_pp), 1)
})

test_that("canonical association rejects ambiguous or inconsistent inputs", {
    y <- c(0, 3, 1, 2, 0); b <- association.fixture(y)
    expect_error(gfcor(y, y, b, b), "Choose support.stage")
    expect_error(gfcor(y + 1, y, b, b, support.stage = "raw"), "construction values")
    other <- b; other$graph.input$vertex.id <- rev(other$graph.input$vertex.id)
    expect_error(gfcor(y, y, b, other, support.stage = "raw"), "vertex order")
    other <- b; other$method <- "superlevel_merge_tree"
    expect_error(gfassoc.membership(other, "raw"), "trajectory_flow")
    other <- b; other$direction <- "max"
    expect_error(gfassoc.membership(other, "raw"), "direction")
    other <- b; other$status <- "failed"
    expect_error(gfassoc.membership(other, "raw"), "successful")
    other <- b; other$parameters$simplify.params$expansion$enabled <- TRUE
    expect_error(gfassoc.membership(other, "retained"), "support.filter")
    expect_s3_class(gfassoc.membership(other, "raw"), "gfassoc_membership")
    expect_error(gfcor(y, y, b, b, vertex.mass = rep(0, 5), support.stage = "raw"), "positive finite")
    expect_error(gfcor(y, y, b, b, vertex.mass = c(NA, 1:4), support.stage = "raw"), "finite")
    m <- gfassoc.membership(b, "raw")
    expect_error(gfassoc.polarity(rev(y), m), "construction values")
    expect_error(gfassoc.overlap(m, gfassoc.membership(b, "retained")), "support.stage")
})

test_that("flat and uncovered vertices remain explicitly invalid", {
    y <- rep(1, 5); b <- association.fixture(y)
    x <- gfcor(y, y, b, b, support.stage = "raw")
    expect_false(any(x$vertex$is_valid))
    expect_equal(x$global$n_invalid, 5)
    y <- c(0, 3, 1, 2, 0)
    b <- association.fixture(y, simplify.params = list(support.filter = list(enabled = TRUE, min.basin.size = 100)))
    m <- gfassoc.membership(b, "retained")
    expect_false(any(gfassoc.polarity(y, m)$is_valid))
})

test_that("retained global association excludes uncovered mass explicitly", {
    y <- c(0,3,1,2,0)
    b <- association.fixture(y, simplify.params=list(support.filter=list(enabled=TRUE,min.basin.size=3)))
    x <- gfcor(y,y,b,b,vertex.mass=1:5,support.stage="retained")
    expect_identical(which(x$vertex$is_valid),2:4)
    p <- c(2*2/(2+1e-10)-1,-1,2/(1+1e-10)-1)
    expect_equal(x$global$A_pol,sum((2:4)*p^2)/9,tolerance=1e-12)
    expect_equal(x$overlap$total_mass,1)
    expect_equal(sum(x$overlap$O_mm),.6)
    expect_output(print(x),"3/5 valid vertices \\(normalized mass 0.6000\\)")
    empty <- association.fixture(y,simplify.params=list(support.filter=list(enabled=TRUE,min.basin.size=100)))
    none <- gfcor(y,y,empty,empty,support.stage="retained")
    expect_equal(none$global$n_invalid,5)
    expect_output(print(none),"reported zeros are placeholders")
})
