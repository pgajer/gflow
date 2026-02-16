source("~/current_projects/gflow/inst/examples/sigmoid_fns.R")

## Generate the 3x3 comparison figure
plot.sigmoid.comparison(
    alpha.values = c(0.125, 0.25, 0.5, 1, 2, 3, 5, 8, 12),
    main.title = "Comparison of Sigmoid Transformations"
)


pdf("sigmoid_comparison.pdf", width = 10, height = 10)
plot.sigmoid.comparison(
    alpha.values = c(0.5, 1, 2, 3, 5, 8, 12, 20, 50),
    main.title = "Comparison of Sigmoid Transformations"
)
dev.off()

## Alternative: PNG for screen viewing
png("sigmoid_comparison.png", width = 900, height = 900, res = 100)
plot.sigmoid.comparison(
    alpha.values = c(0.5, 1, 2, 3, 5, 8, 12, 20, 50),
    main.title = "Comparison of Sigmoid Transformations"
)
dev.off()
