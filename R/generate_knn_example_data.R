
# Function to plot the generated points using base R graphics
mixed.points.plot <- function(points.df,
                            main = "Mixed Point Distribution",
                            point.size = 1.5,
                            point.colors = c(gaussian = "blue", line = "red")) {

    # Set up the plotting area with some margin space for labels
    par(mar = c(4, 4, 3, 1))  # Bottom, left, top, right margins

    # Calculate plot limits with some padding
    x.range <- range(points.df$x)
    y.range <- range(points.df$y)
    x.padding <- diff(x.range) * 0.1
    y.padding <- diff(y.range) * 0.1

    # Create empty plot with proper limits
    plot(NA, NA,
         xlim = x.range + c(-x.padding, x.padding),
         ylim = y.range + c(-y.padding, y.padding),
         xlab = "X coordinate",
         ylab = "Y coordinate",
         main = main,
         type = "n")  # 'n' means no plotting yet

    # Add a grid for better readability
    grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")

    # Plot points by type
    # First plot gaussian points
    gaussian.points <- points.df$type == "gaussian"
    points(points.df$x[gaussian.points],
           points.df$y[gaussian.points],
           pch = 19,  # Filled circles
           cex = point.size,
           col = point.colors["gaussian"])

    # Then plot line points
    line.points <- points.df$type == "line"
    points(points.df$x[line.points],
           points.df$y[line.points],
           pch = 19,  # Filled circles
           cex = point.size,
           col = point.colors["line"])

    # Add a legend
    legend("topright",
           legend = c("Gaussian cluster", "Line-aligned points"),
           pch = 19,
           col = unname(point.colors),
           bty = "n",  # No box around legend
           cex = 0.8)

    # Reset graphical parameters to default
    par(mar = c(5, 4, 4, 2) + 0.1)
}

# Function to generate a mixture of points from:
# 1) A 2D Gaussian centered at origin
# 2) Points along the positive part of y-axis with Gaussian noise
generate.mixed.points <- function(
    n.points,              # Total number of points to generate
    gaussian.prop = 0.3,   # Proportion of points from the central Gaussian
    gaussian.sd = 0.5,     # Standard deviation for the central Gaussian
    line.sd = 0.1,        # Standard deviation for noise around the line
    y.min = 0,            # Minimum y-coordinate for line points
    y.max = 3             # Maximum y-coordinate for line points
) {
    # Input validation
    if (n.points < 1) stop("n.points must be positive")
    if (gaussian.prop < 0 || gaussian.prop > 1) stop("gaussian.prop must be between 0 and 1")
    if (gaussian.sd <= 0 || line.sd <= 0) stop("Standard deviations must be positive")
    if (y.min >= y.max) stop("y.min must be less than y.max")

    # Calculate number of points for each component
    n.gaussian <- round(n.points * gaussian.prop)
    n.line <- n.points - n.gaussian

    # Generate points from the central Gaussian
    gaussian.points <- matrix(rnorm(2 * n.gaussian, 0, gaussian.sd),
                            ncol = 2,
                            dimnames = list(NULL, c("x", "y")))

    # Generate points along the positive part of y-axis
    # First generate evenly spaced y-coordinates within the specified range
    y.coords <- seq(y.min, y.max, length.out = n.line)

    # Add Gaussian noise to create tubular neighborhood
    line.points <- cbind(
        rnorm(n.line, 0, line.sd),    # x-coordinates with noise
        y.coords                       # y-coordinates along the line
    )
    colnames(line.points) <- c("x", "y")

    # Combine both sets of points
    all.points <- rbind(gaussian.points, line.points)

    # Convert to data frame and add point type labels
    points.df <- data.frame(
        x = all.points[, "x"],
        y = all.points[, "y"],
        type = factor(c(rep("gaussian", n.gaussian),
                       rep("line", n.line)))
    )

    return(points.df)
}
