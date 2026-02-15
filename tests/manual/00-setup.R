# Common test data generator
generate.test.data.1d <- function(n.pts = 50, noise_sd = 1, seed = 123) {
  set.seed(seed)

  gm <- generate.1d.gaussian.mixture(
    n.points = n.pts,
    x.knot = c(0, 10),
    y.knot = c(10, 2.5),
    sd.knot = 1.5,
    x.offset = 3
  )

  X <- cbind(gm$x)
  y.smooth <- gm$y.true
  eps <- rnorm(n.pts, 0, noise_sd)
  y <- y.smooth + eps

  list(
    X = X,
    y = y,
    y.smooth = y.smooth,
    gm = gm
  )
}

# Diagnostic plotting function
plot.test.data <- function(test.data) {
  par(mfrow = c(1, 2))
  plot(test.data$X[,1], test.data$y.smooth,
       main = "True Function", xlab = "x", ylab = "y")
  plot(test.data$X[,1], test.data$y,
       main = "Noisy Observations", xlab = "x", ylab = "y")
  par(mfrow = c(1, 1))
}

#' Generate a Timestamp
#'
#' This function returns a timestamp in the format <year>_<month>_<day>_<24h><min><sec>.
#'
#' @return A string representing the current timestamp.
#'
#' @examples
#' timestamp <- generate.timestamp()
#' print(timestamp)
#'
#' @export
generate.timestamp <- function() {
  # Get the current date and time
  current.time <- Sys.time()

  # Format the date and time components
  year <- format(current.time, "%Y")
  month <- format(current.time, "%m")
  day <- format(current.time, "%d")
  hour <- format(current.time, "%H")
  minute <- format(current.time, "%M")
  second <- format(current.time, "%S")

  # Combine components into the desired format
  timestamp <- paste0(year, "_", month, "_", day, "_", hour, minute, second)

  return(timestamp)
}

cat("Setup complete. Use generate.test.data.1d() to create test data.\n")
