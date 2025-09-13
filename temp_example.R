### Name: winsorize.zscore
### Title: Winsorized Z-score normalization
### Aliases: winsorize.zscore

### ** Examples

# Example with random data
set.seed(123)
example.data <- matrix(rnorm(100, 5, 2), ncol=5)
# Add some outliers
example.data[1,1] <- 25
example.data[2,3] <- -15
normalized.data <- winsorize.zscore(example.data)



