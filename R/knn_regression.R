## library(FNN)
## library(caret) # For cross-validation

## # Function to compute the mean of y over the extended kNN set N[i]
## knn_regression <- function(X, y, k) {
##   # Get indices of the k nearest neighbors for each point
##   knn_indices <- get.knn(X, k = k)$nn.index

##   # Initialize predictions vector
##   predictions <- numeric(length(y))

##   # For each observation
##   for (i in 1:length(y)) {
##     # Get the extended neighborhood (kNN plus the point itself)
##     neighborhood <- c(i, knn_indices[i,])

##     # Calculate proportion of 1's in the neighborhood
##     predictions[i] <- mean(y[neighborhood])
##   }

##   return(predictions)
## }

## # Function to find optimal k using cross-validation
## find_optimal_k <- function(X, y, kmin, kmax, nfolds = 10) {
##   # Set up cross-validation
##   set.seed(123)
##   folds <- create.folds(y, k = nfolds)

##   # Initialize error matrix
##   cv_errors <- matrix(0, nrow = kmax - kmin + 1, ncol = nfolds)

##   # For each fold
##   for (fold in 1:nfolds) {
##     # Split data
##     train_indices <- unlist(folds[-fold])
##     test_indices <- folds[[fold]]

##     X_train <- X[train_indices, ]
##     y_train <- y[train_indices]
##     X_test <- X[test_indices, ]
##     y_test <- y[test_indices]

##     # Try different k values
##     for (k in kmin:kmax) {
##       idx <- k - kmin + 1

##       # Train kNN model and get predictions
##       pred <- knn_regression(X_train, y_train, k)

##       # Apply model to test data
##       test_knn_indices <- get.knnx(data = X_train, query = X_test, k = k)$nn.index
##       test_predictions <- numeric(length(test_indices))

##       for (i in 1:length(test_indices)) {
##         # Get predictions based on training set neighbors
##         neighborhood <- test_knn_indices[i,]
##         test_predictions[i] <- mean(y_train[neighborhood])
##       }

##       # Calculate error (using mean squared error)
##       cv_errors[idx, fold] <- mean((y_test - test_predictions)^2)
##     }
##   }

##   # Average errors across folds
##   mean_cv_errors <- rowMeans(cv_errors)

##   # Find k with minimum error
##   optimal_k <- kmin + which.min(mean_cv_errors) - 1

##   return(list(
##     optimal_k = optimal_k,
##     cv_errors = mean_cv_errors,
##     k_values = kmin:kmax
##   ))
## }

## # Example usage:
## # result <- find_optimal_k(X, y, kmin = 1, kmax = 20)
## # optimal_k <- result$optimal_k
## # plot(result$k_values, result$cv_errors, type = "b", xlab = "k", ylab = "CV Error")


## ###
## ### Using the caret package's built-in cross-validation functionality:
## ###
## ## library(caret)
## ## library(FNN)

## # Custom kNN model function
## knn_model <- function(X, y, k) {
##   # Compute distances
##   model <- list(
##     X = X,
##     y = y,
##     k = k
##   )
##   return(model)
## }

## # Prediction function for custom kNN model
## predict_knn <- function(model, newdata) {
##   X_train <- model$X
##   y_train <- model$y
##   k <- model$k

##   # Get indices of k nearest neighbors
##   nn_indices <- get.knnx(data = X_train, query = newdata, k = k)$nn.index

##   # Compute predictions
##   predictions <- numeric(nrow(newdata))
##   for (i in 1:nrow(newdata)) {
##     predictions[i] <- mean(y_train[nn_indices[i,]])
##   }

##   return(predictions)
## }

## # Set up training control for CV
## ctrl <- trainControl(
##   method = "cv",
##   number = 10,
##   summaryFunction = defaultSummary
## )

## # Train and tune the model
## set.seed(123)
## knn_cv <- train(
##   x = X,
##   y = y,
##   method = "knn",
##   tuneGrid = data.frame(k = kmin:kmax),
##   trControl = ctrl,
##   metric = "RMSE"
## )

## Get the best k
##optimal_k <- knn_cv$bestTune$k


###
###  Using kknn package
###
## install.packages("kknn")
## library(kknn)

# Cross-validation with kknn
## set.seed(123)
## knn_cv <- train.kknn(factor(y) ~ ., data = data.frame(y = y, X),
##                       kmax = kmax,
##                       distance = 2,  # Euclidean distance
##                       kernel = "rectangular", # Equal weighting
##                       scale = TRUE)

## Get the best k
## optimal_k <- knn_cv$best.parameters$k
