#' Two-Factor Analysis for Contingency Tables with Proportion Tests
#'
#' @description
#' Performs comprehensive analysis of contingency tables for two categorical variables,
#' including proportion tests and optional LaTeX output. Supports both two-sample
#' proportion comparisons between groups and one-sample tests against expected proportions.
#'
#' @param y factor, the grouping variable (e.g., treatment groups or community state types)
#' @param x factor, the binary outcome variable (e.g., success/failure, case/control)
#' @param out.dir character, path to the output directory where results will be saved
#' @param latex.file character or NULL, optional path to LaTeX output file. If NULL,
#'   defaults to \code{file.path(out.dir, paste0(label, ".tex"))}
#' @param label character or NA, label for output files and LaTeX tables. If NA,
#'   defaults to "fp" for file naming
#' @param do.prop.test logical, whether to perform proportion tests (default: TRUE)
#' @param expected.prop numeric or NULL, if provided (must be between 0 and 1),
#'   performs one-sample proportion test against this expected value. If NULL,
#'   performs two-sample tests between groups
#'
#' @return A list with the following components:
#' \describe{
#'   \item{f}{contingency table of frequencies}
#'   \item{p}{matrix of row percentages}
#'   \item{pval}{numeric vector of p-values from proportion tests}
#'   \item{all.f}{formatted results matrix combining percentages, counts, and p-values}
#' }
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Validates inputs and creates output directory if needed
#'   \item Constructs contingency table and calculates row proportions
#'   \item Performs proportion tests:
#'     \itemize{
#'       \item If \code{expected.prop} is provided: one-sample tests comparing
#'         each group's proportion to the expected value
#'       \item If \code{expected.prop} is NULL: two-sample tests comparing
#'         each group against all other groups combined
#'     }
#'   \item Formats results combining proportions, counts, and p-values
#'   \item Saves outputs as RDA, CSV, and optionally LaTeX files
#' }
#'
#' @examples
#' \dontrun{
#' # Create example data
#' set.seed(123)
#' treatment <- factor(rep(c("A", "B", "C"), each = 100))
#' outcome <- factor(rbinom(300, 1, rep(c(0.3, 0.5, 0.4), each = 100)),
#'                   labels = c("Failure", "Success"))
#'
#' # Two-sample proportion tests
#' result <- two.factor.analysis(
#'   y = treatment,
#'   x = outcome,
#'   out.dir = tempdir(),
#'   label = "treatment_analysis"
#' )
#'
#' # One-sample tests against expected proportion of 0.4
#' result2 <- two.factor.analysis(
#'   y = treatment,
#'   x = outcome,
#'   out.dir = tempdir(),
#'   label = "treatment_vs_expected",
#'   expected.prop = 0.4
#' )
#' }
#'
#' @export
#' @importFrom stats prop.test
#' @importFrom utils write.csv
two.factor.analysis <- function(y,
                                x,
                                out.dir,
                                latex.file = NULL,
                                label = NA,
                                do.prop.test = TRUE,
                                expected.prop = NULL) {
    # Input validation
    if (!is.factor(x)) stop("'x' must be a factor")
    if (!is.factor(y)) stop("'y' must be a factor")
    if (is.null(out.dir) || !is.character(out.dir)) {
        stop("'out.dir' must be a character string specifying the output directory")
    }
    if (!is.null(expected.prop)) {
        if (!is.numeric(expected.prop) || length(expected.prop) != 1) {
            stop("'expected.prop' must be a single numeric value")
        }
        if (expected.prop < 0 || expected.prop > 1) {
            stop("'expected.prop' must be between 0 and 1")
        }
    }

    # Create output directory if needed
    if (!dir.exists(out.dir)) {
        dir.create(out.dir, recursive = TRUE)
    }

    # Get factor levels
    x.levels <- levels(x)
    y.levels <- levels(y)

    if (length(x.levels) != 2) {
        warning("'x' should be a binary factor. Results may be unexpected for factors with more than 2 levels.")
    }

    # Create contingency table
    f <- table(x, y)
    p <- 100 * prop.table(f, 2)  # Column proportions (proportion within each y level)

    # Initialize p-values
    pval <- rep(NA_real_, length(y.levels))
    names(pval) <- y.levels

    # Perform proportion tests if requested
    if (do.prop.test) {
        for (i in seq_along(y.levels)) {
            group_data <- y == y.levels[i]

            if (sum(group_data) > 0) {  # Only test if we have data
                if (!is.null(expected.prop)) {
                    # One-sample test against expected proportion
                    # Count successes (second level of x, typically "Success" or "1")
                    success_count <- sum(x[group_data] == x.levels[2])
                    total_count <- sum(group_data)

                    tryCatch({
                        test_result <- prop.test(
                            x = success_count,
                            n = total_count,
                            p = expected.prop,
                            correct = TRUE
                        )
                        pval[i] <- test_result$p.value
                    }, error = function(e) {
                        warning(sprintf(
                            "Could not perform test for level '%s': %s",
                            y.levels[i], e$message
                        ))
                    })
                } else {
                    # Two-sample test: this group vs all others combined
                    # Create binary indicator for current group
                    group_indicator <- as.integer(y == y.levels[i])

                    # Create contingency table for this comparison
                    comparison_table <- table(x, group_indicator)

                    if (ncol(comparison_table) == 2 && all(colSums(comparison_table) > 0)) {
                        tryCatch({
                            # Test proportions of x.levels[2] in each group
                            test_result <- prop.test(
                                x = comparison_table[2, ],  # Counts of x.levels[2]
                                n = colSums(comparison_table)  # Total counts
                            )
                            pval[i] <- test_result$p.value
                        }, error = function(e) {
                            warning(sprintf(
                                "Could not perform test for level '%s': %s",
                                y.levels[i], e$message
                            ))
                        })
                    }
                }
            }
        }
    }

    # Format results table
    all.f <- matrix(
        "",
        nrow = length(x.levels) + 1,
        ncol = length(y.levels),
        dimnames = list(
            c(x.levels, "p-value"),
            y.levels
        )
    )

    # Fill in percentages and counts
    for (j in seq_along(y.levels)) {
        for (i in seq_along(x.levels)) {
            all.f[i, j] <- sprintf("%.1f%% (%d)", p[i, j], f[i, j])
        }
    }

    # Add p-values
    all.f[nrow(all.f), ] <- ifelse(
        is.na(pval),
        "NA",
        formatC(pval, digits = 3, format = "f")
    )

    # Save results
    save_results(all.f, f, p, pval, out.dir, label, latex.file)

    # Return results
    invisible(list(
        f = f,
        p = p,
        pval = pval,
        all.f = all.f
    ))
}


#' Save Analysis Results in Multiple Formats
#'
#' @description
#' Internal function to save analysis results as RDA, CSV, and optionally LaTeX files.
#'
#' @param all.f formatted results matrix
#' @param f frequency table
#' @param p proportion table
#' @param pval p-values vector
#' @param out.dir output directory
#' @param label file label
#' @param latex_file LaTeX output file path
#'
#' @keywords internal
save_results <- function(all.f, f, p, pval, out.dir, label, latex_file) {
    base_name <- if (!is.na(label)) label else "fp"

    # Save RDA file
    save(f, p, pval, all.f,
         file = file.path(out.dir, paste0(base_name, ".rda")))

    # Save CSV file
    write.csv(all.f,
              file = file.path(out.dir, paste0(base_name, ".csv")),
              quote = FALSE)

    # Handle LaTeX output
    if (is.null(latex_file)) {
        latex_file <- file.path(out.dir, paste0(base_name, ".tex"))
    }

    # Create caption from label if available
    caption <- if (!is.na(label)) {
        paste("Results of the proportion test for", label)
    } else {
        "Results of the proportion test"
    }

    create.latex.table(all.f, latex_file, label, caption)
}


#' Create LaTeX Table from Matrix or Data Frame
#'
#' @description
#' Generates a LaTeX table from a matrix or data frame with proper formatting
#' for inclusion in LaTeX documents.
#'
#' @param data matrix or data.frame to convert to LaTeX format
#' @param file character, path to the output LaTeX file
#' @param label character or NA, LaTeX label for cross-referencing. If not NA,
#'   will be appended with ":tbl" suffix
#' @param caption character, table caption (default: "")
#'
#' @return Invisibly returns the LaTeX content as a character vector
#'
#' @details
#' The function creates a properly formatted LaTeX table with:
#' \itemize{
#'   \item Centered alignment
#'   \item Column headers from the data's column names
#'   \item Row names as the first column
#'   \item Optional caption and label for cross-referencing
#'   \item Proper escaping of special LaTeX characters
#' }
#'
#' @examples
#' \dontrun{
#' # Create example data
#' mat <- matrix(1:12, nrow = 3)
#' colnames(mat) <- paste0("Group", 1:4)
#' rownames(mat) <- paste0("Category", 1:3)
#'
#' # Generate LaTeX table
#' create.latex.table(
#'   data = mat,
#'   file = "output.tex",
#'   label = "results",
#'   caption = "Example results table"
#' )
#' }
#'
#' @export
create.latex.table <- function(data, file, label = NA, caption = "") {
    if (!is.matrix(data) && !is.data.frame(data)) {
        stop("'data' must be a matrix or data frame")
    }

    # Convert to matrix if data.frame
    if (is.data.frame(data)) {
        data <- as.matrix(data)
    }

    # Escape special LaTeX characters
    escape_latex <- function(x) {
        x <- gsub("\\\\", "\\\\textbackslash{}", x, fixed = TRUE)
        x <- gsub("&", "\\\\&", x, fixed = TRUE)
        x <- gsub("%", "\\\\%", x, fixed = TRUE)
        x <- gsub("\\$", "\\\\$", x, fixed = TRUE)
        x <- gsub("#", "\\\\#", x, fixed = TRUE)
        x <- gsub("_", "\\\\_", x, fixed = TRUE)
        x <- gsub("\\{", "\\\\{", x, fixed = TRUE)
        x <- gsub("\\}", "\\\\}", x, fixed = TRUE)
        x <- gsub("~", "\\\\textasciitilde{}", x, fixed = TRUE)
        x <- gsub("\\^", "\\\\textasciicircum{}", x, fixed = TRUE)
        x
    }

    # Start LaTeX table
    latex_content <- c(
        "\\begin{center}",
        paste0("\\begin{tabular}{l", paste(rep("r", ncol(data)), collapse = ""), "}"),
        "\\hline\\hline"
    )

    # Add column headers
    headers <- escape_latex(colnames(data))
    header_row <- paste0(
        "& ",  # Empty cell for row names column
        paste(headers, collapse = " & "),
        " \\\\"
    )
    latex_content <- c(latex_content, header_row, "\\hline")

    # Add data rows
    for (i in seq_len(nrow(data))) {
        row_name <- escape_latex(rownames(data)[i])
        row_values <- escape_latex(data[i, ])
        row <- paste0(
            row_name, " & ",
            paste(row_values, collapse = " & "),
            " \\\\"
        )
        latex_content <- c(latex_content, row)
    }

    # Close table
    latex_content <- c(
        latex_content,
        "\\hline",
        "\\end{tabular}"
    )

    # Add caption if provided
    if (nzchar(caption)) {
        latex_content <- c(
            latex_content,
            paste0("\\captionof{table}{", escape_latex(caption), "}")
        )
    }

    # Add label if provided
    if (!is.na(label)) {
        latex_content <- c(
            latex_content,
            paste0("\\label{", label, ":tbl}")
        )
    }

    latex_content <- c(latex_content, "\\end{center}")

    # Write to file
    writeLines(latex_content, file)

    invisible(latex_content)
}


#' Analyze Categorical Proportions with Mixed Effects Support
#'
#' @description
#' Performs comprehensive analysis of proportions across categorical groups,
#' with optional support for repeated measures using mixed effects models.
#' Calculates confidence intervals, relative proportions, and performs
#' appropriate statistical tests.
#'
#' @param x factor or coercible to factor, the categorical grouping variable
#'   (e.g., community state types)
#' @param y binary variable (factor, logical, or numeric 0/1) representing
#'   the outcome of interest
#' @param subj.ids optional factor or coercible to factor, subject identifiers
#'   for repeated measures analysis. If provided, mixed effects models will
#'   be fitted
#' @param pos.label character, label for the positive/case level of the
#'   binary outcome (default: "sPTB")
#' @param neg.label character, label for the negative/control level of the
#'   binary outcome (default: "TB")
#' @param digits integer, number of decimal places for numerical outputs
#'   (default: 2)
#'
#' @return A list containing:
#' \describe{
#'   \item{prop.test.mat}{matrix with columns for counts, proportions,
#'     relative proportions, confidence intervals, and p-values. If
#'     \code{subj.ids} is provided, includes mixed effects estimates}
#'   \item{contingency}{contingency table of x vs y}
#'   \item{overall.prop}{overall proportion of positive outcomes}
#'   \item{latex.caption}{formatted LaTeX caption describing the analysis}
#'   \item{model}{if mixed effects analysis was performed, the fitted
#'     glmer model object; NULL otherwise}
#' }
#'
#' @details
#' The function performs different analyses based on whether subject IDs are provided:
#'
#' \strong{Without subject IDs:}
#' \itemize{
#'   \item Calculates proportions for each group
#'   \item Performs one-sample proportion tests against the overall proportion
#'   \item Computes 95\% confidence intervals
#'   \item Calculates relative proportions compared to overall rate
#' }
#'
#' \strong{With subject IDs (mixed effects):}
#' \itemize{
#'   \item Fits a generalized linear mixed model with random intercepts for subjects
#'   \item Handles convergence issues by trying multiple optimizers
#'   \item Computes model-based confidence intervals and p-values
#'   \item Falls back to simple proportions if model fitting fails
#' }
#'
#' @examples
#' \dontrun{
#' # Simple analysis without repeated measures
#' set.seed(123)
#' cst <- factor(sample(c("I", "II", "III", "IV"), 200, replace = TRUE))
#' outcome <- rbinom(200, 1, c(0.1, 0.3, 0.2, 0.4)[as.numeric(cst)])
#'
#' result <- analyze.categorical.proportions(
#'   x = cst,
#'   y = outcome,
#'   pos.label = "Case",
#'   neg.label = "Control"
#' )
#'
#' # Analysis with repeated measures
#' subjects <- factor(rep(1:50, each = 4))
#' cst_repeated <- factor(sample(c("I", "II", "III", "IV"), 200, replace = TRUE))
#' outcome_repeated <- rbinom(200, 1, 0.3)
#'
#' result_mixed <- analyze.categorical.proportions(
#'   x = cst_repeated,
#'   y = outcome_repeated,
#'   subj.ids = subjects
#' )
#' }
#'
#' @note
#' For mixed effects models, the function requires the \pkg{lme4} package.
#' If convergence issues occur, the function will attempt multiple optimizers
#' before falling back to simple proportion calculations.
#'
#' @seealso \code{\link[lme4]{glmer}} for mixed effects model details
#'
#' @export
#' @importFrom stats prop.test plogis pnorm
analyze.categorical.proportions <- function(x, y,
                                            subj.ids = NULL,
                                            pos.label = "sPTB",
                                            neg.label = "TB",
                                            digits = 2) {
    # Input validation and conversion
    if (!is.factor(x)) {
        x <- factor(x)
    }

    # Convert y to appropriate binary format
    if (!is.factor(y)) {
        if (is.logical(y)) {
            y <- as.integer(y)
        }
        unique_vals <- unique(y[!is.na(y)])
        if (length(unique_vals) != 2) {
            stop("'y' must be a binary variable")
        }
        if (!all(unique_vals %in% c(0, 1))) {
            # Map to 0/1
            y <- as.integer(factor(y)) - 1
        }
    } else {
        # Convert factor to 0/1
        if (length(levels(y)) != 2) {
            stop("'y' must be a binary factor with exactly 2 levels")
        }
        y <- as.integer(y) - 1
    }

    # Create factor with specified labels
    y <- factor(y,
                levels = c(1, 0),
                labels = c(pos.label, neg.label))

    # Calculate overall proportion
    expected.prop <- mean(y == pos.label, na.rm = TRUE)

    # Get factor levels
    x.levels <- levels(x)

    # Create contingency table
    f <- table(x, y)
    p <- 100 * prop.table(f, 1)  # Row proportions

    # Initialize results matrix
    col_names <- c(
        paste0("n(", pos.label, ")"),
        paste0("n(", neg.label, ")"),
        paste0("%(", pos.label, ")"),
        paste0("Relative %"),
        "95% CI Lower",
        "95% CI Upper",
        "p-value"
    )

    if (!is.null(subj.ids)) {
        col_names <- c(
            col_names,
            "Mixed Effect Est.",
            "Mixed Effect SE",
            "Mixed Effect p-value"
        )
    }

    prop.test.mat <- matrix(
        NA,
        nrow = length(x.levels),
        ncol = length(col_names),
        dimnames = list(x.levels, col_names)
    )

    # Fill in counts
    prop.test.mat[, 1] <- f[, 1]  # Positive counts
    prop.test.mat[, 2] <- f[, 2]  # Negative counts

    if (is.null(subj.ids)) {
        # Simple proportion analysis
        for (i in seq_along(x.levels)) {
            group_data <- x == x.levels[i]
            success_count <- sum(y[group_data] == pos.label)
            total_count <- sum(group_data)

            if (total_count > 0) {
                # Perform proportion test
                test_result <- prop.test(
                    x = success_count,
                    n = total_count,
                    p = expected.prop,
                    correct = TRUE
                )

                p.pos <- 100 * test_result$estimate
                rel.p.pos <- test_result$estimate / expected.prop
                lower.CI <- 100 * test_result$conf.int[1]
                upper.CI <- 100 * test_result$conf.int[2]
                pval <- test_result$p.value

                prop.test.mat[i, 3:7] <- c(
                    format(round(p.pos, digits), nsmall = digits),
                    format(round(rel.p.pos, digits), nsmall = digits),
                    format(round(lower.CI, digits), nsmall = digits),
                    format(round(upper.CI, digits), nsmall = digits),
                    format(round(pval, digits + 1), nsmall = digits + 1)
                )
            } else {
                prop.test.mat[i, 3:7] <- c("0.00", "0.00", "0.00", "NA", "NA")
            }
        }

        model <- NULL

    } else {
        # Mixed effects analysis
        if (!requireNamespace("lme4", quietly = TRUE)) {
            stop("Package 'lme4' is required for mixed effects analysis")
        }

        # Ensure subject IDs are factor
        if (!is.factor(subj.ids)) {
            subj.ids <- factor(subj.ids)
        }

        # Create data frame for analysis
        df <- data.frame(
            outcome = as.integer(y == pos.label),
            group = x,
            subject = subj.ids
        )

        # Remove any rows with missing data
        df <- df[complete.cases(df), ]

        # Try fitting model with different optimizers
        model <- NULL
        optimizers <- list(
            list(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)),
            list(optimizer = "Nelder_Mead", optCtrl = list(maxfun = 100000)),
            list(optimizer = "nloptwrap", optCtrl = list(maxeval = 100000))
        )

        for (opt in optimizers) {
            control <- do.call(lme4::glmerControl, opt)

            model_attempt <- tryCatch({
                lme4::glmer(
                    outcome ~ group + (1|subject),
                    family = binomial,
                    data = df,
                    control = control
                )
            }, error = function(e) NULL, warning = function(w) NULL)

            if (!is.null(model_attempt)) {
                model <- model_attempt
                break
            }
        }

        if (!is.null(model)) {
            # Extract model information
            fixed_effects <- lme4::fixef(model)
            vcov_mat <- tryCatch(
                vcov(model),
                error = function(e) diag(length(fixed_effects))
            )

            # Process each group
            for (i in seq_along(x.levels)) {
                group_data <- x == x.levels[i]
                success_count <- sum(y[group_data] == pos.label)
                total_count <- sum(group_data)

                if (total_count > 0) {
                    # Basic proportions
                    p.pos <- 100 * (success_count / total_count)
                    rel.p.pos <- (success_count / total_count) / expected.prop

                    # Mixed effects estimates
                    if (i == 1) {
                        # Reference level
                        est <- fixed_effects[1]
                        se <- sqrt(vcov_mat[1, 1])
                    } else {
                        # Other levels
                        param_name <- paste0("group", x.levels[i])
                        if (param_name %in% names(fixed_effects)) {
                            est <- fixed_effects[1] + fixed_effects[param_name]
                            idx <- which(names(fixed_effects) == param_name)
                            se <- sqrt(
                                vcov_mat[1, 1] +
                                vcov_mat[idx, idx] +
                                2 * vcov_mat[1, idx]
                            )
                        } else {
                            est <- NA
                            se <- NA
                        }
                    }

                    if (!is.na(est) && !is.na(se)) {
                        # Convert to probability scale
                        prob <- plogis(est)
                        ci <- plogis(est + c(-1.96, 1.96) * se)
                        pval <- 2 * (1 - pnorm(abs(est / se)))

                        prop.test.mat[i, 3:7] <- c(
                            format(round(p.pos, digits), nsmall = digits),
                            format(round(rel.p.pos, digits), nsmall = digits),
                            format(round(100 * ci[1], digits), nsmall = digits),
                            format(round(100 * ci[2], digits), nsmall = digits),
                            format(round(pval, digits + 1), nsmall = digits + 1)
                        )

                        prop.test.mat[i, 8:10] <- c(
                            format(round(est, digits), nsmall = digits),
                            format(round(se, digits), nsmall = digits),
                            format(round(pval, digits + 1), nsmall = digits + 1)
                        )
                    }
                }
            }
        } else {
            warning("Mixed effects model failed to converge. Returning simple proportions only.")
            # Fill in simple proportions
            for (i in seq_along(x.levels)) {
                group_data <- x == x.levels[i]
                success_count <- sum(y[group_data] == pos.label)
                total_count <- sum(group_data)

                if (total_count > 0) {
                    p.pos <- 100 * (success_count / total_count)
                    rel.p.pos <- (success_count / total_count) / expected.prop

                    prop.test.mat[i, 3:4] <- c(
                        format(round(p.pos, digits), nsmall = digits),
                        format(round(rel.p.pos, digits), nsmall = digits)
                    )
                }
            }
        }
    }

    # Generate LaTeX caption
    latex.caption <- if (is.null(subj.ids)) {
        sprintf(
            paste0(
                "Analysis of %s rates across groups. For each group, the table shows ",
                "the number of %s and %s cases, the percentage of %s, relative risk compared ",
                "to the study-wide %s rate (%.1f\\%%), 95\\%% confidence intervals, and p-values."
            ),
            pos.label, pos.label, neg.label, pos.label, pos.label,
            100 * expected.prop
        )
    } else {
        sprintf(
            paste0(
                "Mixed effects analysis of %s rates across groups, accounting for repeated measures. ",
                "Shows counts, proportions, and model-based estimates with standard errors. ",
                "Overall %s rate: %.1f\\%%."
            ),
            pos.label, pos.label, 100 * expected.prop
        )
    }

    # Return results
    list(
        prop.test.mat = prop.test.mat,
        contingency = f,
        overall.prop = expected.prop,
        latex.caption = latex.caption,
        model = model
    )
}


#' Summarize and Test CST Patterns with Binary Outcomes
#'
#' @description
#' Analyzes the relationship between Community State Types (CSTs) and binary
#' outcomes in longitudinal data where subjects have multiple CST measurements
#' but a single, constant outcome. Provides two methods for summarizing
#' repeated CST measurements before testing.
#'
#' @param x factor or character vector of CST measurements, converted to
#'   factor if necessary
#' @param y binary outcome vector that must be constant within each subject.
#'   Can be logical, numeric (0/1), or factor
#' @param subj.ids vector of subject identifiers corresponding to CST
#'   measurements. Must be same length as \code{x}
#' @param pos.label character, label for positive/case outcome
#'   (default: "sPTB")
#' @param neg.label character, label for negative/control outcome
#'   (default: "TB")
#' @param summary.method character, method for summarizing repeated CSTs.
#'   Must be either "most_frequent" (use each subject's modal CST) or
#'   "proportions" (calculate time spent in each CST). Default: "most_frequent"
#'
#' @return A list whose structure depends on \code{summary.method}:
#'
#' If \code{summary.method = "most_frequent"}:
#' \describe{
#'   \item{prop.test.mat}{matrix of proportion test results for modal CSTs}
#'   \item{contingency}{contingency table of modal CSTs vs outcome}
#'   \item{overall.prop}{overall proportion of positive outcome}
#'   \item{latex.caption}{LaTeX caption for results table}
#' }
#'
#' If \code{summary.method = "proportions"}:
#' \describe{
#'   \item{proportions}{matrix where rows are subjects and columns are CSTs,
#'     values are percentage of measurements in each CST}
#'   \item{tests}{list of Wilcoxon test results comparing CST proportions
#'     between outcome groups}
#'   \item{summary}{data frame with median proportions by outcome group
#'     and p-values for each CST}
#' }
#'
#' @details
#' This function addresses the challenge of analyzing longitudinal microbiome
#' data where the outcome of interest (e.g., preterm birth) is measured once
#' per subject but microbiome composition (CST) is measured multiple times.
#'
#' The "most_frequent" method identifies each subject's predominant CST and
#' performs standard proportion tests. The "proportions" method calculates
#' how much time each subject spends in each CST and uses Wilcoxon tests
#' to compare these proportions between outcome groups.
#'
#' @examples
#' \dontrun{
#' # Simulated longitudinal CST data
#' set.seed(123)
#' n_subjects <- 50
#' n_timepoints <- 5
#'
#' # Generate subject IDs and CSTs
#' subj_ids <- rep(1:n_subjects, each = n_timepoints)
#' csts <- sample(c("I", "II", "III", "IV"), n_subjects * n_timepoints,
#'                replace = TRUE, prob = c(0.4, 0.3, 0.2, 0.1))
#'
#' # Binary outcome (constant per subject)
#' outcomes <- rep(rbinom(n_subjects, 1, 0.3), each = n_timepoints)
#'
#' # Analyze using most frequent CST
#' result_mode <- summarize.and.test.cst(
#'   x = csts,
#'   y = outcomes,
#'   subj.ids = subj_ids,
#'   summary.method = "most_frequent"
#' )
#'
#' # Analyze using CST proportions
#' result_prop <- summarize.and.test.cst(
#'   x = csts,
#'   y = outcomes,
#'   subj.ids = subj_ids,
#'   summary.method = "proportions"
#' )
#'
#' print(result_prop$summary)
#' }
#'
#' @note
#' This function assumes that the outcome is constant within subjects.
#' If outcomes vary within subjects, consider alternative approaches or
#' time-to-event analyses.
#'
#' @seealso
#' \code{\link{analyze.weighted.cst}} for time-weighted CST analysis
#'
#' @export
#' @importFrom stats wilcox.test median
summarize.and.test.cst <- function(x, y, subj.ids,
                                   pos.label = "sPTB",
                                   neg.label = "TB",
                                   summary.method = c("most_frequent", "proportions")) {

    # Validate inputs
    summary.method <- match.arg(summary.method)

    if (length(x) != length(y) || length(x) != length(subj.ids)) {
        stop("x, y, and subj.ids must have the same length")
    }

    # Convert to factors if needed
    if (!is.factor(x)) x <- factor(x)

    # Check that outcome is constant within subjects
    outcome_check <- tapply(y, subj.ids, function(vals) length(unique(vals)))
    if (any(outcome_check > 1)) {
        stop("Outcome (y) must be constant within each subject")
    }

    # Create subject-level summary
    unique_subjects <- unique(subj.ids)
    subj_summary <- data.frame(
        subject = unique_subjects,
        outcome = y[match(unique_subjects, subj.ids)]
    )

    # Convert outcome to binary factor
    if (!is.factor(subj_summary$outcome)) {
        subj_summary$outcome <- factor(
            as.integer(subj_summary$outcome),
            levels = c(1, 0),
            labels = c(pos.label, neg.label)
        )
    }

    if (summary.method == "proportions") {
        # Calculate proportion of time in each CST for each subject
        cst_levels <- levels(x)
        n_subjects <- length(unique_subjects)
        n_csts <- length(cst_levels)

        cst_props <- matrix(
            0,
            nrow = n_subjects,
            ncol = n_csts,
            dimnames = list(unique_subjects, cst_levels)
        )

        # Calculate proportions
        for (i in seq_along(unique_subjects)) {
            sid <- unique_subjects[i]
            subject_csts <- x[subj.ids == sid]
            tab <- table(subject_csts)
            props <- tab / length(subject_csts)
            cst_props[as.character(sid), names(props)] <- props
        }

        # Perform Wilcoxon tests for each CST
        tests <- list()
        summary_results <- data.frame(
            CST = character(n_csts),
            median_case = numeric(n_csts),
            median_control = numeric(n_csts),
            p_value = numeric(n_csts),
            stringsAsFactors = FALSE
        )

        for (i in seq_along(cst_levels)) {
            cst <- cst_levels[i]
            props <- cst_props[, cst]

            # Separate by outcome
            case_props <- props[subj_summary$outcome == pos.label]
            control_props <- props[subj_summary$outcome == neg.label]

            # Perform test
            if (length(case_props) > 0 && length(control_props) > 0) {
                test_result <- wilcox.test(case_props, control_props)
                tests[[cst]] <- test_result

                summary_results[i, ] <- data.frame(
                    CST = cst,
                    median_case = median(case_props * 100),
                    median_control = median(control_props * 100),
                    p_value = test_result$p.value
                )
            } else {
                tests[[cst]] <- NULL
                summary_results[i, ] <- data.frame(
                    CST = cst,
                    median_case = ifelse(length(case_props) > 0,
                                         median(case_props * 100), NA),
                    median_control = ifelse(length(control_props) > 0,
                                            median(control_props * 100), NA),
                    p_value = NA
                )
            }
        }

        # Format results
        summary_results$median_case <- sprintf("%.1f%%", summary_results$median_case)
        summary_results$median_control <- sprintf("%.1f%%", summary_results$median_control)
        summary_results$p_value[!is.na(summary_results$p_value)] <-
            sprintf("%.3f", summary_results$p_value[!is.na(summary_results$p_value)])
        summary_results$p_value[is.na(summary_results$p_value)] <- "NA"

        return(list(
            proportions = cst_props * 100,
            tests = tests,
            summary = summary_results
        ))

    } else {  # most_frequent method
        # Find most frequent CST for each subject
        modal_cst <- character(length(unique_subjects))

        for (i in seq_along(unique_subjects)) {
            sid <- unique_subjects[i]
            subject_csts <- x[subj.ids == sid]
            cst_table <- table(subject_csts)
            # In case of tie, take first in alphabetical order
            modal_cst[i] <- names(cst_table)[which.max(cst_table)]
        }

        # Add to subject summary
        subj_summary$modal_cst <- factor(modal_cst, levels = levels(x))

        # Use analyze.categorical.proportions for the analysis
        result <- analyze.categorical.proportions(
            x = subj_summary$modal_cst,
            y = subj_summary$outcome,
            pos.label = pos.label,
            neg.label = neg.label
        )

        # Modify the caption to reflect the summary method
        result$latex.caption <- gsub(
            "Analysis of",
            "Analysis of modal CST and",
            result$latex.caption
        )

        return(result)
    }
}


#' Time-Weighted CST Analysis for Binary Outcomes
#'
#' @description
#' Performs time-weighted analysis of Community State Types (CSTs) in
#' longitudinal studies where sampling intervals may be irregular. Weights
#' each CST observation by its duration to accurately represent the time
#' spent in different microbial states.
#'
#' @param x vector of CST measurements (factor or character), must be
#'   ordered chronologically within each subject
#' @param y binary outcome vector, must be constant within each subject
#' @param subj.ids vector of subject identifiers, same length as \code{x}
#' @param time.points optional numeric vector of time points for each
#'   measurement. If NULL, assumes equal intervals. Must be chronologically
#'   ordered within each subject
#' @param pos.label character, label for positive outcome (default: "sPTB")
#' @param neg.label character, label for negative outcome (default: "TB")
#'
#' @return A list containing:
#' \describe{
#'   \item{weighted_proportions}{matrix where rows are subjects and columns
#'     are CSTs. Values represent the proportion of time spent in each CST}
#'   \item{summary}{data frame with columns:
#'     \itemize{
#'       \item CST: Community State Type identifier
#'       \item median_case: Median time proportion in CST for positive outcome group
#'       \item median_control: Median time proportion in CST for negative outcome group
#'       \item p_value: P-value from Wilcoxon rank-sum test
#'     }}
#'   \item{subject_data}{data frame containing subject IDs and outcomes}
#' }
#'
#' @details
#' The time-weighting algorithm:
#' \enumerate{
#'   \item For each subject, calculates the duration of each CST period
#'   \item When \code{time.points} is provided:
#'     \itemize{
#'       \item First observation: duration = \code{(time[2] - time[1]) / 2}
#'       \item Middle observations: duration = \code{(time[i+1] - time[i-1]) / 2}
#'       \item Last observation: duration = \code{(time[n] - time[n-1]) / 2}
#'     }
#'   \item When \code{time.points} is NULL, each observation gets equal weight
#'   \item Sums time by CST and divides by total time to get proportions
#'   \item Compares proportions between outcome groups using Wilcoxon tests
#' }
#'
#' This approach is particularly useful when:
#' \itemize{
#'   \item Sampling is irregular (e.g., clinical visits at varying intervals)
#'   \item Some CST states are transient and might be missed with regular sampling
#'   \item The duration of CST states is biologically meaningful
#' }
#'
#' @examples
#' \dontrun{
#' # Example with gestational age time points
#' set.seed(123)
#'
#' # Generate data for 20 subjects
#' subjects <- rep(1:20, each = 5)
#' weeks <- rep(c(12, 16, 20, 28, 32), 20)  # Gestational weeks
#'
#' # Simulate CSTs with some subjects more stable than others
#' csts <- character(100)
#' for (i in 1:20) {
#'   if (i <= 10) {  # Stable subjects
#'     csts[(i-1)*5 + 1:5] <- sample(c("I", "III"), 5, replace = TRUE, prob = c(0.8, 0.2))
#'   } else {  # Variable subjects
#'     csts[(i-1)*5 + 1:5] <- sample(c("I", "III", "IV"), 5, replace = TRUE)
#'   }
#' }
#'
#' # Binary outcome
#' outcomes <- rep(c(rep(0, 10), rep(1, 10)), each = 5)
#'
#' # Analyze with time weighting
#' result <- analyze.weighted.cst(
#'   x = csts,
#'   y = outcomes,
#'   subj.ids = subjects,
#'   time.points = weeks
#' )
#'
#' # View results
#' print(result$summary)
#'
#' # Visualize distributions
#' library(ggplot2)
#' prop_df <- as.data.frame(result$weighted_proportions)
#' prop_df$outcome <- result$subject_data$outcome
#' prop_df$subject <- rownames(prop_df)
#'
#' # Plot CST I proportions by outcome
#' ggplot(prop_df, aes(x = outcome, y = I)) +
#'   geom_boxplot() +
#'   geom_point(position = position_jitter(width = 0.1)) +
#'   labs(y = "Proportion of time in CST I",
#'        title = "Time-weighted CST proportions by outcome")
#' }
#'
#' @note
#' \itemize{
#'   \item Time points must be in the same units throughout the dataset
#'   \item The function assumes measurements are ordered chronologically within subjects
#'   \item For regular sampling, the weighted and unweighted approaches will give similar results
#' }
#'
#' @references
#' Gajer P, et al. (2012). Temporal dynamics of the human vaginal microbiota.
#' Science Translational Medicine, 4(132), 132ra52.
#'
#' @seealso
#' \code{\link{summarize.and.test.cst}} for unweighted analysis approaches
#'
#' @export
#' @importFrom stats wilcox.test median
analyze.weighted.cst <- function(x,
                                 y,
                                 subj.ids,
                                 time.points = NULL,
                                 pos.label = "sPTB",
                                 neg.label = "TB") {

    # Input validation
    if (length(x) != length(y) || length(x) != length(subj.ids)) {
        stop("x, y, and subj.ids must have the same length")
    }

    if (!is.null(time.points)) {
        if (length(time.points) != length(x)) {
            stop("time.points must have the same length as x")
        }
        if (!is.numeric(time.points)) {
            stop("time.points must be numeric")
        }
    }

    # Convert x to factor if needed
    if (!is.factor(x)) x <- factor(x)

    # Check that outcome is constant within subjects
    outcome_check <- tapply(y, subj.ids, function(vals) length(unique(vals)))
    if (any(outcome_check > 1)) {
        stop("Outcome (y) must be constant within each subject")
    }

    # Get unique subjects and their outcomes
    unique_subjects <- unique(subj.ids)
    n_subjects <- length(unique_subjects)

    subj_data <- data.frame(
        subject = unique_subjects,
        outcome = factor(
            y[match(unique_subjects, subj.ids)],
            levels = c(1, 0),
            labels = c(pos.label, neg.label)
        )
    )

    # Get CST levels
    cst_levels <- levels(x)
    n_csts <- length(cst_levels)

    # Initialize weighted proportions matrix
    weighted_props <- matrix(
        0,
        nrow = n_subjects,
        ncol = n_csts,
        dimnames = list(as.character(unique_subjects), cst_levels)
    )

    # Calculate weighted proportions for each subject
    for (i in seq_len(n_subjects)) {
        sid <- unique_subjects[i]
        subj_idx <- which(subj.ids == sid)

        subj_csts <- x[subj_idx]
        n_obs <- length(subj_csts)

        if (is.null(time.points)) {
            # Equal weights if no time points provided
            weights <- rep(1, n_obs)
        } else {
            # Calculate time-based weights
            subj_times <- time.points[subj_idx]

            # Ensure times are sorted
            if (any(diff(subj_times) < 0)) {
                warning(paste("Time points not in order for subject", sid,
                              "- sorting automatically"))
                time_order <- order(subj_times)
                subj_times <- subj_times[time_order]
                subj_csts <- subj_csts[time_order]
            }

            # Calculate durations
            if (n_obs == 1) {
                weights <- 1
            } else {
                weights <- numeric(n_obs)
                # First observation
                weights[1] <- (subj_times[2] - subj_times[1]) / 2
                # Middle observations
                if (n_obs > 2) {
                    for (j in 2:(n_obs - 1)) {
                        weights[j] <- (subj_times[j + 1] - subj_times[j - 1]) / 2
                    }
                }
                # Last observation
                weights[n_obs] <- (subj_times[n_obs] - subj_times[n_obs - 1]) / 2
            }
        }

        # Calculate weighted proportions
        total_weight <- sum(weights)
        for (cst in cst_levels) {
            cst_weights <- weights[subj_csts == cst]
            weighted_props[i, cst] <- sum(cst_weights) / total_weight
        }
    }

    # Perform statistical tests
    summary_results <- data.frame(
        CST = cst_levels,
        median_case = numeric(n_csts),
        median_control = numeric(n_csts),
        p_value = numeric(n_csts),
        stringsAsFactors = FALSE
    )

    for (i in seq_along(cst_levels)) {
        cst <- cst_levels[i]
        props <- weighted_props[, cst]

        # Separate by outcome
        case_idx <- subj_data$outcome == pos.label
        control_idx <- subj_data$outcome == neg.label

        case_props <- props[case_idx]
        control_props <- props[control_idx]

        # Calculate medians
        median_case <- if (length(case_props) > 0) median(case_props) else NA
        median_control <- if (length(control_props) > 0) median(control_props) else NA

        # Perform Wilcoxon test
        if (length(case_props) > 0 && length(control_props) > 0) {
            test_result <- wilcox.test(case_props, control_props)
            p_val <- test_result$p.value
        } else {
            p_val <- NA
        }

        summary_results[i, ] <- data.frame(
            CST = cst,
            median_case = median_case,
            median_control = median_control,
            p_value = p_val
        )
    }

    # Format summary table
    summary_results$median_case <- sprintf("%.1f%%", 100 * summary_results$median_case)
    summary_results$median_control <- sprintf("%.1f%%", 100 * summary_results$median_control)
    summary_results$p_value[!is.na(summary_results$p_value)] <-
        sprintf("%.3f", summary_results$p_value[!is.na(summary_results$p_value)])
    summary_results$p_value[is.na(summary_results$p_value)] <- "NA"

    # Return results
    list(
        weighted_proportions = weighted_props,
        summary = summary_results,
        subject_data = subj_data
    )
}
