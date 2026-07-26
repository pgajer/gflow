#' Extract Most Abundant ASVs with Taxonomy Information
#'
#' Extracts the n most abundant Amplicon Sequence Variants (ASVs) from a sample
#' or group of samples, optionally including their taxonomic classification.
#'
#' @param id Character string or numeric index identifying the sample in the
#'   abundance matrix S. If character, must match a row name in S.
#' @param S Numeric matrix of ASV abundances with samples in rows and ASVs in
#'   columns. Row names should contain sample identifiers and column names
#'   should contain ASV identifiers.
#' @param bm.tx Optional named vector containing taxonomic classifications for
#'   ASVs. Names must match column names in S. If NULL, only abundances are
#'   returned.
#' @param n.prof Integer specifying the number of most abundant ASVs to report.
#'   Default is 5.
#' @param k.neighbors Integer specifying the number of nearest neighbors to
#'   include in abundance calculations. Default is 1 (only target sample).
#'   If greater than 1, returns mean abundances across the k nearest neighbors
#'   of the target sample (including the target itself).
#' @param verbose Logical indicating whether to print the resulting profile
#'   matrix. Default is FALSE.
#'
#' @return A matrix with one column containing abundance values (rounded to 2
#'   significant figures). If bm.tx is provided, additional columns contain
#'   taxonomic information. Row names are ASV identifiers.
#'
#' @details When k.neighbors > 1, the function uses k-nearest neighbors based on
#'   Euclidean distance in the abundance space to compute mean abundances. This
#'   can help smooth profiles in sparse data.
#'
#' @importFrom utils head
#' @importFrom FNN get.knn
#'
#' @examples
#' # Create example data
#' S <- matrix(runif(100), nrow = 10, ncol = 10)
#' rownames(S) <- paste0("Sample", 1:10)
#' colnames(S) <- paste0("ASV", 1:10)
#'
#' # Extract profile for first sample
#' prof <- prof.fn("Sample1", S, n.prof = 3)
#'
#' # With taxonomy
#' taxonomy <- paste0("Taxon", 1:10)
#' names(taxonomy) <- colnames(S)
#' prof_with_tax <- prof.fn(1, S, bm.tx = taxonomy, n.prof = 3, verbose = TRUE)
#'
#' @export
prof.fn <- function(id,
                    S,
                    bm.tx = NULL,
                    n.prof = 5,
                    k.neighbors = 1,
                    verbose = FALSE) {
    # Input validation
    if (!is.matrix(S) && !is.data.frame(S)) {
        stop("S must be a matrix or data.frame")
    }
    S <- as.matrix(S)
    if (!is.numeric(S)) {
        stop("S must contain numeric values")
    }
    if (k.neighbors < 1) {
        stop("k.neighbors must be at least 1")
    }
    if (n.prof < 1) {
        stop("n.prof must be at least 1")
    }

    # Handle id parameter - convert to numeric index if character
    if (is.character(id)) {
        id_idx <- which(rownames(S) == id)
        if (length(id_idx) == 0) {
            stop("ID '", id, "' not found in row names of S")
        }
        id_idx <- id_idx[1]  # In case of duplicates
    } else if (is.numeric(id)) {
        if (id < 1 || id > nrow(S)) {
            stop("Numeric id must be between 1 and ", nrow(S))
        }
        id_idx <- as.integer(id)
    } else {
        stop("id must be character or numeric")
    }

    # Check k.neighbors doesn't exceed available samples
    if (k.neighbors > nrow(S)) {
        warning("k.neighbors exceeds number of samples. Using all ", nrow(S), " samples.")
        k.neighbors <- nrow(S)
    }

    if (k.neighbors > 1) {
        # Check if FNN package is available
        if (!requireNamespace("FNN", quietly = TRUE)) {
            stop("Package 'FNN' is required for k-nearest neighbor functionality")
        }

        # Get k nearest neighbors using FNN::get.knn
        nn <- FNN::get.knn(S, k = k.neighbors - 1)  # -1 because we'll include the point itself
        nn.i <- nn$nn.index

        # Get indices for the target point: itself plus its k-1 nearest neighbors
        nn_indices <- c(id_idx, nn.i[id_idx, ])

        # Calculate mean abundances across selected points
        x <- colMeans(S[nn_indices, , drop = FALSE])
    } else {
        # Single sample case
        x <- as.numeric(S[id_idx, ])
    }

    names(x) <- colnames(S)

    # Keep only non-zero abundances
    x <- x[x > 0]

    if (length(x) == 0) {
        warning("No non-zero abundances found")
        return(NULL)
    }

    # Get top n.prof abundances
    n_to_select <- min(length(x), n.prof)
    p <- sort(x, decreasing = TRUE)[seq_len(n_to_select)]

    # Prepare result
    if (!is.null(bm.tx)) {
        # Check that ASV names are in taxonomy
        missing_taxa <- setdiff(names(p), names(bm.tx))
        if (length(missing_taxa) > 0) {
            warning("Taxonomy information missing for ", length(missing_taxa),
                   " ASVs: ", paste(utils::head(missing_taxa, 3), collapse = ", "),
                   if (length(missing_taxa) > 3) "..." else "")
        }

        sp <- bm.tx[names(p)]
        sp[is.na(sp)] <- "Unknown"  # Handle missing taxonomy

        result <- cbind(Taxonomy = sp, Abundance = signif(p, digits = 2))
    } else {
        result <- cbind(Abundance = signif(p, digits = 2))
    }

    if (verbose) {
        print(result)
    }

    return(result)
}

#' Standardize String for File Names
#'
#' Removes special characters from a string and replaces spaces with underscores
#' to create valid file names or standardized identifiers.
#'
#' @param str Character string to be standardized. Typically a name that needs
#'   to be converted to a valid file name or identifier.
#'
#' @return Character string with special characters removed or replaced:
#'   \itemize{
#'     \item Commas, spaces, and forward slashes are replaced with underscores
#'     \item Parentheses, colons, asterisks, apostrophes, and plus signs are removed
#'   }
#'
#' @details This function is useful for converting names (such as metabolite
#'   names, sample identifiers, or feature names) into standardized formats
#'   suitable for use as file names or database keys. The function performs
#'   literal character replacement, not regular expression matching.
#'
#' @examples
#' # Standardize a metabolite name
#' standardize.string("Glucose (6-phosphate)")
#' # Returns: "Glucose_6-phosphate"
#'
#' # Standardize a complex name
#' standardize.string("Compound A/B + C's metabolite*")
#' # Returns: "Compound_A_B__Cs_metabolite"
#'
#' # Use for creating file names
#' sample_name <- "Sample 1: Day 3 (replicate)"
#' file_name <- paste0(standardize.string(sample_name), ".csv")
#' # Results in: "Sample_1_Day_3_replicate.csv"
#'
#' @noRd
standardize.string <- function(str) {
    if (!is.character(str)) {
        stop("Input must be a character string")
    }

    if (length(str) != 1) {
        stop("Input must be a single character string")
    }

    standardized.str <- str

    # Replace characters with underscores
    standardized.str <- gsub(",", "_", standardized.str, fixed = TRUE)
    standardized.str <- gsub(" ", "_", standardized.str, fixed = TRUE)
    standardized.str <- gsub("/", "_", standardized.str, fixed = TRUE)

    # Remove characters completely
    standardized.str <- gsub("(", "", standardized.str, fixed = TRUE)
    standardized.str <- gsub(")", "", standardized.str, fixed = TRUE)
    standardized.str <- gsub(":", "", standardized.str, fixed = TRUE)
    standardized.str <- gsub("*", "", standardized.str, fixed = TRUE)
    standardized.str <- gsub("'", "", standardized.str, fixed = TRUE)
    standardized.str <- gsub("+", "", standardized.str, fixed = TRUE)

    return(standardized.str)
}
