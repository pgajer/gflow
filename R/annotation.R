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

#' Display and Cluster Taxon-Specific Embeddings
#'
#' Reads taxon-specific 3D embedding data and performs or loads clustering
#' analysis. Optionally reorders clusters by size and creates 3D visualization.
#'
#' @param taxon Character string for the taxon name used in file names.
#' @param data.dir Directory containing data files (default: current dir if NULL/"").
#' @param cltr.from.scratch Logical; compute clustering or load existing (default TRUE).
#' @param min.pts Numeric vector of HDBSCAN minPts values (default seq(5, 50, 5)).
#' @param reorder.cltrs Logical; reorder clusters by decreasing size (default TRUE).
#' @param show.plot Logical; show 3D plot if possible (default TRUE).
#' @param n.cores Integer; parallel cores for your internal functions (default 10).
#' @param verbose Logical; print progress (default TRUE).
#'
#' @return (Invisibly) a list with taxon, cltr.ext, cltr, X, pacmap.file, cltr.ext.file
#' @importFrom utils read.csv write.csv
#' @export
show.tx <- function(taxon,
                    data.dir = NULL,
                    cltr.from.scratch = TRUE,
                    min.pts = seq(5, 50, by = 5),
                    reorder.cltrs = TRUE,
                    show.plot = TRUE,
                    n.cores = 10,
                    verbose = TRUE) {

  ## ---- Input validation ----------------------------------------------------
  if (missing(taxon) || is.null(taxon) || !is.character(taxon) || length(taxon) != 1L) {
    stop("'taxon' must be a non-empty character string of length 1")
  }

  # Normalize data.dir and precompute file paths
  if (is.null(data.dir) || identical(data.dir, "")) {
    data.dir <- "."
  }

  if (!dir.exists(data.dir)) {
    if (isTRUE(cltr.from.scratch)) {
      # Friendly behavior: create it if we're going to write outputs
      dir.create(data.dir, recursive = TRUE, showWarnings = FALSE)
      if (!dir.exists(data.dir)) {
        stop("Could not create directory '", data.dir, "'")
      }
    } else {
      stop("Directory '", data.dir, "' does not exist")
    }
  }

  pacmap.file   <- file.path(data.dir, sprintf("%s_pacmap.csv",          taxon))
  cltr.file     <- file.path(data.dir, sprintf("%s_pacmap_cltr.csv",     taxon))
  cltr.ext.file <- file.path(data.dir, sprintf("%s_pacmap_cltr_ext.csv", taxon))

  ## ---- Read embedding ------------------------------------------------------
  if (!file.exists(pacmap.file)) {
    stop("PaCMAP file not found: ", pacmap.file)
  }

  X <- utils::read.csv(pacmap.file, row.names = 1, check.names = FALSE)
  if (!is.data.frame(X) && !is.matrix(X)) {
    stop("Invalid embedding data format; expected data.frame or matrix")
  }

  X <- as.matrix(X)
  storage.mode(X) <- "double"
  n.samples <- nrow(X)
  if (n.samples < 2L) {
    stop("Insufficient samples for clustering (n = ", n.samples, ")")
  }

  max.soft.K <- min(20L, n.samples - 1L)

  ## ---- Clustering (from scratch or load) ----------------------------------
  if (isTRUE(cltr.from.scratch)) {
    if (verbose) message("Performing clustering from scratch...")

    if (length(min.pts) > 1L) {
      # keep min.pts < n.samples
      keep <- min.pts < n.samples
      min.pts <- as.integer(min.pts[keep])
      if (!length(min.pts)) {
        stop("All 'min.pts' values exceed number of samples")
      }

      r <- hdbscan.cltr(
        X,
        min.pts = min.pts,
        soft.K = max.soft.K,
        n.cores = n.cores,
        verbose = verbose,
        method  = "dunn"
      )
      cltr <- r$dunn.cltr
    } else {
      if (min.pts >= n.samples) {
        stop("min.pts (", min.pts, ") must be less than number of samples (", n.samples, ")")
      }
      cltr <- dbscan::hdbscan(X, minPts = as.integer(min.pts))$cluster
    }

    cltr <- as.integer(cltr)

    # Write main clustering (no row names; we read back and +1 later)
    utils::write.csv(
      data.frame(cluster = as.integer(cltr) - 1L),
      file = cltr.file,
      row.names = FALSE
    )

    # Extended clustering (keep row names so we can align later if needed)
    cltr.ext <- kNN.cltr.imputation(X, cltr, K = max.soft.K)
    cltr.ext <- as.integer(cltr.ext)
    names(cltr.ext) <- rownames(X)

    utils::write.csv(
      data.frame(cluster_extended = as.integer(cltr.ext)),
      file = cltr.ext.file,
      row.names = TRUE
    )

  } else {
    if (verbose) message("Loading existing clustering...")

    if (!file.exists(cltr.file)) {
      stop("Clustering file not found: ", cltr.file)
    }

    cltr_data <- utils::read.csv(cltr.file, check.names = FALSE)
    cltr <- as.integer(cltr_data[[1L]]) + 1L

    if (!file.exists(cltr.ext.file)) {
      if (verbose) message("Creating extended clustering (file missing)...")
      cltr.ext <- kNN.cltr.imputation(X, cltr, K = max.soft.K)
      cltr.ext <- as.integer(cltr.ext)
      names(cltr.ext) <- rownames(X)

      utils::write.csv(
        data.frame(cluster_extended = as.integer(cltr.ext)),
        file = cltr.ext.file,
        row.names = TRUE
      )
    } else {
      cltr.ext <- utils::read.csv(cltr.ext.file, row.names = 1, check.names = FALSE)[, 1]
      cltr.ext <- as.integer(cltr.ext)
      # keep names if present; otherwise align to X
      if (is.null(names(cltr.ext))) names(cltr.ext) <- rownames(X)
    }
  }

  ## ---- Optional reordering by size ----------------------------------------
  if (isTRUE(reorder.cltrs)) {
    if (verbose) message("Reordering clusters by size...")

    # Main clusters (skip 0 as "noise", if that's your convention)
    if (length(table(cltr[cltr != 0L])) > 1L) {
      cltr_new <- cltr
      cltr_new[cltr_new != 0L] <- clusters.reorder(cltr_new[cltr_new != 0L])
      if (!identical(cltr_new, cltr)) {
        cltr <- as.integer(cltr_new)
        utils::write.csv(
          data.frame(cluster = as.integer(cltr) - 1L),
          file = cltr.file,
          row.names = FALSE
        )
      }
    }

    # Extended clusters
    if (length(table(cltr.ext[cltr.ext != 0L])) > 1L) {
      cltr.ext.new <- cltr.ext
      cltr.ext.new[cltr.ext.new != 0L] <- clusters.reorder(cltr.ext.new[cltr.ext.new != 0L])
      if (!identical(cltr.ext.new, cltr.ext)) {
        cltr.ext <- as.integer(cltr.ext.new)
        utils::write.csv(
          data.frame(cluster_extended = as.integer(cltr.ext)),
          file = cltr.ext.file,
          row.names = TRUE
        )
      }
    }
  }

  ## ---- Visualization (delegated to plot3D.cltrs) --------------------------
  if (isTRUE(show.plot)) {
    if (!requireNamespace("rgl", quietly = TRUE)) {
      warning("Package 'rgl' is not installed; skipping 3D visualization.")
    } else {
      if (!is.null(cltr.ext) && length(table(cltr.ext)) > 1L) {
        plot3D.cltrs(
          X, cltr.ext,
          legend.title = taxon,
          radius = NA,
          show.cltr.labels = TRUE,
          sort.legend.labs.by.freq = TRUE
        )
      } else {
        plot3D.cltrs(X, cltr = NULL)
      }
    }
  }

  invisible(list(
    taxon = taxon,
    cltr.ext = as.integer(cltr.ext),
    cltr     = as.integer(cltr),
    X = X,
    pacmap.file   = pacmap.file,
    cltr.ext.file = cltr.ext.file
  ))
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
#' @export
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
