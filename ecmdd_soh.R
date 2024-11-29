#!/usr/bin/env Rscript

#' @title Evidential c-Medoids Clustering with DTW Distance
#' @description Implementation of Evidential c-Medoids clustering algorithm with Dynamic Time Warping distance
#' @author Eya Laffet

# Load required packages
suppressPackageStartupMessages({
  library(evclust)
  library(cluster)
  library(dtw)
})

#' Create logging function with standardized levels
#' @param level Character string indicating log level
#' @param message Message to log
#' @param log_level Current logging threshold
#' @return None
create_logger <- function() {
  log_levels <- c("DEBUG" = 1, "INFO" = 2, "WARN" = 3, "ERROR" = 4)
  
  function(level, message, log_level = "INFO") {
    level <- toupper(level)
    log_level <- toupper(log_level)
    
    if (!level %in% names(log_levels)) {
      warning("Invalid log level. Defaulting to INFO")
      level <- "INFO"
    }
    
    if (log_levels[[level]] >= log_levels[[log_level]]) {
      timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
      cat(sprintf("[%s][%s] %s\n", timestamp, level, message))
    }
  }
}

#' Initialize unique medoids for clustering
#' @param n_samples Number of samples
#' @param n_clusters Number of clusters
#' @param dist_mat Distance matrix
#' @param logger Logging function
#' @return Vector of unique medoid indices
initialize_medoids <- function(n_samples, n_clusters, dist_mat, logger) {
  # Use max-min initialization strategy to ensure diverse medoids
  medoids <- numeric(n_clusters)
  
  # Select first medoid randomly
  medoids[1] <- sample.int(n_samples, 1)
  
  # Select remaining medoids using max-min distance
  for (i in 2:n_clusters) {
    # Calculate minimum distances to existing medoids
    min_distances <- apply(dist_mat[, medoids[1:(i-1)], drop = FALSE], 1, min)
    
    # Select point with maximum minimum distance
    medoids[i] <- which.max(min_distances)
  }
  
  logger("DEBUG", sprintf("Initialized medoids: %s", paste(medoids, collapse=", ")))
  return(medoids)
}

#' Find medoids for singleton clusters
#' @param x Input data matrix
#' @param m Membership matrix
#' @param dist_mat Distance matrix
#' @param focal_elements Focal elements matrix
#' @param log_level Logging level
#' @return List containing medoid indices and centroids
find_singleton_medoids <- function(x, m, dist_mat, focal_elements, log_level = "INFO") {
  logger <- create_logger()
  
  # Input validation
  if (!is.matrix(x) || !is.matrix(m) || !is.matrix(dist_mat) || !is.matrix(focal_elements)) {
    stop("All inputs must be matrices")
  }
  
  logger("DEBUG", "Starting singleton medoids identification", log_level)
  
  # Identify singleton classes
  singleton_indices <- which(rowSums(focal_elements[-1, ]) == 1)
  logger("DEBUG", sprintf("Found %d singleton classes", length(singleton_indices)), log_level)
  
  # Filter membership matrix for singleton classes
  m_singleton <- m[, singleton_indices, drop = FALSE]
  
  # Compute weighted distances
  q <- dist_mat %*% m_singleton
  
  # Find medoids ensuring uniqueness
  medoid_indices <- numeric(ncol(q))
  used_indices <- logical(nrow(x))
  
  for (j in 1:ncol(q)) {
    # Get available indices (not yet used as medoids)
    available_indices <- which(!used_indices)
    
    if (length(available_indices) == 0) {
      logger("WARN", "No more unique indices available for medoids", log_level)
      break
    }
    
    # Find best available medoid
    best_idx <- available_indices[which.min(q[available_indices, j])]
    medoid_indices[j] <- best_idx
    used_indices[best_idx] <- TRUE
  }
  
  # Extract centroids
  centroids <- x[medoid_indices, , drop = FALSE]
  
  logger("DEBUG", sprintf("Selected medoids: %s", paste(medoid_indices, collapse=", ")), log_level)
  
  return(list(
    medoid_indices = medoid_indices,
    centroids = centroids
  ))
}

#' Evidential c-Medoids clustering algorithm
#' @param x Input data
#' @param c Number of clusters
#' @param type Type of focal elements ('full' or 'singleton')
#' @param dist_mat Distance matrix
#' @param alpha Mass function parameter
#' @param beta Fuzzifier parameter
#' @param delta Outlier detection parameter
#' @param epsilon Convergence threshold
#' @param display Show iteration progress
#' @param gamma Distance weighting parameter
#' @param eta Prototype selection parameter
#' @param max_trials Number of random initializations
#' @param max_iterations Maximum iterations per trial
#' @param log_level Logging level
#' @return Clustering results object
evidential_cmedoids <- function(x, c, type = 'full', dist_mat = NULL, 
                              alpha = 1, beta = 1.5, delta = 9, epsilon = 1e-3,
                              display = FALSE, gamma = 0.5, eta = 1, 
                              max_trials = 10, max_iterations = 100,
                              log_level = "INFO") {
  
  logger <- create_logger()
  
  # Input validation
  if (!is.matrix(x)) x <- as.matrix(x)
  if (is.null(dist_mat)) stop("Distance matrix is required")
  if (c < 2) stop("Number of clusters must be at least 2")
  
  # Initialize algorithm parameters
  delta_squared <- delta^2
  focal_elements <- makeF(c = c, type = type)
  n_focal <- nrow(focal_elements)
  cardinalities <- rowSums(focal_elements[2:n_focal, ])
  n_samples <- nrow(x)
  
  # Initialize best solution trackers
  best_criterion <- Inf
  best_memberships <- NULL
  best_medoids <- NULL
  
  logger("INFO", sprintf("Starting clustering with %d samples into %d clusters", n_samples, c), log_level)
  
  # Main clustering loop
  for (trial in 1:max_trials) {
    # Initialize medoids using max-min strategy
    medoids <- initialize_medoids(n_samples, c, dist_mat, logger)
    prototypes <- numeric(n_focal - 1)
    
    criterion_old <- Inf
    converged <- FALSE
    iteration <- 0
    
    while (!converged && iteration < max_iterations) {
      iteration <- iteration + 1
      
      # Update prototypes
      for (j in 1:(n_focal - 1)) {
        focal_j <- focal_elements[j + 1, ]
        medoids_j <- medoids[which(focal_j != 0)]
        
        if (sum(focal_j) == 1) {
          prototypes[j] <- medoids_j
        } else {
          distances <- numeric(n_samples)
          for (i in 1:n_samples) {
            dists_to_medoids <- dist_mat[i, medoids_j]
            mean_dist <- mean(dists_to_medoids)
            variance <- (1 / cardinalities[j]) * sum((dists_to_medoids - mean_dist)^2)
            distances[i] <- variance + eta * mean_dist
          }
          prototypes[j] <- which.min(distances)
        }
      }
      
      # Calculate distances to prototypes
      D <- matrix(0, n_samples, n_focal - 1)
      for (j in 1:(n_focal - 1)) {
        focal_j <- focal_elements[j + 1, ]
        medoids_j <- medoids[which(focal_j != 0)]
        for (i in 1:n_samples) {
          dists <- dist_mat[i, medoids_j]
          D[i, j] <- max((dist_mat[i, prototypes[j]] + (gamma/cardinalities[j]) * sum(dists)) / (1 + gamma), 1e-10)
        }
      }
      
      # Update memberships
      m <- matrix(0, n_samples, n_focal - 1)
      for (i in 1:n_samples) {
        for (j in 1:(n_focal - 1)) {
          denominator <- sum((D[i, j]/D[i, ])^(1/(beta - 1)) * 
                           (cardinalities[j]/cardinalities)^(alpha/(beta - 1))) +
                           (cardinalities[j]^alpha * D[i, j]/delta_squared)^(1/(beta - 1))
          m[i, j] <- 1 / max(denominator, 1e-10)
        }
      }
      
      # Update medoids ensuring uniqueness
      medoid_result <- find_singleton_medoids(x, m^beta, dist_mat, focal_elements, log_level)
      medoids <- medoid_result$medoid_indices
      
      # Calculate criterion
      empty_mass <- pmax(0, 1 - rowSums(m))
      criterion <- sum((m^beta) * D * matrix(cardinalities^alpha, n_samples, n_focal - 1, byrow = TRUE)) +
                  delta_squared * sum(empty_mass^beta)
      
      # Check convergence
      converged <- abs(criterion - criterion_old) <= epsilon
      criterion_old <- criterion
      
      if (display) {
        logger("INFO", sprintf("Trial %d, Iteration %d: Criterion = %.6f", 
                             trial, iteration, criterion), log_level)
      }
    }
    
    # Update best solution
    if (criterion < best_criterion) {
      best_criterion <- criterion
      best_memberships <- m
      best_medoids <- medoids
      logger("INFO", sprintf("New best solution found in trial %d with criterion %.6f", 
                           trial, best_criterion), log_level)
    }
  }
  
  # Prepare final results
  if (is.null(best_memberships)) {
    stop("No valid solution found")
  }
  
  final_memberships <- cbind(1 - rowSums(best_memberships), best_memberships)
  
  # Create result object
  result <- extractMass(final_memberships, focal_elements, 
                       g = best_medoids, 
                       method = "ecmdd",
                       crit = best_criterion,
                       param = list(alpha = alpha, beta = beta, delta = delta))
  
  return(result)
}

#' Process CSV data for time series clustering
#' @param file_path Path to CSV file
#' @return List of numeric vectors
process_csv_data <- function(file_path) {
  if (!file.exists(file_path)) {
    stop("CSV file not found")
  }
  
  # Read CSV without headers, skip first row
  data <- read.csv(file_path, header = FALSE, skip = 1)
  
  # Convert columns to numeric vectors and remove NAs
  column_lists <- lapply(data, function(x) {
    clean_data <- as.numeric(na.omit(x))
    if (length(clean_data) == 0) {
      stop("No valid numeric data found in column")
    }
    clean_data
  })
  
  unname(column_lists)
}

#' Create DTW distance matrix
#' @param series List of time series
#' @param window_size DTW window size (0 for no window)
#' @param trace Enable progress tracking
#' @return Distance matrix
create_dtw_matrix <- function(series, window_size = 0, trace = FALSE) {
  if (!is.list(series) || length(series) < 2) {
    stop("Input must be a list of at least two time series")
  }
  
  n_series <- length(series)
  dist_matrix <- matrix(0, nrow = n_series, ncol = n_series)
  
  if (trace) message("Computing DTW distance matrix...")
  
  for (i in 1:(n_series - 1)) {
    for (j in (i + 1):n_series) {
      if (trace) message(sprintf("Computing distance between series %d and %d", i, j))
      
      # Compute DTW distance
      if (window_size == 0) {
        dist <- dtw::dtw(series[[i]], series[[j]], 
                        step.pattern = dtw::symmetric2)$normalizedDistance
      } else {
        effective_window <- min(window_size, 
                              abs(length(series[[i]]) - length(series[[j]])))
        dist <- dtw::dtw(series[[i]], series[[j]],
                        step.pattern = dtw::symmetric2,
                        window.type = "sakoechiba",
                        window.size = effective_window)$normalizedDistance
      }
      
      dist_matrix[i, j] <- dist
      dist_matrix[j, i] <- dist
    }
  }
  
  if (trace) message("Distance matrix computation complete")
  return(dist_matrix)
}

# Example usage
main <- function() {
  # Process data
  column_lists <- process_csv_data('SOH.csv')
  print("Data dimensions:")
  print(sapply(column_lists, length))
  
  # Create distance matrix
  dist_matrix <- create_dtw_matrix(column_lists, trace = TRUE)
  
  # Perform clustering
  results <- evidential_cmedoids(
    x = column_lists,
    c = 3,
    dist_mat = dist_matrix,
    alpha = 1.5,
    beta = 1.5,
    max_trials = 100
  )
  
  # Calculate silhouette score
  sil <- silhouette(x = results$y.bel, dmatrix = dist_matrix)
  
  # Compute metrics
  metrics <- list(
    average_silhouette = mean(sil[, 3]),
    per_cluster_silhouette = tapply(sil[, 3], sil[, 1], mean)
  )
  
  # Print results
  print("Clustering Results:")
  print(results$y.bel)
  print("\nSilhouette Analysis:")
  print(metrics)
}

if (sys.nframe() == 0) {
  main()
}
