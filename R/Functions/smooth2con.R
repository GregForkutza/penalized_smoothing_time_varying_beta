construct_tridiagonal_matrix <- function(alpha, time_steps) {
  
  matrix <- matrix(0, nrow = time_steps, ncol = time_steps)
  diag(matrix) <- 1
  matrix[row(matrix) == col(matrix) + 1] <- -alpha
  matrix[row(matrix) == col(matrix) - 1] <- -alpha
  
  return(matrix)
}


construct_alpha_matrix <- function(alpha, time_steps) {
  # Create an n x n matrix filled with 0's
  matrix <- matrix(0, nrow = time_steps, ncol = time_steps)
  
  # Populate the matrix
  for (row in 1:time_steps) {
    for (col in 1:time_steps) {
      if (row >= col) {
        # Set the value based on the difference between row and col
        matrix[row, col] <- alpha^(row - col)
      }
    }
  }
  return(matrix)
}

scale_matrices <- function(smooth_output) {
  X <- smooth_output$X
  P <- smooth_output$P
  
  # Calculate the sum of squares for each column in X
  m <- colSums(X^2)
  
  # Scale X by the column-wise sum of squares
  X_scaled <- X %*% diag(1/sqrt(m))
  
  # Scale P upward
  P_scaled <- diag(sqrt(m)) %*% P %*% diag(sqrt(m))
  
  return(list(X = X_scaled, P = P_scaled))
}


#' Construct Smoothing Matrices for Epidemiological Models
#'
#' This function constructs smoothing matrices for epidemiological models,
#' supporting various smoothing methods including thin plate splines, cubic regression
#' splines, B-splines, P-splines, adaptive smoothing, Gaussian processes, and an
#' Ornstein-Uhlenbeck-like process. It leverages `mgcv::smoothCon` for most smoothing
#' types but also provides custom handling for an Ornstein-Uhlenbeck-like process.
#' Both scaled and unscaled versions of the basis and penalty matrices are returned.
#'
#' @param smooth A character string specifying the type of smooth to be constructed.
#'        Supported types include "tp" for thin plate splines, "cr" for cubic regression
#'        splines, "bs" for B-splines, "ps" for P-splines, "ad" for adaptive smoothing,
#'        "gp" for Gaussian processes, and "ou" for an Ornstein-Uhlenbeck-like process.
#'        For "gp", `cov_fun` must be provided. For "ou", `alpha` must be specified and
#'        must be between 0 and 1.
#' @param num_variables The number of variables (basis functions) to be used in the
#'        smoothing process, adjusted for identifiability constraints as necessary.
#' @param data A data frame containing the time variable for which smoothing is to be
#'        applied. The data frame should have one column named 'time'.
#' @param cov_fun (Optional) A covariance function specification for the Gaussian process
#'        smoothing. Required if `smooth` is "gp".
#' @param alpha (Optional) A decay parameter used specifically for the "ou" smoothing type
#'        to control the rate of exponential decay. Required and must be between 0 and 1 if
#'        `smooth` is "ou".
#'
#' @return A list containing the following elements:
#'         - `X`: The scaled basis matrix for the smoothing.
#'         - `P`: The scaled penalty matrix associated with the smoothing method.
#'         - `C`: A matrix of constraints applied to the basis functions for identifiability.
#'         - `X_unscaled`: The original unscaled basis matrix.
#'         - `P_unscaled`: The original unscaled penalty matrix.
#'         
#'
#' @examples
#' data <- data.frame(time = 1:100)
#' smooth_params <- smooth2con(smooth = "tp", num_variables = 10, data = data)
#' smooth_params_gp <- smooth2con(smooth = "gp", num_variables = 10, data = data, cov_fun = 2)
#' smooth_params_ou <- smooth2con(smooth = "ou", num_variables = 10, data = data, alpha = 0.5)
#'
#' @export
#'
#' @importFrom mgcv smoothCon
smooth2con <- function(smooth, num_variables, data, cov_fun = NULL, alpha = NULL, absorb.cons = TRUE) {
  dd <- data.frame(time = 1:nrow(data))
  n <- num_variables - 1  # Adjusted for identifiability constraints
  # Validate inputs
  if (smooth == "gp" && is.null(cov_fun)) {
    stop("cov_fun must be provided for 'gp' smooth type")
  }
  if (smooth %in% c("tp", "cr", "bs", "ps", "ts", "cc", "cp", "ds")) {
    # Standard smoothing methods
    s <- smoothCon(object = s(time, bs = smooth, k = num_variables),
                   absorb.cons = absorb.cons, data = dd, knots = NULL)
  } else if (smooth == "ad") {
    s <- smoothCon(object = s(time, bs = smooth, k = num_variables, m = 2, xt=list(bs="cr")),
                            absorb.cons = absorb.cons, data = dd, knots = NULL)
  } else if (smooth == "gp") {
    # Gaussian process smoothing
    s <- smoothCon(object = s(time, bs = smooth, k = num_variables, m = cov_fun),
                   absorb.cons = absorb.cons, data = dd, knots = NULL)
  } else {
    stop("Invalid smooth type specified")
  }
  # For standard smoothing methods and Gaussian processes
  X_unscaled <- s[[1]]$X
  P_unscaled <- s[[1]]$S[[1]]
  C <- s[[1]]$C  
  # Apply scaling
  scaled_matrices <- scale_matrices(list(X = X_unscaled, P = P_unscaled))
    epsilon <- 10^-10
    P_reg <- scaled_matrices$P + epsilon
    P_logdet <- log(det(P_reg))
  # Return both scaled and unscaled matrices
  return(list(
    X = scaled_matrices$X,
    P = scaled_matrices$P, 
    C = C,
    X_unscaled = X_unscaled,
    P_unscaled = P_unscaled,
    P_logdet = P_logdet
  )
  )
}