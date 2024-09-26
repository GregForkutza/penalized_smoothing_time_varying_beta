library(knitr)
library(kableExtra)
compute_aic <- function(simulator) {
  # Compute trace of hat matrix for smoother
  X <- simulator$get$initial("X")
  lambda <- mp_tmb_coef(simulator) %>%
    filter(mat == "log_smooth_sd") %>%
    select(estimate) %>%
    pull() %>%
    exp()
  
  P <- simulator$get$initial("P")
  #A <- X %*% solve(t(X) %*% X + lambda * P) %*% t(X)
  A <- X %*% solve(t(X) %*% X + lambda * P) %*% t(X) %*%  X
  k <- sum(diag(A)) # estimated degrees of freedom 
  # Extract maximized value of the likelihood function for the model 
  L <- simulator$optimization_history$latest()$objective
  
  # Compute AIC
  aic <- 2 * k - 2 * log(L)
  print(log(L))
  print(lambda)
  print(k)
  print(sum(diag(X %*% solve(t(X) %*% X +  P) %*% t(X))))
  print(dim(X))
  print(dim(P))
  return(aic)
}

create_aic_table <- function(simulators, knots, save_to_file = FALSE, output_file = NULL) {
  if (length(knots) != length(simulators)) {
    stop("The length of the 'knots' vector must match the number of simulator sets.")
  }
  
  aic_values <- list()
  for (i in seq_along(simulators)) {
    simulator_set <- simulators[[i]]
    knot_value <- knots[i]
    for (simulator_name in names(simulator_set)) {
      simulator <- simulator_set[[simulator_name]]
 
      aic_value <- compute_aic(simulator)
      if (!simulator_name %in% names(aic_values)) {
        aic_values[[simulator_name]] <- list()
      }
      aic_values[[simulator_name]][[paste0("k = ", knot_value)]] <- aic_value
    }
  }
  
  aic_table <- do.call(rbind, lapply(aic_values, function(x) {
    as.data.frame(t(x))
  }))
  
  aic_table <- cbind(SmoothType = rownames(aic_table), aic_table)
  rownames(aic_table) <- NULL
  
  if (save_to_file && !is.null(output_file)) {
    saveRDS(aic_table, file = paste0(output_file, ".rds"))
  }
  
  return(aic_table)
}












