MSE_trajectory <- function(simulator, incidence, aggregate = FALSE) {
  if (aggregate) {
    # Call the function to get the dataframe
    mp_final_result <- mp_final(simulator)
    
    # Extract the vector where matrix is 'wk_incidence'
    wk_incidence_vector <- mp_final_result[mp_final_result$matrix == "wk_incidence", "value"]
    
    # Convert to a numeric vector
    predicted_incidence <- as.numeric(wk_incidence_vector)
  } else {
    # Get the trajectory and reshape it
    trajectory <- mp_trajectory_sd(simulator) %>%
      select(-row, -col) %>%
      pivot_wider(names_from = matrix, values_from = c(value, sd), names_sep = "_")
    
    # Extract predicted incidence values
    predicted_incidence <- trajectory$value_infection
  }
  
  # Compute the Mean Squared Error
  mse <- sum((incidence - predicted_incidence)^2) / length(incidence)
  
  return(mse)
}

# Function to compute the MSE table
compute_mse_table <- function(simulator_list, incidence, aggregate = FALSE) {
  basis_types <- names(simulator_list)
  results <- numeric(length(basis_types))
  names(results) <- basis_types
  
  for (i in seq_along(basis_types)) {
    # Run the MSE calculation
    mse <- MSE_trajectory(simulator_list[[basis_types[i]]], incidence, aggregate)
    results[i] <- mse
  }
  
  # Convert results to a data frame for kable
  results_df <- as.data.frame(results)
  
  # Rename the column to 'MSE'
  colnames(results_df) <- "MSE"
  
  # Display the table using kable
  return(results_df)
}



