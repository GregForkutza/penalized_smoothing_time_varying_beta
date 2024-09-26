create_simulated_data <- function(sir_spec, smooth_parms_initial, compartmental_parms_initial, step_length, sim_sd) {
  
  # Convert lists to environments
  e_smooth_parms <- list2env(smooth_parms_initial, envir = environment())
  e_compartmental_parms <- list2env(compartmental_parms_initial, envir = environment())
  
  # Manually copy each parameter into the function's environment
  for (name in names(smooth_parms_initial)) {
    assign(name, e_smooth_parms[[name]], envir = environment())
  }
  for (name in names(compartmental_parms_initial)) {
    assign(name, e_compartmental_parms[[name]], envir = environment())
  }

  # Build model to simulate data
  simulator_initial <- construct_simulator(sir_spec, smooth_parms_initial, compartmental_parms_initial, step_length = step_length, eq = TRUE)
  
  # Simulate data and extract incidence trajectory
  sim_data <- generate_simulated_data(simulator_initial, sim_sd, step_length)
  incidence <- sim_data$P_noise
  
  results_initial <- extract_simulation_results(simulator = simulator_initial, conf.levels = c(0.5, 0.95), aggregate = FALSE)
  
  return(list(simulator_initial = simulator_initial, sim_data = sim_data, incidence = incidence))
}

calibrate_multiple_smoothers <- function(smooth_list, smooth_parms_fit, compartmental_parms_fit, model_parms, sir_spec, step_length, eq, observed_vector, aggregate = FALSE) {
  
  # Convert lists to environments
  e_smooth_parms_fit <- list2env(smooth_parms_fit, envir = environment())
  e_compartmental_parms_fit <- list2env(compartmental_parms_fit, envir = environment())
  
  # Manually copy each parameter into the function's environment
  for (name in names(smooth_parms_fit)) {
    assign(name, e_smooth_parms_fit[[name]], envir = environment())
  }
  for (name in names(compartmental_parms_fit)) {
    assign(name, e_compartmental_parms_fit[[name]], envir = environment())
  }
  
  # Store calibrated simulator objects for each smooth type in a nested list
  calibrated_simulators <- list()
  
  # Iterate over each smooth type in the provided list
  for (smooth_type in smooth_list) {
    
    # Adjust smooth model for current smooth type
    smooth_parms_fit <- list(
      smooth = smooth_type,
      num_variables = num_variables,
      data = data,
      cov_fun = cov_fun
    )

    # Reconstruct the simulator using the new smooth type
    # uncomment out the line below if you are using construct_simulator_agg instead of constuct_simulator_sim
   #simulator_fit <- construct_simulator(sir_spec, smooth_parms_fit, compartmental_parms_fit, step_length = step_length, eq = eq, aggregate = aggregate)
    simulator_fit <- construct_simulator(sir_spec, smooth_parms_fit, compartmental_parms_fit, step_length = step_length, eq = eq)
    
    # Update the simulator with the observed data
    update_likelihood_expressions(
      simulator = simulator_fit,
      observed_vector = observed_vector, 
      step_length = step_length,
      model_params = model_params
    )
    # Try to calibrate the model, skip to next smooth type if an error occurs
    fit <- tryCatch({
      simulator_fit$optimize$nlminb()
      calibrated_simulators[[smooth_type]] <- simulator_fit
    }, error = function(e) {
      message("An error occurred during optimization for smooth type ", smooth_type, ": ", e$message)
      NULL  # Return NULL to indicate failure and continue
    })
  }
  
  return(calibrated_simulators)
}

