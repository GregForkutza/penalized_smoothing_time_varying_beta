generate_simulated_data <- function(simulator, sim_sd = 1, step_length) {
  if(step_length != 1) {
    obs_mat <- simulator$report(.phases= "during")
    
    D <- get_parms_groups_sums(data, step_size = step_length)
    
    P_obs <- obs_mat %>%
      filter(matrix == "infection") %>%
      pull(value)
    
    wk_inc_obs <- tapply(P_obs, D$wk, sum) %>% as.vector()
    
    # Assuming wk_inc_obs is a vector and sim_sd is defined
    #wk_inc_obs_noisy <- wk_inc_obs * exp(rnorm(length(wk_inc_obs), mean = 0, sd = sim_sd))
    wk_inc_obs_noisy <- rnorm(n = length(wk_inc_obs), mean = wk_inc_obs, sd = sim_sd)
    
    
    simulator$default[["P_true"]] <- wk_inc_obs 
    simulator$default[["P_noise"]] <- wk_inc_obs_noisy
    return(list(P_true = wk_inc_obs , P_noise =  wk_inc_obs_noisy))
    
  } else {
    obs_mat <- simulator$report(.phases= "during")
    # Extract true values
    obs_I <- obs_mat |> 
      filter(matrix == "I") |> 
      pull(value)
    obs_P <- obs_mat |> 
      filter(matrix == "infection") |> 
      pull(value)
    
    # Generate values with added noise
    I_obs <- obs_mat |> 
      filter(matrix == "I") |> 
      mutate(across(value, ~ rnorm(n(), ., sd = sim_sd))) |> 
      pull(value)

    #P_obs <- obs_mat |> 
     # filter(matrix == "infection") |> 
     # mutate(across(value, ~ rnorm(n(), ., sd = sim_sd))) |> 
      #pull(value)
    
    P_obs <- obs_mat %>%
      filter(matrix == "infection") %>%
      mutate(
        noisy_value = value * exp(rnorm(n(), mean = 0, sd = sim_sd))  # Apply multiplicative noise
      ) %>%
      pull(noisy_value)

    
    # Store the generated data in the simulator object
    #simulator$default[["I_true"]] <- obs_I
    #simulator$default[["I_noise"]] <- I_obs
    #simulator$default[["P_true"]] <- obs_P
    #simulator$default[["P_noise"]] <- P_obs
  
  
  # Optionally, return the generated data if needed outside of the simulator object
  return(list(I_true = obs_I, I_noise = I_obs, P_true = obs_P, P_noise = P_obs))
  }
}