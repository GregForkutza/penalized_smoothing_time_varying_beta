extract_simulation_results <- function(simulator, conf.levels, aggregate = FALSE) {
  if (aggregate == FALSE) {
    extract_trajectory <- function(simulator, conf.levels) {
      
      trajectories <- lapply(conf.levels, function(conf.level) {
        mp_trajectory_sd(simulator, conf.int = TRUE, conf.level = conf.level) %>%
          select(-row, -col) %>%
          pivot_wider(names_from = matrix,  
                      values_from = c(value, sd, conf.low, conf.high),  
                      names_sep = "_")
      })
      
      # Take the first trajectory as the base and rename CI columns
      base_traj <- trajectories[[1]]
      first_conf_level_str <- as.character(conf.levels[1])  
      names(base_traj) <- names(base_traj) %>%
        gsub("conf\\.low_", paste0("conf", first_conf_level_str, "low_"), .) %>%
        gsub("conf\\.high_", paste0("conf", first_conf_level_str, "high_"), .)
      
      # Process additional confidence levels
      for (i in 2:length(trajectories)) {
        traj <- trajectories[[i]]
        conf_level_str <- as.character(conf.levels[i])
        new_names <- names(traj) %>%
          gsub("conf\\.low_", paste0("conf", conf_level_str, "low_"), .) %>%
          gsub("conf\\.high_", paste0("conf", conf_level_str, "high_"), .)
        names(traj) <- new_names
        
        base_traj <- base_traj %>%
          left_join(traj %>% select(time, starts_with("conf")), by = "time")
      }
      
      return(base_traj)
    }
    
    df <- extract_trajectory(simulator, conf.levels)
    return(df)
  } else {
    # Extract trajectory from the .phases = "during" of simulator.
    
    trajectory <- mp_trajectory_sd(simulator) %>%
      select(-row, -col) %>%
      pivot_wider(names_from = matrix,  
                  values_from = c(value, sd),
                  names_sep = "_")
    
    # Define indicator matrix to aggregate fitted values.
    B <- create_indicator_matrix(nrow(simulator$get$initial("X")), step_length)
    
    # Define function to aggregate standard deviation
    aggregate_sd <- function(col_name, B, tmp, V, step_length) {
      base_col_name <- gsub("^sd_", "", col_name)
      index <- with(tmp, which(matrix == base_col_name & !(time %in% c(0, max(time)))))
      index_cov <- V[index, index]
   
      # Calculate the aggregated standard deviation
      index_sd <- ((t(B) %*% index_cov %*% B) / step_length) |>
        diag() |>
        sqrt()
      
      
      return(index_sd)
    }
    
    # Define function to build data lists of aggregated values and sd. 
    aggregate_and_combine <- function(trajectory, B) {
      # Extract column names of value and standard deviation
      value_cols <- grep("^value_", names(trajectory), value = TRUE)
      sd_cols <- grep("^sd_", names(trajectory), value = TRUE)
      
      value_columns <- list()
      
      # Agrregate value columns 
      for (col_name in value_cols) {
        new_column <- NULL
        
        if (col_name == "value_infection") {
          wk_obs <- mp_final(simulator) %>%  
            filter(matrix == "wk_incidence") %>%
            pull(value)
          new_column <- as.matrix(wk_obs)
        } else {
          column_vector <- as.matrix(trajectory[[col_name]])
          #new_column <- (t(B) %*% column_vector) / step_length
          new_column <- (t(B) %*% column_vector) 
        }
        
        value_columns[[col_name]] <- as.vector(new_column)  
      }
      
      # Agrregate sd columns 
      sdr <- simulator$sdreport()
      V <- sdr$cov
      # Get trajectory output for all phases to use as meta data to extract correct elements of V.
      tmp <- simulator$report( .phases=c("before","during","after"), .sort = FALSE)
      sd_columns <- list()
      
      for (col_name in sd_cols) {
        new_column <- NULL
        column_vector <- as.matrix(trajectory[[col_name]])  
        if(col_name == "sd_infection") {
          new_column <- aggregate_sd(col_name, B, tmp, V, step_length = 1/step_length)  
        } else {
          new_column <- aggregate_sd(col_name, B, tmp, V, step_length = 1)
        }
        sd_columns[[col_name]] <- as.vector(new_column)  
      }
      
      
      df_value_columns <- as.data.frame(value_columns)
      df_sd_columns <- as.data.frame(sd_columns)
      aggregated_df <- cbind(df_value_columns, df_sd_columns)
      
      
      return(aggregated_df) 
    }
    df <- aggregate_and_combine(trajectory,B)
    compute_aggregated_CIs <- function(df, conf.levels) {
      value_cols <- grep("^value_", names(df), value = TRUE)
      sd_cols <- grep("^sd_", names(df), value = TRUE)
      
      output_df <- df[, c(value_cols, sd_cols)]
      for (conf.level in conf.levels) {
        z <- abs(qnorm((1-conf.level)/2))
        
        for (i in seq_along(value_cols)) {
          value_col <- value_cols[i]
          sd_col <- sd_cols[i]
          
          ci_lower <- df[[value_col]] - z * df[[sd_col]]
          ci_upper <- df[[value_col]] + z * df[[sd_col]]
          
          lower_col_name <- paste0(sub("^value_", "", value_col), "_", conf.level, "_lower")
          upper_col_name <- paste0(sub("^value_", "", value_col), "_", conf.level, "_upper")
          
          output_df[[lower_col_name]] <- ci_lower
          output_df[[upper_col_name]] <- ci_upper
        }
      }
      return(output_df)
    }
    
    df <- compute_aggregated_CIs(df, conf.levels)  
    
    return(df)
  }
}

