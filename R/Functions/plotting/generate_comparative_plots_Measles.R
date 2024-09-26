library(grid)
library(gridExtra)
library(ggplot2)

create_plots <- function(results, observed_vector,smooth_type, simulator) {
  # Choose the correct x-axis label based on aggregation
  x_label <-  "Day"
  
  # Determine the Date variable based on aggregation
  Date <- results$time
  # Case 3: aggregate = FALSE, results_true is not NULL
  data_combined_incidence <- data.frame(
    Date = Date,
    Incidence = observed_vector,
    I_hat = results$value_infection,
    I_hat_50_lower = results$conf0.5low_infection,
    I_hat_50_upper = results$conf0.5high_infection,
    I_hat_95_lower = results$conf0.95low_infection,
    I_hat_95_upper = results$conf0.95high_infection
  )
  
  p1 <- ggplot(data_combined_incidence, aes(x = Date)) +
    geom_line(aes(y = Incidence), colour = "red") +
    geom_line(aes(y = I_hat), colour = "blue") +
    geom_ribbon(aes(ymin = I_hat_50_lower, ymax = I_hat_50_upper), alpha = 0.4, fill = "lightblue") +
    geom_ribbon(aes(ymin = I_hat_95_lower, ymax = I_hat_95_upper), alpha = 0.2, fill = "blue") +
    labs(x = x_label, y = "Incidence") +
    ggtitle(smooth_type) +
    theme_minimal()
  
  data_combined_beta <- data.frame(
    Date = Date,
    Beta_hat = exp(results$value_theta),
    Beta_hat_50_lower = exp(results$conf0.5low_theta),
    Beta_hat_50_upper = exp(results$conf0.5high_theta),
    Beta_hat_95_lower = exp(results$conf0.95low_theta),
    Beta_hat_95_upper = exp(results$conf0.95high_theta)
  )
  
  p2 <- ggplot(data_combined_beta, aes(x = Date)) +
    geom_line(aes(y = Beta_hat), color = "green") +
    geom_ribbon(aes(ymin = Beta_hat_50_lower, ymax = Beta_hat_50_upper), alpha = 0.4, fill = "lightgreen") +
    geom_ribbon(aes(ymin = Beta_hat_95_lower, ymax = Beta_hat_95_upper), alpha = 0.2, fill = "green") +
    labs(x = x_label, y = expression(beta)) +
    ggtitle(smooth_type) +
    theme_minimal()
  
  
  data_combined_R_t <- data.frame(
    Date = Date,
    R_t = exp(results$value_R_t),
    R_t_50_lower = exp(results$conf0.5low_R_t),
    R_t_50_upper = exp(results$conf0.5high_R_t),
    R_t_95_lower = exp(results$conf0.95low_R_t),
    R_t_95_upper = exp(results$conf0.95high_R_t)
  )
  
  p3 <- ggplot(data_combined_R_t, aes(x = Date)) +
    geom_line(aes(y = R_t), color = "purple") +
    geom_ribbon(aes(ymin = R_t_50_lower, ymax = R_t_50_upper), alpha = 0.4, fill = "#D7B9D5") +
    geom_ribbon(aes(ymin = R_t_95_lower, ymax = R_t_95_upper), alpha = 0.2, fill = "purple") +
    labs(x = x_label, y = (expression(R_t))) +
    ggtitle(smooth_type) +
    theme_minimal()
  
  return(list(plot_incidence = p1, plot_beta = p2, plot_R_t = p3))
}


generate_comparative_plots <- function(simulators, aggregate = FALSE, conf.levels = c(0.5, 0.95), data_name) {
  all_plots <- list()  # Initialize empty list to store all plots
  
  # Identify the calibrated model types for this initial smooth
  calibrated_model_types <- names(simulators)
  
  # Iterate over each calibrated model type
  for (smooth_type in calibrated_model_types) {
    #Extract results for calibrated model
    results_fitted <- extract_simulation_results(simulator = simulators[[smooth_type]], conf.levels = c(0.5, 0.95), aggregate = FALSE)
    plots <- create_plots(results = results_fitted, observed_vector = incidence, smooth_type, simulators[[smooth_type]])
    plot_list_name <- paste(data_name, smooth_type, "plots", sep = "_")
    all_plots[[plot_list_name]] <- plots
  }
  return(all_plots)
}


all_plots <- function(all_comparative_plots, output_dir = NULL, file_name = "comparative_plot") {
  
  # Number of lists of plots passed in
  num_plot_lists <- length(all_comparative_plots)
  
  # Identifiers for each of the three plots
  plot_identifiers <- c("incidence", "transmission", "reproduction")
  
  # Create an empty list to store the final combined plots for each smooth type
  final_combined_plots <- list()
  
  # Iterate over each smooth type (assuming one smooth type in each list)
  for (smooth_type_index in 1:length(all_comparative_plots[[1]])) {
    
    # Initialize an empty list for storing the combined plots for this smooth type
    smooth_type_plots <- list()
    
    # Iterate over each variable
    for (variable_index in 1:length(plot_identifiers)) {
      
      # Initialize an empty list for storing the plots for this row
      row_plots <- list()
      
      # Iterate over each list of plots
      for (plot_list_index in 1:num_plot_lists) {
        
        # Extract the plot for the current variable and smooth type
        row_plots[[paste("Plot", plot_list_index)]] <- all_comparative_plots[[plot_list_index]][[smooth_type_index]][[variable_index]]
      }
      
      # Combine the row plots into a single row
      smooth_type_plots[[paste("Row", variable_index)]] <- do.call(arrangeGrob, c(row_plots, ncol = num_plot_lists))
    }
    
    # Combine all rows into a single pane for the current smooth type
    final_combined_plots[[paste("SmoothType", smooth_type_index)]] <- arrangeGrob(grobs = smooth_type_plots, ncol = 1)
  }
  
  # Display all combined plots
  for (smooth_type in names(final_combined_plots)) {
    grid.newpage()
    grid.draw(final_combined_plots[[smooth_type]])
    
    # Save the combined plot as a .png file if output_dir is provided
    if (!is.null(output_dir)) {
      file_path <- file.path(output_dir, paste0(file_name, "_", smooth_type, ".png"))
      ggsave(filename = file_path, plot = final_combined_plots[[smooth_type]], width = 15, height = 12)
    }
  }
  
  return(final_combined_plots)
}
