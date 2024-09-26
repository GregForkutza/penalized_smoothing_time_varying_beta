library(grid)
library(gridExtra)
library(ggplot2)
library(tidyr)

create_plots <- function(results, results_true, observed_vector, smooth_type, simulator) {

  # Choose the correct x-axis label based on aggregation
  Date <- seq(1:length(results$value_infection))
  B <- create_indicator_matrix(nrow(simulator$get$initial("X")), step_length)
  Beta_hat_true = exp((t(B) %*% results_true$value_theta))[,1]
  R_t_true = exp((t(B) %*% results_true$value_R_t))[,1]
  
  # Data for Incidence Plot
  data_combined_incidence <- data.frame(
    Date = Date,
    Incidence = observed_vector,
    I_hat = results$value_infection,
    I_hat_50_lower = results$infection_0.5_lower,
    I_hat_50_upper = results$infection_0.5_upper,
    I_hat_95_lower = results$infection_0.95_lower,
    I_hat_95_upper = results$infection_0.95_upper
  )
  
  p1 <- ggplot(data_combined_incidence, aes(x = Date)) +
    geom_line(aes(y = Incidence), colour = "red") +
    geom_line(aes(y = I_hat), colour = "blue") +
    geom_ribbon(aes(ymin = I_hat_50_lower, ymax = I_hat_50_upper), alpha = 0.4, fill = "lightblue") +
    geom_ribbon(aes(ymin = I_hat_95_lower, ymax = I_hat_95_upper), alpha = 0.2, fill = "blue") +
    labs(x = NULL, y = NULL) +  # Remove axis labels
    theme_minimal() 
  
  # Plot for R_t with 50% and 95% CIs
  data_combined_beta <- data.frame(
    Date = Date,
    Beta_hat = exp(results$value_theta),
    Beta_hat_true = Beta_hat_true,
    Beta_50_lower = exp(results$theta_0.5_lower),
    Beta_50_upper = exp(results$theta_0.5_upper),
    Beta_95_lower = exp(results$theta_0.95_lower),
    Beta_95_upper = exp(results$theta_0.95_upper)
  )
  
  p2 <- ggplot(data_combined_beta, aes(x = Date)) +
    geom_line(aes(y = Beta_hat_true), color = "pink") +
    geom_line(aes(y = Beta_hat), color = "green") +
    geom_ribbon(aes(ymin = Beta_50_lower, ymax = Beta_50_upper), alpha = 0.4, fill = "lightgreen") +
    geom_ribbon(aes(ymin = Beta_95_lower, ymax = Beta_95_upper), alpha = 0.2, fill = "green") +
    labs(x = NULL, y = NULL) +  # Remove axis labels
    theme_minimal() 
  
  # Data for R_t Plot
  data_combined_R_t <- data.frame(
    Date = Date,
    R_t = exp(results$value_R_t),
    R_t_true = R_t_true,
    R_t_50_lower = exp(results$R_t_0.5_lower),
    R_t_50_upper = exp(results$R_t_0.5_upper),
    R_t_95_lower = exp(results$R_t_0.95_lower),
    R_t_95_upper = exp(results$R_t_0.95_upper)
  )
  
  p3 <- ggplot(data_combined_R_t, aes(x = Date)) +
    geom_line(aes(y = R_t), color = "purple") +
    geom_line(aes(y = R_t_true), color = "black") +
    geom_ribbon(aes(ymin = R_t_50_lower, ymax = R_t_50_upper), alpha = 0.4, fill = "#D7B9D5") +
    geom_ribbon(aes(ymin = R_t_95_lower, ymax = R_t_95_upper), alpha = 0.2, fill = "purple") +
    labs(x = NULL, y = NULL) +  # Remove axis labels
    theme_minimal() 
  
  # Data for Basis Plot
  plot4_data <- as.data.frame(simulator$default$X)
  plot4_data$Time <- seq(1:nrow(simulator$default$X))
  plot4_data_long <- gather(plot4_data, key = "Basis", value = "Value", -Time)
  p4 <- ggplot(plot4_data_long, aes(x = Time, y = Value, color = Basis)) +
    labs(x = NULL, y = NULL) + 
    geom_line() +
    theme(legend.position = "none") 
  
  return(list(plot_incidence = p1, plot_beta = p2, plot_R_t = p3, plot_basis = p4))
}

generate_comparative_plots <- function(simulators, aggregate = FALSE, conf.levels = c(0.5, 0.95), knots_vector) {

  all_plots <- list()  # Initialize empty list to store all plots
  calibrated_model_types <- names(simulators$calibrated_simulator)
  #sim_initial <- simulators$simulator_initial[[initial_model_type]]

  #results <- extract_simulation_results(simulator = sim_initial, conf.levels = conf.levels, aggregate = FALSE)
  results <- simulators$initial_results
  for (smooth_type in calibrated_model_types) {
 
    results_fitted <- extract_simulation_results(simulator = simulators$calibrated_simulator[[smooth_type]], conf.levels = c(0.5, 0.95), aggregate = aggregate)
    plots <- create_plots(results = results_fitted, results_true = results, observed_vector = simulators$incidence, smooth_type, simulators$calibrated_simulator[[smooth_type]])
    plot_list_name <- paste( initial_smooth,smooth_type, "plots", sep = "_")
    all_plots[[plot_list_name]] <- plots
  }
  return(all_plots)
}

extract_smooth_types <- function(all_comparative_plots) {

  all_names <- unlist(lapply(all_comparative_plots, names))
  smooth_types <- unique(sub("^[^_]+_([^_]+)_.*$", "\\1", all_names))
  
  return(smooth_types)
}

all_plots <- function(all_comparative_plots, output_dir = NULL, file_name = "comparative_plot", knots_vector, x_label = "Day") {
  
  smooth_types <- extract_smooth_types(all_comparative_plots)
  num_smooth_types <- length(smooth_types)
  num_plot_lists <- length(all_comparative_plots)
  plot_identifiers <- c("plot_incidence", "plot_beta", "plot_R_t", "plot_basis")
  
  y_labels <- list(
    plot_incidence = "Incidence",
    plot_beta = expression(beta[t]),
    plot_R_t = expression(R[t])
  )
  
  final_combined_plots <- list()
  
  for (i in seq_along(plot_identifiers)) {
    variable_plots <- list()
    
    for (j in seq_along(smooth_types)) {
      smooth_type_name <- smooth_types[j]
      row_plots <- list()
      
      for (k in seq_along(all_comparative_plots)) {
        smooth_type_plots <- all_comparative_plots[[k]]
        matching_plots <- smooth_type_plots[grepl(paste0("_", smooth_type_name, "_plots$"), names(smooth_type_plots))]
        
        if (length(matching_plots) > 0) {
          plot_list <- matching_plots[[1]]
          row_plots[[paste("Plot", k)]] <- plot_list[[plot_identifiers[i]]]
        }
      }
      
      row_combined <- arrangeGrob(grobs = row_plots, ncol = num_plot_lists)
      
      row_label <- textGrob(smooth_type_name, gp = gpar(fontsize = 15))
      
      row_plots_with_label <- arrangeGrob(
        arrangeGrob(row_label, nrow = 1),
        row_combined,
        ncol = 2,
        widths = unit.c(unit(2, "lines"), unit(1, "npc") - unit(2, "lines"))
      )
      
      variable_plots[[paste("Row", j)]] <- row_plots_with_label
    }
    
    x_grob <- textGrob(x_label, gp = gpar(fontsize = 15))
    y_grob <- textGrob(y_labels[[plot_identifiers[i]]], gp = gpar(fontsize = 15), rot = 90)
    
    combined_rows <- arrangeGrob(grobs = variable_plots, ncol = 1, left = y_grob, bottom = x_grob)
    
    col_labels <- lapply(knots_vector, function(knot) textGrob(paste("k =", knot), gp = gpar(fontsize = 15)))
    col_combined <- arrangeGrob(grobs = col_labels, ncol = length(knots_vector))
    
    final_combined_plots[[plot_identifiers[i]]] <- arrangeGrob(
      col_combined,
      combined_rows,
      nrow = 2,
      heights = unit.c(unit(1, "lines"), unit(1, "npc") - unit(1, "lines"))
    )
  }
  
  for (identifier in plot_identifiers) {
    grid.newpage()
    grid.draw(final_combined_plots[[identifier]])
    
    if (!is.null(output_dir)) {
      file_path <- file.path(output_dir, paste0(file_name, "_", identifier, ".png"))
      ggsave(filename = file_path, plot = final_combined_plots[[identifier]], width = 5, height = 6)
    }
  }
  
  return(final_combined_plots)
}
