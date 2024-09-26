# Load necessary libraries
library(macpan2)
library(tidyverse)
library(mgcv)
library(patchwork)
library(gridExtra)
library(grid)

source("R/Functions/construct_simulator/construct_simulator_sim.R")
source("R/Functions/smooth2con.R")
source("R/Functions/update_likelihood/update_likelihood_sim.R")
source("R/Functions/calibrate_multiple_smoothers/calibrate_multiple_smoothers.R")
source("R/Functions/generate_simulated_data.R")
source("R/Functions/extract_simulation_results.R")
source("R/Functions/plotting/generate_comparative_plots_sim.R")
source("R/Functions/AIC_multplie_sim.R")
source("R/Functions/MSE.R")

# Model parameters   
# Data parameters
set.seed(1L)
time_steps <- 200
step_length <- 1

# Smoothing parameters
data <- data.frame(time = 1:time_steps)
initial_smooth <- "cr" # Define basis used to generate data. 
num_variables <- 20 # Define number of knots used to generate data. 
num_variables_fit = 20 # Define number of knots used to fit data. 
smooth_list <-c("gp","tp", "cr", "bs") # mgcv basis used to fit data. 

# Kernel Parameters (Gaussian Process Only)
kernel = -2 # Matern with kappa=1.5
range = 2 # Determines how the spatial correlation changes with distance.
cov_fun = c(kernel,range,1)


# Starting Compartmental Parameters
sim_sd = 0.2 # Add noise to simulated data
N = 4000000
b_sd= 1 # sd for initial values for smoothing coefficients. 
beta = 1
gamma = 1/7 # Recovery rate.
phi = 0.001 # Waning immunity rate.


# Place parms into lists 
smooth_parms_initial = list(smooth = initial_smooth,
                            num_variables = num_variables,
                            data = data,
                            cov_fun = cov_fun
)

smooth_parms_fit = list(smooth = smooth,
                        num_variables = num_variables_fit,
                        data = data,
                        cov_fun = cov_fun
)

compartmental_parms_initial = list(beta = beta,
                                   gamma = gamma,
                                   phi = phi,
                                   N = N,
                                   b_sd = b_sd
)       


compartmental_parms_fit  = list(beta = beta,
                                gamma = gamma,
                                phi = phi,
                                N = N,
                                b_sd = b_sd
)    

# Starting values for parms to be optimized. 
model_params <- list(
  log_I_sd = 1, # sd for predicted incidence.
  log_smooth_sd = 1, # sd for smoothing coefficient.
  b_0 = 1 # Intercept for smooth. 
)

# Load Model Spec
sir_spec <- mp_tmb_library("starter_models", "sir_waning", package = "macpan2")

# Create simulated data
sim_data_results <- create_simulated_data(sir_spec, smooth_parms_initial, compartmental_parms_initial, step_length, sim_sd)
simulator_initial <- sim_data_results$simulator_initial
sim_data <- sim_data_results$sim_data
incidence <- sim_data_results$incidence
results_initial <- extract_simulation_results(simulator = simulator_initial, conf.levels = c(0.5, 0.95), aggregate = FALSE)

# Calibrate multiple smoothers across smoothing basis
calibrated_smoothers <- calibrate_multiple_smoothers(smooth_list, smooth_parms_fit,
                                                     compartmental_parms_fit, model_parms, sir_spec, step_length, eq=TRUE, incidence)

# Combine results into a list 
return_list <- list(
  simulator_initial = list(smooth_initial = simulator_initial, sim_data = sim_data),
  calibrated_simulator = calibrated_smoothers,
  incidence = incidence
)

# Create plots 

  # This function creates panes of plots that vary pairs of smoothing basis and estimated quantity. 
  plots <- generate_comparative_plots(return_list)
  
  # Extract the plots from plots
  plot_incidence_gp <- plots$smooth_initial_gp_plots$plot_incidence
  plot_beta_gp <- plots$smooth_initial_gp_plots$plot_beta
  plot_R_t_gp <- plots$smooth_initial_gp_plots$plot_R_t
  
  plot_incidence_tp <- plots$smooth_initial_tp_plots$plot_incidence
  plot_beta_tp <- plots$smooth_initial_tp_plots$plot_beta
  plot_R_t_tp <- plots$smooth_initial_tp_plots$plot_R_t
  
  plot_incidence_cr <- plots$smooth_initial_cr_plots$plot_incidence
  plot_beta_cr <- plots$smooth_initial_cr_plots$plot_beta
  plot_R_t_cr <- plots$smooth_initial_cr_plots$plot_R_t
  
  # Assuming you have the bs plots as well (replace these with actual plot objects)
  plot_incidence_bs <- plots$smooth_initial_bs_plots$plot_incidence
  plot_beta_bs <- plots$smooth_initial_bs_plots$plot_beta
  plot_R_t_bs <- plots$smooth_initial_bs_plots$plot_R_t
  
  # Create column labels
  col_labels <- textGrob(c("Incidence", expression(beta[t]), expression(R[t])), 
                         x = unit(c(0.17, 0.5, 0.83), "npc"), y = unit(0.5, "npc"), 
                         gp = gpar(fontsize = 15, fontface = "bold"))
  
  # Create row labels
  row_labels <- c("gp", "tp", "cr", "bs")
  row_label_grobs <- lapply(row_labels, function(label) {
    textGrob(label, rot = 0, gp = gpar(fontsize = 15, fontface = "bold"))
  })
  
  # Arrange the individual plots
  plots_gp <- arrangeGrob(plot_incidence_gp, plot_beta_gp, plot_R_t_gp, ncol = 3)
  plots_tp <- arrangeGrob(plot_incidence_tp, plot_beta_tp, plot_R_t_tp, ncol = 3)
  plots_cr <- arrangeGrob(plot_incidence_cr, plot_beta_cr, plot_R_t_cr, ncol = 3)
  plots_bs <- arrangeGrob(plot_incidence_bs, plot_beta_bs, plot_R_t_bs, ncol = 3)
  
  # Arrange the combined plots with row labels
  combined_plots <- arrangeGrob(
    row_label_grobs[[1]], plots_gp,
    row_label_grobs[[2]], plots_tp,
    row_label_grobs[[3]], plots_cr,
    row_label_grobs[[4]], plots_bs,
    ncol = 2,
    widths = c(0.1, 0.9)
  )
  
  # Final layout with column labels
  final_plot <- arrangeGrob(
    col_labels, combined_plots,
    nrow = 2,
    heights = c(0.1, 0.9)
  )
  
  # Draw the final plot
  grid.newpage()
  grid.draw(final_plot)

# Save the final plot
# ggsave("sim_agg_combined_cr.png", plot = final_plot, width = 5, height = 6)

  # Compute conditional AIC and create table. 
all_simulators <- list(s3 = calibrated_smoothers)
table <- create_aic_table(all_simulators, knots = c(20), save_to_file = FALSE, output_file =output_file )

# Compute MSE of trajectory against noisy data and create table
compute_mse_table(simulator_list = calibrated_smoothers, incidence = incidence, aggregate = FALSE)
