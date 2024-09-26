# Load libraries
library(macpan2)
library(tidyverse)
library(mgcv)
library(outbreaks)
library(readxl)
library(patchwork)

# Load Scripts 
source("R/Functions/construct_simulator/construct_simulator_ScarletandMeasles.R")
source("R/Functions/calibrate_multiple_smoothers/calibrate_multiple_smoothers.R")
source("R/Functions/smooth2con.R")
source("R/Functions/update_likelihood/update_likelihood_expressions_Scarlett.R")
source("R/Functions/extract_simulation_results.R")
source("R/Functions/plotting/generate_comparative_plots_IrelandandScarlet.R")
source("R/Functions/AIC_multplie_sim.R")
source("R/Functions/MSE.R")

set.seed(1L)

#______________________Load and Prep Data_______________________________________
scarlet_fever_ontario <- readRDS(file = "Data/scarlet_fever_ontario.RDS")
# get data between 1929-08-01 and 1930-10-01
observed_data = (scarlet_fever_ontario
                 ## select the variables to be modelled -- a time-series of case reports.
                 |> select(period_end_date, cases_this_period)
                 
                 ## change the column headings so that they match the columns
                 ## in the simulated trajectories.
                 |> mutate(matrix = "reports")
                 |> rename(value = cases_this_period)
                 
                 ## create a `time` column with the time-step IDs that will correspond
                 ## to the time-steps in the simulation. this column heading also 
                 ## must match the column with the time-steps in the simulated trajectories
                 |> mutate(time = seq_along(period_end_date))
)
incidence = observed_data$value

#____________________________Model Preparation__________________________________

# Define time scale
step_length = 1


# Smoothing parameters
time_steps <- length(incidence) * step_length
data <- data.frame(time = 1:time_steps)
num_variables = 9
smooth_list <-c("gp", "tp", "cr", "bs")

# Kernel Parameters (Gaussian Process Only)
kernel = -2 # Matern with kappa= 1.5
range = 2 # Determines how the spatial correlation changes with distance.
cov_fun <- c(kernel,range,1)

# Starting Compartmental Parameters
N = 3392636
I_0 = 50
mean_log_I_0 = log(I_0)
sd_log_I_0 =  0.1
beta = 1
gamma = 1/10 # Recovery rate
mean_log_gamma <- log(gamma)
sd_log_gamma  = 0.1



# Place parms into lists 
smooth_parms = list(smooth = smooth,
                     num_variables = num_variables,
                     data = data,
                     cov_fun = cov_fun
)

compartmental_parms = list(gamma = gamma,
                           beta = beta,
                           mean_log_gamma = mean_log_gamma,
                           sd_log_gamma = sd_log_gamma,
                           mean_log_I_0 = mean_log_I_0,
                           sd_log_I_0 = sd_log_I_0,  
                           N = N,
                           I_0 = I_0
)  

# Starting values for parms to be estiamted
model_params <- list(
  log_I_sd = 1, # sd for predicted incidence
  log_smooth_sd = 1, # sd for smoothing coefficient 
  b_0 = 1, # intercept for smooth
  log_I_0 = 1,  
  log_gamma = 1 
)

# Load Model  Spec
sir_spec <- mp_tmb_library("starter_models", "sir", package = "macpan2")


# Calibrate using basis in smooth_list
calibrated_smoothers <- calibrate_multiple_smoothers( smooth_list, smooth_parms,
                                                      compartmental_parms,
                                                      model_parms,
                                                      sir_spec, step_length,
                                                      eq=FALSE, incidence)
#______________________________ Plot ------------------------------------------#
x_label <- "Week"
start_date <- as.Date("1929-08-03")
time_steps <- length(incidence)
dates_vector <- seq.Date(from = start_date, by = "week", length.out = time_steps)


plots <- generate_comparative_plots(calibrated_smoothers, data_name = "p3",
                                    aggregate = FALSE, Date = dates_vector,
                                    observed_vector = incidence)


plot_incidence_gp <- plots$p3_gp_plots$plot_incidence
plot_beta_gp <- plots$p3_gp_plots$plot_beta
plot_R_t_gp <- plots$p3_gp_plots$plot_R_t

plot_incidence_tp <- plots$p3_tp_plots$plot_incidence
plot_beta_tp <-  plots$p3_tp_plots$plot_beta
plot_R_t_tp <- plots$p3_tp_plots$plot_R_t

plot_incidence_cr <- plots$p3_cr_plots$plot_incidence
plot_beta_cr <- plots$p3_cr_plots$plot_beta
plot_R_t_cr <- plots$p3_cr_plots$plot_R_t

plot_incidence_bs <- plots$p3_bs_plots$plot_incidence
plot_beta_bs <- plots$p3_bs_plots$plot_beta
plot_R_t_bs <- plots$p3_bs_plots$plot_R_t

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


# Compute conditional AIC and create table 
all_simulators <- list( s3 = calibrated_smoothers)
table <- create_aic_table(all_simulators, knots = c(9), save_to_file = FALSE,
                          output_file = "Tables/Scarlet_combined_k9")

# Compute MSE of predicted incidence and data. Create Table. 
compute_mse_table(simulator_list = calibrated_smoothers, incidence = incidence,
                  aggregate = FALSE)

# Combine MSE and AIC Tables 
  mse_data <- data.frame(
    `Basis Type` = c("gp", "tp", "cr", "bs"),
    MSE = c(1283.815, 1191.876, 1234.858, 1329.298),
    check.names = FALSE
  )
  
  delta_cAIC_data <- data.frame(
    `Basis Type` = c("gp", "tp", "cr", "bs"),
    cAIC = c(-12.60388, -12.17104, -13.32448, -12.76761),
    check.names = FALSE
  )
  
  # Merge the data frames
  combined_df <- merge(mse_data, delta_cAIC_data, by = "Basis Type",
                       check.names = FALSE)
  
  # Subtract the minimum value from each element in the 'cAIC' column
  min_cAIC <- min(combined_df$cAIC)
  min_MSE <- min(combined_df$MSE)
  combined_df$cAIC <- combined_df$cAIC - min_cAIC
  combined_df$MSE <- combined_df$MSE - min_MSE
  

