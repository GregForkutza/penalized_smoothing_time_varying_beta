# Load necessary libraries
library(macpan2)
library(tidyverse)
library(mgcv)
library(outbreaks)
library(readxl)
library(patchwork)

source("R/Functions/construct_simulator/construct_simulator_Ireland.R")
source("R/Functions/calibrate_multiple_smoothers/calibrate_multiple_smoothers.R")
source("R/Functions/smooth2con.R")
source("R/Functions/update_likelihood/update_likelihood_expressions_Ireland.R")
source("R/Functions/extract_simulation_results.R")
source("R/Functions/plotting/generate_comparative_plots_IrelandandScarlet.R")
source("R/Functions/AIC_multplie_sim.R")
source("R/Functions/MSE.R")





set.seed(1L)

# Load Data and aggreate to weekly scale. 
Data <- read_excel("Data/pcbi.1010206.s009.xlsx")
incidence <- Data$ReportedCases
step_length = 7
time_steps <- length(incidence)
data <- data.frame(time = 1:time_steps)
D <- get_parms_groups_sums(data, step_size = step_length)
wk_inc_obs <- tapply(incidence, D$wk, sum) %>% as.vector()


# Smoothing parameters
data <- data.frame(time = 1:time_steps)
num_variables = 7
smooth_list <-c("gp", "tp", "cr")

# Kernel Parameters (Gaussian Process Only)
kernel = -2
range =1
cov_fun = c(kernel, range,1)

# Starting Compartmental Parameters
beta = 1 # Example value for beta
N = 5000000
I_0 = 10
mean_log_I_0= log(I_0)
sd_log_I_0 =  0.1
gamma = 1/6 # Recovery rate
mean_log_gamma <- log(gamma)
sd_log_gamma  <- 0.1
b_sd = 1

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

# Starting values for optimized parameters
model_params <- list(
  log_I_sd = 1, 
  log_smooth_sd = 1,  
  b_0 = 1, # intercept for smooth
  log_I_0 = 1,  # initial value for prevalence
  log_gamma = 1
)

# Load Model Spec
sir_spec <- mp_tmb_library("starter_models", "sir", package = "macpan2")

calibrated_smoothers <- calibrate_multiple_smoothers(smooth_list, smooth_parms, 
                                                     compartmental_parms, 
                                                     model_parms, sir_spec,
                                                     step_length, eq=FALSE, 
                                                     observed_vector = wk_inc_obs
                                                     )

# Create Plots
x_label <- "Day"
start_date <- as.Date("2020-02-20")
time_steps <- length(wk_inc_obs)
dates_vector <- seq.Date(from = start_date, by = "week", length.out = time_steps)

# Generate plots
plots <- generate_comparative_plots(calibrated_smoothers, data_name = "p3", aggregate = TRUE, Date = dates_vector, observed_vector = wk_inc_obs)

# Extract the plots from plots
plot_incidence_gp <- plots$p3_gp_plots$plot_incidence
plot_beta_gp <- plots$p3_gp_plots$plot_beta
plot_R_t_gp <- plots$p3_gp_plots$plot_R_t

plot_incidence_tp <- plots$p3_tp_plots$plot_incidence
plot_beta_tp <-  plots$p3_tp_plots$plot_beta
plot_R_t_tp <- plots$p3_tp_plots$plot_R_t

plot_incidence_cr <- plots$p3_cr_plots$plot_incidence
plot_beta_cr <- plots$p3_cr_plots$plot_beta
plot_R_t_cr <- plots$p3_cr_plots$plot_R_t

# Create column labels
col_labels <- textGrob(c("Incidence", expression(beta[t]), expression(R[t])), 
                       x = unit(c(0.17, 0.5, 0.83), "npc"), y = unit(0.5, "npc"), 
                       gp = gpar(fontsize = 15, fontface = "bold"))

# Create row labels
row_labels <- c("gp", "tp", "cr")
row_label_grobs <- lapply(row_labels, function(label) {
  textGrob(label, rot = 0, gp = gpar(fontsize = 15, fontface = "bold"))
})

# Arrange the individual plots
plots_gp <- arrangeGrob(plot_incidence_gp, plot_beta_gp, plot_R_t_gp, ncol = 3)
plots_tp <- arrangeGrob(plot_incidence_tp, plot_beta_tp, plot_R_t_tp, ncol = 3)
plots_cr <- arrangeGrob(plot_incidence_cr, plot_beta_cr, plot_R_t_cr, ncol = 3)


# Arrange the combined plots with row labels
combined_plots <- arrangeGrob(
  row_label_grobs[[1]], plots_gp,
  row_label_grobs[[2]], plots_tp,
  row_label_grobs[[3]], plots_cr,
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

# Compute Conditional AIC and create Table. 
all_simulators <- list(s1 = calibrated_smoothers)
output_file = "Tables/Ireland_agg_k(7)_bsd1_beta1_gamma6_sd01"
table <- create_aic_table(all_simulators, knots = c(7),
                          save_to_file = FALSE, output_file =output_file )

# Compute MSE of predicted incidence and data. 
compute_mse_table(simulator_list = calibrated_smoothers,
                  incidence = wk_inc_obs, aggregate = TRUE)

# Combine AIC and MSE table
  mse_data <- data.frame(
    `Basis Type` = c("gp", "tp", "cr"),
    MSE = c(46647.80, 106161.65, 86354.08),
    check.names = FALSE
  )
  
  delta_cAIC_data <- data.frame(
    `Basis Type` = c("gp", "tp", "cr"),
    cAIC = c(-8.62613, -9.727719, -9.865612),
    check.names = FALSE
  )
  
  # Merge the data frames
  combined_df <- merge(mse_data, delta_cAIC_data,
                       by = "Basis Type", check.names = FALSE)
  
  # Subtract the minimum value from each element in the 'cAIC' column
  min_cAIC <- min(combined_df$cAIC)
  min_MSE <- min(combined_df$MSE)
  combined_df$cAIC <- combined_df$cAIC - min_cAIC
  combined_df$MSE <- combined_df$MSE - min_MSE

