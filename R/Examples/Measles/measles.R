# Load libraries
library(macpan2)
library(tidyverse)
library(mgcv)
library(outbreaks)
library(readxl)
library(patchwork)

# Load scripts 
source("R/Functions/construct_simulator/construct_simulator_ScarletandMeasles.R")
source("R/Functions/calibrate_multiple_smoothers/calibrate_multiple_smoothers.R")
source("R/Functions/smooth2con.R")
source("R/Functions/update_likelihood/update_likelihood_expressions_Scarlett.R")
source("R/Functions/extract_simulation_results.R")
source("R/Functions/plotting/generate_comparative_plots_Measles.R")
source("R/Functions/AIC_multplie_sim.R")
set.seed(1L)

#______________________Load and Prep Data_______________________________________
measles_london <- read_csv(file = "Data/meas_uk__lon_1944-94_wk.csv")[-(1:6),]
measles1944to1946 = measles_london[2:nrow(measles_london),]
incidence = measles1944to1946$...4 %>% as.numeric() 
incidence = incidence[1:2140] # Exclude last 10 years of data. 
#____________________________Model Preperation__________________________________
step_length = 1

# Smoothing parameters
time_steps <- length(incidence) * step_length
data <- data.frame(time = 1:time_steps)
num_variables = 30
smooth_list <- "gp"

# Kernel Parameters (Gaussian Process Only)
kernel = -2 # Matern with kappa= 1.5
range = 20 # Determines how the spatial correlation changes with distance.
cov_fun = c(kernel,range,1)

# Starting Compartmental Parameters
N = 6000000
I_0 = 250
mean_log_I_0= log(I_0)
sd_log_I_0 =  0.1
beta = 1
gamma = 1/8 # Recovery rate
mean_log_gamma <- log(gamma)
sd_log_gamma  <- 0.1



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

# Starting values for estimatd parameters
model_params <- list(
  log_I_sd = 1, 
  log_smooth_sd = 1, 
  b_0 = 1, # intercept for smooth
  log_I_0 = 1,  # initial value for prevalence
  log_gamma = 1
)

# Load Model Spec
sir_spec <- mp_tmb_library("starter_models", "sir", package = "macpan2")

# Calibrate multiple smoothers
calibrated_smoothers <- calibrate_multiple_smoothers( smooth_list, smooth_parms,
                                                      compartmental_parms, 
                                                      model_parms,
                                                      sir_spec, step_length, 
                                                      eq=FALSE, incidence)


# Create plots
plots <- generate_comparative_plots(calibrated_smoothers,
                                     data_name = "Scarlet_k6")

all_comparative_plots <- list(p1 =plots)
all_plots(all_comparative_plots, output_dir = "Plots",
          file_name = "UK_Measles_plot_test")




