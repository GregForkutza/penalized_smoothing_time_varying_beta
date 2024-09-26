# There are two degrees of freedom here. 
# 1. simulated data vs data. This involves 
# 2. agg vs non agg. 


# For 1. It should only matter for calibrate_multiple_smoothers and plotting. 
# For 2/ Its relavent for constuct_simulator, Update_likelihood, calibrate_multiple_smoother, plotting, 

# Construct Simulator 
source("R/Functions/Base/construct_simulator.R") # Measles and Scarlet
source("R/Functions/Base/Simulations/construct_simulator_sim.R") # Sim 
source("R/Functions/Base/Simulations/Aggregated/construct_simulator_sim_agg.R") # Sim Agg
source("R/Functions/Base/RealExamples/construct_simulator_agg.R") # Ireland Covid 

# Update Likelihood
source("R/Functions/Base/Simulations/update_likelihood_sim.R") # Sim
source("R/Functions/Base/Simulations/Aggregated/update_likelihood_expressions_sim_agg.R") # sim Agg
source("R/Functions/Base/update_likelihood_expressions.R") # Scarlett
source("R/Functions/Base/RealExamples/update_likelihood_expressions_agg.R") # Ireland Covid 
source("R/Functions/Base/update_likelihood_expressions.R") # Measles 

# Smooth2Con
source("R/Functions/Base/smooth2con.R") # All

# Generate Sim data
source("R/Functions/Base/Simulations/generate_simulated_data.R") # Sim and Sim Agg

# Calibrate Multiple Smoothers
source("R/Functions/Base/Simulations/Aggregated/calibrate_multiple_smoothers_agg.R") # Sim Agg
source("R/Functions/Base/Simulations/factored_calibrate_multiple_smoothers.R") # Sim, Scarlett, Ireland, Measles 

# Extract Sim results
source("R/Functions/Base/extract_simulation_results.R") # All

# Plotting
source("R/Functions/Plotting/Simulations/generate_comparative_plots_fixed.R") # Sim, 
source("R/Functions/Base/Simulations/Aggregated/generate_comparative_plots_sim_agg.R") # Sim Agg, 
source("R/Functions/Plotting/RealExamples/generate_comparative_plots_agg.R") # Ireland, Scarlet
source("R/Functions/Base/RealExamples/generate_comparative_plots_singelsmooth_stackedparms.R") # Measles

# AIC
source("R/Functions/Base/Simulations/AIC_multplie_sim.R") # All

# MSE
source("R/Functions/Base/MSE.R") # All
