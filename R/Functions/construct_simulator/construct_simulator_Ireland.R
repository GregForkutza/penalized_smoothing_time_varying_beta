create_indicator_matrix <- function(length_of_vector, period_days) {
  # Calculate the number of periods needed to cover the length_of_vector
  num_periods <- ceiling(length_of_vector / period_days)
  
  # Initialize the indicator matrix
  indicator_matrix <- matrix(0, nrow = length_of_vector, ncol = num_periods)
  
  for (period in 1:num_periods) {
    # Calculate start and end indices for each period
    start_index <- (period - 1) * period_days + 1
    end_index <- min(period * period_days, length_of_vector)
    
    # Calculate the length of the current period
    current_period_length <- end_index - start_index + 1
    
    # Assign 1/current_period_length to each day in the period to average over the correct number of days
    indicator_matrix[start_index:end_index, period] <- 1 / current_period_length
  }
  
  return(indicator_matrix)
}

get_parms_groups_sums <- function(data, step_size = 7, offset = FALSE) {
  # Set offset to start counting days at 1 or 0
  index_offset <- ifelse(offset, 1, 0)
  
  # Determine the number of steps (e.g., days) in the data
  num_steps <- nrow(data)
  
  # Calculate week indices, adjusting for offset
  wk <- as.integer(((1:num_steps) - 1) %/% step_size + index_offset)
  
  # Calculate the total number of weeks, adjusting for offset
  nweek <- max(wk) + (1 - index_offset)
  
  return(list(wk = wk, nweek = nweek))
}

#' Construct a Smooth Simulator for Epidemiological Models
#'
#' Constructs a simulator object for epidemiological modeling with specified 
#' smoothing parameters. This function integrates the construction of basis and 
#' penalty matrices, initialization of model parameters, and sets up expressions 
#' for dynamic model simulation.
#'
#' @param model The epidemiological model from the `mp_tmb_library` to be used 
#'        in the simulation.
#' @param smooth A character string specifying the type of smoothing to be applied.
#'        Supported types include "tp" for thin plate splines, "cr" for cubic regression 
#'        splines, "bs" for B-splines, "ps" for P-splines, and "ad" for adaptive smoothing.
#' @param num_variables The number of variables (basis functions) to be used in the 
#'        smoothing process.
#' @param data A data frame containing the covariate values for which the model 
#'        simulation is to be conducted. It should include a 'time' column.
#' @param cov_fun (Optional) A specification for the covariance function in cases 
#'        where Gaussian process smoothing is applied. Not used for other types of smoothing.
#' @param b_sd Standard deviation for the normal distribution from which the 
#'        starting values for basis coefficients are drawn.
#' @param alpha The decay parameter used specifically for Ornstein-Uhlenbeck-like 
#'        smoothing. It controls the rate of exponential decay.
#' @param beta The transmission rate parameter in the epidemiological model.
#' @param gamma The recovery rate parameter in the epidemiological model.
#' @param phi The waning immunity rate parameter in the epidemiological model.
#' @param N The total population size.
#'
#' @return An object of class `mp_simulator` representing the configured epidemiological 
#'         model simulator, ready for simulation.
#'
#' @examples
#' data <- data.frame(time = 1:100)
#' simulator <- construct_smooth_simulator(model = "sir_waning",
#'                                         smooth = "tp",
#'                                         num_variables = 10,
#'                                         data = data,
#'                                         alpha = 0.1,
#'                                         beta = 1,
#'                                         gamma = 0.3,
#'                                         phi = 0.01,
#'                                         N = 1000)
#'
#' @export
#' @importFrom stats rnorm log
construct_simulator <- function(spec,
                                smooth_parms,
                                compartmental_parms,
                                step_length = 1,
                                prevalence = FALSE,
                                eq = FALSE,
                                aggregate = FALSE
) {
  # Convert lists to environments
  e_smooth_parms <- list2env(smooth_parms, envir = environment())
  e_compartmental_parms <- list2env(compartmental_parms, envir = environment())
  
  # Manually copy each parameter into the function's environment
  for (name in names(smooth_parms)) {
    assign(name, e_smooth_parms[[name]], envir = environment())
  }
  for (name in names(compartmental_parms)) {
    assign(name, e_compartmental_parms[[name]], envir = environment())
  }

  if(eq) {
    spec = mp_tmb_insert(spec,
                         phase = "before",
                         expressions = list(S ~ (gamma * N) / beta),
                         at = 1
    )
    
    spec = mp_tmb_insert(spec,
                         phase = "before",
                         expressions = list(I ~ (phi * N * (beta - gamma)) / (beta * (phi + gamma))),
                         at = 1
    )
    
    spec = mp_tmb_insert(spec,
                         phase = "before",
                         expressions = list( R ~ (gamma * N * (beta - gamma)) / (beta * (phi + gamma))),
                         at = 1
    )
  } else {
    # Initialize compartment starting values
    spec = mp_tmb_insert(spec,
                         phase = "before",
                         expressions = list(S ~ N - I_0),
                         at = 1
    )
    
    spec = mp_tmb_insert(spec,
                         phase = "before",
                         expressions = list(I ~ I_0),
                         at = 1
    )
    
    spec = mp_tmb_insert(spec,
                         phase = "before",
                         expressions = list(R ~ 0),
                         at = 1
    )
  }
  
  if(step_length != 1) {
    # Insert step-specific expression for aggregation
    
    # Define the output type 
    output_type <-  "incidence"
    observation_phase <-  "wk_incidence"
    group_sum_fun <- "rbind_time(infection)"
    
    # Add expression for aggregating trajectory before optimization
    D <- get_parms_groups_sums(data, step_size = step_length)
    spec = mp_tmb_insert(spec,
                         at = Inf,
                         integers = list(wk = D$wk),
                         phase = "during",
                         default = c(list(nweek = D$nweek)
                         )
    )
    
    expression_text <- paste0(observation_phase, " ~ group_sums(", group_sum_fun, ", wk, rep(0, nweek))")
    expression_formula <- as.formula(expression_text)
    
    spec = mp_tmb_insert(spec,
                         phase = "after",
                         expressions = list(expression_formula),
                         must_save = observation_phase,
                         at = 1
    )
    
  } 
  
  
  # Add expressions to record effective and reproductive number
  spec = mp_tmb_insert(spec,
                       phase = "before",
                       expressions = list(eta ~ b_0 + (X %*% b)),
                       at = Inf
  )
  
  spec = mp_tmb_insert(spec,
                       phase = "during",
                       expressions = list(theta ~ eta[time_step(1)]),
                       at = 1
                       
  )
  
  spec = mp_tmb_insert(spec,
                       phase = "during",
                       expressions = list(beta ~ exp(eta[time_step(1)])),
                       at = 2
                       
  )
  
  spec = mp_tmb_insert(spec,
                       phase = "during",
                       expressions = list(R_t ~ (log(beta) - log(gamma) + log(S) - log(N))),
                       at = 3
  )
  
  # Define the output type based on prevalence
  output_type <- "incidence"
  
  
  # Extract the number of time steps from the data
  time_steps <- nrow(data)
  
  # Construct basis and penalty matrices
  smooth_output <- smooth2con(smooth, num_variables, data, cov_fun, alpha)
  # Generate starting values for basis coefficients 
  b <- rnorm(ncol(smooth_output$X-1), sd = b_sd)
  
  
  # Construct the base list of outputs
  outputs_base <- c("S", "I", "R", "infection", "beta","theta", "R_t")
  
  # Conditionally append observation_phase to outputs if step_length is not equal to 1
  outputs <- if (step_length != 1) c(outputs_base, observation_phase) else outputs_base
  # Initialize the simulator with the conditionally constructed outputs
  simulator <- mp_simulator(
    spec,
    time_steps = time_steps,
    outputs = outputs,
    default = c(list(gamma = gamma,
                     beta = beta,
                     N = N,
                     I_0 = I_0,
                     X = smooth_output$X,
                     P = smooth_output$P,
                     P_logdet = smooth_output$P_logdet,
                     b = b,
                     b_0 = log(beta),
                     k = num_variables,
                     mean_log_gamma = mean_log_gamma,
                     sd_log_gamma = sd_log_gamma,
                     mean_log_I_0 = mean_log_I_0,
                     sd_log_I_0 = sd_log_I_0,
                     nweek = D$nweek,
                     wk = D$wk
    )
    )
  )
  
  simulator$default <- list(
    gamma = gamma,
    beta = beta, 
    N = N, 
    I_0 = I_0,
    X = smooth_output$X,
    P = smooth_output$P,
    P_logdet = smooth_output$P_logdet,
    b = b,
    b_0 = log(beta),
    step_length = step_length
  )
  
  return(simulator)
}
