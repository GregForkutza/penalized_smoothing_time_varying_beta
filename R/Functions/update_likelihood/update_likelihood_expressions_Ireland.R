#' Update Likelihood Expressions in Simulator
#'
#' This function updates a given simulator object with new likelihood expressions
#' based on observed data. It allows for the incorporation of both default and 
#' custom statistical parameters for smooth and non-linear least squares (NLL) 
#' components. It's designed to facilitate the integration of simulated or real 
#' observational data into the simulator's model fitting process.
#'
#' @param simulator An object of class `mp_simulator` representing the 
#'        epidemiological model simulator to be updated.
#' @param data A data frame or list containing the observational data to be 
#'        used in the likelihood expressions. The data should be indexed by the 
#'        name provided in the `name` parameter.
#' @param name A character string specifying the column name or list key in `data`
#'        that contains the observational data to be used in the likelihood calculations.
#' @param custom_params (Optional) A list containing custom parameter frames for 
#'        `smooth_params` and/or `nll_params` to override the default parameter 
#'        settings. Each should be a data frame matching the structure required 
#'        for smooth and NLL parameters, respectively. If `NULL`, default 
#'        parameters are used.
#'
#' @details
#' The function first checks for custom parameters provided by the user; if none
#' are found, it falls back to default parameter settings. It then extracts the 
#' observational data from the `data` argument based on the `name` provided, and 
#' updates the simulator object with new matrices and expressions for the 
#' likelihood calculations. This update includes adding noise to the observational 
#' data, setting up matrices for observed values, standard deviation of the 
#' observations, and the smoothness parameter, and inserting the new expressions 
#' for calculating the log likelihood.
#'
#' @return 
#' The function does not explicitly return a value but modifies the `simulator`
#' object in place by adding new matrices and updating its expressions for
#' likelihood calculations.
#'
#' @examples
#' # Assuming `simulator` is a pre-existing mp_simulator object
#' # and `data` contains observational data named "I_obs"
#' update_likelihood_expressions(simulator, data, "I_obs")
#'
#' # With custom parameters
#' custom_params <- list(
#'   smooth_params = data.frame(mat = c("log_I_sd"), row = c(0), col = c(0), default = c(-1)),
#'   nll_params = data.frame(mat = "b", row = 0:3, col = rep(0, 4), default = rep(0.5, 4))
#' )
#' update_likelihood_expressions(simulator, data, "I_obs", custom_params = custom_params)
#'
#' @export
#' @importFrom dplyr filter mutate across pull
#' @importFrom stats dnorm rbind
update_likelihood_expressions <- function(
    simulator,
    observed_vector,
    model_params,
    step_length,
    prevalence = FALSE
) {
  # This function accepts a named list and constructs the default parameters matrix
  read_model_params <- function(params_list) {
    
    params_df <- as.data.frame(unlist(params_list))
    names(params_df) <- c("default")
    params_df$mat <- rownames(params_df)
    params_df$row <- 0
    params_df$col <- 0
    params_df <- params_df[, c("mat", "row", "col", "default")]
    
    return(params_df)
  }
  # Format list of model parameters into data frame.
  nll_params <- read_model_params(model_params)
  # Function to construct the likelihood expression
  construct_likelihood_expression <- function(observation_type) {
    paste(
      "-sum(dnorm(I_obs, ", observation_type, ", I_sd)) + ",
      "-dnorm(log_gamma, mean_log_gamma, sd_log_gamma)",
      " +k/2*log(2 *3.141593)",
      " -P_logdet + ",
      " -dnorm(log_I_0, mean_log_I_0, sd_log_I_0) +",
      "2 * log(smooth_sd) + ",
      "((t(b) %*% P %*% b) / smooth_sd^2)",
      sep = ""
    )
  }

  # Construct matrix of smoothing coefficients. 
  num_variables <- length(simulator$get$initial("b")) + 1
  smooth_params <- data.frame(
    mat = "b",
    row = 0:(num_variables - 2),
    col = 0,
    default = 0.1
  )
  
  # Add matrices to simulator
  simulator$add$matrices(I_obs = observed_vector,
                         log_lik = empty_matrix,
                         I_sd = 1,
                         smooth_sd = 1
  )
  
  # Construct likelihood expression
  if (step_length == 1) {
    observation_type <- if(prevalence) "rbind_time(I)" else "rbind_time(infection)"
  } else {
    observation_type <- if(prevalence) "wk_prevalence" else "wk_incidence"
  }
  likelihood_expr <- construct_likelihood_expression(observation_type)
  
  # Insert expression into simulator
  simulator$insert$expressions(
    reformulate(termlabels = likelihood_expr, response = "log_lik"),
    .phase = "after",
    .at = Inf
  )
  
  # Add transformations and update frames
  simulator$replace$obj_fn(~ log_lik)
  simulator$add$transformations(Log("I_sd"))
  simulator$add$transformations(Log("smooth_sd"))
  simulator$add$transformations(Log("gamma"))
  simulator$add$transformations(Log("I_0"))
  simulator$replace$params_frame(nll_params)
  simulator$replace$random_frame(smooth_params)
}