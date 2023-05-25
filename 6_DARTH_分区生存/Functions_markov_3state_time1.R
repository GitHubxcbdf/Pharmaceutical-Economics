#------------------------------------------------------------------------------#
####                         Decision Model                                 ####
#------------------------------------------------------------------------------#
#' 3-State Decision Model
#'
#' \code{decision_model} implements the 3-state health-sick-dead decision model.
#'
#' @param l_params_all List with all parameters of decision model
#' @param verbose Logical variable to indicate print out of messages
#' @return The transition probability array and the cohort trace matrix.
#' @export
decision_model <- function(l_params_all, verbose = FALSE) {
  with(as.list(l_params_all), {
    ########################### Process model inputs ###########################
    ## Model states
    v_names_states  <- c("Healthy", "Sick", "Dead")  # state names
    n_states        <- length(v_names_states)        # number of health states 
    
    ## Cycle names
    v_names_cycles  <- paste("cycle", 0:n_cycles)
    
    ## Probability of dying when healthy (age-dependent) - this is now a sequence of numbers
    
    # All starting healthy
    v_m_init <- c("Healthy" = 1, "Sick" = 0, "Dead" = 0)  
    
    ###################### Construct state-transition models ###################
    ### Initialize cohort trace for SoC 


    ## Store the cohort traces in a list 
    l_m_M <- list(A   =  m_M_trtA,
                  B   =  m_M_trtB)
    names(l_m_M) <- v_names_str
    
    ########################################## RETURN OUTPUT  ##########################################
    out <- list(l_m_M = l_m_M)
    
    return(out)
  }
  )
}

#------------------------------------------------------------------------------#
####              Calculate cost-effectiveness outcomes                     ####
#------------------------------------------------------------------------------#
#' Calculate cost-effectiveness outcomes
#'
#' \code{calculate_ce_out} calculates costs and effects for a given vector of parameters using a simulation model.
#' @param l_params_all List with all parameters of decision model
#' @param n_wtp Willingness-to-pay threshold to compute net monetary benefits (
#' NMB)
#' @return A dataframe with discounted costs, effectiveness and NMB.
#' @export
calculate_ce_out <- function(l_params_all, n_wtp = 100000){ # User defined
  with(as.list(l_params_all), {
    ### Run decision model to get transition dynamics array
    l_m_M <- decision_model(l_params_all = l_params_all)$l_m_M
    
    ### State rewards
    ## Scale by the cycle length 
    # Vector of state utilities under strategy SoC
    # Vector of state utilities under treatment A
    v_u_trtA   <- c(H  = u_H/12, 
                    S  = u_S/12, 
                    D  = u_D/12) * cycle_length
    # Vector of state costs under treatment A
    v_c_trtA   <- c(H  = c_H + c_trtA, 
                    S  = c_S, 
                    D  = c_D) * cycle_length
    # Vector of state utilities under treatment B
    v_u_trtB   <- c(H  = u_H/12, 
                    S  = u_S/12, 
                    D  = u_D/12) * cycle_length
    # Vector of state costs under treatment B
    v_c_trtB   <- c(H  = c_H + c_trtB, 
                    S  = c_S, 
                    D  = c_D) * cycle_length
    
    ## Store state rewards 
    # Store the vectors of state utilities for each strategy in a list 
    l_u   <- list(A  = v_u_trtA,
                  B  = v_u_trtB)
    # Store the vectors of state cost for each strategy in a list 
    l_c   <- list(A  = v_c_trtA,
                  B  = v_c_trtB)
    
    # assign strategy names to matching items in the lists
    names(l_u) <- names(l_c) <- v_names_str
    
    # Create empty vectors to store total utilities and costs 
    v_tot_qaly <- v_tot_cost <- vector(mode = "numeric", length = n_str)
    names(v_tot_qaly) <- names(v_tot_cost) <- v_names_str
    
    ## Loop through each strategy and calculate total utilities and costs 
    for (i in 1:n_str) {
      v_u_str <- l_u[[i]]   # select the vector of state utilities for the i-th strategy
      v_c_str <- l_c[[i]]   # select the vector of state costs for the i-th strategy
      
      ### Expected QALYs and costs per cycle 
      ## Vector of QALYs and Costs
      # Apply state rewards 
      v_qaly_str <- l_m_M[[i]] %*% v_u_str # sum the utilities of all states for each cycle
      v_cost_str <- l_m_M[[i]] %*% v_c_str # sum the costs of all states for each cycle
      
      ### Discounted total expected QALYs and Costs per strategy and apply within-cycle correction if applicable
      # QALYs
      v_tot_qaly[i] <- t(v_qaly_str) %*% (v_dwe * v_wcc)
      # Costs
      v_tot_cost[i] <- t(v_cost_str) %*% (v_dwc * v_wcc)
    }
    
    ## Vector with discounted net monetary benefits (NMB)
    v_nmb <- v_tot_qaly * n_wtp - v_tot_cost
    
    ## data.frame with discounted costs, effectiveness and NMB
    df_ce <- data.frame(Strategy = v_names_str,
                        Cost     = v_tot_cost,
                        Effect   = v_tot_qaly,
                        NMB      = v_nmb)
    
    return(df_ce)
  }
  )
}


#------------------------------------------------------------------------------#
####             Generate a PSA input parameter dataset                     ####
#------------------------------------------------------------------------------#
#' Generate parameter sets for the probabilistic sensitivity analysis (PSA)
#'
#' \code{generate_psa_params} generates a PSA dataset of the parameters of the 
#' cost-effectiveness analysis.
#' @param n_sim Number of parameter sets for the PSA dataset
#' @param seed Seed for the random number generation
#' @return A data.frame with a PSA dataset of he parameters of the 
#' cost-effectiveness analysis
#' @export
generate_psa_params <- function(n_sim = 1000, seed = 071818){
  set.seed(seed) # set a seed to be able to reproduce the same results
  df_psa <- data.frame(
    # Costs
    c_H       = rgamma(n_sim, shape = 16, scale = 25),        # cost of one cycle in healthy state
    c_S       = rgamma(n_sim, shape = 100, scale = 10),       # cost of one cycle in sick state
    c_D       = 0,                                            # cost of one cycle in dead state
    c_trtA    = 800,                                          # cost of treatment A (per cycle) in healthy state
    c_trtB    = 1500,                                         # cost of treatment B (per cycle) in healthy state
    
    # Utilities
    u_H       = rbeta(n_sim, shape1 =  1.5, shape2 = 0.0015), # utility when healthy 
    u_S       = rbeta(n_sim, shape1 = 49.5, shape2 = 49.5),   # utility when sick
    u_D       = 0                                             # utility when dead
  )
  return(df_psa)
}
