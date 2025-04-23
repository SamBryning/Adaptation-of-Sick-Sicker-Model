# ==============
# Adapted from code written by Robert Smith & Paul Schneider
# ==============

# Set-up
library(truncnorm)  # load the package trunc-norm
library(faux)

rm(list=ls())

f_wrapper <- function(
  
  #================================================================
  #                         Shiny inputs
  #================================================================
  n_sim   = 100,    
  WTP = 20000,
  model = "Sick Model",
  include_discrepancy = "No",
  discrepancy_spec = "Non-informative",
  NID_pwcorr = 0,
  NID_sd = 0.05, 
  ID1_1 = 3, 
  ID1_2 = 15, 
  ID2_1 = 100,
  ID2_2 = 15,
  ID3_1 = -0.01,
  ID3_2 = 0.01,
  ID4_1 = 100,
  ID4_2 = 20
  ){
  
  
  #================================================================
  #                         Non-shiny inputs
  #================================================================
  n_age_init = 25   # age at baseline default is 25
  n_age_max  = 110  # maximum age of follow up default is 110
  d_r     = 0.035  # discount rate for costs & QALYS (NICE 3.5%)
  c_Trt   = 2000 # cost of treatment
  
  n_t <- n_age_max - n_age_init # time horizon, number of cycles
  v_n <- c("H", "S1", "S2", "D") # the health states of the model:
  n_states <- length(v_n) # number of health states 
  
 
  
  
  #==========================================================
  #                   PSA INPUT FUNCTION
  #==========================================================
  
# Function that generates random sample for PSA, note dependant on wrapper input.
gen_psa <- function(n_sim = 1000){
   
  df_psa <- data.frame(
    
    # Transition probabilities (per cycle)
    p_HS1   = rbeta(n_sim, 30, 170),        # probability to become sick when healthy
    p_S1H   = rbeta(n_sim, 60, 60),        # probability to become healthy when sick
    p_S1S2  = rbeta(n_sim, 84, 716),        # probability to become sicker when sick
    
    p_HD    = rbeta(n_sim, 10, 1990)      ,  # probability to die when healthy
    hr_S1   = rlnorm(n_sim, log(3),  0.01),  # rate ratio of death in S1 vs healthy
    hr_S2   = rlnorm(n_sim, log(10), 0.02),  # rate ratio of death in S2 vs healthy 
    
    # Cost vectors with length n_sim
    c_H   = rgamma(n_sim, shape = 100, scale = 20)    , # cost of remaining one cycle in state H
    c_S1  = rgamma(n_sim, shape = 177.8, scale = 22.5), # cost of remaining one cycle in state S1
    c_S2  = rgamma(n_sim, shape = 225, scale = 66.7)  , # cost of remaining one cycle in state S2
    c_D   = 0                                         , # cost of being in the death state
    c_Trt = c_Trt,                                # cost of treatment (per cycle)
  
    
    # Utility vectors with length n_sim 
    u_H   = rtruncnorm(n_sim, mean =  0.8, sd = 0.01, b = 1), # utility when healthy
    u_S1  = rtruncnorm(n_sim, mean = 0.6, sd = 0.02, b = 1), # utility when sick
    u_S2  = rtruncnorm(n_sim, mean = 0.45, sd = 0.03, b = 1), # utility when sicker
    u_D   = 0                                               , # utility when dead
    u_Trt = rtruncnorm(n_sim, mean = 0.75, sd = 0.02, b = 1),  # utility when being treated
  
  # # Non-informative discrepancy terms
  NID = cbind(rnorm_multi(n_sim,
                    mu = c(0,0,0,0),
                    sd = c(as.numeric(NID_sd)*(etruncnorm(mean = 0.8, sd = 0.01, b = 1)), as.numeric(NID_sd)*mean(rgamma(n_sim, shape = 100, scale = 20)), 
                           as.numeric(NID_sd)*(etruncnorm(mean = 0.75, sd = 0.02, b = 1)), as.numeric(NID_sd)*c_Trt),
                    r = NID_pwcorr,
                    varnames = c("utility", "cost", "utility_trt", "cost_trt"))),
  
  # Informative discrepancy
  ID_utility = rbeta(n_sim, ID1_1, ID1_2),
  ID_cost = rgamma(n_sim, shape = ID2_1, scale = ID2_2),
  ID_utility_trt = rtruncnorm(n_sim, mean = ID3_1, sd = ID3_2, a=-0.2, b=0.2), # rbeta(n_sim, 2, 100))
  ID_cost_trt = mean(rgamma(n_sim, shape = ID4_1, scale = ID4_2))
  
   )
  return(df_psa)
}

#==============================================================================
#                                    MODEL FUNCTION
#==============================================================================

f_MM_sicksicker <- function(params) {
  with(as.list(params), {
    
    # compute internal parameters as a function of external parameter
    r_HD    = - log(1 - p_HD) # rate of death in healthy
    r_S1D   = hr_S1 * r_HD 	  # rate of death in sick
    r_S2D   = hr_S2 * r_HD  	# rate of death in sicker
    p_S1D   = 1 - exp(-r_S1D) # probability to die in sick
    p_S2D   = 1 - exp(-r_S2D) # probability to die in sicker
    
    v_dwe <- v_dwc <- 1 / (1 + d_r) ^ (0:n_t) # calculate discount weight for each cycle based on discount rate d_r
    
    # create transition probability matrix for NO treatment
    m_P <- reactive ({
      m_P <- matrix(0,
                  nrow = n_states, ncol = n_states,
                  dimnames = list(v_n, v_n))
    # fill in the transition probability array
    ### From Healthy
    m_P["H", "H"]  <- 1 - (p_HS1 + p_HD)
    m_P["H", "S1"] <- p_HS1
    m_P["H", "D"]  <- p_HD
    
    ### From Sick
    m_P["S1", "H"]  <- if(model == "Sick-Sicker Model"){
        p_S1H
      }
     else{
       1-p_S1D
     }
    
    m_P["S1", "S1"] <- if(model == "Sick-Sicker Model") {
      1 - (p_S1H + p_S1S2 + p_S1D)
    }
      else{
        0 # In Sick Model illness is temporary and patients transition back to healthy after one model cycle
      }
    
    m_P["S1", "S2"] <- if(model == "Sick-Sicker Model") {
        p_S1S2
      }
      else{
         0
        }

    m_P["S1", "D"]  <- p_S1D
    
    ### From Sicker (in Sick Sicker Model only)
    m_P["S2", "S2"] <- if(model == "Sick-Sicker Model") {
        1 - p_S2D
      }
      else{
        0
        }
   
    m_P["S2", "D"]  <- if(model == "Sick-Sicker Model") {
       p_S2D
      }
      else{
         0
      }

    ### From Dead
    m_P["D", "D"] <- 1

   return(m_P)
  })
   # print(m_P())
    
    
    m_TR <- reactive({
      
    m_TR <- matrix(NA, nrow = n_t + 1 , ncol = n_states, 
                   dimnames = list(0:n_t, v_n))     # create Markov trace (n_t + 1 because R doesn't understand  Cycle 0)
    
    m_TR[1, ] <- c(1, 0, 0, 0)                      # initialize Markov trace
   # print(m_TR)
    ############# PROCESS ###########################################

    for (t in 1:n_t){                              # throughout the number of cycles
      m_TR[t + 1, ] <- m_TR[t, ] %*%  m_P()        # estimate the Markov trace for cycle the next cycle (t + 1)
     on.exit(return(m_TR))
    }
    })

    ############ OUTPUT  ###########################################
    # create vectors of utility and costs for each state
    v_u_trt <- reactive({
      if(model == "Sick Model" & include_discrepancy == 'Yes' & discrepancy_spec == 'Non-informative'){
        v_u_trt <- c(min(c(u_H+NID.utility,1)), min(c(u_Trt+NID.utility_trt,1)), u_S2, u_D)
      }
      else if(model == "Sick Model" & include_discrepancy == 'Yes' & discrepancy_spec == 'Informative'){
        v_u_trt <- c(min(c(u_H-ID_utility,1)), min(c(u_Trt+ID_utility_trt,1)), u_S2, u_D)
      }
      else {
        v_u_trt <- c(u_H, u_Trt, u_S2, u_D)
      }
    })
    
    v_u_no_trt <- reactive({
      if(model == "Sick Model" & include_discrepancy == 'Yes' & discrepancy_spec == 'Non-informative'){
        v_u_no_trt <- c(min(c(u_H+NID.utility,1)), u_S1, u_S2, u_D)
      }
      else if(model == "Sick Model" & include_discrepancy == 'Yes' & discrepancy_spec == 'Informative'){
        v_u_no_trt <- c(min(c(u_H-ID_utility,1)), u_S1, u_S2, u_D)
      }
      else{
        v_u_no_trt <- c(u_H, u_S1, u_S2, u_D)
      }
    })
    
    v_c_trt <- reactive({
      if(model == "Sick Model" & include_discrepancy == 'Yes'& discrepancy_spec == 'Non-informative'){
        v_c_trt <- c(max(c(c_H + NID.cost,0)), c_S1 + max(c(c_Trt+NID.cost_trt,0)), c_S2 + c_Trt, c_D)
      }
      else if (model == "Sick Model" & include_discrepancy == 'Yes' & discrepancy_spec == 'Informative'){
        v_c_trt <- c(max(c(c_H + ID_cost,0)), c_S1 + max(c(c_Trt+ID_cost_trt,0)), c_S2 + c_Trt, c_D)
      }
        else{
        v_c_trt <- c(c_H, c_S1 + c_Trt, c_S2 + c_Trt, c_D)
      }
    })
    
    v_c_no_trt <- reactive({
      if(model == "Sick Model" & include_discrepancy == 'Yes' & discrepancy_spec == 'Non-informative'){
        v_c_no_trt <- c(max(c(c_H + NID.cost,0)), c_S1, c_S2 + c_Trt, c_D)
      }
      else if (model == "Sick Model" & include_discrepancy == 'Yes' & discrepancy_spec == 'Informative'){
        v_c_no_trt <- c(max(c(c_H + ID_cost,0)), c_S1, c_S2 + c_Trt, c_D)
      }
      else{
        v_c_no_trt <- c(c_H, c_S1, c_S2, c_D)
      }
    })

    # estimate mean QALYs and costs 
    v_E_no_trt <- reactive({
      v_E_no_trt <- m_TR() %*% v_u_no_trt()
    })
    
    v_E_trt <- reactive({
        v_E_trt <- m_TR() %*% v_u_trt()
    })
    
    v_C_no_trt <- reactive({
      v_C_no_trt <-m_TR() %*% v_c_no_trt()
    })
    
    v_C_trt <- reactive({
      v_C_trt <- m_TR() %*% v_c_trt()
    })
    
    te_no_trt <- reactive({
      t(v_E_no_trt()) %*% v_dwe  # 1x31 %*% 31x1 -> 1x1
    })
    
    te_trt <- reactive({
      te_trt <- t(v_E_trt()) %*% v_dwe
    })
    
    tc_no_trt <- reactive({
      t(v_C_no_trt()) %*% v_dwc
    })
    
    tc_trt <- reactive({
      t(v_C_trt()) %*% v_dwc
    })
    results <- reactive({
      c("NID_utility" = NID.utility,
        "NID_cost" = NID.cost,
        "NID_utility_trt" = NID.utility_trt,
        "NID_cost_trt" = NID.cost_trt,
        "Cost_NoTrt" = tc_no_trt(),
        "Cost_Trt"   = tc_trt(), 
        "QALY_NoTrt" = te_no_trt(), 
        "QALY_Trt"   = te_trt(),
        "ICER"       = (tc_trt() - tc_no_trt())/(te_trt() - te_no_trt()),
        "NMB_NoTrt"  = (te_no_trt()*WTP) - tc_no_trt(),
        "NMB_Trt"    = (te_trt()*WTP) - tc_trt(),
        "pCE"       = if((te_no_trt()*WTP) - tc_no_trt()<(te_trt()*WTP) - tc_trt()){
          1
        }
        else{
          0
        }
    )
    })
    return(results())
  }
  )
}
  
  #==============================================================================
  #                   Probabilistic Sensitivity Analysis
  #==============================================================================
  
## Draw random sample for PSA
df_psa <- gen_psa(n_sim = n_sim)
  

## Initialize matrix of results/outcomes
df_out <- matrix(NaN, nrow = n_sim, ncol = 12)
colnames(df_out) <- c("NID_utility","NID_cost",
                      "NID_utility_trt","NID_cost_trt",
                      "Cost_NoTrt", "Cost_Trt",
                      "QALY_NoTrt", "QALY_Trt",
                      "ICER", "NMB_NoTrt", "NMB_Trt", "pCE")

## Run PSA
withProgress(message = "Running PSA", value = 0,{  # display the progress of the simulation
for(i in 1:n_sim){
  df_out[i,] <- f_MM_sicksicker(df_psa[i,])
  cat('\r', paste(round(i/n_sim * 100), "% done", sep = " "))
  incProgress(1/n_sim, detail = paste("Iteration", i))
  Sys.sleep(.1)
  }
  })

df_out <- as.data.frame(df_out)

return(df_out)


}
