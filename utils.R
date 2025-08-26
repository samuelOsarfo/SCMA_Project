## -------------------------------------------------------------------
## 1.  Positive mediators  (R = 1)  –– li1_1
## -------------------------------------------------------------------
# par_vec : initial parameter vector  [gamma, theta, beta_k, alpha_k, tau_k, zeta_k, b_k,
#                                     omega0_k1, omega1_k1, omega2_k1 (vector),
#                                     omega0_k2, omega1_k2, omega2_k2 (vector),
#                                     phi]


li1_1_vec <- function(par_vec,
                      t_vec, m_k_vec, x_vec,
                      conf_matrix, d_vec) {
  
  n <- length(t_vec)              
 
  
  vapply(seq_len(n), function(i)
    li1_1(par_vec,
          ti        = t_vec[i],
          m_ik      = m_k_vec[i],
          xi        = x_vec[i],
          conf_vec  = if (is.null(conf_matrix)) NULL else conf_matrix[i, , drop=FALSE],
          di        = d_vec[i]),
    numeric(1L))                            # FUN.VALUE: one number back
}



## -------------------------------------------------------------------
## 2.  True-zero mediators (R = 0) –– li0_2
## -------------------------------------------------------------------


li0_2_vec <- function(par_vec,
                      t_vec, x_vec,
                      conf_matrix, d_vec) {
  
  n <- length(t_vec)
  vapply(seq_len(n), function(i)
    li0_2(par_vec,
          ti        = t_vec[i],
          xi        = x_vec[i],
          conf_vec  = if (is.null(conf_matrix)) NULL else conf_matrix[i, , drop=FALSE],
          di        = d_vec[i]),
    numeric(1L))
}


## -------------------------------------------------------------------
## 3.  Undetected-zero mediators (R = 0) –– li1_2
## -------------------------------------------------------------------


li1_2_vec <- function(par_vec,
                      t_vec, x_vec,
                      conf_matrix, d_vec, s_vec) {
  
  n <- length(t_vec)
  vapply(seq_len(n), function(i)
    li1_2(par_vec,
          ti        = t_vec[i],
          xi        = x_vec[i],
          conf_vec  = if (is.null(conf_matrix)) NULL else conf_matrix[i, , drop=FALSE],
          di        = d_vec[i],
          si        = s_vec[i]),
    numeric(1L))
}




unpack_par <-function(par_vec, n_conf=0){
  list(
    gamma = par_vec[1],
    theta = if(n_conf>0)  par_vec[2:(1 + n_conf)] else numeric(0),
    beta_k = par_vec[2 + n_conf],
    alpha_k = par_vec[3 + n_conf],
    tau_k = par_vec[4 + n_conf],
    zeta_k = par_vec[5 + n_conf],
    b_k = par_vec[6 + n_conf],
    
    # zero-inflation parameters (Equation 15)
    omega0_k1 = par_vec[7 + n_conf],
    omega1_k1 = par_vec[8 + n_conf],
    omega2_k1 = if(n_conf >0) par_vec[9 + n_conf + (0:(n_conf-1))] else  numeric(0),
    
    # non-zero mean parameters (Equation 16)
    omega0_k2 = par_vec[9 + 2*n_conf],
    omega1_k2 = par_vec[10 + 2*n_conf ],
    omega2_k2= if(n_conf>0) par_vec[(11 + 2*n_conf):(10+ 3*n_conf)] else numeric(0),
    phi = par_vec[11 + 3*n_conf]
  )
}


# non-zero uniform values generator over an interval (-a, b) where b>0.
non_zero_unif <- function(n, min, max){
  sapply(1:n, function (x){
    val <-0
    while(val==0) val <-runif(1, min, max)
    return(val)
  })
}  



################################################################################################
#####################likelihood function for non-zero observed mediators########################
################################################################################################

#Inputs
# par_vec : initial parameter vector  [gamma, theta, beta_k, alpha_k, tau_k, zeta_k, b_k,
#                                     omega0_k1, omega1_k1, omega2_k1 (vector),
#                                     omega0_k2, omega1_k2, omega2_k2 (vector),
#                                     phi]
# m_ik : observed relative abundance of taxon k for subject i
# ti : time-to-event value
# xi : exposure/treatment value
# conf_vec :  vector of confounder values
# di: event indicator
# si: sequencing depth
#
#output
# # log_lik : Log-likelihood contribution for subject i and taxon k



li1_1 <-function(par_vec, ti, m_ik, xi, conf_vec=NULL,di){
  
  n_conf <- length(conf_vec)
  p <- unpack_par(par_vec, n_conf)
  
  ###  calculate mediator parameters
  # zero-inflation probability (logistic) :: (Equations 15-16)
  logit_delta <- p$omega0_k1 + p$omega1_k1 * xi
  if(n_conf > 0) logit_delta <- logit_delta + sum(p$omega2_k1 * conf_vec)
  delta_ik <- plogis(logit_delta)
  
  
  # non-zero mean (beta regression)
  logit_mu <- p$omega0_k2 + p$omega1_k2 * xi
  if(n_conf > 0) logit_mu <- logit_mu + sum(p$omega2_k2 * conf_vec)
  mu_ik <- plogis(logit_mu)
  
  
  
  # Accelerated Failure Time model (Equation 11)
  yi <- log(ti)
  outcome_model <- p$gamma * xi + p$beta_k * m_ik + p$alpha_k * 1 + p$tau_k * xi * m_ik + p$zeta_k * xi * 1
  if(n_conf > 0) outcome_model <- outcome_model + sum(p$theta * conf_vec)
  
  epsilon_ik <- (yi - outcome_model) / p$b_k
  
  
  
  # survival component of the AFT log-normal
  if (di == 1) {
    # observed event
    log_surv <- -log(p$b_k) + dnorm(epsilon_ik, log = TRUE)
  } else {
    # censored
    log_surv <- pnorm(epsilon_ik, lower.tail = FALSE, log.p = TRUE)
  }
  
  # mediator component :: Zero-inflated Beta
  # Beta density (Equation 13)
  a <- mu_ik * p$phi
  b <- (1 - mu_ik) * p$phi
  log_beta_density <- (a - 1)*log(m_ik) + (b - 1)*log(1 - m_ik) - lbeta(a, b)
  
  # Account for zero-inflation (Equation 14)
  log_med <- log(1 - delta_ik) + log_beta_density
  
  
  #Total log-likelihood = Survival + Mediator
  log_lik <- log_surv + log_med
  log_lik
  
  
}




################################################################################################
#####################likelihood function for truly-zero observed mediators######################
################################################################################################

li0_2 <- function(par_vec, ti, xi, conf_vec=NULL,di) {
  # get parameters
  n_conf <- length(conf_vec)
  p <- unpack_par(par_vec, n_conf)
  
  
  
  
  # calculate Delta_ik (prob of true zero)
  logit_delta <- p$omega0_k1 + p$omega1_k1 * xi
  if(n_conf > 0) logit_delta <- logit_delta + sum(p$omega2_k1 * conf_vec)
  delta_ik <- plogis(logit_delta)
  
  
  
  # true zero component (c_i = 0)
  yi <- log(ti)
  outcome_model <- p$gamma * xi
  if(n_conf > 0) outcome_model <- outcome_model + sum(p$theta * conf_vec)
  
  epsilon_ik <- (yi - outcome_model) / p$b_k
  
  
  #survival componetn
  if (di == 1) {
    log_surv <- -log(p$b_k) + dnorm(epsilon_ik, log = TRUE)
  } else {
    log_surv <- pnorm(epsilon_ik, lower.tail = FALSE, log.p = TRUE)
  }
  
  # mediator component: log(P(true zero)) = log(Δ_ik)
  log_med <- log(delta_ik)
  
  
  # total log-likelihood
  log_lik <- log_surv + log_med
  log_lik
}






################################################################################################
#####################likelihood function for undetected-zero observed mediators#################
################################################################################################

li1_2 <- function(par_vec, ti, xi, conf_vec=NULL, di, si) {
  # get parameters
  n_conf <- length(conf_vec)
  p <- unpack_par(par_vec, n_conf)
  
  #  calculate mediator parameters (plogis for expit)
  logit_delta <- p$omega0_k1 + p$omega1_k1 * xi
  if(n_conf > 0) logit_delta <- logit_delta + sum(p$omega2_k1 * conf_vec)
  delta_ik <- plogis(logit_delta)
  
  logit_mu <- p$omega0_k2 + p$omega1_k2 * xi
  if(n_conf > 0) logit_mu <- logit_mu + sum(p$omega2_k2 * conf_vec)
  mu_ik <- plogis(logit_mu)
  
  
  
  # integrand function for m in (0, 1/S_i)
  integrand <- function(m) {
    # AFT outcome model
    outcome_model <- p$gamma * xi + p$beta_k * m + p$alpha_k * 1 +
      p$tau_k * xi * m + p$zeta_k * xi * 1
    if(n_conf > 0) outcome_model <- outcome_model + sum(p$theta * conf_vec)
    
    
    
    yi <- log(ti)
    epsilon_ik <- (yi - outcome_model) / p$b_k
    
    # survival component
    if (di == 1) {
      # observed event
      surv_part <- (1/p$b_k) * dnorm(epsilon_ik)
    } else {
      # censored
      surv_part <- pnorm(epsilon_ik, lower.tail = FALSE)
    }
    
    
    surv_part
  }
  
  # numerical integration
  upper_lim <- 1 / si  # is it divided by m*si ?
  int_result <- tryCatch({
    integrate(integrand, lower = 0, upper = upper_lim,
              rel.tol = 1e-6, subdivisions = 1000)$value
  },
  error = function(e) {
    warning("Integration failed for si=", si, ": ", e$message)
    1e-100  # fallback value
  }
  )
  
  # Beta density for mediator
  beta_part <- beta(mu_ik * p$phi, (1 - mu_ik) * p$phi)
  
  # likelihood component
  joint_density <- (1 - delta_ik) *(1/beta_part) *int_result
  
  log(max(joint_density, 1e-100))
}






### ----  E–STEP  ------------------------------------------------------------


#How do I think of the data? Let's get all pieces for the E-step
# Inputs
# par_vec : vector of initial parameters
# x : vector of exposures
# m_k : vector of values for mediator M_k
# conf_mat : matrix of confounders
# t : vector
# d : vector
# s : vector
# Output:
#


#-------------------------------------------------------------------------------
# E-STEP: Compute posterior weights
#-------------------------------------------------------------------------------
Estep <- function(par_vec, x, m_k, d, t, s, conf_mat=NULL) {
  n_conf <- if (is.null(conf_mat)) 0 else ncol(conf_mat)
  p <- unpack_par(par_vec, n_conf)
  
  
  n <- length(m_k)
  g1 <- which(m_k != 0)
  g2 <- which(m_k == 0)
  
  
  # Compute pi (P(Ci=0) using logistic model
  logit_pi <- p$omega0_k1 + p$omega1_k1 * x
  if (n_conf > 0) logit_pi <- logit_pi + conf_mat %*% p$omega2_k1
  pi_vec <- plogis(logit_pi)
  
  
  eta0 <- numeric(n)
  eta0[g1] <- 0
  
  
  if(length(g2) > 0) {
    # compute log-likelihood components
    l2i0 <- li0_2_vec(par_vec, t[g2], x[g2], if (!is.null(conf_mat)) conf_mat[g2, , drop = FALSE] else NULL, d[g2])
    l2i1 <- li1_2_vec(par_vec, t[g2], x[g2], if (!is.null(conf_mat)) conf_mat[g2, , drop = FALSE] else NULL, d[g2], s[g2])
    
    # compute posterior weights
    numerator <- pi_vec[g2] * exp(l2i0)
    denominator <- numerator + (1 - pi_vec[g2]) * exp(l2i1)
    eta0[g2] <- numerator / denominator
  }
  
  eta0  # return P(Ci=0|data) for all subjects
}





###### Q- function############
## I read that Optim() function search for a minimum: using negative -likelihood function.
Q_func <- function(par_vec, x, m_k, d, t, s, conf_mat=NULL, eta0) {
  n_conf <- if (is.null(conf_mat)) 0 else ncol(conf_mat)
  p <- unpack_par(par_vec, n_conf)
  
  
  n <- length(m_k)
  # group-1 indices
  g1 <- which(m_k !=0)
  
  # group-2 indices
  g2  <- which(m_k == 0)
  
  
  
  ## pi's given theta
  logit_pi <- p$omega0_k1 + p$omega1_k1 * x
  if (n_conf > 0) logit_pi <- logit_pi + conf_mat %*% p$omega2_k1
  pi_vec <- plogis(logit_pi)
  
  
  ## group-1 contribution to the likelihood function
  ll_g1 <-  log(1 - pi_vec[g1]) + li1_1_vec(par_vec, t[g1], m_k[g1], x[g1],  if (!is.null(conf_mat)) conf_mat[g1, , drop = FALSE] else NULL, d[g1])
  
  
  ## group-2 contribution to the likelihood function
  l2i0 <- li0_2_vec(par_vec, t[g2], x[g2],  if (!is.null(conf_mat)) conf_mat[g2, , drop = FALSE] else NULL, d[g2])
  l2i1 <- li1_2_vec(par_vec, t[g2], x[g2],  if (!is.null(conf_mat)) conf_mat[g2, , drop = FALSE] else NULL, d[g2], s[g2])
  
  
  # get posterior weights
  eta0_g2 <- eta0[g2]
  eta1_g2 <- 1 - eta0_g2
  
  
  
  ll_g2 <- eta0_g2 * (log(pi_vec[g2]) + l2i0) + eta1_g2 * (log(1 - pi_vec[g2]) + l2i1)
  
  
  -(sum(ll_g1) + sum(ll_g2))         # negative Q for optim()
}


#-------------------------------------------------------------------------------
#  Observed log-lokelihood value
#-------------------------------------------------------------------------------
observed_loglik <- function(par_vec, t_vec, m_k_vec, x_vec, d_vec, s_vec, conf_matrix = NULL) {
  # Pre-calculate log-likelihoods for all subjects for each possible case
  # Note: For a given subject, only one of these will be used based on m_k_vec[i]
  
  # Likelihood for observed positive mediators (m_k > 0)
  loglik_pos <- li1_1_vec(par_vec, t_vec, m_k_vec, x_vec, conf_matrix, d_vec)
  
  # Likelihood for true-zero mediators (contribution if C_i=0)
  loglik_truezero <- li0_2_vec(par_vec, t_vec, x_vec, conf_matrix, d_vec)
  
  # Likelihood for undetected-zero mediators (contribution if C_i=1)
  loglik_undetected <- li1_2_vec(par_vec, t_vec, x_vec, conf_matrix, d_vec, s_vec)
  
  
  # Initialize the total log-likelihood
  total_loglik <- 0
  n <- length(t_vec)
  
  for (i in 1:n) {
    if (m_k_vec[i] > 0) {
      # Case 2a: Observed positive -> definitely C_i=1
      total_loglik <- total_loglik + loglik_pos[i]
      
    } else if (m_k_vec[i] == 0) {
      # Could be either true zero (C_i=0) or undetected (C_i=1)
      # We need to sum the likelihoods from both latent classes
      # log( P(C=0)*f(data|C=0) + P(C=1)*f(data|C=1) )
      #    = log( exp(loglik_truezero[i]) + exp(loglik_undetected[i]) )
      
      # Use log-sum-exp for numerical stability
      max_val <- max(loglik_truezero[i], loglik_undetected[i])
      total_loglik <- total_loglik + max_val + log(exp(loglik_truezero[i] - max_val) + exp(loglik_undetected[i] - max_val))
    }
  }
  
  return(total_loglik)
}



#-------------------------------------------------------------------------------
# MAIN EM ALGORITHM
#-------------------------------------------------------------------------------
EM_algorithm <- function(par_init, x, m_k, d, t, s, conf_mat=NULL,  tol = 1e-6, max_iter = 100) {
  # Initialize
  par_current <- par_init
  iter <- 0
  hess_final   <- NULL
  converged <- FALSE
  prev_loglik <- -Inf
  
  
  while(!converged && iter < max_iter) {
    iter <- iter + 1
    
    # E-step: Compute posterior weights
    eta0 <- Estep(par_current, x, m_k, d, t, s, conf_mat)
    
    # M-step: Maximize Q-function
    opt_result <- optim(
      par = par_current,
      fn = Q_func,
      method = "L-BFGS-B",
      x = x,
      m_k = m_k,
      d = d,
      t = t,
      s = s,
      conf_mat = conf_mat,
      eta0 = eta0, # Use FIXED weights from E-step
      hessian = TRUE
    ) #we can add
    
    # Update parameters
    par_current <- opt_result$par
    
    current_obs_loglik <- observed_loglik(par_current, t, m_k, x, d, s, conf_mat)   
    hess_final  <- opt_result$hessian
    
  
    # Check convergence (relative log-likelihood change)
    if(iter > 1 && abs(current_obs_loglik - prev_obs_loglik) < tol) {
      converged <- TRUE
    }
    prev_obs_loglik <- current_obs_loglik
    
    # Optional: Print progress
    cat(sprintf("Iteration %d: log-likelihood = %.4f\n", iter, current_obs_loglik))
  }
  
  # Return results
  list(
    par = par_current,
    loglik = current_obs_loglik,
    iterations = iter,
    hessian    = hess_final,
    converged = converged
  )
}


#-------------------------------------
#  Get STD_Errors for Parameter Estimates
#------------------------------------


#EM_result <- EM_algorithm(par_vec, x, m_k, d, t, s)


# # After EM converges:
# EM_stderror <-function(em_result){
#   
#   hessian <- em_result$hessian  # Extract Hessian (this is already -d²Q/dΘ²)
#   
#   
#   cov_matrix <- ginv(hessian)   # Invert to get covariance matrix
#   
#   sqrt(diag(cov_matrix))
#   
# }
# 
# 
# stdError <-EM_stderror(EM_result)

#-------------------------------------
#  NIE Estimation
#-------------------------------------

expit<- function(x){1/(1 + exp(-x))}

NIE_func <- function(par, x1 = 0, x2 = 1) {
  
  beta_k <- par[2]   # beta
  tau_k  <- par[4]   # tau
  alpha_k <- par[3]  # alpha
  zeta_k <- par[5]   # zeta
  omega0_1 <- par[7] # omega_0_1
  omega1_1 <- par[8] # omega_1_1
  omega0_2 <- par[9] # omega_0_2
  omega1_2 <- par[10] # omega_1_2

  
 
  
  #----------------------------------------------------------------------
  # Calculate Components for NIE_1 (Equation 28) :: There are no covariates for now
  #----------------------------------------------------------------------
  
  NIE_1 <- (beta_k + tau_k*x2)*( (expit(omega0_2 + omega1_2*x2)* (1- expit(omega0_1 + omega1_1*x2))) - (expit(omega0_2 + omega1_2*x1)* (1- expit(omega0_1 + omega1_1*x1))) )
  
  NIE_11 <-(beta_k + tau_k*x2)*( (expit(omega0_2 + omega1_2*x2) - expit(omega0_2 + omega1_2*x1)) - (expit(omega0_2 + omega1_2*x2)*expit(omega0_1 + omega1_1*x2) -expit(omega0_2 + omega1_2*x1)*expit(omega0_1 + omega1_1*x1)) )
  #----------------------------------------------------------------------
  # Calculate NIE_2 (Equation 29) :: There are no covariates for now
  #----------------------------------------------------------------------
  # Effect through change in presence/absence probability
  NIE_2 <- (alpha_k + zeta_k * x2) * (expit(omega0_2 + omega1_2*x1) - expit(omega0_2 +omega1_2*x2))
  
  #----------------------------------------------------------------------
  # Total Conditional NIE (Sum of both pathways)
  #----------------------------------------------------------------------
  NIE_total <- NIE_1 + NIE_2
  print(NIE_1)
  print(NIE_11)
  print(NIE_2)
  return(NIE_total)
}


#NIE <- NIE_func(EM_result$par, x1 = 0, x2 = 1)

#-------------------------------
# Std_Error for NIE From the Delta Method
#------------------------------




# m_k <-as.vector(out$dat$M[,10])
# par_vec <-as.vector(unname(out$par_mat[10,]))
# EM_result <- EM_algorithm(par_vec, x, m_k, d, t, s)
# stdError <-EM_stderror(EM_result)
# NIE <- NIE_func(EM_result$par, x1 = 0, x2 = 1)
# se.nie <-SE_nie(EM_result$par, EM_result$hessian, NIE_func, x1 = 0, x2 = 1, alpha = 0.05)
# se.nie







SE_nie <- function(par_est, nie_function, x1 = 0, x2 = 1, alpha = 0.05,cov_mat) {
  # Step b: Compute the Gradient (G) of the NIE function
  #         with respect to the parameters, at the estimated values.
  gradient <- grad(func = nie_function, x = par_est, x1 = x1, x2 = x2)
  
  # Step c: Compute the Covariance Matrix and then the Variance of the NIE
  #         The Hessian from the M-step is -d²Q/dΘ², which is the observed Fisher information.
  var_nie <- t(gradient) %*% cov_mat %*% gradient # G' * Cov(Θ) * G
  se_nie <- sqrt(var_nie) # Standard error of the NIE
  
  # Step d: Compute the Point Estimate, Confidence Interval, and P-value
  nie_estimate <- nie_function(par_est, x1 = x1, x2 = x2) # Point estimate
  
  # Wald test statistic
  z_value <- nie_estimate / se_nie
  
  # Two-sided p-value
  p_value <- 2 * pnorm(-abs(z_value))
  
  # Critical value for (1-alpha/2) confidence (e.g., 1.96 for 95% CI)
  z_critical <- qnorm(1 - alpha/2)
  
  # Confidence Interval
  ci_lower <- nie_estimate - z_critical * se_nie
  ci_upper <- nie_estimate + z_critical * se_nie
  
  # Step e: Organize and Return Results
  results <- data.frame(
    Estimate = nie_estimate,
    SE = as.numeric(se_nie), # Convert 1x1 matrix to numeric
    CI_lower = ci_lower,
    CI_upper = ci_upper,
    z_value = z_value,
    p_value = p_value
  )
  
  
  # Return the results invisibly (so they can be stored in a variable)
  return(results)
}

