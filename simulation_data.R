#---------------------------------
#  Generate Mediator Matrix
#---------------------------------


gen_meds <- function(X,
                     omega_0_1, omega_1_1,     # zero part: logit P(M=0)
                     omega_0_2, omega_1_2,     # positive part: logit E[M | M>0]
                     phi  # precision : scalar
                     ) { 
  n<- length(X)
  k<- length(omega_0_1)
  phi <- rep(phi, k)
  
  
  # ---- allocate ----
  delta <- matrix(0, n, k)   # P(M=0)
  mu    <- matrix(0, n, k)   # E[M | M>0]
  Ipos  <- matrix(0L, n, k)         # presence indicator: 1 if present, 0 if zero
  M     <- matrix(0, n, k)          # mediator values in [0,1]
  
  
  # ---- parameter matrices (k x 2), then transpose to 2 x k ----
  params_1 <- cbind(omega_0_1, omega_1_1) # zero part
  params_2 <- cbind(omega_0_2, omega_1_2) # positive part
  
  
  
  # ---- linear predictors ----
  Xmat  <- cbind(1, X)              # n x 2 (intercept + X)
  eta_1 <- Xmat %*% t(params_1)     # n x k
  eta_2 <- Xmat %*% t(params_2)     # n x k
  delta <- plogis(eta_1)            # P(M=0) in (0,1)
  mu    <- plogis(eta_2)            # mean in (0,1) for M|M>0
  
  
  # ---- simulate presence ----
  Ipos[] <- rbinom(n * k, size = 1, prob = 1 - delta)  # 1 = present
  
  # ---- simulate positive abundances ----
  for (j in seq_len(k)) {
    idx <- which(Ipos[, j] == 1)
    if (length(idx)) {
      shape1 <- mu[idx, j] * phi[j]
      shape2 <- (1 - mu[idx, j]) * phi[j]
      M[idx, j] <- rbeta(length(idx), shape1 = shape1, shape2 = shape2)
    }
  }
  
  # return matrices and vectors
  list(M = M, Ipos = Ipos, mu = mu, delta = delta)
}









#--------------------------
# Generate T
#--------------------------
gen_T <- function(X,
                  M,
                  Ipos,
                  gamma,
                  beta, 
                  alpha,
                  tau = NULL, 
                  zeta= NULL,
                  b) {
  
  n <- nrow(M); k <- ncol(M)
  if (is.null(tau))  tau  <- rep(0, k)
  if (is.null(zeta)) zeta <- rep(0, k)
  
  
  
  LT <- gamma * X +  M%*% beta + Ipos %*% alpha + (M* as.vector(X)) %*% tau + (Ipos*as.vector(X)) %*% zeta
  
  
  log_T  <- LT + b* rnorm(n)
  T_true <- exp(log_T)
  
  
  list(T_true = T_true, log_T = log_T)
}



#---------------------------
# Generate censoring status
#---------------------------

gen_censoring_status <-function(T_true, censoring_prop){
  
  # Generate censoring times to achieve desired censoring proportion
  # A common way is to generate from an exponential distribution (Is this correct?)
  C_time <- matrix(rexp(n, rate = 1/quantile(T_true, 1-censoring_prop)), nrow = n, ncol = 1)
  
  obs_time <- pmin(T_true, C_time)
  
  d <- matrix((T_true <= C_time) + 0, nrow = n, ncol = 1) # Event indicator
  
  return(d)
}


#-----------------------------------------------------
# Generating Sequencing Depth :: Using the Log-normal
#-----------------------------------------------------

# Draw depths with a log-normal; specify median and IQR (Q1,Q3) for interpretability.
sdepth_lognorm <- function(n, median, q1 = NULL, q3 = NULL, sdlog = NULL) {

    if (!is.null(sdlog)) {
    mu <- log(median)              # median = exp(mu)
    sigma <- sdlog
  } else {
    stopifnot(!is.null(q1), !is.null(q3), q1 > 0, q3 > q1)
    # For log-normal: Q3/Q1 = exp(1.349 * sigma)
    sigma <- log(q3 / q1) / 1.349
    # median = exp(mu) = sqrt(Q1*Q3)
    mu <- log(sqrt(q1 * q3))
  }
  pmax(1L, as.integer(round(rlnorm(n, meanlog = mu, sdlog = sigma))))
}





#-------------------------------------
# Generating Complete Simulation Data
#-------------------------------------


gen_data <- function(
    X,
    # gen_meds parameters
    omega_0_1, omega_1_1,        # zero part: logit P(M=0)
    omega_0_2, omega_1_2,        # positive part: logit E[M|M>0]
    phi,                         # precision (scalar)
    # gen_T parameters
    gamma,                       # direct X effect
    beta,                        # length k (abundance effects)
    alpha,                       # length k (presence effects)
    tau   = NULL,                # length k (X×abundance), default 0
    zeta  = NULL,                # length k (X×presence),  default 0
    b,                   # AFT error scale
    # censoring
    censoring_prop = 0.15
) {
  

  
  n <- length(X)
  k <- length(omega_0_1)
  
  
  if (is.null(tau))  tau  <- rep(0, k)
  if (is.null(zeta)) zeta <- rep(0, k)

  
  # --- mediators ---
  med <- gen_meds(X = X, omega_0_1 = omega_0_1, omega_1_1 = omega_1_1, omega_0_2 = omega_0_2, omega_1_2 = omega_1_2,
    phi = phi)
  
  M    <- med$M
  Ipos <- med$Ipos
  mu   <- med$mu
  delta<- med$delta
  
  
  # ---  Times-event (AFT) ---
  T_out <- gen_T(
    X = X,
    M = M,
    Ipos = Ipos,
    gamma = gamma,
    beta = beta, 
    alpha = alpha,
    tau = tau, 
    zeta = zeta,
    b = b)
  
  T_true <- as.numeric(T_out$T_true)
  log_T  <- as.numeric(T_out$log_T)
  
  # --- censoring  ---
  d <- gen_censoring_status(T_true, censoring_prop)
  
  
  # --- Sequencing depth-------
  # Examples
  #S_16S     <- sdepth_lognorm(n = 300, median = 20000, q1 = 12000, q3 = 40000)
  #S_shotgun <- rdepth_lognorm(n = 300, median = 5e6,  q1 = 2e6,   q3 = 1e7,    seed = 1)
  
  s <- sdepth_lognorm(n, median = 20000, q1 = 12000, q3 = 40000)
  
  # return data parts
  list(
    # sizes
    n = n, k = k,
    
    #exposure
    X=X,
    
    # mediators
    M = M, Ipos = Ipos, mu = mu, delta = delta,
    
    # event times
    T_true = T_true, log_T = log_T,
    
    # censoring
    d = d,
    
    # sequencing depth
    s=s
  )
}



