
rm(list = ls())

require(doParallel)
require(foreach)
require(MASS)
require(numDeriv)
require(MCMCpack)
require(truncdist)


num_cores <- 11

cl <- parallel::makeCluster(num_cores)   # start workers
doParallel::registerDoParallel(cl)  # register with foreach

on.exit({ stopCluster(cl); registerDoSEQ() }, add = TRUE)


foreach::getDoParWorkers()             # check number of workers



seednum <- 2025 # set the number for set.seed()

num_sim <- 50 # the total number of independent simulation runs

#my_path <- "/home/sosarfo/hdjmt/" # the directory in the server where the source code is located

#my_path <- 'C://Users//fsosa//OneDrive//Masters_Course_Materials//Thesis//Code//Sim 2//'

source("utils.R")
source("simulation_data.R")



## Parameter initialization
# get command line arguments
args = commandArgs(trailingOnly = TRUE)
n = as.numeric(args[1])
k = as.numeric(args[2])
s11 =as.numeric(args[3])
#s12 =as.numeric(args[4])
#s13 =as.numeric(args[5])
#b = as.numeric(args[6])
censoring_prop = as.numeric(args[7])



n<-200                # sample size
k <- 800        # total number of mediators
s11 <-5              # number of true mediators under scenario 1
#s12 <-3              # number of true mediators under scenario 2
#s13 <-4              # number of true mediators under scenario 1
b <- 0.8             # scale parameter for log-weibull distribution (in AFT model)
censoring_prop <- round(runif(1, 0.1, 0.25),2)     # censored proportion




my_path <-"C:/Users/fsosa/OneDrive/Masters_Course_Materials/Thesis/Papers/Mediation Analysis_With Survival/SCMA_Project/"

file_name = sprintf("%ssim_res_%s_%s_%s_%s_%s_%s_%s.RData", my_path, n, k, s11, s12, s13, b, censoring_prop)



data_step_func <- function(n, k, s11, b, censoring_prop, seednum){
  
  set.seed(seednum)

  #k<-800
  # Baseling Parameters for ZIB Mediator model for m_k=0
  omega_0_1 <- rep(0, k) #intercept
  omega_1_1 <- rep(0, k)  #X coefficient
  
  
  
  # Baseline Parameters for ZIB Mediator model for m_k>0
  omega_0_2 <- rep(0, k) # intercept for abundance part
  omega_1_2 <- rep(0, k)  # effect of X on log-odds of mean abundance
  
  
  
  # Baseline Parameters for AFT Model
  gamma <- 0.5             # Direct effect of X on log(T)
  beta <- rep(0, k)     # effect of mediator M_k on log(T)
  alpha <- rep(0, k)    # effect of presence I(M_k>0) on log(T)
  tau <-rep(0, k)       # interaction coeff. for X*m_k
  zeta <- rep(0, k)     # interaction coeff.  for X*presence
  
  
  
  # -------------------------------
  # Choosing Parameters to Test Indirect Effect
  # -------------------------------
  alpha[0:s11]     <-  -1.5
  zeta[0:s11]      <- -1.5
  
  
  beta[0:s11] <-  2
  tau[0:s11] <-   2
  
  #set.seed(1)
  omega_0_1[0:k] <- qlogis(0.45) #seem to work max{k<-6500}
  omega_1_1[0:k] <- qlogis(0.55)-qlogis(0.45)

  
  # omega_0_1[0:k] <- qlogis(0.75)
  # omega_1_1[0:k] <- qlogis(0.55)-qlogis(0.75)
  # 
  
  # omega_0_2[0:k] <- qlogis(0.02)
  # omega_1_2[0:k] <- qlogis(0.03)-qlogis(0.02)
  
  omega_0_2[0:k] <- qlogis(c(rep(0.01, round(k*0.05)), rep(0.001, round(k*0.20)), rep(0.00015, round(k*0.25)), rep(0.00015, round(k*0.5))))
  omega_1_2[0:k] <- qlogis(c(rep(0.02, round(k*0.05)), rep(0.002, round(k*0.20)), rep(0.0002, round(k*0.25)), rep(0.0002, round(k*0.5)))) - omega_0_2
  
  
  # omega_0_2[0:k] <- qlogis(0.00015)
  # omega_1_2[0:k] <- qlogis(0.00025)-qlogis(0.00015)
  
  # omega_0_2[0:k] <- qlogis(0.00025)
  # omega_1_2[0:k] <- qlogis(0.00015)-qlogis(0.00025)
  
  # omega_0_2[0:k] <- qlogis(0.00015)
  # omega_1_2[0:k] <- sample(c(qlogis(0.00025)-qlogis(0.00015), -qlogis(0.00025)+qlogis(0.00015)), k,replace = T )

  
  
  
  X <- matrix(rbinom(n, 1, 0.5), nrow = n, ncol = 1)
  #phi <- 5  # dispersion parameter for beta distribution
  
  
  dat <- gen_data(
    X,
    # gen_meds parameters
    omega_0_1, omega_1_1,        # zero part: logit P(M=0)
    omega_0_2, omega_1_2,        # positive part: logit E[M|M>0]
    # gen_T parameters
    gamma,                       # direct X effect
    beta,                        # length k (abundance effects)
    alpha,                       # length k (presence effects)
    tau,                         # length k (X×abundance), default 0
    zeta,                        # length k (X×presence),  default 0
    b,                           # AFT error scale
    # censoring
    censoring_prop
  ) 
  
  #parameter matrix
  par_mat <-data.frame(gamma=rep(gamma, k), beta, alpha, tau, zeta, b=rep( b,k), omega_0_1, omega_1_1, omega_0_2, omega_1_2)
  
  par_mat <-as.matrix(par_mat)
  
  list(dat=dat, 
       par_mat=par_mat,                       
       censoring_prop=censoring_prop
       )
}


## Add another function to call the EM_algorithm, get the standard errors, get the NIE estimates from the Delta method
## using parallel compute.



all_res <- list()

all_res <- foreach(i = 1:num_sim, .packages = c( "foreach", "doParallel")) %dopar% {
  tryCatch(data_step_func(n, k, s11, b, censoring_prop,seednum + i), error=function(e) NA)
}


## testing out some outs
out <-all_res[[1]]

m_k <-as.vector(out$dat$M[,2])
x<-as.vector(out$dat$X)
d<-as.vector(out$dat$d)
t<-as.vector(out$dat$T_true)
s<-as.vector(out$dat$s)
par_vec <-as.vector(unname(out$par_mat[2,]))
EM_result <- EM_algorithm(par_vec, x, m_k, d, t, s)


# 2. Compute the OBSERVED-DATA Hessian at the final parameter estimates
H_observed <- hessian(func = observed_loglik, x=EM_result$par,t_vec = t, m_k_vec = m_k, x_vec = x, d_vec = d, s_vec = s)

# 3. Invert the observed Hessian to get the CORRECT covariance matrix
#    Note: The observed Hessian is -d²L/dΘ², so we invert -H_observed to get Cov(Θ)
#cov_matrix <- ginv(-H_observed)

cov_matrix <- ginv(EM_result$hessian)

NIE <- NIE_func(EM_result$par, x1 = 0, x2 = 1)
se.nie <-SE_nie(EM_result$par, NIE_func, x1 = 0, x2 = 1, alpha = 0.05, cov_matrix)
se.nie





#save(all_res, file = file_name)

if (exists("cl") && inherits(cl, "cluster")) parallel::stopCluster(cl)

