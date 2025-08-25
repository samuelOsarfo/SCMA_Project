
rm(list = ls())

require(doParallel)
require(foreach)
require(MASS)


num_cores <- 11

cl <- parallel::makeCluster(num_cores)   # start workers
doParallel::registerDoParallel(cl)  # register with foreach

on.exit({ stopCluster(cl); registerDoSEQ() }, add = TRUE)


foreach::getDoParWorkers()             # check number of workers



seednum <- 2025 # set the number for set.seed()

num_sim <- 5 # the total number of independent simulation runs

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
s12 =as.numeric(args[4])
s13 =as.numeric(args[5])
b = as.numeric(args[6])
censoring_prop = as.numeric(args[7])



n<-25                # sample size
k <- 30              # total number of mediators
s11 <-3              # number of true mediators under scenario 1
s12 <-3              # number of true mediators under scenario 2
s13 <-4              # number of true mediators under scenario 1
b <- 0.8             # scale parameter for log-weibull distribution (in AFT model)
censoring_prop <- round(runif(1, 0.1, 0.25),2)     # censored proportion



# 
# # -------------------------------
# # Choosing Parameters to Establish The 3 Scenarios Below
# # -------------------------------
# 
# ## X -> M_k (relative abundance) ->T
# ra_indices <-1:3
# omega_1_2[ra_indices]<- non_zero_unif(n=3, min=-2.5, max=2.5)
# beta[ra_indices]<-non_zero_unif(n=3, -0.7, 0.7)
# tau[ra_indices]<- runif(3, -0.2, 0.2)
# omega_0_1[ra_indices] <- runif(3, -1.5, 0.5) 
# omega_0_2[ra_indices] <- runif(3, -1.5, -0.5)
# 
# 
# ## X -> I(M_k >0) (presence)-> T
# pres_indices <-4:6
# omega_1_1[pres_indices]<-non_zero_unif(3, -2.5, -2)
# alpha[pres_indices]<- non_zero_unif(3, -0.7, 0.7)
# zeta[pres_indices] <-runif(3, -0.3, 0.3)
# omega_0_1[pres_indices] <- runif(3, 2, 3)
# 
# 
# ## X affects both abundance and presence the both affect T 
# both_indices <-7:10
# omega_1_1[both_indices] <- runif(4, -2.5, -1.5 )
# omega_1_2[both_indices] <- runif(4, -1.5, 1.5)
# alpha[both_indices] <- non_zero_unif(4, -0.4, 0.8)
# beta[both_indices] <- non_zero_unif(4, -0.4, 0.4)
# zeta[both_indices]  <- runif(4, -0.3, 0.3)
# tau[both_indices] <- runif(4, -0.4, 0.4)
# omega_0_1[both_indices] <- runif(4, 0.5, 1.5) 
# omega_0_2[both_indices] <- runif(4, -1, 0.5)
# 


my_path <-"C:/Users/fsosa/OneDrive/Masters_Course_Materials/Thesis/Papers/Mediation Analysis_With Survival/SCMA_Project/"

file_name = sprintf("%ssim_res_%s_%s_%s_%s_%s_%s_%s.RData", my_path, n, k, s11, s12, s13, b, censoring_prop)



data_step_func <- function(n, k, s11, s12, s13, b, censoring_prop, seednum){
  
  set.seed(seednum)
  s <- s11 + s12 + s13  # total number of true mediators
  
  
  # Baseling Parameters for ZIB Mediator model for m_k=0
  baseline_prevalence <- rbeta(k, shape1 = 1, shape2 = 10)
  omega_0_1 <- log( (1 - baseline_prevalence) / baseline_prevalence ) #intercept
  omega_1_1 <- rep(0, k)  #X coefficient
  
  
  
  # Baseline Parameters for ZIB Mediator model for m_k>0
  baseline_abundance <- rbeta(k, shape1 = 1, shape2 = 5)
  omega_0_2 <- log( baseline_abundance / (1 - baseline_abundance) ) # intercept for abundance part
  omega_1_2 <- rep(0, k)  # effect of X on log-odds of mean abundance
  
  
  
  # Baseline Parameters for AFT Model
  gamma <- 0.5             # Direct effect of X on log(T)
  beta <- rep(0, k)     # effect of mediator M_k on log(T)
  alpha <- rep(0, k)    # effect of presence I(M_k>0) on log(T)
  tau <-rep(0, k)       # interaction coeff. for X*m_k
  zeta <- rep(0, k)     # interaction coeff.  for X*presence
  
  
  
  # -------------------------------
  # Choosing Parameters to Establish The 3 Scenarios Below
  # -------------------------------
  
  ## X -> M_k (relative abundance) ->T
  ra_indices <-1:s11
  omega_1_2[ra_indices]<- non_zero_unif(s11, min=-2.5, max=2.5)
  beta[ra_indices]<-non_zero_unif(s11, -0.7, 0.7)
  tau[ra_indices]<- runif(s11, -0.2, 0.2)
  omega_0_1[ra_indices] <- non_zero_unif(s11, -1.5, 0.5) 
  omega_0_2[ra_indices] <- runif(s11, -1.5, -0.5)
  
  
  ## X -> I(M_k >0) (presence)-> T
  pres_indices <-(s11+1):(s11+s12)
  omega_1_1[pres_indices]<-non_zero_unif(s12, -2.5, -2)
  alpha[pres_indices]<- non_zero_unif(s12, -0.7, 0.7)
  zeta[pres_indices] <-runif(s12, -0.3, 0.3)
  omega_0_1[pres_indices] <- runif(s12, 2, 3)
  
  
  ## X affects both abundance and presence the both affect T 
  both_indices <-(s11+s12+1): s
  omega_1_1[both_indices] <- runif(s13, -2.5, -1.5 )
  omega_1_2[both_indices] <- runif(s13, -1.5, 1.5)
  alpha[both_indices] <- non_zero_unif(s13, -0.4, 0.8)
  beta[both_indices] <- non_zero_unif(s13, -0.4, 0.4)
  zeta[both_indices]  <- runif(s13, -0.3, 0.3)
  tau[both_indices] <- runif(s13, -0.4, 0.4)
  omega_0_1[both_indices] <- runif(s13, 0.5, 1.5) 
  omega_0_2[both_indices] <- non_zero_unif(s13, -1, 0.5)
  
  
  
  X <- matrix(rbinom(n, 1, 0.6), nrow = n, ncol = 1)
  phi <- 10   # dispersion parameter for beta distribution
  
  
  dat <- gen_data(
    X,
    # gen_meds parameters
    omega_0_1, omega_1_1,        # zero part: logit P(M=0)
    omega_0_2, omega_1_2,        # positive part: logit E[M|M>0]
    phi,                         # precision (scalar)
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
  par_mat <-data.frame(gamma=rep(gamma, k), beta, alpha, tau, zeta, b=rep( b,k), omega_0_1, omega_1_1, omega_0_2, omega_1_2, phi=rep( phi, k))
  
  par_mat <-as.matrix(par_mat)
  
  list(dat=dat, 
       par_mat=par_mat,                       
       censoring_prop=censoring_prop
       )
}


## Add another to call the EM_algorithm, get the standard errors, get the NIE estimates from the Delta method




all_res <- list()

all_res <- foreach(i = 1:num_sim, .packages = c( "foreach", "doParallel")) %dopar% {
  tryCatch(data_step_func(n, k, s11, s12, s13, b, censoring_prop,seednum + i), error=function(e) NA)
}




save(all_res, file = file_name)

if (exists("cl") && inherits(cl, "cluster")) parallel::stopCluster(cl)

