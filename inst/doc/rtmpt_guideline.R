## ----eval=TRUE----------------------------------------------------------------
eqn = "
# CORE MODEL
## EQN SYNTAX

## tree ; cat ;          mpt
      0 ;   0 ;           Do
      0 ;   0 ;     (1-Do)*g
      0 ;   1 ; (1-Do)*(1-g)

      1 ;   3 ;           Dn
      1 ;   2 ;     (1-Dn)*g
      1 ;   3 ; (1-Dn)*(1-g)
"

## ----eval=TRUE----------------------------------------------------------------
mdl = "
# CORE MODEL
## MDL SYNTAX

### targets
Do+(1-Do)*g
(1-Do)*(1-g)

### lure
(1-Dn)*g
Dn+(1-Dn)*(1-g)
"

## ----eval=TRUE----------------------------------------------------------------
library(rtmpt)
# MAKING AN ERTMPT MODEL
#   using the MDL syntax:
e2HTM <- to_ertmpt_model(mdl_file = mdl)
#   using the EQN syntax:
e2HTM <- to_ertmpt_model(eqn_file = eqn)
e2HTM

# MAKING A DRTMPT MODEL
#   using the MDL syntax:
d2HTM <- to_drtmpt_model(mdl_file = mdl)
#   using the EQN syntax:
d2HTM <- to_drtmpt_model(eqn_file = eqn)
d2HTM
# 

## ----eval=TRUE----------------------------------------------------------------
theta2const(model = e2HTM, names = "g", constants = 0.5)

## ----eval=TRUE----------------------------------------------------------------
theta2theta(model = e2HTM, names = c("Do", "Dn"), keep_consts = FALSE)

## ----eval=TRUE----------------------------------------------------------------
tau2zero(model = e2HTM, names = "g", outcomes = "minus", values = 0)

## ----eval=TRUE----------------------------------------------------------------
tau2tau(model = e2HTM, names = c("Do", "Dn"), outcome = "minus", keep_zeros = FALSE)

## ----eval=TRUE----------------------------------------------------------------
a2const(model = d2HTM, names = "g", constants = 1)

## ----eval=TRUE----------------------------------------------------------------
a2a(model = d2HTM, names = c("Do", "Dn"), keep_consts = FALSE)

## ----eval=TRUE----------------------------------------------------------------
nu2const(model = d2HTM, names = "g", constants = 1)

## ----eval=TRUE----------------------------------------------------------------
nu2nu(model = d2HTM, names = c("Do", "Dn"), keep_consts = FALSE)

## ----eval=TRUE----------------------------------------------------------------
omega2const(model = d2HTM, names = "g", constants = 0.5)

## ----eval=TRUE----------------------------------------------------------------
omega2omega(model = d2HTM, names = c("Do", "Dn"), keep_consts = FALSE)

## ----eval=TRUE----------------------------------------------------------------
delta2delta(model = e2HTM, trees = c(0,1), categories = c(1,3), mappings = c(1,1))

## ----eval=TRUE----------------------------------------------------------------
set.seed(2021)
raw_data <- data.frame(tree = rep(1:2, 8), rt = round(1000*(runif(16)+.3)), 
                       group = rep(1:2, each=8), subj = rep(1:4, each = 4), cat = rep(1:4, 4))
raw_data
data <- to_ertmpt_data(raw_data = raw_data, model = e2HTM)
data

## ----eval=FALSE---------------------------------------------------------------
#  ## do not run
#  ertmpt_out <- fit_ertmpt(model = e2HTM, data = data)
#  ## end not run

## ----eval=FALSE---------------------------------------------------------------
#  ## do not run
#  ertmpt_out <- fit_ertmpt(model = restr_2HTM,
#                           data = data,
#                           n.chains = 4,
#                           n.iter = 5000,
#                           n.burnin = 200,
#                           n.thin = 1,
#                           Rhat_max = 1.05,
#  					               Irep = 1000,
#                           prior_params = NULL,
#                           indices = FALSE,
#                           save_log_lik = FALSE)
#  ## end not run

## ----eval=FALSE---------------------------------------------------------------
#  ## do not run
#  drtmpt_out <- fit_drtmpt(model = d2HTM, data = data)
#  ## end not run

## ----eval=FALSE---------------------------------------------------------------
#  ## do not run
#  drtmpt_out <- fit_drtmpt(model = restr_2HTM,
#                           data = data,
#                           n.chains = 4,
#                           n.iter = 1000,
#                           n.phase1 = 1000,
#                           n.phase1 = 2000,
#                           n.thin = 1,
#                           Rhat_max = 1.05,
#  					               Irep = 1000,
#                           prior_params = NULL,
#                           flags = NULL,
#  					               control = NULL)
#  ## end not run

## ----eval=FALSE---------------------------------------------------------------
#  ertmpt_out$diags$R_hat
#  coda::traceplot(ertmpt_out$samples[, 1:9])
#  summary(ertmpt_out)

## ----eval=FALSE---------------------------------------------------------------
#  drtmpt_out$diags$R_hat
#  coda::traceplot(drtmpt_out$samples[, 1:9])
#  summary(drtmpt_out)

## ----eval=FALSE---------------------------------------------------------------
#  eqn = "
#  # CORE MODEL
#  ## tree ; cat ;         path
#        0 ;   0 ;           Do
#        0 ;   0 ;     (1-Do)*g
#        0 ;   1 ; (1-Do)*(1-g)
#  
#        1 ;   3 ;           Dn
#        1 ;   2 ;     (1-Dn)*g
#        1 ;   3 ; (1-Dn)*(1-g)
#  "

## ----eval=FALSE---------------------------------------------------------------
#  mdl = "
#  # CORE MODEL
#  ## targets:
#  Do+(1-Do)*g
#  (1-Do)*(1-g)
#  
#  ## lures
#  (1-Dn)*g
#  Dn+(1-Dn)*(1-g)
#  "

## ----eval=TRUE----------------------------------------------------------------
# randomly drawn group-level mean values
mdl_2HTM <- "
# targets
do+(1-do)*g     ; 0
(1-do)*(1-g)    ; 1

# lures
(1-dn)*g        ; 0
dn+(1-dn)*(1-g) ; 1

# do: detect old; dn: detect new; g: guess
"

model <- to_ertmpt_model(mdl_file = mdl_2HTM)

# random group-level parameters
params <- list(mean_of_mu_alpha = 0, 
               var_of_mu_alpha = 1,         # delete this line to fix mean_of_mu_alpha to 0
               mean_of_exp_mu_beta = 10, 
               var_of_exp_mu_beta = 10,     # delete this line to fix mean_of_exp_mu_beta to 10
               mean_of_mu_gamma = 0.5, 
               var_of_mu_gamma = 0.0025,    # delete this line to fix mean_of_mu_gamma to 0.5
               mean_of_omega_sqr = 0.005, 
               var_of_omega_sqr = 0.000025, # delete this line to fix mean_of_omega_sqr to 0.005
               df_of_sigma_sqr = 10, 
               sf_of_scale_matrix_SIGMA = 0.1, 
               sf_of_scale_matrix_GAMMA = 0.01, 
               prec_epsilon = 10,
               add_df_to_invWish = 5)

sim_dat <- sim_ertmpt_data(model, seed = 123, n.subj = 40, n.trials = 30, params = params)
head(sim_dat$data_frame)


## ----eval=TRUE----------------------------------------------------------------
R <- 2000
rand_rankmat <- matrix(data = sample(0:99, R*393, replace=TRUE), nrow = R, ncol = 393)

## pearsons' chi-square for testing uniformity
x <- apply(rand_rankmat[1:R,], 2, table)
expect <- R/100 # 100 = number of bins/cells (0:99)
pearson <- apply(X = x, MARGIN = 2, FUN = function(x) {sum((x-expect)^2/expect)})
z95 <- qchisq(0.95, 99) # 99 = degrees of freedom
sum(pearson>z95) / length(pearson)

