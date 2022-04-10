## ---- eval=TRUE---------------------------------------------------------------
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

## ---- eval=TRUE---------------------------------------------------------------
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

## ---- eval=TRUE---------------------------------------------------------------
library(rtmpt)
# using the MDL syntax:
TwoHTM <- to_rtmpt_model(mdl_file = mdl)
# using the EQN syntax:
TwoHTM <- to_rtmpt_model(eqn_file = eqn)
TwoHTM
# 

## ---- eval=TRUE---------------------------------------------------------------
theta2const(model = TwoHTM, names = "g", constants = 0.5)

## ---- eval=TRUE---------------------------------------------------------------
theta2theta(model = TwoHTM, names = c("Do", "Dn"), keep_consts = FALSE)

## ---- eval=TRUE---------------------------------------------------------------
tau2zero(model = TwoHTM, names = "g", outcomes = "minus", values = 0)

## ---- eval=TRUE---------------------------------------------------------------
tau2tau(model = TwoHTM, names = c("Do", "Dn"), keep_zeros = FALSE)

## ---- eval=TRUE---------------------------------------------------------------
delta2delta(model = TwoHTM, trees = c(0,1), categories = c(1,3), mappings = c(1,1))

## ---- eval=TRUE---------------------------------------------------------------
set.seed(2021)
raw_data <- data.frame(tree = rep(1:2, 8), rt = round(1000*(runif(16)+.3)), 
                       group = rep(1:2, each=8), subj = rep(1:4, each = 4), cat = rep(1:4, 4))
raw_data
data <- to_rtmpt_data(raw_data = raw_data, model = TwoHTM)
data

## ---- eval=FALSE--------------------------------------------------------------
#  ## do not run
#  rtmpt_out <- fit_rtmpt(model = TwoHTM, data = data)
#  ## end not run

## ---- eval=FALSE--------------------------------------------------------------
#  ## do not run
#  rtmpt_out <- fit_rtmpt(model = restr_2HTM,
#                         data = data,
#                         n.chains = 4,
#                         n.iter = 5000,
#                         n.burnin = 200,
#                         n.thin = 1,
#                         Rhat_max = 1.05,
#  					             Irep = 1000,
#                         prior_params = NULL,
#                         indices = FALSE,
#                         save_log_lik = FALSE)
#  ## end not run

## ---- eval=FALSE--------------------------------------------------------------
#  rtmpt_out$diags$R_hat
#  coda::traceplot(rtmpt_out$samples[, 1:9])
#  summary(rtmpt_out)

## ---- eval=FALSE--------------------------------------------------------------
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

## ---- eval=FALSE--------------------------------------------------------------
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

## ---- eval=TRUE---------------------------------------------------------------
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

model <- to_rtmpt_model(mdl_file = mdl_2HTM)

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

sim_dat <- sim_rtmpt_data(model, seed = 123, n.subj = 40, n.trials = 30, params = params)
head(sim_dat$data_frame)


## ---- eval=TRUE---------------------------------------------------------------
R <- 2000
rand_rankmat <- matrix(data = sample(0:99, R*393, replace=TRUE), nrow = R, ncol = 393)

## pearsons' chi-square for testing uniformity
x <- apply(rand_rankmat[1:R,], 2, table)
expect <- R/100 # 100 = number of bins/cells (0:99)
pearson <- apply(X = x, MARGIN = 2, FUN = function(x) {sum((x-expect)^2/expect)})
z95 <- qchisq(0.95, 99) # 99 = degrees of freedom
sum(pearson>z95) / length(pearson)

