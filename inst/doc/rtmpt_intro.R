## ---- eval=TRUE----------------------------------------------------------
eqn = "
# CORE MODEL
## tree ; cat ;          mpt
      0 ;   0 ;           Do
      0 ;   0 ;     (1-Do)*g
      0 ;   1 ; (1-Do)*(1-g)

      1 ;   3 ;           Dn
      1 ;   2 ;     (1-Dn)*g
      1 ;   3 ; (1-Dn)*(1-g)

# OPTIONAL RESTRICTIONS / SPECIFICATIONS
const_prob: g=0.5

##   Tree ; Cat ; Resp
resp:   0 ;   0 ; 0
resp:   0 ;   1 ; 1
resp:   1 ;   2 ; 0
resp:   1 ;   3 ; 1
"

## ---- eval=TRUE----------------------------------------------------------
library(rtmpt)
restr_2HTM <- to_rtmpt_model(eqn_file = eqn)
restr_2HTM

## ---- eval=TRUE----------------------------------------------------------
crazy_model <- restr_2HTM
# change parameters of the model:
crazy_model <- set_params(model = crazy_model, parameter = "probs", 
                            names = "g", values = NA)  # g will now be estimated again
crazy_model <- set_params(model = crazy_model, parameter = "tau_minus",
                            names = "Dn", values = 0) # suppress process-completion time for a 
                                                      # failure to detect a new item as "new"
# change responses:
crazy_model <- set_resps(model = crazy_model, tree = 0, categories = 1, values = 0)
crazy_model <- set_resps(model = crazy_model, tree = 1, categories = 3, values = 0)
crazy_model

## ---- eval=TRUE----------------------------------------------------------
set.seed(2018)
raw_data <- data.frame(tree = rep(1:2, 8), rt = round(1000*(runif(16)+.3)), 
                       group = rep(1:2, each=8), subj = rep(1:4, each = 4), cat = rep(1:4, 4))
raw_data
data <- to_rtmpt_data(raw_data = raw_data, model = restr_2HTM)
data

## ---- eval=FALSE---------------------------------------------------------
#  ## do not run
#  rtmpt_out <- fit_rtmpt(model = restr_2HTM, data = data)
#  rtmpt_out$diags$R_hat
#  ## end not run

## ---- eval=FALSE---------------------------------------------------------
#  ## do not run
#  rtmpt_out <- fit_rtmpt(model = restr_2HTM,
#                         data = data,
#                         n.chains = 4,
#                         n.iter = 5000,
#                         n.burnin = 200,
#                         n.thin = 1,
#                         Rhat_max = 1.05,
#  					   Irep = 1000,
#                         prior_params = NULL,
#                         indices = FALSE,
#                         save_log_lik = FALSE)
#  ## end not run

## ---- eval=FALSE---------------------------------------------------------
#  eqn = "
#  # CORE MODEL
#  ## tree ; cat ;          mpt
#        0 ;   0 ;           Do
#        0 ;   0 ;     (1-Do)*g
#        0 ;   1 ; (1-Do)*(1-g)
#  
#        1 ;   3 ;           Dn
#        1 ;   2 ;     (1-Dn)*g
#        1 ;   3 ; (1-Dn)*(1-g)
#  
#  # OPTIONAL RESTRICTIONS / SPECIFICATIONS
#  const_prob: g=0.5, Dn = .5
#  suppress_tau: Dn-, Do+
#  
#  ##   Tree ; Cat ; Resp
#  resp:   0 ;   0 ; 0
#  resp:   0 ;   1 ; 1
#  resp:   1 ;   2 ; 0
#  resp:   1 ;   3 ; 1
#  "

## ---- eval=FALSE---------------------------------------------------------
#  eqn = "
#  # CORE MODEL
#  ## tree ; cat ;          mpt
#        0 ;   0 ;           Do
#        0 ;   0 ;     (1-Do)*g
#        0 ;   1 ; (1-Do)*(1-g)
#  
#        1 ;   3 ;           Dn
#        1 ;   2 ;     (1-Dn)*g
#        1 ;   3 ; (1-Dn)*(1-g)
#  "

## ---- eval=FALSE---------------------------------------------------------
#  mdl = "
#  # CORE MODEL
#  ## targets:
#  Do+(1-Do)*g
#  (1-Do)*(1-g)
#  
#  ## lure
#  (1-Dn)*g
#  Dn+(1-Dn)*(1-g)
#  "

## ---- eval=TRUE----------------------------------------------------------
mdl = "
# CORE MODEL
## MDL         ; RESP

### targets
Do+(1-Do)*g    ;    0
(1-Do)*(1-g)   ;    1

### lure
(1-Dn)*g       ;    0
Dn+(1-Dn)*(1-g);    1

# OPTIONAL RESTRICTIONS
const_prob: g=0.5, Dn = .5
suppress_tau: Dn-, Do+
"
to_rtmpt_model(mdl_file = mdl)

## ---- eval=TRUE----------------------------------------------------------
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


## ---- eval=TRUE----------------------------------------------------------
R <- 2000
rand_rankmat <- matrix(data = sample(0:99, R*393, replace=TRUE), nrow = R, ncol = 393)

## pearsons' chi-square for testing uniformity
x <- apply(rand_rankmat[1:R,], 2, table)
expect <- R/100 # 100 = number of bins/cells (0:99)
pearson <- apply(X = x, MARGIN = 2, FUN = function(x) {sum((x-expect)^2/expect)})
z95 <- qchisq(0.95, 99) # 99 = degrees of freedom
sum(pearson>z95) / length(pearson)

