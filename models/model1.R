library(mvmstan)

# read the model specifications
m <- read_mvm("./model1.txt")

# plot the non-specification model without the variable assignment
plot_mvm(m, show_prob = T, show_code = F, fname_dot="./model1.dot", start_xdot = T)

# plot the partial specification model with the variable assignment
plot_mvm(m, show_prob = T, show_code = T, fname_dot="./model1.dot", start_xdot = F)



### Declare variables and log-likelihood functions to be used in the stan model

# independent variables
iv_vars = c(iv_ev='int<lower=0,upper=1>', iv_yes='int<lower=0,upper=1>')

# dependent variables
dv_vars = c(response_yes='int<lower=0,upper=1>')

# parameters to be estimated
par_vars = c(R='real', C='real', G='real',
             R_diff='real', C_diff='real', G_diff='real')

# relationship between variables assigned in the model and the log-likelihood of the DVs
logLik = c(response='bernoulli_log(response_yes, pY )')




### Generate Stan code for the model

# simple model (no random effects)
transforms = c(R = 'inv_logit( R + R_diff)', 
               C = 'inv_logit( C + C_diff)',
               G = 'inv_logit( C + G_diff)'
)
mvm_generate_code(m, iv_vars=iv_vars, par_vars=par_vars, dv_vars=dv_vars, logLik=logLik, raneff=transforms, file="./model1_simple.stan")

# now let's add random effects (all by-subject)
raneff = c(R = 'inv_logit( R + R_diff + subj)', 
           C = 'inv_logit( C + C_diff + subj)',
           G = 'inv_logit( G + G_diff + subj)'
)

# model with random effects
mvm_generate_code(m, iv_vars=iv_vars, par_vars=par_vars,  dv_vars=dv_vars, logLik=logLik, raneff=raneff, file="./model1_ranef.stan")



### Compile Stan models

library(rstan)

m_simple <- stan_model(model_name="nonspec_simple", model_code = read_file("./model1_simple.stan"))
m_ranef <- stan_model(model_name="nonspec_ranef", model_code = read_file("./model1_ranef.stan"))


### Optimize Stan models (in order to get all parameter names, and just for fun)

res_simple <- rstan::optimizing(m_simple, data=d_lst)
res_simple

res_ranef <- rstan::optimizing(m_ranef, data=d_lst)
res_ranef


### Fit Stan models

# select parameters of interest
pars = c(names(res_simple$par), names(res_ranef$par[grep("_sd", names(res_ranef$par))]))
pars

fit_ranef <- sampling(m_ranef, data=d_lst, verbose=T, pars=pars, chains=1, iter=500)
