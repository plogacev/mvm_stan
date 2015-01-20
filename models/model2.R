detach("package:mvmstan", unload = T)
library(mvmstan)

# read the model specifications
m <- read_mvm("./model2.txt")

# plot the non-specification model without the variable assignment
plot_mvm(m, show_prob = T, show_code = F, fname_dot="./model2.dot", start_xdot = T)

# plot the partial specification model with the variable assignment
plot_mvm(m, show_prob = T, show_code = T, fname_dot="./model2_complete.dot", start_xdot = T)



### Declare variables and log-likelihood functions to be used in the stan model

# independent variables
iv_vars = c(iv_ev='int<lower=0,upper=1>', iv_yes='int<lower=0,upper=1>', 
            iv_high_interference='int<lower=0,upper=1>')

# dependent variables
dv_vars = c(response_yes='int<lower=0,upper=1>', RT='int<lower=0>')

# parameters to be estimated
par_vars = c(R='real<lower=0,upper=1>', C='real<lower=0,upper=1>', G='real<lower=0,upper=1>',
             R_diff='real', C_diff='real', G_diff='real',
             muB='real<lower=0>', sigmaB='real<lower=0>', muR='real<lower=0>', sigmaR='real<lower=0>')

# relationship between variables assigned in the model and the log-likelihood of the DVs
logLik = c(response='bernoulli_log(response_yes, pY )', response='lognormal_log(RT, mu, sigma)')


### Generate Stan code for the model

# simple model (no random effects)
# NOTE: modifications of path probabilities have to be handled via the code originally intended for random effects # TODO: Think of a cleaner solution
transforms = c(R = 'inv_logit( R + R_diff)', C = 'inv_logit( C + C_diff)', G = 'inv_logit( G + G_diff)')
mvm_generate_code(m, iv_vars=iv_vars, par_vars=par_vars, dv_vars=dv_vars, logLik=logLik, raneff=transforms, file="./model2_simple.stan")


# now let's add random effects (all by-subject)
transforms2_ranef = c(R_diff = 'R_diff + subj', C_diff = 'C + C_diff', G_diff = 'G + G_diff + subj',
                      R = 'inv_logit( R + R_diff + subj)', C = 'inv_logit( C + C_diff + subj)', G = 'inv_logit( G + G_diff + subj)',
                      mu='exp(log(mu) + subj)', sigma='exp(log(sigma) + subj)', muR='exp(log(muR) + subj)', sigmaR='exp(log(sigmaR) + subj)'
)


# model with random effects
mvm_generate_code(m, iv_vars=iv_vars, par_vars=par_vars,  dv_vars=dv_vars, logLik=logLik, raneff=transforms2_ranef, file="./model2_ranef.stan")



### Compile Stan models

library(rstan)

m_simple <- stan_model(model_name="nonspec_simple", model_code = read_file("./model2_simple.stan"))
m_ranef <- stan_model(model_name="nonspec_ranef", model_code = read_file("./model2_ranef.stan"))


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
