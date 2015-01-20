
source("./mpt_models.R")
dir <- "."

# read the model
m <- read_mpt(sprintf("./model1.txt", dir))

# plot model without the variable assignment
plot_mpt(m, show_prob = T, show_code = F, fname_dot="./model1.dot", start_xdot = F)

# plot model with the variable assignment
plot_mpt(m, show_prob = T, show_code = T, fname_dot="./model1.dot", start_xdot = F)



# independent variables
iv_vars = c(iv_amb='int<lower=0,upper=1>')

# dependent variables
dv_vars = c(response='int<lower=0,upper=1>', reading_time='real<lower=0>')

# parameters to be estimated
par_vars = c(p_yes_m1='real<lower=0,upper=1>', p_yes_m2='real<lower=0,upper=1>',
             p_m1_c1='real<lower=0,upper=1>', p_m1_c2='real<lower=0,upper=1>',
             mu11='real<lower=0>', mu12='real<lower=0>', sigma='real<lower=0>'
             )

# relationship between variables assigned in the model and the log-likelihood of the DVs
logLik = c(response='bernoulli_log(response, pY)', reading_time='lognormal_log(reading_time, mu1, sigma)')

# parameters for which random effects will be specified
# NOTE: This is not really implemented yet
ranEff = c(mu11='mu11 + ritem_mu11[item] + rsubj_mu11[subj]', mu12='mu12 + ritem_mu11[item] + rsubj_mu11[subj]')
  
yield = model_generate_code(m, iv_vars=iv_vars, par_vars=par_vars, 
                            dv_vars=dv_vars, logLik=logLik, ranEff=ranEff)

cat(yield)
cat( yield, file=sprintf("%s/model1_auto.stan", dir))




library(rstan)

sm <- stan_model(model_name="sm1", model_code = read_file("./model1_auto.stan"))

generate_data <- function(p_yes, mu1, mu2, mix, n) {
  mix = round(n*mix)
  response = c(runif(mix[1])<p_yes[1], runif(n-mix[1])<p_yes[2],
               runif(mix[2])<p_yes[1], runif(n-mix[2])<p_yes[2])
  reading_time = c(rlnorm(mix[1], meanlog=mu1[1]), rlnorm(n-mix[1], meanlog=mu1[2]), 
                   rlnorm(mix[2], meanlog=mu1[1]), rlnorm(n-mix[2], meanlog=mu1[2]))
  response_time = c(rlnorm(mix[1], meanlog=mu2[1]), rlnorm(n-mix[1], meanlog=mu2[2]), 
                    rlnorm(mix[2], meanlog=mu2[1]), rlnorm(n-mix[2], meanlog=mu2[2]))
  class = c(rep(1,mix[1]), rep(2,n-mix[1]), 
            rep(1,mix[2]), rep(2,n-mix[2]))
  data.frame(response = response, reading_time = reading_time, response_time = response_time,
             mixclass=class, iv_amb= rep(0:1, each=n))
}

d <- generate_data(p_yes=c(.2, .85), mu1=c(6, 6.1), mu2=c(2, 2.2), mix=c(.2, .8), n=1000)
d_lst <- c(c(n_obs=nrow(d), n_subj=1, n_item=1), as.list(d))


res3 <- rstan::optimizing(sm, data=d_lst)

res3


tapply(d$response, list(d$mixclass, d$iv_amb), mean)

tapply(d$reading_time, list(d$mixclass, d$iv_amb), mean)

tapply(d$mixclass==1, list(d$iv_amb), mean)


###########################################################################
# 
#  fit <- sampling(sm, data=d_lst, chains=0)
#  fn <- function(par) log_prob(fit, par)
#  gr <- function(par) grad_log_prob(fit, par)
#  
#  p <- c(0, 0, 0, 0, 0, 0, 0)
#  fn(p)
#  gr(p)
#  
#  p2logodds <- function(p) log(p/(1-p))
#  logodds2p <- function(lodds) exp(lodds)/(1+exp(lodds))
#  
#  res <- optim(p, fn, gr, method="Nelder-Mead", control=list(fnscale=-1, maxit=10^3))
#  res
#  
# logodds2p(res$par)
# 
# res2 <- optim(p, fn, gr, method="BFGS", control=list(fnscale=-1), hessian=TRUE)
# 
# logodds2p(res2$par)
# 
# sapply(seq(0,1,.05), function(p1) sapply(seq(0,1,.05), function(p2) fn(c(p1,2))))

