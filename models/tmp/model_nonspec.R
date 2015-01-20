rm(list=ls())
detach("package:mvmstan", unload = TRUE) 
library(mvmstan)

# read the model
m <- read_mvm(sprintf("./model_nonspec.txt", dir))

# plot model without the variable assignment
plot_mvm(m, show_prob = T, show_code = F, fname_dot="./model_nonspec.dot", start_xdot = T)

# plot model with the variable assignment
plot_mvm(m, show_prob = T, show_code = T, fname_dot="./model_nonspec.dot", start_xdot = F)

  
# TODO: The checks fail if a data column is included which is not used in the model. This should be a warning, not a error.
# independent variables
iv_vars = c(iv_cond='int<lower=0,upper=2>',iv_questN1='int<lower=0,upper=1>',crit_region_cnt='int<lower=1>')

# dependent variables
dv_vars = c(response_yes='int<lower=0,upper=1>', reading_time='real<lower=0>', response_RT='real<lower=0>')

# parameters to be estimated
par_vars = c(p_uspec='real<lower=0,upper=1>', p_att_n1='real<lower=0,upper=1>',
             p_retrieval_fail='real<lower=0,upper=1>', p_guess_yes='real<lower=0,upper=1>',
             rate='real<lower=0>', # TODO: 'rate' was missing here, and the code was still generated. Fix this behavior.
             alpha0='real<lower=0>', alpha1='real<lower=0>', alpha2='real<lower=0>',
             beta0='real<lower=0>', beta1='real<lower=0>'
             )

# relationship between variables assigned in the model and the log-likelihood of the DVs
logLik = c(response='bernoulli_log(response_yes, pY)', reading_time='gamma_log(reading_time, reading_shape, rate)', 
           reading_time='gamma_log(response_RT, response_shape, rate)')

# parameters for which random effects will be specified
# NOTE: This is not really implemented yet
ranEff = c(p_uspec = 'inv_logit( logit(p_uspec) + subj)', p_att_n1 = 'inv_logit( logit(p_att_n1) + subj + item)')


mvm_generate_code(m, iv_vars=iv_vars, par_vars=par_vars,  dv_vars=dv_vars, logLik=logLik, raneff=ranEff,
                  file="./model_nonspec_auto.stan")



# load Swets' data
source("~/Documents/MyPapers/models_of_underspecification/Data_Swets_et_al/RScriptSwetsData.R")
d <- read.swets.data("~/Documents/MyPapers/models_of_underspecification/Data_Swets_et_al/Data_All.csv")

d = subset(d, qtype%in%c('RC questions'))
d = d[with(d, order(subj,trial,item,rid)),]
d.critical <- subset(d, rid=="reflexive") 
d.critical2 <- subset(d, rid=="region9") 
d.critical$RT = d.critical$RT + d.critical2$RT

# find participants with more than 35% correct answers in the unambiguous conditions
subject.acc = with(d.critical, tapply(response.yes==correct.response.yes, list(subj, attachment), mean ))
nmap(asc(d.critical$subj),  subject.acc[,'N1'], as.double) -> d.critical$acc.N1
nmap(asc(d.critical$subj),  subject.acc[,'N2'], as.double) -> d.critical$acc.N2

acc.min = .5
d.qrc = subset(d.critical, acc.N1 > acc.min & acc.N2 > acc.min)
d.excluded = unique( subset(d.critical, !(acc.N1 > acc.min & acc.N2 > acc.min))[,c('subj','qtype','acc.N1','acc.N2')] )

excluded.subjects.rc.N1 = unique(subset(d.excluded, acc.N1 > acc.min)$subj)
excluded.subjects.rc.N2 = unique(subset(d.excluded, acc.N2 > acc.min)$subj)

#d.qrc$attachment = ordered(d.qrc$attachment, levels=c("N1","N2","ambiguous"))
d.qrc$attachment.n2 = with(d.qrc, ifelse(attachment=="N2", 1, 0))
d.qrc$response.correct = NA
d.qrc$response.correct[d.qrc$attachment!="ambiguous"] = with(subset(d.qrc, attachment!="ambiguous"), ifelse(response.yes==correct.response.yes, 1, 0))
d.qrc$response.yes.str = nmap(d.qrc$response.yes, c('0'="`no'", '1'="`yes'"))
d.qrc$response.correct.str = nmap(asc(d.qrc$response.correct), c('1'='correct response', '0'='incorrect response'))
d.qrc$crit.region.cnt = d.qrc$pc.cnt + 1

d.qrc = subset(d.qrc, resp.RT < 15000 )




d <- subset(d.qrc, qtype=="RC questions")[,c('item','subj','trial')]
d$iv_questN1 <- as.integer(d.qrc$questionNP=="N1")
d$crit_region_cnt <- d.qrc$crit.region.cnt
d$item <- as.integer(as.factor(d.qrc$item))
d$subj <- as.integer(as.factor(d.qrc$subj))
d$response_yes <- d.qrc$response.yes
d$response_RT <- d.qrc$resp.RT
d$reading_time <- d.qrc$RT
d$iv_cond <- ifelse(d.qrc$attachment=="ambiguous", 0, ifelse(d.qrc$attachment=="N1", 1, 2))

d_lst <- c(c(n_obs=nrow(d), n_subj=length(unique(d$subj)), n_item=length(unique(d$item))), as.list(d))



library(rstan)

sm <- stan_model(model_name="sm1", model_code = read_file("./model_nonspec_auto.stan"))
sm <- stan_model(model_name="sm1", model_code = read_file("./model_nonspec_multi.stan"))

init = list(p_uspec=0.1, p_att_n1=0.5, p_retrieval_fail=0.3, p_guess_yes=0.5, rate=.1, alpha0=2, alpha1=2, alpha2=2, beta0=3, beta1=3)
res3 <- rstan::optimizing(sm, data=d_lst, init=init) # , algorithm="Newton"
res3

res3 <- rstan::optimizing(sm, data=d_lst)#, algorithm="Newton")
head(res3$par, 10)

tail(res3$par, 10)


plot(function(x) dgamma(x, shape =12.507e+00 , rate= 2.3890e-03), xlim=c(0,5000))

tapply(d$response_yes, list(d$iv_cond, d$iv_questN1), mean)

tapply(d$reading_time, list(d$mixclass, d$iv_amb), mean)

tapply(d$mixclass==1, list(d$iv_amb), mean)


###########################################################################
# 
  fit <- sampling(sm, data=d_lst, chains=0)
  fn <- function(par) log_prob(fit, par)
  gr <- function(par) grad_log_prob(fit, par)
  
  p <- c(0, 0, 0, 0, 0, 0)
  fn(p)
  gr(p)
  
  p2logodds <- function(p) log(p/(1-p))
  logodds2p <- function(lodds) exp(lodds)/(1+exp(lodds))
  
  res <- optim(p, fn, gr, method="SANN", control=list(fnscale=-1, maxit=10^3))
  res
#  
# logodds2p(res$par)
# 
# res2 <- optim(p, fn, gr, method="BFGS", control=list(fnscale=-1), hessian=TRUE)
# 
# logodds2p(res2$par)
# 
# sapply(seq(0,1,.05), function(p1) sapply(seq(0,1,.05), function(p2) fn(c(p1,2))))

