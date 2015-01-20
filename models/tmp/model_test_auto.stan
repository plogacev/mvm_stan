data
{
	int<lower=1> n_obs;
	int<lower=0,upper=1> iv_amb[n_obs];
	int<lower=0,upper=1> response[n_obs];
}
parameters
{
	real<lower=0,upper=1> p_yes;
	real<lower=0,upper=1> p_err;
}
model
{
	int i;
	real prob_trial_type;
	real logProb_path;
	real pY;
	real pX;
	int i_obs;

print("lp 1",  get_lp());

//	for(i_obs in 1:n_obs)
	{
	i_obs <- 1;
		if(iv_amb[i_obs] == 0 && p_err == 0)
		{
			real logLik[1];
			
			// path: start-NU-CR
			logProb_path <- 0;
			pY <- 0.7;
			logLik[1] <- logProb_path + bernoulli_log(response[i_obs], pY);
			print("ll 1 ", log_sum_exp(logLik));
			increment_log_prob(log_sum_exp(logLik));
		};
		if(iv_amb[i_obs] == 0 && p_err == 1)
		{
			real logLik[1];
			
			// path: start-NU-guess
			logProb_path <- 0;
			pY <- p_yes; 
			logLik[1] <- logProb_path + bernoulli_log(response[i_obs], pY);
			print("ll 2 ", log_sum_exp(logLik));
			increment_log_prob(log_sum_exp(logLik));
		};
		if(iv_amb[i_obs] == 1)
		{
			real logLik[2];
			
			// path: start-U-guess
			logProb_path <- log((0.5));
			pY <- p_yes;
			logLik[1] <- logProb_path + bernoulli_log(response[i_obs], pY);
			
			// path: start-U-guess2
			logProb_path <- log((0.5));
			pY <- p_yes/2;
			logLik[2] <- logProb_path + bernoulli_log(response[i_obs], pY);
			print("ll 3 ", log_sum_exp(logLik));
	print("lp 2a ",  get_lp());
			increment_log_prob(log_sum_exp(logLik));
	print("lp 2b ",  get_lp());
		};
	};
	print("lp 3 ",  get_lp());
}