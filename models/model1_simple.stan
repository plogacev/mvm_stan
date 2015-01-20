
data
{
	int<lower=1> n_obs;
	int<lower=0,upper=1> iv_ev[n_obs];
	int<lower=0,upper=1> iv_yes[n_obs];
	int<lower=0,upper=1> response_yes[n_obs];
}

parameters
{
	real<lower=0,upper=1> R;
	real<lower=0,upper=1> C;
	real<lower=0,upper=1> G;
}

model
{
	int i;
	real prob_trial_type;
	real logProb_path;
	real pY;
	
	for(i_obs in 1:n_obs)
	{
		
		if(iv_ev[i_obs] == 0)
		{
			real logLik[2];
			
			// path: start-MV-NR;
			logProb_path <- log((1 - R));
			pY <- G;
			logLik[1] <- logProb_path + bernoulli_log(response_yes[i_obs], pY);
			
			// path: start-MV-R_MV;
			logProb_path <- log((R));
			pY <- (iv_yes[i_obs] == 1);
			logLik[2] <- logProb_path + bernoulli_log(response_yes[i_obs], pY);
			
			increment_log_prob(log_sum_exp(logLik));
		};
		
		if(iv_ev[i_obs] == 1)
		{
			real logLik[3];
			
			// path: start-EV-NR;
			logProb_path <- log((1 - R));
			pY <- G;
			logLik[1] <- logProb_path + bernoulli_log(response_yes[i_obs], pY);
			
			// path: start-EV-R_EV-CR_EV;
			logProb_path <- log((R) * (C));
			pY <- (iv_yes[i_obs] == 1);
			logLik[2] <- logProb_path + bernoulli_log(response_yes[i_obs], pY);
			
			// path: start-EV-R_EV-IR_EV;
			logProb_path <- log((R) * (1 - C));
			pY <- (iv_yes[i_obs] == 0);
			logLik[3] <- logProb_path + bernoulli_log(response_yes[i_obs], pY);
			
			increment_log_prob(log_sum_exp(logLik));
		};
	};
}