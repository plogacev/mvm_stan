
data
{
	int<lower=1> n_obs;
	int<lower=0,upper=1> iv_ev[n_obs];
	int<lower=0,upper=1> iv_yes[n_obs];
	int<lower=0,upper=1> iv_high_interference[n_obs];
	int<lower=0,upper=1> response_yes[n_obs];
}

parameters
{
	real R_low;
	real C_low;
	real G_low;
	real R_diff;
	real C_diff;
	real G_diff;
}

model
{
	int i;
	real prob_trial_type;
	real logProb_path;
	real R;
	real C;
	real G;
	real pY;
	
	for(i_obs in 1:n_obs)
	{
		
		if(iv_ev[i_obs] == 0)
		{
			real logLik[2];
			
			// path: start0-start-MV-NR;
			R <- R_low + R_diff * iv_high_interference[i_obs];
			C <- C_low + C_diff * iv_high_interference[i_obs];
			G <- G_low + G_diff * iv_high_interference[i_obs];
			pY <- G;
			logProb_path <- log((1 - R));
			logLik[1] <- logProb_path + bernoulli_log(response_yes[i_obs], pY);
			
			// path: start0-start-MV-R_MV;
			R <- R_low + R_diff * iv_high_interference[i_obs];
			C <- C_low + C_diff * iv_high_interference[i_obs];
			G <- G_low + G_diff * iv_high_interference[i_obs];
			pY <- (iv_yes[i_obs] == 1);
			logProb_path <- log((R));
			logLik[2] <- logProb_path + bernoulli_log(response_yes[i_obs], pY);
			
			increment_log_prob(log_sum_exp(logLik));
		};
		
		if(iv_ev[i_obs] == 1)
		{
			real logLik[3];
			
			// path: start0-start-EV-NR;
			R <- R_low + R_diff * iv_high_interference[i_obs];
			C <- C_low + C_diff * iv_high_interference[i_obs];
			G <- G_low + G_diff * iv_high_interference[i_obs];
			pY <- G;
			logProb_path <- log((1 - R));
			logLik[1] <- logProb_path + bernoulli_log(response_yes[i_obs], pY);
			
			// path: start0-start-EV-R_EV-CR_EV;
			R <- R_low + R_diff * iv_high_interference[i_obs];
			C <- C_low + C_diff * iv_high_interference[i_obs];
			G <- G_low + G_diff * iv_high_interference[i_obs];
			pY <- (iv_yes[i_obs] == 1);
			logProb_path <- log((R) * (C));
			logLik[2] <- logProb_path + bernoulli_log(response_yes[i_obs], pY);
			
			// path: start0-start-EV-R_EV-IR_EV;
			R <- R_low + R_diff * iv_high_interference[i_obs];
			C <- C_low + C_diff * iv_high_interference[i_obs];
			G <- G_low + G_diff * iv_high_interference[i_obs];
			pY <- (iv_yes[i_obs] == 0);
			logProb_path <- log((R) * (1 - C));
			logLik[3] <- logProb_path + bernoulli_log(response_yes[i_obs], pY);
			
			increment_log_prob(log_sum_exp(logLik));
		};
	};
}