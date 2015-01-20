
data
{
	int<lower=1> n_obs;
	int<lower=1> n_subj;
	int<lower=0,upper=1> iv_ev[n_obs];
	int<lower=0,upper=1> iv_yes[n_obs];
	int<lower=0,upper=1> response_yes[n_obs];
	int<lower=1,upper=n_subj> subj[n_obs];
}

parameters
{
	real R;
	real C;
	real G;
	real R_diff;
	real C_diff;
	real G_diff;
	vector[n_subj] R_subj;
	real<lower=0> R_subj_sd;
	vector[n_subj] C_subj;
	real<lower=0> C_subj_sd;
	vector[n_subj] G_subj;
	real<lower=0> G_subj_sd;
}

model
{
	int i;
	real prob_trial_type;
	real logProb_path;
	real pY;
	real cur_R;
	real cur_C;
	real cur_G;
	int cur_subj;
	
	for(i_obs in 1:n_obs)
	{
		cur_subj <- subj[i_obs];
		cur_R <- inv_logit(R + R_diff + R_subj[cur_subj]);
		cur_C <- inv_logit(C + C_diff + C_subj[cur_subj]);
		cur_G <- inv_logit(G + G_diff + G_subj[cur_subj]);
		
		if(iv_ev[i_obs] == 0)
		{
			real logLik[2];
			
			// path: start-MV-NR;
			logProb_path <- log((1 - cur_R));
			pY <- cur_G;
			logLik[1] <- logProb_path + bernoulli_log(response_yes[i_obs], pY);
			
			// path: start-MV-R_MV;
			logProb_path <- log((cur_R));
			pY <- (iv_yes[i_obs] == 1);
			logLik[2] <- logProb_path + bernoulli_log(response_yes[i_obs], pY);
			
			increment_log_prob(log_sum_exp(logLik));
		};
		
		if(iv_ev[i_obs] == 1)
		{
			real logLik[3];
			
			// path: start-EV-NR;
			logProb_path <- log((1 - cur_R));
			pY <- cur_G;
			logLik[1] <- logProb_path + bernoulli_log(response_yes[i_obs], pY);
			
			// path: start-EV-R_EV-CR_EV;
			logProb_path <- log((cur_R) * (cur_C));
			pY <- (iv_yes[i_obs] == 1);
			logLik[2] <- logProb_path + bernoulli_log(response_yes[i_obs], pY);
			
			// path: start-EV-R_EV-IR_EV;
			logProb_path <- log((cur_R) * (1 - cur_C));
			pY <- (iv_yes[i_obs] == 0);
			logLik[3] <- logProb_path + bernoulli_log(response_yes[i_obs], pY);
			
			increment_log_prob(log_sum_exp(logLik));
		};
	};
	R_subj ~ normal(0, R_subj_sd);
	C_subj ~ normal(0, C_subj_sd);
	G_subj ~ normal(0, G_subj_sd);
}