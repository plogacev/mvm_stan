
data
{
	int<lower=1> n_obs;
	int<lower=1> n_subj;
	int<lower=0,upper=1> iv_ev[n_obs];
	int<lower=0,upper=1> iv_yes[n_obs];
	int<lower=0,upper=1> iv_high_interference[n_obs];
	int<lower=0,upper=1> response_yes[n_obs];
	int<lower=1,upper=n_subj> subj[n_obs];
}

parameters
{
	real R_low;
	real C_low;
	real G_low;
	real R_diff;
	real C_diff;
	real G_diff;
	vector[n_subj] R_low_subj;
	real<lower=0> R_low_subj_sd;
	vector[n_subj] C_low_subj;
	real<lower=0> C_low_subj_sd;
	vector[n_subj] G_low_subj;
	real<lower=0> G_low_subj_sd;
	vector[n_subj] R_diff_subj;
	real<lower=0> R_diff_subj_sd;
	vector[n_subj] C_diff_subj;
	real<lower=0> C_diff_subj_sd;
	vector[n_subj] G_diff_subj;
	real<lower=0> G_diff_subj_sd;
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
	real cur_R_low;
	real cur_C_low;
	real cur_G_low;
	real cur_R_diff;
	real cur_C_diff;
	real cur_G_diff;
	int cur_subj;
	
	for(i_obs in 1:n_obs)
	{
		cur_subj <- subj[i_obs];
		cur_R_low <- inv_logit(R_low + R_low_subj[cur_subj]);
		cur_C_low <- inv_logit(C_low + C_low_subj[cur_subj]);
		cur_G_low <- inv_logit(G_low + G_low_subj[cur_subj]);
		cur_R_diff <- inv_logit(R_diff + R_diff_subj[cur_subj]);
		cur_C_diff <- inv_logit(C_diff + C_diff_subj[cur_subj]);
		cur_G_diff <- inv_logit(G_diff + G_diff_subj[cur_subj]);
		
		if(iv_ev[i_obs] == 0)
		{
			real logLik[2];
			
			// path: start0-start-MV-NR;
			R <- cur_R_low + cur_R_diff * iv_high_interference[i_obs];
			C <- cur_C_low + cur_C_diff * iv_high_interference[i_obs];
			G <- cur_G_low + cur_G_diff * iv_high_interference[i_obs];
			pY <- G;
			logProb_path <- log((1 - R));
			logLik[1] <- logProb_path + bernoulli_log(response_yes[i_obs], pY);
			
			// path: start0-start-MV-R_MV;
			R <- cur_R_low + cur_R_diff * iv_high_interference[i_obs];
			C <- cur_C_low + cur_C_diff * iv_high_interference[i_obs];
			G <- cur_G_low + cur_G_diff * iv_high_interference[i_obs];
			pY <- (iv_yes[i_obs] == 1);
			logProb_path <- log((R));
			logLik[2] <- logProb_path + bernoulli_log(response_yes[i_obs], pY);
			
			increment_log_prob(log_sum_exp(logLik));
		};
		
		if(iv_ev[i_obs] == 1)
		{
			real logLik[3];
			
			// path: start0-start-EV-NR;
			R <- cur_R_low + cur_R_diff * iv_high_interference[i_obs];
			C <- cur_C_low + cur_C_diff * iv_high_interference[i_obs];
			G <- cur_G_low + cur_G_diff * iv_high_interference[i_obs];
			pY <- G;
			logProb_path <- log((1 - R));
			logLik[1] <- logProb_path + bernoulli_log(response_yes[i_obs], pY);
			
			// path: start0-start-EV-R_EV-CR_EV;
			R <- cur_R_low + cur_R_diff * iv_high_interference[i_obs];
			C <- cur_C_low + cur_C_diff * iv_high_interference[i_obs];
			G <- cur_G_low + cur_G_diff * iv_high_interference[i_obs];
			pY <- (iv_yes[i_obs] == 1);
			logProb_path <- log((R) * (C));
			logLik[2] <- logProb_path + bernoulli_log(response_yes[i_obs], pY);
			
			// path: start0-start-EV-R_EV-IR_EV;
			R <- cur_R_low + cur_R_diff * iv_high_interference[i_obs];
			C <- cur_C_low + cur_C_diff * iv_high_interference[i_obs];
			G <- cur_G_low + cur_G_diff * iv_high_interference[i_obs];
			pY <- (iv_yes[i_obs] == 0);
			logProb_path <- log((R) * (1 - C));
			logLik[3] <- logProb_path + bernoulli_log(response_yes[i_obs], pY);
			
			increment_log_prob(log_sum_exp(logLik));
		};
	};
	R_low_subj ~ normal(0, R_low_subj_sd);
	C_low_subj ~ normal(0, C_low_subj_sd);
	G_low_subj ~ normal(0, G_low_subj_sd);
	R_diff_subj ~ normal(0, R_diff_subj_sd);
	C_diff_subj ~ normal(0, C_diff_subj_sd);
	G_diff_subj ~ normal(0, G_diff_subj_sd);
}