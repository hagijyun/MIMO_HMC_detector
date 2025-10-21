/*
 * Stanを活用したAnnealed Lagevin法
 *
 *
 *	stan_Lang1a.stan
 */


functions {
	// 硬判定関数
	vector Q(vector z, 				int n, int COMPONENT_NUM, vector symb_real, vector BREAKS) {
		vector[n] x;
	
		for (j in 1:n){
					 if (z[j] < BREAKS[1]								){	x[j] = symb_real[						 1]; }
			else if (BREAKS[COMPONENT_NUM-1] <= z[j]){	x[j] = symb_real[COMPONENT_NUM]; }
			else{
				for (ct in 2:(COMPONENT_NUM-1)){ if ((BREAKS[ct-1] <= z[j]) && (z[j] < BREAKS[ct])){ x[j] = symb_real[ct]; break; } }
			}
		}
	
		return x;
	}

	// 対数スケーリング関数
	vector log_scaling(vector log_vector){
		return(log_vector - max(log_vector));
	}
}

data {
  int<lower=1>   				 	n;            		// 送信アンテナ本数*2
  int<lower=1>   				 	m;           		 	// 受信アンテナ本数*2
	int											min_n_m;					// n, mの最小値
  vector[m]		  				  y;            		// 観測値
	matrix[m, n]						H;								// チャネル行列
	vector[min_n_m]					s;								// チャネル行列の特異値分解（Σの対角要素）
	matrix[m, min_n_m]			U;								// チャネル行列の特異値分解（U）
	matrix[min_n_m, n]			V_t;							// チャネル行列の特異値分解（t(V)）
	real										sigma_w2;					// 熱雑音電力
	real										sigma_0;					// 熱雑音電圧
	int											L;								// 熱浴段階
	vector[L]								sigma;						// 混合正規分布の標準偏差
	int											T;								// 更新ステップ数
	real										eps;							// ε
	real										tau;							// τ

	int<lower=1>   			 		COMPONENT_NUM;		// I/Q軸上の変調多値数
	vector[COMPONENT_NUM]		symb_real;				// シンボル位置
	vector[COMPONENT_NUM-1]	BREAKS;						// シンボル硬判定のための中点
}

transformed data {
	vector[min_n_m]					eta;
	eta = U' * y;
}

parameters {
	real<lower=0>						dummy;						// ダミーのパラメータ
}

model {																			// 空行→尤度に何の制約もない
}

generated quantities {
	vector[min_n_m]					chi;
	vector[min_n_m]					chi_prev;
	vector[n]								x;

	vector[min_n_m]					Lambda;
	vector[min_n_m]					w;

	vector[min_n_m]					log_pos;
	vector[min_n_m]					log_lik;
	vector[n]								log_pri_x;
	vector[COMPONENT_NUM]		exp_part;

	/* チャンピオン値と対数MLコスト */
	vector[n] 							x_champ;
	real 										log_ML_cost_champ;


	for (j in 1:min_n_m){
		chi_prev[j] = uniform_rng(-1, 1);
	}
	for (l in 1:L){
		for (j in 1:min_n_m){
			if (sigma[l] * s[j] <= sigma_0){
				Lambda[j] = eps * sigma[l]^2 / sigma[L]^2 * (1					- (sigma[l]/sigma_0	)^2	* s[j]^2);
			}else{
				Lambda[j] = eps							 / sigma[L]^2 * (sigma[l]^2	- (sigma_0/s[j]			)^2					);
			}
		}

		for (t in 1:T){
			for (j in 1:min_n_m){
				w[j] = normal_rng(0, 1);
			}

			log_lik = s ./ fabs(rep_vector(sigma_0^2, min_n_m) - square(sigma[l] * s)) .* (eta - s .* chi_prev);

			x = V_t' * chi_prev;
			for (j in 1:n){
				exp_part = exp(log_scaling(-square(x[j] - symb_real) / (2 * sigma[l]^2)));
				log_pri_x[j] =  (1/sum(exp_part) * sum(symb_real .* exp_part) - x[j]) / sigma[l]^2;
			}

			for (j in 1:min_n_m){
				if (s[j] == 0){
					log_pos[j] = 					 (V_t * log_pri_x)[j];
				} else if (sigma_0 >= sigma[l] * s[j]){
					log_pos[j] = (log_lik + V_t * log_pri_x)[j];
				} else if (sigma_0 <  sigma[l] * s[j]){
					log_pos[j] = log_lik[j]										 ;
				}
			}

			chi = chi_prev + Lambda .* log_pos + sqrt(2*Lambda*tau) .* w;
			chi_prev = chi;				// 次のループに備える
		}
	}


	/* チャンピオン値 */
	x_champ = Q(V_t' * chi,			n, COMPONENT_NUM, symb_real, BREAKS);
	log_ML_cost_champ = -dot_self(y - H * x_champ) / (2 * sigma_w2);		// ノルム（実数+虚数）/雑音電力（実数*2）
}
