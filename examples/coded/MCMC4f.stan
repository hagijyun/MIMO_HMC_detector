// MCMC4f.stan
// uncoded HMC用

data{
  int<lower=1>   			 	n;            		// 送信アンテナ本数*2
  int<lower=1>   			 	m;           		 	// 受信アンテナ本数*2
  vector[m]		  			  y;           			// 観測値
  matrix[m, n]    		  H;           		 	// チャネル行列
	real<lower=0>					sigma_w2;					// 雑音電力
//vector<lower=0>[m]		sigma_w;					// 雑音電力

	int<lower=1>   			 		COMPONENT_NUM;	// I/Q軸上の変調多値数
	vector[COMPONENT_NUM]		symb_real;			// シンボル位置
	vector[COMPONENT_NUM-1]	BREAKS;					// シンボル硬判定のための中点

	real<lower=0>						sigma;					// 標準偏差

  vector[n]		  			  regular_mu;  		 	// 人為正則化済の尤度：平均ベクトル
  matrix[n, n]    		  regular_OMEGA;		// 人為正則化済の尤度：精度行列

  matrix[COMPONENT_NUM, n]   log_q;				// 各アンテナにおけるシンボルの先見（対数）加重

	real<lower=0>						nu;							// t分布の自由度
}

parameters{
  vector[n]							u;								// シンボル
}

model{
  // 尤度の部分
	target += -0.5 * quad_form(regular_OMEGA/sigma_w2, u - regular_mu);

  // 事前分布の部分
	for (i in 1:n){
	  target += log_sum_exp(log_q[, i] - 0.5 * (nu + 1) * log1p(1 / nu * square((u[i] - symb_real) / sigma)));							// 混合t分布：nu = 大で正規分布と同等
	}
}

generated quantities{
	vector[n] u2;
	real ML_cost;

	// 先に硬判定
	for (j in 1:n){
				 if (u[j] < BREAKS[1]								){	u2[j] = symb_real[						1]; }
		else if (BREAKS[COMPONENT_NUM-1] <= u[j]){	u2[j] = symb_real[COMPONENT_NUM]; }
		else{
			for (ct in 2:(COMPONENT_NUM-1)){ if ((BREAKS[ct-1] <= u[j]) && (u[j] < BREAKS[ct])){ u2[j] = symb_real[ct]; break; } }
		}
	}

	// そのMLコスト
	ML_cost = - dot_self(y - H*u2) / (2 * sigma_w2);
}
