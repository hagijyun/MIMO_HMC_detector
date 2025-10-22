/*
 * Stanを活用したMHGD法
 *
 *
 *	stan_MHGD1a.stan
 *
 *		m = nに限定
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

	// 複素ベクトルを展開した実ベクトルから、再び複素ベクトルの要素毎の絶対値を求める関数
  vector complex_abs(vector r, int n) {
    return sqrt(square(r[1:n]) + square(r[(n+1):(2*n)]));
  }
}

data {
  int<lower=1>   				 	n;            		// 送信アンテナ本数*2
  int<lower=1>   				 	m;           		 	// 受信アンテナ本数*2
  vector[m]		  				  y;            		// 観測値
  matrix[m, n]   				  H;  		      		// チャネル行列
	real<lower=0>						sigma_w2;					// 雑音電力

	int<lower=1>   			 		COMPONENT_NUM;		// I/Q軸上の変調多値数
	vector[COMPONENT_NUM]		symb_real;				// シンボル位置
	vector[COMPONENT_NUM-1]	BREAKS;						// シンボル硬判定のための中点

	real<lower=0>						d_qam;						// QAM格子間の最小距離
	vector[n]		  			 		z0;			 	  			// MMSE解
	matrix[n, m]						tau_Mp_H_dash;		// τMpH'
	matrix[n, n]						Mc;								// Mc
	int<lower=1>   			 		Ns;   						// 繰り返し回数の上限
}

parameters {
	real<lower=0>						sigma;						// ダミーのパラメータ
}

model {																			// 空行→尤度に何の制約もない
}

generated quantities {
	vector[n]	z;
	vector[n]	v;

	vector[n]	x;
	vector[n]	x_prev;
	real			l;
	real			l_prev;
	vector[m]	r;
	vector[m]	r_prev;

	real			log_Pacc;
	real			log_Puni;

	/* 対数MLコストとチャンピオン値 */
	real 			log_ML_cost;
	vector[n] x_champ;
	real 			log_ML_cost_champ;


	// Step 0
	z = z0;
	x_prev = Q(z,			n, COMPONENT_NUM, symb_real, BREAKS);
	r_prev = y - H * x_prev; l_prev = dot_self(r_prev);

	/* チャンピオン値の初期化 */
	x_champ = x_prev;
	log_ML_cost_champ = -dot_self(y - H * x_prev) / (2 * sigma_w2);		// ノルム（実数+虚数）/雑音電力（実数*2）

	for (i in 1:Ns){
		// Step 1
		for (j in 1:n){
			v[j] = normal_rng(0, sqrt(1.0/2.0));
		}

		// Step 2
		z = (x_prev + tau_Mp_H_dash * r_prev) + fmax(d_qam, max(complex_abs(r_prev, m/2))/sqrt(n)) * Mc * v;
		x = Q(z,			n, COMPONENT_NUM, symb_real, BREAKS);

		// Step 3
		r = y - H * x; l = dot_self(r);

		// Step 4
		log_Pacc = (l_prev - l) / 2;
		log_Puni = log(uniform_rng(0, 1));

		// Step 5
		if (log_Pacc >= log_Puni){
			; // 提案値は既にx, l, rに設定されているのでそれをそのまま使う
		}else{
			x = x_prev; l = l_prev; r = r_prev;
		}

		// Step 6
		/* 対数MLコストが上がっていたら、チャンピオン値を更新 */
		log_ML_cost = -dot_self(y - H * x) / (2 * sigma_w2);
		if (log_ML_cost_champ <= log_ML_cost){
			x_champ = x;
			log_ML_cost_champ = log_ML_cost;
		}

		/* 次のループに備える */
		x_prev = x; l_prev = l; r_prev = r;
	}
}
