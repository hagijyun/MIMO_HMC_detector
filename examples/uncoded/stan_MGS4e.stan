/*
 * Stanを活用したMixed Gibbs法
 *
 *
 *	stan_MGS4e.stan
 */


functions {
	// 対数スケーリング関数
	vector log_scaling(vector log_vector){
		return(log_vector - max(log_vector));
	}
}

data {
  int<lower=1>   			 	n;            // 送信アンテナ本数*2
  int<lower=1>   			 	m;            // 受信アンテナ本数*2
  vector[m]		  			  y;            // 観測値
  matrix[m, n]   			  H;  		      // チャネル行列
	real<lower=0>					sigma_w2;			// 雑音電力

	int<lower=1>   			 		COMPONENT_NUM;	// I/Q軸上の変調多値数
	vector[COMPONENT_NUM]		symb_real;			// シンボル位置

	vector[n]		  			  MMSE;			 	  // 未使用

	int<lower=1>   			 	MGS_it_MAX;   // 繰り返し回数の上限
}

transformed data {
	int MGS_K = COMPONENT_NUM;					// I/Q軸上の変調多値数
	real MGS_q = 1.0 / n;								// マルコフ連鎖の混合比
}

parameters {
	real<lower=0>					sigma;				// ダミーのパラメータ
}

model {																// 空行→尤度に何の制約もない
}

generated quantities {
	// 宣言と初期化
	/// κ
	real MGS_kappa;

	/// シンボルベクトルと対数MLコスト
	vector[n] MGS_u_hat;
	real MGS_gamma;

	/// シンボルベクトルと対数MLコストのチャンピオン
	vector[n] MGS_u_hat_champ;
	real MGS_gamma_champ;

	/// シンボル候補に対する重み（対数・線形）
	vector[MGS_K] MGS_log_w;
	vector[MGS_K] MGS_w;

	/// 高速化対応
	vector[m] Hu;
	vector[m] Hu_tmp;

	/// 64QAM対応
	int symb_num[n];
	int prev_symb_num;


	// 初期化
	for (MGS_n in 1:n){									// ランダムに初期化
		 symb_num[MGS_n] = categorical_rng(rep_vector(1.0/MGS_K, MGS_K));
		MGS_u_hat[MGS_n] = symb_real[symb_num[MGS_n]];
	}
	Hu = H * MGS_u_hat;
	MGS_gamma = -dot_self(y - Hu) / (2 * sigma_w2);		// ノルム（実数+虚数）/雑音電力（実数*2）
	MGS_u_hat_champ = MGS_u_hat;
	MGS_gamma_champ = MGS_gamma;

	// 繰り返し
	for (MGS_it in 2:MGS_it_MAX){				// iterを変える
		/// シンボルベクトルの更新
		for (MGS_n in 1:n){
			//// 0～1の一様乱数を取得し、マルコフ連鎖の混合のための指標に活用する
			MGS_kappa = uniform_rng(0, 1);														// MGSの場合
			Hu	-= col(H, MGS_n) * MGS_u_hat[MGS_n];
			//// 候補シンボルに対する重みの算出
			if (MGS_q < MGS_kappa){				// 逆温度 = 1の場合
				for (ct in 1:MGS_K){																		// あるアンテナの受信シンボルが送信シンボル（±）であった場合の、それぞれの対数コスト関数
					MGS_u_hat[MGS_n] = symb_real[ct];
					Hu_tmp = Hu + col(H, MGS_n)	 * MGS_u_hat[MGS_n];
					MGS_log_w[ct] = -dot_self(y - Hu_tmp) / (2 * sigma_w2);
				}
				MGS_w = exp(log_scaling(MGS_log_w));										// 比率は変えず、扱いやすくなるようにスケーリングを施す
				for (ct in 1:MGS_K){
					if (is_inf(MGS_w[ct]) || is_nan(MGS_w[ct])){ MGS_w[ct] = 1e+30; }		// 値が大きすぎた場合は、倍精度数値の最大値以下の大きめの値で置換する
				}
			}else{												// 逆温度 = ∞
				MGS_w = rep_vector(0, MGS_K);	// 一旦0クリア
				if (MGS_K == 8){	// 64QAMなら
					prev_symb_num = symb_num[MGS_n];
					MGS_w[prev_symb_num] = 1.0;
							 if (prev_symb_num == 1){															  MGS_w[prev_symb_num+1] = 1.0;	}
					else if (prev_symb_num == 8){	MGS_w[prev_symb_num-1] = 1.0;															  }
					else												{	MGS_w[prev_symb_num-1] = 1.0;	MGS_w[prev_symb_num+1] = 1.0;	}
				}else{
					for (ct in 1:MGS_K){																	// あるアンテナの受信シンボルが送信シンボル（±）であった場合の、それぞれの対数コスト関数
						MGS_w[ct] = uniform_rng(0, 1);
					}
				}
			}

			//// シンボル候補をサンプリング→規格化された対数コストに応じてカテゴリカル分布から整数を選び、その整数に応じてシンボルを選ぶ
			 symb_num[MGS_n] = categorical_rng(MGS_w/sum(MGS_w));
			MGS_u_hat[MGS_n] = symb_real[symb_num[MGS_n]];
			Hu += col(H, MGS_n)	* MGS_u_hat[MGS_n];
		}

		//// 対数MLコストが上がっていたら、チャンピオン値を更新
		MGS_gamma = -dot_self(y - Hu) / (2 * sigma_w2);
		if (MGS_gamma_champ <= MGS_gamma){
			MGS_u_hat_champ  = MGS_u_hat;
			MGS_gamma_champ  = MGS_gamma;
		}
	}
}
