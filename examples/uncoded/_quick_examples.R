# uncoded HMC MIMO detector
#
# sample based on MIMO5_para.Rmd
#




# 前処理


## ●デバグ用にStanの結果を保存するか？
DEBUG_OUT			<- FALSE			 # TRUE		# FALSE																						# ●デバグ用にStanの結果を保存する場合★

## ●Matlabで作成したデータシンボルや、更に試行結果を利用するか？
SAVED_SYMBOL	<- FALSE																																		# ●誤り訂正復号用に準備したビットデータを使う（誤り訂正なしの場合はFALSE）★
SAVED_TRIAL		<- FALSE																																		# ●Hやwを完全に再現できるため、提案法の繰り返し復号や方式別の個別復号が可能、読み込みメモリサイズを考慮してある程度小分、並列実行★

## ファイルから設定値を読み込む
  source("LDPC_setting.R")																														# ● LDPC符号の場合★

## バッチ実行時の引数を設定（args[1]：並列化ブロック数、args[2]：提案手法の繰り返し復号回数、args[3]：提案手法の繰り返し復号が初発 or 次発以降）
args <- commandArgs(trailingOnly = TRUE); if (length(args) == 0){ args[1] <- "1"; args[2] <- "0"; args[3] <- "FIRST" }
if (args[3] == "FIRST"){ EXTRINSIC_LLR <- FALSE }else{ EXTRINSIC_LLR <- TRUE }


## 関数定義

### 複素行列の実表現
real.mat <- function(mat){
	return(
		matrix(c(1, 0,  0, 1), 2, 2) %x% (Re(mat))
		+
		matrix(c(0, 1, -1, 0), 2, 2) %x% (Im(mat))
	)
}

### 実行列の複素表現
comp.mat <- function(mat){
	col_len <- ncol(mat)
	row_len <- nrow(mat)

	return(
		matrix(complex(re = mat[1:(row_len/2), 1:(col_len/2)], im = mat[(row_len/2 + 1):row_len, 1:(col_len/2)]), ncol = col_len/2)
	)
}

### 複素ベクトルの実表現
real.vec <- function(vec){
	return(
		c(Re(vec), Im(vec))
	)
}

### 実ベクトルの複素表現
comp.vec <- function(vec){
	len <- length(vec)
	return(
		complex(re = vec[1:(len/2)], im = vec[(len/2 + 1):len])
	)
}

### 行列の平方根
sqrt.mat <- function(mat){
	mat.svd <- svd(mat)
	return(
		mat.svd$u %*% diag(sqrt(mat.svd$d)) %*% t(Conj(mat.svd$v))
	)
}

### logsumexp：線形領域での規格化（入力：規格化前の対数値ベクトル、戻り値：規格化後の対数値ベクトル）
normalize <- function(l){
  # 入力の対数値ベクトルが最大値を取る番号
  max_ind <- which.max(l)

  # スケーリングを施すことでアンダーフローを極力抑止
  return(
    l - l[max_ind] - log1p(sum(exp(l[-max_ind] - l[max_ind])))
  )
}

### 総送信電力（送信アンテナ数を考慮した値）
sigma <- function (n, Es, SNR){
	return(
		n * Es * 10^(-SNR/10)
	)
}

### ●送信シンボルの元データ
if (SAVED_SYMBOL == TRUE){
	MATLAB_symbols <- scan(file = "e:/datasymbols.R", sep = ",")													# MATLAB出力を読み込み
	MATLAB_symbols <- MATLAB_symbols[!is.na(MATLAB_symbols)]															# NAを取り除く
	MATLAB_symbols <- t(matrix(MATLAB_symbols, nrow = 96, ncol = LONG_FRAME))							# アンテナ数 = 96前提、行列形式で1符号語分を格納
}

### 送信・伝送路・受信に関する処理
send.channel.receive <- function (SNR, k_tmp = k){
	## 送信シンボル
	if (SAVED_SYMBOL == TRUE){
		u_ID <- MATLAB_symbols[(k_tmp-1)%%LONG_FRAME + 1, ]																	# 1符号語の中からある送信タイミングの中身（行）を順繰りに（循環的に）抽出
	}else{
		u_ID <- sample(x = rep(x = 0:(COMPONENT_NUM-1), length.out = n), size = n, replace = F) * COMPONENT_NUM + sample(x = rep(x = 0:(COMPONENT_NUM-1), length.out = n), size = n, replace = F) + 1		# I/Q軸で独立にサイコロを振って組み合わせる、送信アンテナ数がQPSK/16QAM/64QAMの時に2/4/8の倍数となる必要あり
	}
	u <- QPSK_SYMB[u_ID]																																	# シンボルのID（整数値）をIQ平面上の実数値に変換（接頭語はQPSKだが他変調次数でもそのまま使える）
	u <- matrix(real.vec(u), ncol = 1)

	## 伝送路
	# 無相関なチャネル要素
	g <- complex(re = rnorm(n = n*m, mean = 0, sd = sqrt(CHANNEL.POW)), im = rnorm(n = n*m, mean = 0, sd = sqrt(CHANNEL.POW)))
	G <- matrix(g, ncol = n, nrow = m)

	# クロネッカモデル
	H <- real.mat(Pr %*% G %*% t(Conj(Pt)))

	## 雑音
	sigma.w2 <- sigma(n = n, Es = Es, SNR = SNR)
	w <- complex(re = rnorm(n = m, mean = 0, sd = sqrt(sigma.w2)), im = rnorm(n = m, mean = 0, sd = sqrt(sigma.w2)))
	w <- matrix(real.vec(w), ncol = 1)

	## 受信信号
	y <- H %*% u + w

	# 戻り値にまとめる
	return(
		list(y = y, H = H, sigma.w2 = sigma.w2, u_ID = u_ID, u = u, w = w)
	)
}

### シンボルの硬判定
hard.descision <- function(u_hat){
	sapply(u_hat, function(x){
			return(
				which.min(abs(x - QPSK_SYMB))
			)
	})
}


## 記号定数

### 各種の値
rho <- 0.5 * 0													# ●アンテナ間相関、0なら無相関★
  MCMC_it_MAX <- 1000										# ●HMCの繰り返し回数★
MR_MAX			<-  1												# 　HMCにおけるMultiple Restartの最大回数（初回含）
if (SAVED_SYMBOL){
	k_MAX <- LONG_FRAME * Turbo_trial			# ●試行回数の最大値（符号化あり向）
}else{
	k_MAX <- 10														# ●試行回数の最大値（符号化なし向）★
}
SNR_BY  <- +1														# ●SNRの刻み値

n <-   96																# ●送信アンテナ数、シンボルの混合比が等確率であればQPSK/16QAM/64QAMの時に2/4/8の倍数となる必要あり、SDを使う場合はn <= mとなる必要あり★
m <-   96																# ●受信アンテナ数★

MGS_iter <- 8 * n * 2^COMPONENT_BITS		# MCMC_it_MAX
Tikhonov_w_coded <- c(15, 62, 230)			# ●提案法の超パラメータ：伝送路符号化あり時の正則化加重
 HS_sigma <- c(3.5, 5.0, 3.0)						# ●提案法の超パラメータ馬蹄尤度を適用する際のσ_wに対する係数		 LDPC符号の場合★

EP_MAX <- 1															# EPの初期値変更試行回数

COMPONENT_NUM <- 2^COMPONENT_BITS				# I/Q軸上の状態数
COMPONENT_INT_SCALING_FACTOR <- c(1/sqrt(2), 1/sqrt(10), 1/sqrt(42))	# I/Q軸上の原数列に対するスケーリング係数（複素平均シンボル電力を1とした時）
SYMB_bits <- 2 * COMPONENT_BITS					# 1シンボルあたりのビット数（I/Qで2倍）
M <- 2^SYMB_bits												# シンボルの変調多値数
Es <- 1/2																# 実平均シンボル電力
CHANNEL.POW <- 1/2											# 実チャネルの平均電力

### QPSK/QAMシンボル（実数版）
QPSK_SYMB_real <- seq(from = -(COMPONENT_NUM - 1), to = COMPONENT_NUM - 1, by = 2) * COMPONENT_INT_SCALING_FACTOR[COMPONENT_BITS]
BREAKS <- (filter(x = QPSK_SYMB_real, filter = c(1, 1)) / 2)[-COMPONENT_NUM]			# 領域判定のための中点の算出
if (length(BREAKS) == 1){ dim(BREAKS) <- 1 }

### QPSK/QAMシンボル（Grayコード）
			if (COMPONENT_BITS == 3){					# MATLABに整合させる
	GRAY_ORDER <- order(c(4, 5, 7, 6, 2, 3, 1, 0, 12, 13, 15, 14, 10, 11, 9, 8, 28, 29, 31, 30, 26, 27, 25, 24, 20, 21, 23, 22, 18, 19, 17, 16, 52, 53, 55, 54, 50, 51, 49, 48, 60, 61, 63, 62, 58, 59, 57, 56, 44, 45, 47, 46, 42, 43, 41, 40, 36, 37, 39, 38, 34, 35, 33, 32))
}else if (COMPONENT_BITS == 2){
	GRAY_ORDER <- order(c(2, 3, 1, 0, 6, 7, 5, 4, 14, 15, 13, 12, 10, 11, 9, 8))
}else if (COMPONENT_BITS == 1){
	GRAY_ORDER <- order(c(1, 0, 3, 2))
}
SYMB_LOCATION <- apply(expand.grid(1:COMPONENT_NUM, 1:COMPONENT_NUM), 1, function(x){ complex(re = QPSK_SYMB_real[x[2]], im = QPSK_SYMB_real[x[1]]) })
QPSK_SYMB <- SYMB_LOCATION[GRAY_ORDER]

### 整数→ビット列の変換（ビット幅：SYMB_bits）
bit <- function(i){
	return(rev(as.integer(intToBits(i)[1:SYMB_bits])))
}

### ビットと整数の変換用行列
bits <- sapply(0:(M-1), function(x){ bit(x) })

### 相関行列
corr_mat_sqrt <- function(size, rho = rho){
	non_diag <- toeplitz(0:(size-1))
	for (i in 1:(size-1)){ non_diag[non_diag == i] <- rho^i }
	return(sqrt.mat(non_diag + diag(size)))
}
Pt <- corr_mat_sqrt(size = n, rho = rho)
Pr <- corr_mat_sqrt(size = m, rho = rho)


## 初期化

### SNRのリスト
SNR_LIST <- seq(from = SNR_MIN, to = SNR_MAX, by = SNR_BY)
names(SNR_LIST) <- sprintf("SNR=%.1f", SNR_LIST)

### 試行回数のリスト
k_LIST <- 1:k_MAX
names(k_LIST) <- sprintf("k=%d", k_LIST)

### 乱数種★
# set.seed(123)
  set.seed(Sys.time())

### ライブラリの読み込み、事前設定
library(ggplot2)			# AWS★
require(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
theme_set(theme_get() + theme(aspect.ratio = 3/4))

require(doParallel)
cl <- makeCluster(multi_processor, outfile = sprintf("./parallell_log.txt"))
registerDoParallel(cl)

require(purrr)

library(R.matlab)		# MATLABファイルの読み書き用

### 自作検出関数の読み込み

#### EP
source("EP_SIC6a.R")

#### MCMC
stan_mod_regular				<- stan_model(file = "MCMC4f.stan")					# 混合t分布+人為正則化
stan_MGS_mod						<- stan_model(file = "stan_MGS4e.stan")			# MGS、GS+高速化、64QAMのリフレッシュ方法を書籍とあわせる
stan_MHGD_mod						<- stan_model(file = "stan_MHGD1a.stan")		# MHGD、m = nの場合に限定
stan_Lang_mod						<- stan_model(file = "stan_Lang1a.stan")		# 熱欲ランジュバン

### 混合t分布における最適なパラメータ（σは送信アンテナ数が多いほど少し小さくなる）
NOT_OVERLOADING <- !(n > m)
MIXED_NORMAL_SIGMA_ADJ <- function(fine_tuning){
	if (fine_tuning == TRUE){	return(0.02 * log2(n/16))	}
	else{											return(0)									}
}
MCMC_FIXED_SIGMA <- function(fine_tuning = NOT_OVERLOADING){
			 if (COMPONENT_BITS == 1){	return((0.3 - MIXED_NORMAL_SIGMA_ADJ(fine_tuning)		 ) / (COMPONENT_NUM/2))	}					#  QPSK
	else if (COMPONENT_BITS == 2){	return((0.3 - MIXED_NORMAL_SIGMA_ADJ(fine_tuning)		 ) / (COMPONENT_NUM/2))	}					# 16QAM
	else if (COMPONENT_BITS == 3){	return((0.3 - MIXED_NORMAL_SIGMA_ADJ(fine_tuning)/1.5) / (COMPONENT_NUM/2))	}					# 64QAM
	else{														return(-1)		}
}

OPTIMAL_t_SIGMA_WEIGHT <- function(){
			 if (COMPONENT_BITS == 1){	return(0.5)	}						#  QPSK
	else if (COMPONENT_BITS == 2){	return(0.5)	}						# 16QAM
	else if (COMPONENT_BITS == 3){	return(0.8)	}						# 64QAM
	else{														return(-1)	}
}

OPTIMAL_nu <- function(){
			 if (COMPONENT_BITS == 1){	return(1.8)		}					#  QPSK
	else if (COMPONENT_BITS == 2){	return(1.8)		}					# 16QAM
	else if (COMPONENT_BITS == 3){	return(2.5)		}					# 64QAM
	else{														return(-1)		}
}


### ●保存済みtrialを読み込む場合
if (SAVED_TRIAL == TRUE){
			 if (COMPONENT_BITS == 1){ trial_file <- sprintf("e:/trail_%s%s.RData", "QPSK" , args[1]) }
	else if (COMPONENT_BITS == 2){ trial_file <- sprintf("e:/trail_%s%s.RData", "16QAM", args[1]) }
	else if (COMPONENT_BITS == 3){ trial_file <- sprintf("e:/trail_%s%s.RData", "64QAM", args[1]) }
	load(file = trial_file)
}

### ●シンボルLLRの設定
if (EXTRINSIC_LLR == TRUE){
	LLR <- readMat(sprintf("e:/LLR%s.mat", args[1]))[[1]]
	source("LLR.R")
}else{
	LLR <-  array(log(1/COMPONENT_NUM), c(COMPONENT_NUM, 2*n, k_MAX, length(SNR_LIST)))																												# 等確率の場合
}




# 本処理

## 送信・伝送路・受信・判定の実行
	res <- foreach(ct = 1:length(SNR_LIST), .packages = c("rstan", 											"foreach"			 )) %dopar% {														# lapply(SNR_LIST, function(snr){	# SNR、並列時は%dopar%、直列時は%do%		SDを使う場合は追加
	snr <- SNR_LIST[ct]
	lapply(k_LIST, function(k){																																																								# 試行回数
	  ### 進捗バーの表示
		if (floor(k %% (k_MAX/100)) == 0){ cat(sprintf("%s: SNR = %.1f, k = %06d / %06d, (blk, iter) = (%s, %s)\n", Sys.time(), snr, k, k_MAX, args[1], args[2])) }		# 試行回数に応じて変わる、デフォルトでは100記録分出力される

		### 送信・伝送路・受信処理
		if (SAVED_TRIAL == TRUE){
			trial <- trial_stored[[ct]][[k]]$trial
		}else{
			trial <- send.channel.receive(SNR = snr, k_tmp = k)																																										# 試行回数を入力に設定#
		}

		### リッジ正則化（精度行列と平均ベクトル）
		svd_H <- svd(trial$H)
		opt <- sapply(simplify = FALSE, USE.NAMES = TRUE, c("HMC"), function(opt_method){
#			regular_lambda2 <- (max(svd_H$d))^2 / Tikhonov_w_coded[COMPONENT_BITS]																																# ●正則化：あり（　　HMC）、コメントアウトでどちらか一つ★
			regular_lambda2 <- 0																																																									# ●正則化：なし（素のHMC）、コメントアウトでどちらか一つ★
			regular_OMEGA <- t(trial$H) %*% trial$H + regular_lambda2 * diag(2*n)
			regular_mu <- as.vector(svd_H$v %*% diag(svd_H$d / ((svd_H$d)^2 + regular_lambda2)) %*% t(svd_H$u) %*% trial$y)

			return(list(regular_OMEGA = regular_OMEGA, regular_mu = regular_mu))
		})

		### 軟判定（MMSE）
		u_hat_MMSE <- solve(t(trial$H) %*% trial$H + trial$sigma.w2/Es * diag(n*2)) %*% t(trial$H) %*% trial$y
		### 判定（MMSE）
		u_hat_ID_MMSE <- 								comp.vec(t(u_hat_MMSE))																																									# 軟値

		ML_cost_champ <- -1e+300
		foreach (EP_ct = 1:EP_MAX) %do% {
			### 軟判定（EP）
			EP_res <- EP_SIC(trial = trial, SIC_MAX = 0, SIC_first_init = NULL)																																		# SICの一番最初の初期値（NULLでMMSE相当）
			u_hat_EP <- EP_res$u_hat_EP_SIC
			### 硬判定（EP）
			u_hat_ID_EP_tmp <- hard.descision(comp.vec(t(u_hat_EP)))

			# MLコストを抽出
			ML_cost <- -crossprod(trial$y - trial$H %*% real.vec(QPSK_SYMB[u_hat_ID_EP_tmp]))

			# MLコストの最大値とその時の結果を更新
			if (ML_cost_champ <= ML_cost){
				ML_cost_champ <- ML_cost
				u_hat_ID_EP <- comp.vec(t(u_hat_EP))																																																# 軟値
			}
		}


if (TRUE){																																																																	# ●素のHMC、どちらか一つ★
		### 軟判定（HMC、正則化+混合t分布）：HSなし, warm=12, 統計量はML値
		fit_stan_q <- purrr::quietly(sampling)(object = stan_mod_regular,
												 data = list(n = n*2, m = m*2, y = as.vector(trial$y), H = trial$H, sigma_w2 = trial$sigma.w2, COMPONENT_NUM = COMPONENT_NUM, symb_real = QPSK_SYMB_real, BREAKS = BREAKS, sigma = MCMC_FIXED_SIGMA()*OPTIMAL_t_SIGMA_WEIGHT(), nu = OPTIMAL_nu(), regular_mu = opt$HMC$regular_mu, regular_OMEGA = opt$HMC$regular_OMEGA, log_q = LLR[ , , k, ct]),
												 chains = MR_MAX * round(MCMC_it_MAX/(n*2)), iter = 				n*2, warmup = 12,
										 		 seed = 123)
		fit_stan_l <- fit_stan_q$result

		# 軟判定シンボルを出力
		u_hat_MCMC_l_chains <- extract(fit_stan_l, permuted = FALSE, inc_warmup = TRUE, par = c("u", "ML_cost"))																# 硬判定シンボルに基づくMLコスト
		old_dim_l <- dim(u_hat_MCMC_l_chains); dim(u_hat_MCMC_l_chains) <- c(old_dim_l[1] * old_dim_l[2], old_dim_l[3])
		champ_chain_l <- which.max(u_hat_MCMC_l_chains[, dim(u_hat_MCMC_l_chains)[2]])
		u_hat_MCMC_l <- u_hat_MCMC_l_chains[champ_chain_l, -dim(u_hat_MCMC_l_chains)[2]]
		### 判定（MCMC）
		u_hat_ID_MCMC_fixed_sigma  <- 							 comp.vec(t(u_hat_MCMC_l))																																	# 軟値

		# デバグ用にStanの結果を保存する場合
		if (DEBUG_OUT){
			saveRDS(list(COMPONENT_BITS = COMPONENT_BITS, SNR = SNR_LIST[ct], k = k, u_true = trial$u, fit_stan_l = fit_stan_l), file = sprintf("e:/DEBUG_fit_stan%s%s_snr%02d_trial%04d_%s.RData", args[1], args[2], SNR_LIST[ct], k, format(Sys.time(), "%Y%m%d%H%M%S")))
		}
}else{
	u_hat_ID_MCMC_fixed_sigma <- NULL
}


if (FALSE){
		### 軟判定（MGS）
		fit_stan_MGS_q <- purrr::quietly(sampling)(object = stan_MGS_mod,
												 data = list(n = n*2, m = m*2, y = as.vector(trial$y), H = trial$H, sigma_w2 = trial$sigma.w2, COMPONENT_NUM = COMPONENT_NUM, symb_real = QPSK_SYMB_real, MMSE = as.vector(u_hat_MMSE), MGS_it_MAX = MGS_iter),
												 chains = round(10*MCMC_it_MAX / MGS_iter), iter = 1, warmup =  0,
										 		 seed = 123)
		fit_stan_MGS   <- fit_stan_MGS_q$result

		u_hat_MGS <- 				rstan::extract(fit_stan_MGS, pars = "MGS_u_hat_champ", permuted = FALSE, inc_warmup = FALSE)
		MGS_gamma_champ <-  rstan::extract(fit_stan_MGS, pars = "MGS_gamma_champ", permuted = FALSE, inc_warmup = FALSE)
		MGS_gamma_champ_max_ind <- which.max(MGS_gamma_champ)
		u_hat_MGS <- matrix(u_hat_MGS[1, MGS_gamma_champ_max_ind, ], ncol = 1, nrow = 2*n)
		### 判定（MGS）
		u_hat_ID_MGS <- 							 comp.vec(t(u_hat_MGS))																																										# 軟値
}else{
	u_hat_ID_MGS <- NULL
}

if (FALSE){
		H_inv <- solve(trial$H); d_qam = 2/2 * COMPONENT_INT_SCALING_FACTOR[COMPONENT_BITS]
		### 軟判定（MHGD）：m = nの場合に限定
		fit_stan_MHGD_q <- purrr::quietly(sampling)(object = stan_MHGD_mod,
												 data = list(n = n*2, m = m*2, y = as.vector(trial$y), H = trial$H, sigma_w2 = trial$sigma.w2, COMPONENT_NUM = COMPONENT_NUM, symb_real = QPSK_SYMB_real, BREAKS = BREAKS, d_qam = d_qam, z0 = as.vector(u_hat_MMSE), tau_Mp_H_dash = 1 * solve(t(trial$H) %*% trial$H + 1*trial$sigma.w2/(d_qam^2) * diag(2*n)) %*% t(trial$H), Mc = H_inv %*% diag(1/sqrt(colSums(H_inv^2))), Ns = n*2 * 8),
												 chains = round(10*MCMC_it_MAX / (n*2 * 8)), iter = 1, warmup =  0,
										 		 seed = 123)
		fit_stan_MHGD   <- fit_stan_MHGD_q$result

		u_hat_MHGD <- 			 rstan::extract(fit_stan_MHGD, pars = "x_champ"					 , permuted = FALSE, inc_warmup = FALSE)
		MHGD_gamma_champ <-  rstan::extract(fit_stan_MHGD, pars = "log_ML_cost_champ", permuted = FALSE, inc_warmup = FALSE)
		MHGD_gamma_champ_max_ind <- which.max(MHGD_gamma_champ)
		u_hat_MHGD <- matrix(u_hat_MHGD[1, MHGD_gamma_champ_max_ind, ], ncol = 1, nrow = 2*n)
		### 判定（MHGD）
		u_hat_ID_MHGD <- 								comp.vec(t(u_hat_MHGD))																																									# 軟値
}else{
	u_hat_ID_MHGD <- NULL
}

if (FALSE){
		### 軟判定（Lang）
		fit_stan_Lang_q <- purrr::quietly(sampling)(object = stan_Lang_mod,
												 data = list(n = n*2, m = m*2, min_n_m = min(n, m) * 2, y = as.vector(trial$y), H = trial$H, s = svd_H$d, U = svd_H$u, V_t = t(svd_H$v), sigma_w2 = trial$sigma.w2, sigma_0 = sqrt(trial$sigma.w2), sigma = exp(seq(log(1), log(0.01), length.out = 20)), L = 20, T = 70, eps = 3.0e-5, tau = 1/2, COMPONENT_NUM = COMPONENT_NUM, symb_real = QPSK_SYMB_real, BREAKS = BREAKS),
												 chains = round(10*MCMC_it_MAX / (20*70)), iter = 1, warmup =  0,
										 		 seed = 123)
		fit_stan_Lang   <- fit_stan_Lang_q$result

		u_hat_Lang <- 			 rstan::extract(fit_stan_Lang, pars = "x_champ"					 , permuted = FALSE, inc_warmup = FALSE)
		Lang_gamma_champ <-  rstan::extract(fit_stan_Lang, pars = "log_ML_cost_champ", permuted = FALSE, inc_warmup = FALSE)
		Lang_gamma_champ_max_ind <- which.max(Lang_gamma_champ)
		u_hat_Lang <- matrix(u_hat_Lang[1, Lang_gamma_champ_max_ind, ], ncol = 1, nrow = 2*n)
		### 判定（Lang）
		u_hat_ID_Lang <- 								comp.vec(t(u_hat_Lang))																																									# 軟値
}else{
	u_hat_ID_Lang <- NULL
}


		### dummyの設定
		dummy <- NULL

		### 戻り値にまとめる
		return(
			list(
				u_ID							= trial$u_ID,																																																			# どの検出方法でも必ず含める「正解」のシンボル

				u_hat_ID_MMSE			= u_hat_ID_MMSE,
				u_hat_ID_EP				= u_hat_ID_EP,

				u_hat_ID_MCMC_fixed_sigma = u_hat_ID_MCMC_fixed_sigma,

				u_hat_ID_MGS			= u_hat_ID_MGS,
				u_hat_ID_MHGD			= u_hat_ID_MHGD,
				u_hat_ID_Lang			= u_hat_ID_Lang,

				dummy							= dummy
			)
		)
	})
}

stopCluster(cl)





# ●硬判定

res_hard <- lapply(res, function(by_snr){																																# SNR
	lapply(by_snr, function(by_trial){																																		# 試行回数
		by_trial$u_hat_ID_MMSE							<- hard.descision(by_trial$u_hat_ID_MMSE)
		by_trial$u_hat_ID_EP								<- hard.descision(by_trial$u_hat_ID_EP)
		by_trial$u_hat_ID_MCMC_fixed_sigma	<- hard.descision(by_trial$u_hat_ID_MCMC_fixed_sigma)
		by_trial$u_hat_ID_turbo							<- hard.descision(by_trial$u_hat_ID_turbo)
		by_trial$u_hat_ID_MGS								<- hard.descision(by_trial$u_hat_ID_MGS)
		by_trial$u_hat_ID_MHGD							<- hard.descision(by_trial$u_hat_ID_MHGD)
		by_trial$u_hat_ID_Lang							<- hard.descision(by_trial$u_hat_ID_Lang)
		by_trial$u_hat_SD										<- hard.descision(by_trial$u_hat_ID_SD	)
#
		return(by_trial)
	})
})


## 誤り率の計算

### 補誤差関数、SISO AWGNでのBER（理論値）を求める関数
erfc <- function(x){ return(2 * pnorm(x*sqrt(2), lower = FALSE)) }
qpsk.awgn.ber <- function(SNR.dB){												# SISO AWGNの理論値
	SNR.True <- 10^(SNR.dB/10) * (m/n)											# 真値に変換、送受信アンテナ数の比率に基づき補正
				if (COMPONENT_BITS == 1){													# QPSK
		return(1/ 2 * erfc(sqrt(SNR.True/ 2)) )
	}else if (COMPONENT_BITS == 2){													# 16QAM
		return(3/ 8 * erfc(sqrt(SNR.True/10)) )
	}else if (COMPONENT_BITS == 3){													# 64QAM
		return(7/24 * erfc(sqrt(SNR.True/42)) )
	}else{
		return(-1)
	}
}

### SISO AWGNでのBER（理論値）
AWGN_BER <- qpsk.awgn.ber(SNR_LIST)

### 平均BERを計算する関数
average_BER <- function(res, method){
	res_err <- sapply(res, function(by_snr){								# SNR
		sapply(by_snr, function(by_trial){										# 試行回数
					 if (method == "MMSE"							){ u_hat_bits <- as.vector(bits[, by_trial$u_hat_ID_MMSE							])	}
			else if (method == "MMSE-SIC"					){ u_hat_bits <- as.vector(bits[, by_trial$u_hat_ID_MMSE_SIC					])	}
			else if (method == "EP"								){ u_hat_bits <- as.vector(bits[, by_trial$u_hat_ID_EP								])	}
			else if (method == "HMC-fixed-sigma"	){ u_hat_bits <- as.vector(bits[, by_trial$u_hat_ID_MCMC_fixed_sigma	])	}
			else if (method == "TURBO"						){ u_hat_bits <- as.vector(bits[, by_trial$u_hat_ID_turbo							])	}
			else if (method == "MGS"							){ u_hat_bits <- as.vector(bits[, by_trial$u_hat_ID_MGS								])	}
			else if (method == "MHGD"							){ u_hat_bits <- as.vector(bits[, by_trial$u_hat_ID_MHGD							])	}
			else if (method == "Lang"							){ u_hat_bits <- as.vector(bits[, by_trial$u_hat_ID_Lang							])	}
#
			else if (method == "SD"								){ u_hat_bits <- as.vector(bits[, by_trial$u_hat_ID_SD								])	}

			u_bits <- as.vector(bits[, by_trial$u_ID])
			return(
				sum(xor(u_hat_bits, u_bits))
			)
		})
	})
	
	return(colMeans(res_err) / (n * SYMB_bits))
}

### 各方式のBER
MMSE_BER								<- average_BER(res = res_hard, method = "MMSE"						)
EP_BER									<- average_BER(res = res_hard, method = "EP"							)
MCMC_FIXED_SIGMA_BER		<- average_BER(res = res_hard, method = "HMC-fixed-sigma"	)
#MGS_BER								<- average_BER(res = res_hard, method = "MGS"							)
#MHGD_BER								<- average_BER(res = res_hard, method = "MHGD"						)
#Lang_BER								<- average_BER(res = res_hard, method = "Lang"						)



## 結果の確認

### 値の一覧
cat("SNR_LIST\n"); print(SNR_LIST)
cat("rho\n"); print(rho)
cat("n\n"); print(n)
cat("m\n"); print(m)
cat("k_MAX\n"); print(k_MAX)
cat("AWGN_BER\n"); print(AWGN_BER)
cat("MMSE_BER\n"); print(MMSE_BER)
cat("EP_BER\n"); print(EP_BER)
cat("MCMC_FIXED_SIGMA_BER\n"); print(MCMC_FIXED_SIGMA_BER)
