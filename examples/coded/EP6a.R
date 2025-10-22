# EP法
# EP6a.R


# 重み付分散
weighted.variance <- function(x, w){
	pmean <- weighted.mean(x = x, w = w)

	return(
		weighted.mean(x = (x - pmean)^2, w = w)
	)
}


# EP法の処理：EP-SICから呼び出す形態に修正
EP_in_SIC <- function(ep_trial, IT_MAX = 10, SIC_init = NULL){
	## 初期化
	### アンテナ数（複素→実表現で2倍）
	n2 <- ncol(ep_trial$H) * 2

	### uの複素演算用の各種行列
	ep_HH <- real.mat(t(Conj(ep_trial$H)) %*% ep_trial$H)
	ep_Hy <- real.vec(t(Conj(ep_trial$H)) %*% ep_trial$y)

	### 停止条件の閾値
	EP.eps <- 1e-4									# .Machine$double.eps		# 1e-4

	### 忘却係数
	beta <- 0.2											# 1.0										# 0.2


	## 初期値の設定
	gamma_pr <- rep(1/Es, n2)
	gamma_mu <- rep(   0, n2); if (!is.null(SIC_init)){ gamma_mu <- SIC_init }

	## 繰り返し
	for (it in 1:IT_MAX){
		### 尤度の考慮 ###
		#### 十分統計量の計算
		SIGMA2_tmp <- solve(1/ep_trial$sigma.w2 * ep_HH + diag(gamma_pr))
		    mu_tmp <- as.vector(SIGMA2_tmp %*% (1/ep_trial$sigma.w2 * ep_Hy + gamma_mu))
		#### 前回との差分が小さくなっていたら終了
		if (it != 1){
			SIGMA2_diff <- abs(SIGMA2_tmp - SIGMA2)
			    mu_diff <- abs(    mu_tmp -     mu)
##	if (any(SIGMA2_diff < EP.eps) ||	any(mu_diff < EP.eps)){
		if (all(SIGMA2_diff < EP.eps) &&	all(mu_diff < EP.eps)){
				break
			}
		}
		#### 値を更新
		SIGMA2 <- SIGMA2_tmp
		    mu <-     mu_tmp

		### cavity ###
		#### 十分統計量の計算
		h2 <- diag(SIGMA2) / (1 - diag(SIGMA2) * gamma_pr)
		t  <- h2 * (mu /diag(SIGMA2) - gamma_mu)

		### 本物分布の考慮 ###
		#### 重みの設定
		weight <- sapply(QPSK_SYMB_real, function(point){ dnorm(x = point, mean = t, sd = sqrt(h2)) })
		#### 十分統計量の計算
		sigma2_p <- sapply(1:(n2), function(ct){ weighted.variance(x = QPSK_SYMB_real, w = weight[ct, ]) })
		    mu_p <- sapply(1:(n2), function(ct){ weighted.mean(    x = QPSK_SYMB_real, w = weight[ct, ]) })
		#### 分散が小さすぎる場合は適度な値に置換する
		sigma2_p <- pmax(5.0e-7, sigma2_p)

		### 事前分布の更新 ###
		#### 十分統計量の計算
		gamma_pr_tmp <- beta * (   1/sigma2_p - 1/h2) + (1 - beta) * gamma_pr
		gamma_mu_tmp <- beta * (mu_p/sigma2_p - t/h2) + (1 - beta) * gamma_mu
		#### 分散が正の部分だけ更新
		ind <- which(0 < gamma_pr_tmp)
		gamma_pr[ind] <- gamma_pr_tmp[ind]
		gamma_mu[ind] <- gamma_mu_tmp[ind]
	}


	# 戻り値にまとめる
	return(list(u_hat_EP = matrix(comp.vec(mu), ncol = 1), SIC_cov = comp.mat(SIGMA2), it = it))
}
