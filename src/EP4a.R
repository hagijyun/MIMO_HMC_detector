# EP法
# EP4a.R


# 重み付分散
weighted.variance <- function(x, w){
	pmean <- weighted.mean(x = x, w = w)

	return(
		weighted.mean(x = (x - pmean)^2, w = w)
	)
}


# EP法の処理
EP <- function(IT_MAX = 10, trial){
	## 初期値の設定
	gamma_pr <- rep(1/Es, 2*n)
	gamma_mu <- rep(   0, 2*n)

	## 繰り返し
	for (it in 1:IT_MAX){
		### 尤度の考慮 ###
		#### 十分統計量の計算
		SIGMA2_tmp <- solve(1/trial$sigma.w2 * t(trial$H) %*% trial$H + diag(gamma_pr))
		    mu_tmp <- as.vector(SIGMA2_tmp %*% (1/trial$sigma.w2 * t(trial$H) %*% trial$y + gamma_mu))
		#### 前回との差分が小さくなっていたら終了
		if (it != 1){
			SIGMA2_diff <- abs(SIGMA2_tmp - SIGMA2)
			    mu_diff <- abs(    mu_tmp -     mu)
#			cat("it = ", it, ": ", min(SIGMA2_diff), min(mu_diff), "\n")
#			if (any(1.5e-8 < SIGMA2_diff & SIGMA2_diff < 1.0e-4) ||
#					any(1.5e-8 <     mu_diff &     mu_diff < 1.0e-4)){ 
			if (any(SIGMA2_diff < .Machine$double.eps) ||	any(mu_diff < .Machine$double.eps)){
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
		sigma2_p <- sapply(1:(2*n), function(ct){ weighted.variance(x = QPSK_SYMB_real, w = weight[ct, ]) })
		    mu_p <- sapply(1:(2*n), function(ct){ weighted.mean(    x = QPSK_SYMB_real, w = weight[ct, ]) })
		#### 分散が小さすぎる場合は適度な値に置換する
		sigma2_p <- pmax(5.0e-7, sigma2_p)

		### 事前分布の更新 ###
		#### 十分統計量の計算
		gamma_pr_tmp <-    1/sigma2_p - 1/h2
		gamma_mu_tmp <- mu_p/sigma2_p - t/h2
		#### 分散が正の部分だけ更新
		ind <- which(0 < gamma_pr_tmp)
		gamma_pr[ind] <- gamma_pr_tmp[ind]
		gamma_mu[ind] <- gamma_mu_tmp[ind]
	}


	# 戻り値にまとめる
	return(list(u_hat_EP = mu, sigma2 = diag(SIGMA2)))
}
