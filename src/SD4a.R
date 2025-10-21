# SD4a.R
#
#
# 超球半径が小さすぎる場合、解が見つからずに終わる可能性は一定存在する
# column pivoting非対応
#
# n <= 16であれば低SNRでも何とか回る、16 < nでは低中SNRは極端に時間がかかる
#		36*36& QPSK		14dB以上が現実的
#		32*32&64QAM	　28dB以上が現実的


# 前処理
require(SD)
#library(Rcpp); sourceCpp("C:/cygwin/home/hagijyun/C/cpp_lib/SD/src/sphere_decoding4a.cpp")


# SDの処理
SD <- function(trial){
	# 前処理
	## 初期化
	qr_H <- qr(trial$H * COMPONENT_INT_SCALING_FACTOR[COMPONENT_BITS]*2)						# 　補正：シンボル間隔が1
	SD_R <- qr.R(qr_H); sgn <- sign(diag(SD_R));	SD_R <- diag(sgn) %*% SD_R				# 　補正：Rの対角項が負にならないようにする
	SD_Q <- qr.Q(qr_H);														SD_Q <- SD_Q %*% diag(sgn)				# 逆補正：Rの対角項が負にならないようにする
#	pivot <- 1:(2*n)		# pivot <- order(diag(SD_R))	 															# 　補正：絶対値の大きいRの対角項から計算するように行を入れ替え（数値演算精度向上）
#	SD_R <- SD_R[pivot, pivot]
#	SD_Q <- t(t(SD_Q)[pivot, ])
#	SD_y <- matrix(t(SD_Q[, 1:(2*n)]) %*% trial$y, ncol = 1)												# n <= mを仮定
	SD_y <- matrix(t(SD_Q)						%*% trial$y, ncol = 1)												# n <= mを仮定
	SD_d2 <- 2 * qchisq(p = 0.99, df = n*2) * trial$sigma.w2												# カイ二乗分布の99%値
#	if (n < m){
#		SD_dash2 <- 0		# crossprod(t(SD_Q[, (2*n+1):(2*m)]) %*% trial$y)
#	}else{
		SD_dash2 <- 0
#	}
	SD_m <- 2 * n

	# 本処理
	ret_tmp <- sphere_decoding(SD_R, SD_y, SD_d2, SD_dash2, SD_m, COMPONENT_NUM)
	u_hat_SD <- ret_tmp$s_champ * COMPONENT_INT_SCALING_FACTOR[COMPONENT_BITS]*2		# 逆補正：シンボル間隔が1
#	u_hat_SD <- u_hat_SD[pivot]																											# 逆補正：絶対値の大きいRの対角項から計算するように行を入れ替え（数値演算精度向上）

	# 後処理
	## 戻り値にまとめる
	return(list(u_hat_SD = matrix(u_hat_SD, ncol = 1)))
}
