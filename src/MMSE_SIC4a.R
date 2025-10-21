# MMSE_SIC4a.R


# MMSE-SICの処理
MMSE_SIC <- function(trial){
	# 前処理
	## 初期化
	SIC_y <- matrix(comp.vec(trial$y), ncol = 1)										# 受信信号
	SIC_H <- comp.mat(trial$H); colnames(SIC_H) <- 1:ncol(SIC_H)		# チャネル
	u_hat_SIC <- rep(NA_complex_, n)																# 戻り値

	# 本処理
	## 干渉（他アンテナからの影響）の除去を繰り返す
	while(TRUE){
		### 共分散（MMSE）
		SIC_cov <- solve(t(Conj(SIC_H)) %*% SIC_H + trial$sigma.w2/Es * diag(ncol(SIC_H)))			# trial$sigma.w2/Esは比率なので複素表現でも不変

		### 軟判定（MMSE）
		SIC_u_hat <- SIC_cov %*% t(Conj(SIC_H)) %*% SIC_y

		### S/Nが最大となるアンテナ番号
		i_SN_max <- which.min(abs(diag(SIC_cov)))

		### 軟判定結果を戻り値へ格納
		u_hat_SIC[as.numeric(names(i_SN_max))] <- SIC_u_hat[i_SN_max, 1]

		### 干渉除去処理
		if (2 <= ncol(SIC_H)){																				### 干渉源がまだある
			#### 硬判定
			SIC_u_hat_ID <- hard.descision(t(SIC_u_hat))

			#### レプリカの作成
			SIC_u_replica <- QPSK_SYMB[SIC_u_hat_ID]

			#### 干渉除去
			SIC_y <- SIC_y - SIC_H[, i_SN_max, drop = FALSE] * SIC_u_replica[i_SN_max]
			SIC_H <- SIC_H[, -i_SN_max, drop = FALSE]
		}else{																												### 干渉源がもうない
			#### 繰り返し処理終了
			break
		}
	}

	# 後処理
	## 戻り値にまとめる
	return(list(u_hat_MMSE_SIC = matrix(real.vec(u_hat_SIC), ncol = 1)))
}
