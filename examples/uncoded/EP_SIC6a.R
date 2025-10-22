# EP_SIC6a.R


# EP処理を行う関数の読み込み
source("EP6a.R")


# EP-SICの処理
EP_SIC <- function(trial, SIC_MAX = n-1, SIC_first_init = NULL){	# SIC_MAX：単発処理に加えてSICを何回施すか、0～n-1
	# 前処理
	## 初期化
	SIC_first <- TRUE																								# SICの一番最初か？
	SIC_y <- matrix(comp.vec(trial$y), ncol = 1)										# 受信信号
	SIC_H <- comp.mat(trial$H); colnames(SIC_H) <- 1:ncol(SIC_H)		# チャネル
	u_hat_SIC <- rep(NA_complex_, n)																# 戻り値

	# 本処理
	## 干渉（他アンテナからの影響）の除去を繰り返す
	while(TRUE){
		### SICの一番最初かどうかの判定
		if (SIC_first == TRUE){	SIC_init <- SIC_first_init; SIC_first <- FALSE	}
		else{										SIC_init <- NULL					 											}

		### 軟判定（EP）
#		EP_res <- EP_in_SIC(ep_trial = list(H = SIC_H, y = SIC_y, sigma.w2 = trial$sigma.w2), IT_MAX =  2, SIC_init = SIC_init)
		EP_res <- EP_in_SIC(ep_trial = list(H = SIC_H, y = SIC_y, sigma.w2 = trial$sigma.w2), IT_MAX = 10, SIC_init = SIC_init)
		SIC_u_hat <- EP_res$u_hat_EP

		### S/Nが最大となるアンテナ番号
		i_SN_max <- which.min(abs(diag(EP_res$SIC_cov)))

		### 軟判定結果を戻り値へ格納
		if (ncol(SIC_H) == n){ u_hat_SIC <- SIC_u_hat[, 1] }
		u_hat_SIC[as.numeric(names(i_SN_max))] <- SIC_u_hat[i_SN_max, 1]

		### 干渉除去処理
		if ((n-SIC_MAX) < ncol(SIC_H)){																### 干渉源がまだある
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
	return(list(u_hat_EP_SIC = matrix(real.vec(u_hat_SIC), ncol = 1), it = EP_res$it))
}
