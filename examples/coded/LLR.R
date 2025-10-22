# LLR.R
#
# ビットLLRをシンボルLLRに変換


# 整数とビット列の変換
bit_LLR <- function(i){
	return(rev(as.integer(intToBits(i)[1:COMPONENT_BITS])))
}

# ビット列の全パタン
bits_LLR <- sapply(0:(COMPONENT_NUM-1), function(x){
	bit_LLR(x)
})
if(is.vector(bits_LLR)){ bits_LLR <- matrix(bits_LLR, nrow = 1) }				# QPSK用

# 1が立つビット位置の組み合わせ（全0のビットパタンを除く）
combination_LLR <- lapply(2:COMPONENT_NUM, function(j){
	which(bits_LLR[, j] == 1)
})


if (length(dim(LLR)) == 3){ LLR <- array(LLR, dim = c(dim(LLR), 1)) } 	# SNRが1種類のみの場合を想定し、該当する4次元目を1で拡張
# 列（アンテナ）毎に組み合わせの和を求め、対数領域で規格化
tmp_LLR <- apply(X = LLR, MARGIN = c(2, 3, 4), FUN = function(GAMMA){
	normalize(c(0, sapply(combination_LLR, function(combi){ sum(GAMMA[combi]) })))
})


# Gray符号の順番に変更
GRAY_ORDER_I <- (c(0, 1, 3, 2, 6, 7, 5, 4) + 1)[1:COMPONENT_NUM]
GRAY_ORDER_Q <- rev(GRAY_ORDER_I)

tmp_LLR[,		 1 :	 n , ,] <- tmp_LLR[GRAY_ORDER_I,		1 : 	n , ,]
tmp_LLR[, (n+1):(2*n), ,] <- tmp_LLR[GRAY_ORDER_Q, (n+1):(2*n), ,]


# 最終出力
LLR <- tmp_LLR
