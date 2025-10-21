# グレイ符号化前処理
#
# 	gray_preprocess.R


# 変調次数(COMPONENT_BITS)の読み込み
source("e:/LDPC_setting.R")

# 変調次数に応じたGRAY_ORDERテーブル選択
			if (COMPONENT_BITS == 3){					# MATLABに整合させる
	GRAY_ORDER <- order(c(4, 5, 7, 6, 2, 3, 1, 0, 12, 13, 15, 14, 10, 11, 9, 8, 28, 29, 31, 30, 26, 27, 25, 24, 20, 21, 23, 22, 18, 19, 17, 16, 52, 53, 55, 54, 50, 51, 49, 48, 60, 61, 63, 62, 58, 59, 57, 56, 44, 45, 47, 46, 42, 43, 41, 40, 36, 37, 39, 38, 34, 35, 33, 32))
}else if (COMPONENT_BITS == 2){
	GRAY_ORDER <- order(c(2, 3, 1, 0, 6, 7, 5, 4, 14, 15, 13, 12, 10, 11, 9, 8))
}else if (COMPONENT_BITS == 1){
	GRAY_ORDER <- order(c(1, 0, 3, 2))
}

# 変換の実行
MATLAB_symbols <- scan(file = "e:/datasymbols.R", sep = ",")													# MATLAB出力を読み込み
MATLAB_symbols <- MATLAB_symbols[!is.na(MATLAB_symbols)]															# NAを取り除く
converted_data <- GRAY_ORDER[MATLAB_symbols]
cat(converted_data-1, file=sprintf("e:/datasymbols_%d.py", COMPONENT_BITS), sep=",")	# pythonで入力するため0始まりで-1
