# LDPC符号用に頻繁に変更する設定一覧
#
# LDPC_setting.R

COMPONENT_BITS		<-  3 #	 3																															# ●I/Q軸上のビット数：QPSK = 1, 16QAM = 2, 64QAM = 3★
SNR_MIN						<- 21	# 12																															# ●SNRの最小値★
SNR_MAX						<- 30	# 17																															# ●SNRの最大値★
Turbo_trial				<- 10																																		# ●LDPC符号として成立するフレームを何回試行？（おおよそ10～50）★


# 3GPP符号
LDPC_Z <- 96																																							# LDPC符号のリフティング数：符号長[ビット]は32 -> 2176, 64 -> 4352, 88 -> 5984, 96 -> 6528
LONG_FRAME_TABLE	<- ceiling(68 * LDPC_Z / (96 * c(2, 4, 6)));														# 1符号長がまたがる送信タイミング数（変調次数により異なる）
# Gallager符号
# LONG_FRAME_TABLE	<- ceiling(c(6144, 6144, 6336)[COMPONENT_BITS] / (96 * c(2, 4, 6)));	# 1符号長がまたがる送信タイミング数（変調次数により異なる）
# IEEE符号
# LONG_FRAME_TABLE	<- ceiling(1944 / (96 * c(2, 4, 6)))																	# 符号長：1944

LONG_FRAME				<- LONG_FRAME_TABLE[COMPONENT_BITS]																			# 変調次数毎に設定
multi_processor		<- min(SNR_MAX - SNR_MIN + 1, 16-2)																			# SNRの並列化数 * スクリプトによる並列化ブロック数 < マシンのコア数		# 16-2
