# ☆MIMOの信号検出プログラムを並列実行させるPowerShellスクリプト
#
#								誤り訂正：あり（初発 + 繰返）
#
# para5.ps1


# 前処理
PowerShell Set-ExecutionPolicy RemoteSigned
chcp 932

# 本処理（並列化ブロック）
Start-Job -ScriptBlock {
	param($j, $i_MAX = 5)

	# 初回
	& 'C:\Program Files\R\R-3.6.1\bin\Rscript.exe' --encoding=UTF-8 'C:\cygwin\home\hagijyun\D\無ア開分（バックアップ+α）\MIMO\VB\conf\MIMO5_para.Rmd' "$j" "0" "FIRST"
#	& 'C:\Program Files\R\R-4.1.2\bin\Rscript.exe' --encoding=UTF-8 'C:\cygwin64\home\hokudai\conf\MIMO5_para.Rmd' "$j" "0" "FIRST"
	Start-Sleep -s 1
	# ●コメントアウトでどれか一つ
#	& 'C:\Program Files\MATLAB\R2024b\bin\matlab.exe' -batch "send_receive2_para($j, 0, {'HMC', 'EP', 'MMSE', 'MGS', 'MHGD', 'Lang'})"		# 素のHMCと既存法全てのパタン
	& 'C:\Program Files\MATLAB\R2024b\bin\matlab.exe' -batch "send_receive2_para($j, 0, {'HMC', 'EP', 'MMSE'})"														# HMC+EP+MMSEのパタン
	Start-Sleep -s 1
	#	HMC?.matとLLR?.matを別名で保存（後の再利用に備える）
	Get-ChildItem "E:\HMC$j.mat" | % { Copy-Item $_.FullName -Destination ("E:\$($_.BaseName)_first$($_.Extension)") }
	Get-ChildItem "E:\LLR$j.mat" | % { Copy-Item $_.FullName -Destination ("E:\$($_.BaseName)_first$($_.Extension)") }
	Start-Sleep -s 1

	# 繰返
  for ($i=1; $i -le $i_MAX; $i++){
		& 'C:\Program Files\R\R-3.6.1\bin\Rscript.exe' --encoding=UTF-8 'C:\cygwin\home\hagijyun\D\無ア開分（バックアップ+α）\MIMO\VB\conf\MIMO5_para.Rmd' "$j" "$i" "not_FIRST"
#		& 'C:\Program Files\R\R-4.1.2\bin\Rscript.exe' --encoding=UTF-8 'C:\cygwin64\home\hokudai\conf\MIMO5_para.Rmd' "$j" "$i" "not_FIRST"
		Start-Sleep -s 1
		& 'C:\Program Files\MATLAB\R2024b\bin\matlab.exe' -batch "send_receive2_para($j, $i, {'HMC'})"
		Start-Sleep -s 1
  }
} -ArgumentList $args




# ①設定変更、並列化ブロック数[$j_MAX]の決定（SNRの並列化数 * スクリプトによる並列化ブロック数 < マシンのコア数）
#		e:/LDPC_setting.R
#		e:/LDPC_setting.m
#
# ②送信データの作成（MATLAB出力のdatasymbols.Rを作業ディレクトリに保存）
# 	★send_receive2.m
#
#	③受信データの作成（提案法の繰り返し復号や方式別に復号処理をする場合でもHやwの再現性を確保、メモリに読み込むバイナリファイルサイズの都合で小分けに分割して並列処理を行う）			MIMO5_data.Rmdの中でTurbo/LDPC符号の切り替えを設定★
#	PowerShell Set-ExecutionPolicy RemoteSigned
#	chcp 932
#	for ($j=1; $j -le [$j_MAX]; $j++){
#		& 'C:\Program Files\R\R-3.6.1\bin\Rscript.exe' --encoding=UTF-8 'C:\cygwin\home\hagijyun\D\無ア開分（バックアップ+α）\MIMO\VB\conf\MIMO5_data.Rmd' "$j" "-1" "FIRST"
#	#	& 'C:\Program Files\R\R-4.1.2\bin\Rscript.exe' --encoding=UTF-8 'C:\cygwin64\home\hokudai\conf\MIMO5_data.Rmd' "$j" "-1" "FIRST"
#		Start-Sleep -s 1
#	}
#
# ④初回+繰り返し復号・誤り訂正
#		.\para5.ps1 1～[$j_MAX]
#
# ⑤結果の平均化
#		mean_numErrs.m
