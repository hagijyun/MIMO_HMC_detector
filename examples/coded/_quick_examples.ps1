# ☆MIMOの信号検出プログラムを並列実行させるPowerShellスクリプト
#
#								誤り訂正：あり（初発 + 繰返）
#
# based on para5.ps1


# 前処理
PowerShell Set-ExecutionPolicy RemoteSigned
chcp 932

# PATHを追加
#	$env:PATH += ";C:\Program Files\R\R-4.1.2\bin"					# Rscript
#	$env:PATH += ";C:\rtools40\usr\bin"     								# Rtools(make など)
#	$env:PATH += ";C:\rtools40\mingw64\bin" 								# Rtools(gcc, g++ など)

& 'Rscript' --encoding=UTF-8 'stan_pre_compilation.R'

# 本処理（並列化ブロック）
Start-Job -ScriptBlock {
	param($j, $i_MAX = 5)

  # PATHを追加
#	$env:PATH += ";C:\Program Files\R\R-4.1.2\bin"					# Rscript

  # 作業ディレクトリを設定
  Set-Location 'e:/'

	# 初回
	& 'Rscript' --encoding=UTF-8 'MIMO5_para.Rmd' "$j" "0" "FIRST"
	Start-Sleep -s 1
	# ●コメントアウトでどれか一つ
#	& 'matlab' -batch "send_receive2_para($j, 0, {'MMSE', 'EP', 'HMC', 'MGS', 'MHGD', 'Lang'})"		# 全てのパタン
	& 'matlab' -batch "send_receive2_para($j, 0, {'MMSE', 'EP', 'HMC'})"													# MMSE+EP+HMCのパタン
	Start-Sleep -s 1

	# 繰返
  for ($i=1; $i -le $i_MAX; $i++){
		& 'Rscript' --encoding=UTF-8 'MIMO5_para.Rmd' "$j" "$i" "not_FIRST"
		Start-Sleep -s 1
		& 'matlab' -batch "send_receive2_para($j, $i, {'HMC'})"
		Start-Sleep -s 1
  }
} -ArgumentList $args
