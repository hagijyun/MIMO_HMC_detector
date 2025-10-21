# OAMP_Net2_LMMSE_QPSK_Ma5.ps1


$mod_types = @("QPSK", "16QAM", "64QAM")
$mod_types = @("64QAM")														# ★変調次数

# 変調方式ごとに異なる SNR の範囲を定義
$snrdb_values = @{
    "QPSK"  = (-1.. 4)		# (0..20)								# 0 から 20 まで 1 刻み
    "16QAM" = ( 6..11)		# (0, 5, 10, 15, 20)		# 任意の値
    "64QAM" = (12..17)		# 21 から 30 まで 1 刻み
}

# 変調方式ごとに異なるファイル名を定義
$inputFiles = @{
    "QPSK"  = "datasymbols_1.py"
    "16QAM" = "datasymbols_2.py"
    "64QAM" = "datasymbols_3.py"
}


foreach ($mod in $mod_types) {
	$inputFile = $inputFiles[$mod]
 	foreach ($snr in $snrdb_values[$mod]) {
		python OAMP_Net2_LMMSE_QPSK_Ma5.py --snrdb $snr --mod_type $mod --input_file $inputFile
	}
}
