# OAMP_Net2_LMMSE_QPSK_Ma5.ps1


$mod_types = @("QPSK", "16QAM", "64QAM")
$mod_types = @("64QAM")														# ���ϒ�����

# �ϒ��������ƂɈقȂ� SNR �͈̔͂��`
$snrdb_values = @{
    "QPSK"  = (-1.. 4)		# (0..20)								# 0 ���� 20 �܂� 1 ����
    "16QAM" = ( 6..11)		# (0, 5, 10, 15, 20)		# �C�ӂ̒l
    "64QAM" = (12..17)		# 21 ���� 30 �܂� 1 ����
}

# �ϒ��������ƂɈقȂ�t�@�C�������`
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
