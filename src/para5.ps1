# ��MIMO�̐M�����o�v���O�����������s������PowerShell�X�N���v�g
#
#								�������F����i���� + �J�ԁj
#
# para5.ps1


# �O����
PowerShell Set-ExecutionPolicy RemoteSigned
chcp 932

# �{�����i���񉻃u���b�N�j
Start-Job -ScriptBlock {
	param($j, $i_MAX = 5)

	# ����
	& 'C:\Program Files\R\R-3.6.1\bin\Rscript.exe' --encoding=UTF-8 'C:\cygwin\home\hagijyun\D\���A�J���i�o�b�N�A�b�v+���j\MIMO\VB\conf\MIMO5_para.Rmd' "$j" "0" "FIRST"
#	& 'C:\Program Files\R\R-4.1.2\bin\Rscript.exe' --encoding=UTF-8 'C:\cygwin64\home\hokudai\conf\MIMO5_para.Rmd' "$j" "0" "FIRST"
	Start-Sleep -s 1
	# ���R�����g�A�E�g�łǂꂩ���
#	& 'C:\Program Files\MATLAB\R2024b\bin\matlab.exe' -batch "send_receive2_para($j, 0, {'HMC', 'EP', 'MMSE', 'MGS', 'MHGD', 'Lang'})"		# �f��HMC�Ɗ����@�S�Ẵp�^��
	& 'C:\Program Files\MATLAB\R2024b\bin\matlab.exe' -batch "send_receive2_para($j, 0, {'HMC', 'EP', 'MMSE'})"														# HMC+EP+MMSE�̃p�^��
	Start-Sleep -s 1
	#	HMC?.mat��LLR?.mat��ʖ��ŕۑ��i��̍ė��p�ɔ�����j
	Get-ChildItem "E:\HMC$j.mat" | % { Copy-Item $_.FullName -Destination ("E:\$($_.BaseName)_first$($_.Extension)") }
	Get-ChildItem "E:\LLR$j.mat" | % { Copy-Item $_.FullName -Destination ("E:\$($_.BaseName)_first$($_.Extension)") }
	Start-Sleep -s 1

	# �J��
  for ($i=1; $i -le $i_MAX; $i++){
		& 'C:\Program Files\R\R-3.6.1\bin\Rscript.exe' --encoding=UTF-8 'C:\cygwin\home\hagijyun\D\���A�J���i�o�b�N�A�b�v+���j\MIMO\VB\conf\MIMO5_para.Rmd' "$j" "$i" "not_FIRST"
#		& 'C:\Program Files\R\R-4.1.2\bin\Rscript.exe' --encoding=UTF-8 'C:\cygwin64\home\hokudai\conf\MIMO5_para.Rmd' "$j" "$i" "not_FIRST"
		Start-Sleep -s 1
		& 'C:\Program Files\MATLAB\R2024b\bin\matlab.exe' -batch "send_receive2_para($j, $i, {'HMC'})"
		Start-Sleep -s 1
  }
} -ArgumentList $args




# �@�ݒ�ύX�A���񉻃u���b�N��[$j_MAX]�̌���iSNR�̕��񉻐� * �X�N���v�g�ɂ����񉻃u���b�N�� < �}�V���̃R�A���j
#		e:/LDPC_setting.R
#		e:/LDPC_setting.m
#
# �A���M�f�[�^�̍쐬�iMATLAB�o�͂�datasymbols.R����ƃf�B���N�g���ɕۑ��j
# 	��send_receive2.m
#
#	�B��M�f�[�^�̍쐬�i��Ė@�̌J��Ԃ�����������ʂɕ�������������ꍇ�ł�H��w�̍Č������m�ہA�������ɓǂݍ��ރo�C�i���t�@�C���T�C�Y�̓s���ŏ������ɕ������ĕ��񏈗����s���j			MIMO5_data.Rmd�̒���Turbo/LDPC�����̐؂�ւ���ݒ聚
#	PowerShell Set-ExecutionPolicy RemoteSigned
#	chcp 932
#	for ($j=1; $j -le [$j_MAX]; $j++){
#		& 'C:\Program Files\R\R-3.6.1\bin\Rscript.exe' --encoding=UTF-8 'C:\cygwin\home\hagijyun\D\���A�J���i�o�b�N�A�b�v+���j\MIMO\VB\conf\MIMO5_data.Rmd' "$j" "-1" "FIRST"
#	#	& 'C:\Program Files\R\R-4.1.2\bin\Rscript.exe' --encoding=UTF-8 'C:\cygwin64\home\hokudai\conf\MIMO5_data.Rmd' "$j" "-1" "FIRST"
#		Start-Sleep -s 1
#	}
#
# �C����+�J��Ԃ������E������
#		.\para5.ps1 1�`[$j_MAX]
#
# �D���ʂ̕��ω�
#		mean_numErrs.m
