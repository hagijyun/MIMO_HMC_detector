# ��MIMO�̐M�����o�v���O�����������s������PowerShell�X�N���v�g
#
#								�������F����i���� + �J�ԁj
#
# based on para5.ps1


# �O����
PowerShell Set-ExecutionPolicy RemoteSigned
chcp 932

# PATH��ǉ�
#	$env:PATH += ";C:\Program Files\R\R-4.1.2\bin"					# Rscript
#	$env:PATH += ";C:\rtools40\usr\bin"     								# Rtools(make �Ȃ�)
#	$env:PATH += ";C:\rtools40\mingw64\bin" 								# Rtools(gcc, g++ �Ȃ�)

& 'Rscript' --encoding=UTF-8 'stan_pre_compilation.R'

# �{�����i���񉻃u���b�N�j
Start-Job -ScriptBlock {
	param($j, $i_MAX = 5)

  # PATH��ǉ�
#	$env:PATH += ";C:\Program Files\R\R-4.1.2\bin"					# Rscript

  # ��ƃf�B���N�g����ݒ�
  Set-Location 'e:/'

	# ����
	& 'Rscript' --encoding=UTF-8 'MIMO5_para.Rmd' "$j" "0" "FIRST"
	Start-Sleep -s 1
	# ���R�����g�A�E�g�łǂꂩ���
#	& 'matlab' -batch "send_receive2_para($j, 0, {'MMSE', 'EP', 'HMC', 'MGS', 'MHGD', 'Lang'})"		# �S�Ẵp�^��
	& 'matlab' -batch "send_receive2_para($j, 0, {'MMSE', 'EP', 'HMC'})"													# MMSE+EP+HMC�̃p�^��
	Start-Sleep -s 1

	# �J��
  for ($i=1; $i -le $i_MAX; $i++){
		& 'Rscript' --encoding=UTF-8 'MIMO5_para.Rmd' "$j" "$i" "not_FIRST"
		Start-Sleep -s 1
		& 'matlab' -batch "send_receive2_para($j, $i, {'HMC'})"
		Start-Sleep -s 1
  }
} -ArgumentList $args
