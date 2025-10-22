# 前提
## Path
作業ディレクトリから以下が実行できるようにパスを通しておく
- Rscript.exe
- Rtools (make, gcc等)
- matlab.exe

## Working directory
以下のファイルでは作業ディレクトリにe:/を想定した記述がある
- _quick_examples.ps1
- MIMO5_para.Rmd
- send_receive2_para.m
- mean_numErrs.m

# Usage
それなりに実行時間かかるが、ログファイルparallell_log.txtで進行状況が確認できる

1. Download all files to the working directory.
2. 3GPP_LDPC_BG1.Z96とsend_receive2_para.mを、MATLABを起動したときのMATLABのデフォルトディレクトリ(WindowsであればC:\Users\ユーザー名\Documents\MATLAB)へ移動
3. WSHを管理者権限で実行し作業ディレクトリに移る
4. ``.\\_quick_examples.ps1 1 5''として実行
5. MATLABでmean_numErrs.mを実行(numErrs_mean0.csvに初発、numErrs_mean1to5.csvに繰り返し回数毎のBERが記録される)
