# 前提
## Path
作業ディレクトリからRscript.exe, matlab.exeが実行できるようにパスを通しておく

## MATLABのスクリプト
send_receive2_para.m and mean_numErrs.m suppose working directory as e:\.


# Usage
それなりに実行時間かかるが、ログファイルparallell_log.txtで進行状況が確認できる

1. Download all files to the working directory.
2. 3GPP_LDPC_BG1.Z96とsend_receive2_para.mを、MATLABを起動したときのMATLABのデフォルトディレクトリ(WindowsであればC:\Users\ユーザー名\Documents\MATLAB)へ移動
3. WSHを管理者権限で実行し作業ディレクトリに移る
4. _quick_examples.ps1を実行
5. MATLABでmean_numErrs.mを実行
