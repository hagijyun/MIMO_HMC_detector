% MatlabでnumErrs??.matの平均化
%
% 	mean_numErrs.m


% 前処理
dir_name = "e:";																			% ●ファイル入出力のディレクトリ名

PARA		= 1;																					% 並列ブロック数
SNR			= 6;																					% SNR数
METHOD	= 3;																					% 初発における復号法の種類（HMC, EP, MMSE, 			MGS, MHGD, Lang）
ITER		= 5;																					% 次発における繰返数


% 初発
first_data = zeros(SNR, METHOD, PARA);

for para = 1:PARA
  load(sprintf("%s/numErrs%d0.mat", dir_name, para));	% ファイルからデータを読み込む
  first_data(:, :, para) = numErrs_out;							  % 読み込んだデータを配列に追加
end

first_mean = mean(first_data, 3);											% 並列試行を通じた平均値を計算
disp(first_mean);
fname = sprintf("%s/numErrs_mean0.csv", dir_name);
writecell({'MMSE','EP','HMC#0'}, fname);
writematrix(first_mean, fname, 'WriteMode', 'append');


% 次発
next_data = zeros(SNR, ITER, PARA);

for para = 1:PARA
	for iter = 1:ITER
    load(sprintf("%s/numErrs%d%d.mat", dir_name, para, iter));	% ファイルからデータを読み込む
    next_data(:, iter, para) = numErrs_out;					  % 読み込んだデータを配列に追加
	end
end

next_mean = mean(next_data, 3);												% 並列試行を通じた平均値を計算
disp(next_mean);
fname = sprintf("%s/numErrs_mean1to%d.csv", dir_name, ITER);
writecell(cellstr(compose("HMC#%d", 1:ITER)), fname);
writematrix(next_mean, fname, 'WriteMode', 'append');
