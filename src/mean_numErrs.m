% MatlabでnumErrs??.matの平均化
%
% 	mean_numErrs.m


% 前処理
PARA		= 2;																					% 並列ブロック数
SNR			= 6;																					% SNR数
METHOD	= 3;																					% 初発における復号法の種類（HMC, EP, MMSE, MGS, MHGD, Lang）
ITER		= 5;																					% 次発における繰返数


% 初発
first_data = zeros(SNR, METHOD, PARA);

for para = 1:PARA
  load(sprintf("e:/numErrs%d0.mat", para));						% ファイルからデータを読み込む
  first_data(:, :, para) = numErrs_out;							  % 読み込んだデータを配列に追加
end

first_mean = mean(first_data, 3);											% 並列試行を通じた平均値を計算
disp(first_mean);
%	writematrix(first_mean, sprintf("e:/numErrs_mean-1.csv"));
	writematrix(first_mean, sprintf("e:/numErrs_mean0.csv" ));


% 次発
next_data = zeros(SNR, ITER, PARA);

for para = 1:PARA
	for iter = 1:ITER
    load(sprintf("e:/numErrs%d%d.mat", para, iter));	% ファイルからデータを読み込む
    next_data(:, iter, para) = numErrs_out;					  % 読み込んだデータを配列に追加
	end
end

next_mean = mean(next_data, 3);												% 並列試行を通じた平均値を計算
disp(next_mean);
writematrix(next_mean, sprintf("e:/numErrs_mean1to%d.csv", ITER));
