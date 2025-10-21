% MATLABでのLDPC符号処理（受信側）
%
% C:\Users\hagijyun\Documents\MATLAB\send_receive2_para.m


function send_receive2_para(blk_ct, iter_ct, method_list)
%
% 前処理
%
rng(10);

dir_name = "e:";																																	% ●ファイル入出力のディレクトリ名
DEBUG_OUT = false;																																% ●デバグ用にHMCの誤りパタンの結果を保存する場合★

% スクリプトによる並列化ブロック数
blk_ct = sprintf("%d", blk_ct);

% 提案方式の繰り返し復号回数
iter_ct = sprintf("%d", iter_ct);

% ファイルから設定値を読み込む（テキストなので評価が必要）
fid = fopen(sprintf('%s/LDPC_setting.m', dir_name), 'rt'); codeFromFile = fread(fid, '*char')'; fclose(fid); eval(codeFromFile);

SYMB_bits = [2, 4, 6];																														% 1シンボル当たりのビット数

% LDPC符号（3GPP）
baseMatQC = load(sprintf("3GPP_LDPC_BG1.Z%d", Z));																% Matlabフォルダにあらかじめ保存しておく
pcmatrix = ldpcQuasiCyclicMatrix(Z, baseMatQC);

% LDPC符号（Gallager）
% load('gallager_ldpc.mat');																											% 系統的形式のパリティチェック行列を読み込む
% pcmatrix = logical(H);																													% スパース行列をlogical型に変換

% LDPC符号（IEEE）
% [cfgLDPCEnc, cfgLDPCDec] = generateConfigLDPC(1/2, 1944);												% IEEE (rates: 1/2, 2/3, 3/4, and 5/6, codeword lengths: 648, 1296, and 1944)

% 受信処理の準備
cfgLDPCDec = ldpcDecoderConfig(pcmatrix);																					% 3GPP, Gallager

numIter = 5;																																			% 内部繰り返し数

long_frame = ceil(cfgLDPCDec.BlockLength ./ (96 * SYMB_bits));										% 1符号長がまたがる送信タイミング数（変調次数により異なる）
padding_len = 96 * SYMB_bits(modOrder) * long_frame(modOrder) - cfgLDPCDec.BlockLength;
blkLength = cfgLDPCDec.NumInformationBits;																				% 情報パケット長[ビット]

% パディング分の読み込み
load(sprintf('%s/padding_bits.mat', dir_name));																		% padding_bitsを読み込み
padding_llr = 1 - 2 * padding_bits;																								% padding_bitsに相当するLLR（0→+1.0, 1→-1.0）

% デインターリーブインデックスを読み込み
load(sprintf('%s/intrIndices.mat', dir_name), '-mat');


%
% 本処理
%
	method_list = string(method_list);
%	method_list = ["HMC", "EP", "MMSE", "MGS", "MHGD", "Lang"];											% すべての場合
%	method_list = ["HMC", "EP", "MMSE"											 ];											% HMC+EP+MMSEの場合
method_MAX = length(method_list);

trial_MAX = long_frame(modOrder) * TURBO_trial_MAX;																% 総試行回数

% SNRの設定
snr_list = SNR_LIST;
snr_MAX = length(snr_list);

% ビットLLRの格納領域
LLR = zeros(long_frame(modOrder) * SYMB_bits(modOrder)*96, trial_MAX/long_frame(modOrder), snr_MAX);

% LLRダンピング用に旧LLRを読み込み
LLR_old = LLR;
if (0 < str2double(iter_ct)) load(sprintf('%s/LLR%s_dumping.mat', dir_name, blk_ct)); end

% 誤り率の格納領域
numErrs = zeros(method_MAX, trial_MAX/long_frame(modOrder), snr_MAX);

% 誤り率の判定のために、元データのビット乱数を変数dataに読み込み
load(sprintf('%s/databits.mat', dir_name), '-mat')

% デバグ用に保存する場合
if (DEBUG_OUT) fid = fopen(sprintf('%s/DEBUG_coded_biterr_ptn%s%s_%s.log', dir_name, blk_ct, iter_ct, datestr(now, 'yyyymmddHHMMSS')), 'a'); end

% 誤り数を計算
for method = 1:method_MAX																													% 手法毎
	% 軟判定シンボルを変数receivedSignalに読み込み
	load(sprintf('%s/%s%s.mat', dir_name, method_list(method), blk_ct))

	for snr = 1:snr_MAX																															% SNR毎
		% 前処理
		snrdB_modified = snr_list(snr) + 10 * log10(cfgLDPCDec.NumInformationBits/cfgLDPCDec.BlockLength);																					% LLRの適切なスケーリングのために、demodSignalに渡す雑音分散に符号化率を考慮
		% 本処理
		for trial = 1:long_frame(modOrder):trial_MAX																	% 符号語に相当する試行回数の塊毎
			% QAM復調
			demodSignal = qamdemod(reshape(receivedSignal(:, trial:(trial + long_frame(modOrder)-1), snr), [], 1), 2^SYMB_bits(modOrder), OutputType = 'llr', UnitAveragePower = true, NoiseVariance = (10.^(-snrdB_modified/10)));
			demodSignal = demodSignal(1:cfgLDPCDec.BlockLength);																																											% 0パディング分を除いて考える
%			demodSignal = deintrlv(demodSignal, intrIndices);																																													% ●ビットデインターリービング（未使用）

			% LDPC復号
			% 硬判定情報
			receivedBits = ldpcDecode(demodSignal, cfgLDPCDec, numIter, DecisionType = 'hard', OutputFormat = 'info');
			% 軟判定情報
			if (method_list(method) == "HMC")
				tmp_LLR = ldpcDecode(demodSignal, cfgLDPCDec, numIter, DecisionType = 'soft', OutputFormat = 'whole'); % - demodSignal;									% 外部LLRとするにはemodSignalを引く
				clipping_level = max(abs(tmp_LLR(:))) / 1;																																															% LLRのクリッピング用、「/ 1」でクリッピングなし
				clipped_LLR = tmp_LLR;															% min(max(tmp_LLR, -clipping_level), clipping_level);
				d = cfgLDPCDec.BlockLength;													% cfgLDPCDec.NumInformationBits;
				LLR(				1:d,																								(trial - 1)/long_frame(modOrder) + 1, snr) = clipped_LLR(1:d);					% intrlv(clipped_LLR, intrIndices);	% clipped_LLR(1:d);		% ●ビットデインターリービングを適用する場合は切り替え（未使用）
				if (padding_len ~= 0)
					LLR((d + 1):(long_frame(modOrder) * SYMB_bits(modOrder)*96),	(trial - 1)/long_frame(modOrder) + 1, snr) = clipping_level/10 * padding_llr;					% パディング分のLLRを追加
				end
			end

			% 誤り率の計算
			[biterr_num, ~, biterr_ptn] = biterr(data, receivedBits);
			numErrs(method, (trial - 1)/long_frame(modOrder) + 1, snr) = biterr_num / blkLength;

			% ●デバグ用に保存する場合
			if (DEBUG_OUT && contains(method_list(method), "HMC"))
				fprintf(fid, '%s\n', sprintf("method:%s, snr:%02d, trial:%04d, err_num:%04d\n%s", method_list(method), snr_list(snr), trial, biterr_num, sprintf("%d", biterr_ptn)));
			end
		end
	end
end
LLR = -LLR;																																																																			% Turbo符号とあわせて符号を逆転

% デバグ用に保存する場合
if (DEBUG_OUT) fclose(fid); end

% ●LLRの後加工（未使用）
% 対数圧縮：0付近をα/T倍に強調・Tより大きな値は減らす
% alpha =  50.0; T =  50.0; LLR = alpha * sign(LLR) .* log(1 + 1/T * abs(LLR));
% alpha =  50.0; T =  50.0; LLR = alpha * tanh(1/T * LLR);
% 指数伸長：0付近をα倍に強調・Tより大きな値はそのまま
% alpha = 3; T = 10; LLR = LLR .* (1 + alpha * exp(-(abs(LLR)/T).^2));

% ダンピング（実質未使用）
beta = 0.0; LLR = beta * LLR_old + (1 - beta) * LLR;
LLR_old = LLR;
save(sprintf('%s/LLR%s_dumping.mat', dir_name, blk_ct), 'LLR_old');


%
% 後処理
%
% numErrsのファイル出力（試行回数を通じた平均誤り率）
numErrs_tmp = mean(numErrs, 2);
numErrs_out = reshape(numErrs_tmp, method_MAX, snr_MAX)';
save(sprintf('%s/numErrs%s%s.mat', dir_name, blk_ct, iter_ct), 'numErrs_out');

% LLRのファイル出力（複素数の実数展開に対応して、IとQを離して配置）
LLR_tmp = reshape(LLR    , SYMB_bits(modOrder)  * 96  , trial_MAX, snr_MAX);
LLR_tmp = reshape(LLR_tmp, SYMB_bits(modOrder)/2, 96*2, trial_MAX, snr_MAX);
LLR_out = [LLR_tmp(:, 1:2:(96*2), :, :), LLR_tmp(:, 2:2:(96*2), :, :)];
save(sprintf('%s/LLR%s.mat', dir_name, blk_ct), 'LLR_out');

end
