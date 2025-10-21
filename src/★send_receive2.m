% MATLABでのLDPC符号処理（送信側）
%
% C:\Users\hagijyun\Documents\MATLAB\★send_receive2.m


%
% 前処理
%
rng(10);

dir_name = "e:";																																	% ●ファイル入出力のディレクトリ名

% ファイルから設定値を読み込む（テキストなので評価が必要）
fid = fopen(sprintf('%s/LDPC_setting.m', dir_name), 'rt'); codeFromFile = fread(fid, '*char')'; fclose(fid); eval(codeFromFile);

SYMB_bits = [2, 4, 6];																														% 1シンボル当たりのビット数

% LDPC符号（3GPP）の下準備
baseMatQC = load(sprintf("3GPP_LDPC_BG1.Z%d", Z));																% Matlabフォルダにあらかじめ保存しておく
pcmatrix = ldpcQuasiCyclicMatrix(Z, baseMatQC);

% LDPC符号（Gallager）
% load('gallager_ldpc.mat');																											% 系統的形式のパリティチェック行列を読み込む
% pcmatrix = logical(H);																													% スパース行列をlogical型に変換

% LDPC符号（IEEE）
% [cfgLDPCEnc, cfgLDPCDec] = generateConfigLDPC(1/2, 1944);												% IEEE (rates: 1/2, 2/3, 3/4, and 5/6, codeword lengths: 648, 1296, and 1944)

% 送信処理の準備
cfgLDPCEnc = ldpcEncoderConfig(pcmatrix);																					% 3GPP, Gallager

long_frame = ceil(cfgLDPCEnc.BlockLength ./ (96 * SYMB_bits));										% 1符号長がまたがる送信タイミング数（変調次数により異なる）
padding_len = 96 * SYMB_bits(modOrder) * long_frame(modOrder) - cfgLDPCEnc.BlockLength;
blkLength = cfgLDPCEnc.NumInformationBits;																				% 情報パケット長[ビット]

intrIndices = randperm(cfgLDPCEnc.BlockLength)';																	% ビット単位のRandom interleaving
save(sprintf('%s/intrIndices.mat', dir_name), 'intrIndices', '-mat');							% インターリーブインデックスをファイル出力


%
% 本処理
%
% ビット乱数生成、誤り率計算のためにファイル出力
data = randi([0 1], blkLength, 1);
save(sprintf('%s/databits.mat', dir_name), 'data', '-mat');

% LDPC符号
encodedData = ldpcEncode(data, cfgLDPCEnc);
% encodedData = intrlv(encodedData, intrIndices);																	% ●ビットインターリービング（未使用）
padding_bits = randi([0 1], padding_len, 1);		% zeros(padding_len, 1);					% パディングビット (padding_len == 0でもエラーなし)
save(sprintf('%s/padding_bits.mat', dir_name), 'padding_bits', '-mat');
encodedData = [encodedData; padding_bits];

% ビットをシンボルに変換、Rで読み込むためにファイル出力
fileID = fopen(sprintf('%s/datasymbols.R', dir_name), 'w');
fprintf(fileID, '%d,', bit2int(encodedData, SYMB_bits(modOrder))+1);							% Rの配列は1始まりなので+1
fclose(fileID);





% ●MIMO5_para.Rmdでの処理
%
% 入力：MATLAB出力を読み込んで送信シンボルに設定
% 出力：軟判定シンボルをファイルに保存（方式名.mat）
