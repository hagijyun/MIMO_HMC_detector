% LDPC符号用に頻繁に変更する設定一覧
%
% LDPC_setting.m

modOrder = 3;								% ●変調次数：QPSK = 1, 16QAM = 2, 64QAM = 3★
SNR_LIST = (12):17;					% ●SNRの範囲★
TURBO_trial_MAX = 50;				% ●LDPC符号として成立するフレームを何回試行？（おおよそ10～50）★


% 3GPP符号の場合
Z = 96;											% 符号長[ビット]：32 -> 2176, 64 -> 4352, 88 -> 5984, 96 -> 6528
% Gallager符号, IEEE符号の場合
% 特になし
