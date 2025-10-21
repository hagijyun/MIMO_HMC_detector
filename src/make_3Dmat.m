% make_3Dmat.m
%

%--- 設定 ---------------------------------------
m = 1;			% ★変調次数

% 各変調方式ごとに対応する SNR リスト（dB）
snrLists = { (-1):4, 6:11, 12:17 };
modTypes = {'QPSK','16QAM','64QAM'};


%--- 各変調方式ごとに処理 --------------------------------
% for m = 1:numel(modTypes)
    mt      = modTypes{m};
    snrList = snrLists{m};
    nSnr    = numel(snrList);

    % (1) 最初のファイルでサイズを取得
    fn0 = sprintf('e:/%s_SNR%ddB_softIQ.mat', mt, snrList(1));
    tmp = load(fn0, 'soft_iq');
    [rows, cols] = size(tmp.soft_iq);

    % (2) receivedSignal を予約 (rows × cols × nSnr)
    receivedSignal = zeros(rows, cols, nSnr);

    % (3) 各 SNR ファイルを読み込んで 3rd 次元に格納
    for k = 1:nSnr
        fn = sprintf('e:/%s_SNR%ddB_softIQ.mat', mt, snrList(k));
        d  = load(fn, 'soft_iq');
        receivedSignal(:,:,k) = d.soft_iq;
    end

    % (4) 方式名をファイル名にして保存
    outFn = sprintf('e:/OAMP1.mat', mt);	% mtは無視される
    save(outFn, 'receivedSignal');
% end
