# OAMP_Net2_LMMSE_QPSK_Ma5.py
#
#		OAMP-Net2 MIMO detector - Extended version
#		Original: https://github.com/hehengtao/OAMP-Net
#		Original authors: H. He, C. Wen, S. Jin, and G. Y. Li
#		Extensions: Support for higher-order modulations and channel coding
#
#		BER・シンボル出力対応、ファイル入力対応
#
# トレーニングモード（計算長い、学習パラメータをファイル保存）
#		python OAMP_Net2_LMMSE_QPSK_Ma5.py --snrdb 10.0 --mod_type QPSK --train
#
# 　　　学習済モード（学習パラメータをファイル入力）★
#		python OAMP_Net2_LMMSE_QPSK_Ma5.py --snrdb 10.0 --mod_type QPSK --input_file datasymbols_1.py
#			0. 前処理
#				a. 学習済みのmatファイルを準備
#				b. LDPC_setting.R, LDPC_setting.mを設定（Turbo_trial, TURBO_trial_MAX = 100）
#				c. ★send_receive2.m
#				d. & 'C:\Program Files\R\R-3.6.1\bin\Rscript.exe' --encoding=UTF-8 'C:\cygwin\home\hagijyun\C\IEEETran\投稿\20250218\OAMP\gray_preprocess.R'
#			1. 本処理
#				a. OAMP_Net2_LMMSE_QPSK_Ma5.py							(batch_size = 12																		# ★LONG_FRAMEで置き換える		32 16 11 for Turbo, 34 17 12 for LDPC	)
#				b. OAMP_Net2_LMMSE_QPSK_Ma5.ps1							($mod_types = @("64QAM")														# ★変調次数																												)
#			2. 後処理
#				a. make_3Dmat.m 														(m = 3;			% ★変調次数)
#				b. & 'C:\Program Files\MATLAB\R2024b\bin\matlab.exe' -batch "send_receive2_para(1, 0, {'OAMP'})"
#


import time
import tensorflow as tf
import numpy as np
import scipy.io as sc 
import os
from scipy.linalg import toeplitz
os.environ['CUDA_VISIBLE_DEVICES'] = '0'
reuse=tf.compat.v1.AUTO_REUSE


import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--train", action = "store_true")
parser.add_argument("--snrdb", type=float, default=10.0)
parser.add_argument("--mod_type", type=str, default="QPSK")
parser.add_argument("--output_file", type=str, default="soft_iq.mat")																												# 指定しない（実質使用されないため）
parser.add_argument("--input_file", type=str, default=None, help="Input file containing symbol indices")
parser.add_argument("--recycle_count", type=int, default=100, help="Number of times to recycle the input symbols")					# 指定しない（デフォルトの100を常用）、Turbo_trialと同じ値
args = parser.parse_args()


#parameters
num_trials = 50						# trainなしの場合の、SER導出の固定の試行回数（50（num_trials）× 100（batch_size）= 5,000回のシンボルベクトル送信）
mod_type = args.mod_type	# 'QPSK', '16QAM', '64QAM'
N = 96											# ●デバグ用に変更	96
M = 96											# ●デバグ用に変更	96
#T=M
snrdb_train=np.array([args.snrdb], dtype=np.float64) # training SNR
snr_train = 10.0 ** (snrdb_train/10.0)  #train_SNR_linear
batch_size = 100
epochs = 100
itermax = 10

# ファイル入力モードなら上書き
if args.input_file:
    num_trials = args.recycle_count				# Turbo_trialと同じ
    batch_size = 12												# ★LONG_FRAMEで置き換える		32 16 11 for Turbo, 34 17 12 for LDPC
    print(f"[input_file] batch_size={batch_size}, num_trials={num_trials}")

train_size = 5000
valid_size = 10000
errsum = 1000  #计算错误总数
rho=0.00	# 0.50
snrdb_test=snrdb_train
snr_test = 10.0 ** (snrdb_test/10.0)  #train_SNR_linear
weight_mat = f"{args.mod_type}_SNR{args.snrdb:g}_Rayleigh MIMO"		# 'Rayleigh MIMO'

Rx = rho ** np.arange(M)																					#★np.logspace( 0, np.log10( rho**(M-1)),M)

Tx = rho ** np.arange(N)																					#★np.logspace( 0, np.log10( rho**(N-1)),N)

sigma2=1/snr_test

# ファイル入力に伴う追加（グローバル変数定義）
input_symbols = None
symbol_index = 0


# ===== MATLAB Gray 互換 ======================
def _gray2bin(arr):
    b = arr.copy()
    for sh in (1, 2, 4, 8, 16):
        b ^= b >> sh
    return b

# index→bit 変換
def _idx2bits(idx, mod_type):
    """
    MATLAB qammod/qamdemod ('gray') と同一のビット並びを返す。
    idx : (batch, N)  シンボル index
    戻り : (batch, N*log2M) ビット列
    """
    M   = len(get_constellation(mod_type)) ** 2   # 総シンボル数
    bps = int(np.log2(M))                         # bits / symbol
    L   = int(np.sqrt(M))                         # 軸ごとのレベル数

    # 実軸・虚軸に分離
    i_idx = idx // L
    q_idx = idx  % L

    i_bin = i_idx																	# _gray2bin(i_idx)
    q_bin = q_idx																	# _gray2bin(q_idx)

    bpi = bps // 2                                # 軸あたりビット数
    i_bits = ((i_bin[..., None] >> np.arange(bpi-1, -1, -1)) & 1)
    q_bits = ((q_bin[..., None] >> np.arange(bpi-1, -1, -1)) & 1)

    bits = np.concatenate([i_bits, q_bits], axis=-1)  # [I|Q]
    return bits.reshape(idx.shape[0], idx.shape[1] * bps)

# インデックスからシンボルへの変換関数の追加
def indices_to_symbols(indices, mod_type):
    constellation = get_constellation(mod_type)
    L = len(constellation)
    i_idx = indices // L
    q_idx = indices % L
    real_part = constellation[i_idx]
    imag_part = constellation[q_idx]
    return real_part + 1j * imag_part
# ===== MATLAB Gray 互換 ======================

# データ読み込みと拡張関数の追加
def load_recycled_symbols(filename, mod_type, recycle_count):
    with open(filename, 'r') as f:
        content = f.read().strip()
        indices = np.array([int(x) for x in content.split(',') if x])
    print(f"Loaded {len(indices)} symbols")
    symbols = indices_to_symbols(indices, mod_type)
    recycled_symbols = np.tile(symbols, recycle_count)
    print(f"Recycled to {len(recycled_symbols)} symbols")
    return recycled_symbols


def Network_ini(theta):
    update=[]
    for var in tf.compat.v1.trainable_variables():
#       update.append(tf.compat.v1.assign(tf.compat.v1.get_default_graph().get_tensor_by_name(var.name),tf.constant(np.reshape(theta[var.name],var.shape))))		# 旧コード
        key = var.name.replace('/', '_').replace(':', '_')        # 例: "gamma:0" → "gamma_0" のように置換
        new_value = tf.constant(np.reshape(theta[key], var.shape)); update.append(var.assign(new_value))																										# 新コード：直接 var に対して assign する
    return update

def Variable(shape):
    gamma=tf.compat.v1.get_variable('gamma', shape=shape, initializer=tf.ones_initializer,dtype='float64')
    return gamma
    
def Variable1(shape):    
    gamma1=tf.compat.v1.get_variable('gamma1', shape=shape, initializer=tf.ones_initializer,dtype='float64')
    return gamma1 
   
def Variable2(shape):    
    corre=tf.compat.v1.get_variable('corre', shape=shape, initializer=tf.ones_initializer,dtype='float64')
    return corre  

def Variable3(shape):    
    expec=tf.compat.v1.get_variable('expec', shape=shape, initializer=tf.zeros_initializer,dtype='float64')
    return expec     
    
def generate_data_iid_test(B,M,N,SNR_dB):  
    global input_symbols, symbol_index				# グローバル変数の追加
    y_real=np.zeros([2*M,B])   
    x_real=np.zeros([2*N,B])
    H_real=np.zeros([B,2*M,2*N])    
#    snrdb_test_inv=10.**(-snrdb_test/10.)
    snr = 10.0 ** (SNR_dB/10.0)  #train_SNR_linear
    sigma2=1/snr
    H_=np.zeros([B,M,N],dtype=complex)      
#    H_ = (np.random.randn(B,M,N)+1j*np.random.randn(B,M,N))
#    H_ = (np.random.randn(B,M,N)+1j*np.random.randn(B,M,N))/np.sqrt(2*M)										#★コメントアウト
    for i in range(B):																																			#★ループ
      H_[i,:,:] = Correlated(N, M, Rx, Tx) / np.sqrt(M)  																		#★以前のスケール(1/√M)に合わせる
    x_=np.zeros([N,B],dtype=complex)
# ファイル入力に伴う追加
    if input_symbols is not None:	    # 入力シンボルがある場合はそれを使用
        symbols_needed = N * B
        if symbol_index + symbols_needed > len(input_symbols):
            symbol_index = 0
        used_symbols = input_symbols[symbol_index:symbol_index+symbols_needed]
        symbol_index += symbols_needed
#       x_ = used_symbols.reshape(N, B)
        x_ = used_symbols.reshape(B, N).T		# データの並びが間違っていたため修正
    else:													    # 従来通りランダム生成
        x_ = generate_symbols(N, B, mod_type)
# ファイル入力に伴う追加
    y_=np.zeros([M,B],dtype=complex)
    w=np.sqrt(1/2)*(np.random.randn(M,B)+1j*np.random.randn(M,B))*np.sqrt(sigma2)
    for i in range(B):
      H=H_[i,:,:]
      y_[:,i]=np.matmul(H,x_[:,i])+w[:,i]
      y_real[:,i]=np.hstack((np.real(y_[:,i]), np.imag(y_[:,i]))).T #stack 要注意加()
      H_real[i,:,:]=np.vstack((np.hstack((np.real(H),-np.imag(H))),np.hstack((np.imag(H),np.real(H)))))
      x_real[:,i]=np.hstack((np.real(x_[:,i]),np.imag(x_[:,i]))).T
    return y_real,H_real,x_real
    
def Correlated(N,M,Rx,Tx):  
    T=toeplitz(Tx)
#    h=np.math.sqrt(1/2)*(np.random.normal(0,1,(M*N,1))+1j*np.random.normal(0,1,(M*N,1)))		#★コメントアウト
    R=toeplitz(Rx)
    Lr = np.linalg.cholesky(R)            																									#★効率向上 M×M
    Lt = np.linalg.cholesky(T)             																									#★効率向上	N×N
    W  = (np.random.randn(M,N)+1j*np.random.randn(M,N))/np.sqrt(2)													#★効率向上
    return Lr @ W @ Lt.T               				    # Tが実対称なら .T でOK（一般なら .conj().T）		#★効率向上
#    C=np.kron(T,R)																																					#★コメントアウト
#    AAA=np.matmul(np.linalg.cholesky(C),h)																									#★コメントアウト
#    A=np.reshape(AAA,(M,N))																																#★コメントアウト
#    return A    																																						#★コメントアウト

# def Modulation(bits):                                        
#     x_re= (2*bits-1)/np.sqrt(2)               
#     return x_re
def get_constellation(mod_type):
    if mod_type == 'QPSK':
        # QPSKの場合、実部・虚部ともに [-1, 1] を正規化
        return np.array([-1, 1]) / np.sqrt(2)
    elif mod_type == '16QAM':
        # 16QAMの場合、実部・虚部ともに [-3, -1, 1, 3] を正規化
        return np.arange(-3, 4, 2) / np.sqrt(10)
    elif mod_type == '64QAM':
        # 64QAMの場合、実部・虚部ともに [-7, -5, -3, -1, 1, 3, 5, 7] を正規化
        return np.arange(-7, 8, 2) / np.sqrt(42)
    else:
        raise ValueError("Unsupported modulation type")

def generate_symbols(N, B, mod_type):
    constellation = get_constellation(mod_type)
    # 各シンボルは、指定されたコンステレーションから無作為に選択
    real_part = np.random.choice(constellation, size=(N, B))
    imag_part = np.random.choice(constellation, size=(N, B))
    x_sym = real_part + 1j * imag_part
    return x_sym

#def shrink_bg_QPSK(r,rvar):
#    b=2*phi(-1/np.sqrt(2),r,rvar)    
##    tp = phi(-1/np.sqrt(2),r,rvar) + phi(1/np.sqrt(2),r,rvar) + eps
#    tp = phi(-1/np.sqrt(2),r,rvar) + phi(1/np.sqrt(2),r,rvar)
#    a=1-tf.compat.v1.div(b,tp)
##    print(b,tp,a)
#    xhat_= a/np.sqrt(2)
#    return (xhat_)    
def shrink_bg_QPSK(r, rvar):
    # 例: r: shape=(batch_size, 2*N, 1) とする場合、rvar を (batch_size, 2*N, 1) に合わせる
    rvar_expanded = tf.tile(tf.reshape(rvar, [-1, 1, 1]), [1, tf.shape(r)[1], 1])
    constellation = get_constellation(mod_type)
    
    likelihoods = []
    for point in constellation:
        likelihood = tf.exp(-tf.square(r - point) / (2 * rvar_expanded)) + 1e-10
        likelihoods.append(likelihood)
    likelihoods = tf.stack(likelihoods, axis=-1)
    weights = likelihoods / tf.reduce_sum(likelihoods, axis=-1, keepdims=True)
    constellation_tensor = tf.constant(constellation, dtype=tf.float64)
    xhat = tf.reduce_sum(weights * constellation_tensor, axis=-1)
    return xhat

def phi(x,r,rvar):
    rvar=tf.tile(tf.expand_dims(tf.expand_dims(rvar, axis=-1), axis=-1), [1,2*N,1])
#    print(r)
    y=tf.maximum(tf.exp(-tf.square(x-r)/(2*rvar)),eps)
#    print(rvar,'yyy')
    return y

def Train_batch(sess):
    _loss = list()
    _ser = list()
    packet=valid_size//batch_size
    for offset in range(packet):
        batch_Y, batch_H, batch_X = generate_data_iid_test(batch_size,M,N,snrdb_train)
        _, b_loss, b_ser = sess.run([optimizer,cost,ser], {Y_: batch_Y.T, A_: batch_H, X_: batch_X.T})
        _loss.append(b_loss)
        _ser.append(b_ser)
    return np.mean(_loss), np.mean(_ser)
    
def Valid_batch(sess):
    _loss = []
    _ser = []
    packet=valid_size//batch_size
    for offset in range(packet):
        batch_Y, batch_H, batch_X = generate_data_iid_test(batch_size,M,N,snrdb_test)
        b_loss, b_ser = sess.run([cost,ser], {Y_: batch_Y.T, A_: batch_H, X_: batch_X.T})
        _loss.append(b_loss)
        _ser.append(b_ser)
    return np.mean(_loss), np.mean(_ser)
    
def Train():
    print("\nTraining ...")     
    with tf.compat.v1.Session() as sess:        
        tf.compat.v1.global_variables_initializer().run()              
        weight_file=weight_mat+'oamp.mat'
        best_valid_loss=0    
        for i in range(epochs): 
            start_time = time.time()          
            train_loss, ser_train = Train_batch(sess)
            valid_loss, ser_valid = Valid_batch(sess)
            time_taken = time.time() - start_time
            print("Epoch %d Valid Loss: %.8f, Valid SER: %.8f, Time Cost: %.2f s"
                  % (i+1,valid_loss, ser_valid, time_taken))
            if i==0 or (i>0 and valid_loss < best_valid_loss):
                best_valid_loss = valid_loss #保存最优的网络loss
        #        print(best_valid_loss)
                Save(weight_file)																			# トレーニングモードでは、後のために学習したパラメータをファイルに保存
        print("\nTraining is finished.") 
        
def Test_batch(sess, ebn0_test):
		_loss = []
		_ser = []
		# BER対応
		all_soft_iq = []
		err_total   = 0
		C  = get_constellation(mod_type)
		constellation = (C[:, None] + 1j*C).flatten()    # shape=(M,)
		bits_sent = 0          													 # 追加: 送信ビット数
		bps = int(np.log2(len(constellation)))  				 # 追加：1シンボル当たりのビット数
		# BER対応

		if args.train:
			sernum = 0
			while sernum < errsum:
				batch_Y, batch_H, batch_X= generate_data_iid_test(batch_size,M,N,snrdb_test)
				b_loss, b_ser, errnum, s_, gamma_, gamma1_, corre_, expec_ = sess.run([cost,ser,err,s, gamma, gamma1, corre, expec], {Y_: batch_Y.T, A_: batch_H, X_: batch_X.T})
				_loss.append(b_loss)
				_ser.append(b_ser)      
				# BER対応
				soft_iq = s_[:, :N] + 1j * s_[:, N:2*N]
				all_soft_iq.append(soft_iq)
				bx = batch_X.T                							# shape (B, 192) = (100, 192)
				tx = bx[:, :N] + 1j * bx[:, N:2*N]   				# shape (100, 96)
				tx_idx = np.argmin(np.abs(tx[..., None]      - constellation)**2, axis=2)
				rx_idx = np.argmin(np.abs(soft_iq[..., None] - constellation)**2, axis=2)
				err_total += np.sum(_idx2bits(tx_idx, mod_type) != _idx2bits(rx_idx, mod_type))
				bits_sent += bps * N * batch_Y.shape[1]			# 追加: 送信ビット数
				# BER対応
				sernum=sernum+errnum
		else:
			for _ in range(num_trials):
				batch_Y, batch_H, batch_X= generate_data_iid_test(batch_size,M,N,snrdb_test)
				b_loss, b_ser, errnum, s_, gamma_, gamma1_, corre_, expec_ = sess.run([cost,ser,err,s, gamma, gamma1, corre, expec], {Y_: batch_Y.T, A_: batch_H, X_: batch_X.T})
				_loss.append(b_loss)
				_ser.append(b_ser)      
				# BER対応
				soft_iq = s_[:, :N] + 1j * s_[:, N:2*N]
				all_soft_iq.append(soft_iq)
				bx = batch_X.T                							# shape (B, 192) = (100, 192)
				tx = bx[:, :N] + 1j * bx[:, N:2*N]   				# shape (100, 96)
				tx_idx = np.argmin(np.abs(tx[..., None]      - constellation)**2, axis=2)
				rx_idx = np.argmin(np.abs(soft_iq[..., None] - constellation)**2, axis=2)
				err_total += np.sum(_idx2bits(tx_idx, mod_type) != _idx2bits(rx_idx, mod_type))
				bits_sent += bps * N * batch_Y.shape[1]			# 追加: 送信ビット数
				# BER対応

#		return np.mean(_loss), np.mean(_ser), s_, gamma_, gamma1_, corre_, expec_
		return (np.mean(_loss), np.mean(_ser), err_total, bits_sent, s_, gamma_, gamma1_, corre_, expec_, all_soft_iq)

def Test(snr_test):
    with tf.compat.v1.Session() as sess: 
        tf.compat.v1.global_variables_initializer().run()
        sess.run(update)
#       loss_test, ser_test,s_, gamma_, gamma1_, corre_, expec_ = Test_batch(sess, snrdb_test)
        loss_test, ser_tmp, err_total, bits_sent, s_, gamma_, gamma1_, corre_, expec_, all_soft_iq = Test_batch(sess, snrdb_test)
				# BER対応
        ber = err_total / bits_sent																																									# BERの分母に戻り値を活用
        sc.savemat(f"{args.mod_type}_SNR{int(args.snrdb)}dB_softIQ.mat", {"soft_iq": np.vstack(all_soft_iq).T})			# ●デバグ時はファイル出力をコメントアウト
				# BER対応
#       print("Eb/N0: %.0f, Testing Loss: %.8f, Testing SER: %.8f"% (snr_test, loss_test, ser_test))
        print("SNR: %d dB  Loss: %.8f  BER: %.8f" % (snr_test, loss_test, ber))
#   return loss_test, ser_test, s_, gamma_, gamma1_, corre_, expec_
    return loss_test, ber, s_, gamma_, gamma1_, corre_, expec_

def Save(weight_file):
    dict_name={}
    for varable in tf.compat.v1.trainable_variables():  
        key = varable.name.replace('/', '_').replace(':', '_')        # 例: "gamma:0" → "gamma_0" のように置換
        dict_name[key]=varable.eval()
    sc.savemat(weight_file, dict_name)     

with tf.compat.v1.Graph().as_default():
    #tensorflow placeholders, the input given to the model in order to train and test the network
    A_ = tf.compat.v1.placeholder(tf.float64,shape=[None,2*M,2*N])
    X_ = tf.compat.v1.placeholder(tf.float64,shape=[None,2*N])
    Y_ = tf.compat.v1.placeholder(tf.float64,shape=[None,2*M])

# ファイル入力に伴う追加
    if args.input_file:
        input_symbols = load_recycled_symbols(args.input_file, mod_type, args.recycle_count)
# ファイル入力に伴う追加

    s=tf.zeros((batch_size,2*N,1),dtype=tf.float64)#给s多加一个维度
    damping=0.1
    tau2=1
    sigma2=10.**(-snrdb_test/10.)
    eps=1e-10
    v2=tf.ones((batch_size,),dtype=tf.float64)
    beta=5e-1
#    I=tf.eye(2*M,batch_shape=[batch_size])
    IM=tf.eye(2*M,batch_shape=[batch_size],dtype='float64')
    IN=tf.eye(2*N,batch_shape=[batch_size],dtype='float64')
    with tf.compat.v1.variable_scope('gamma', reuse=reuse):
        gamma=Variable((itermax,))
    with tf.compat.v1.variable_scope('gamma1', reuse=reuse):
        gamma1=Variable1((itermax,)) 
        
    with tf.compat.v1.variable_scope('corre'):
        corre=Variable2((itermax,)) 
        
    with tf.compat.v1.variable_scope('expec'):
        expec=Variable3((itermax,))                 
#    W=tf.matmul(A_,tf.matrix_inverse(tf.matmul(A_,A_,adjoint_b = True)+beta*I),adjoint_a = True)#LMMSE Matrix     
    for t in range(itermax):
        
        v2M=tf.tile(tf.expand_dims(tf.expand_dims(v2, axis=-1), axis=-1),[1,2*M,2*M])
        
        v2N=tf.tile(tf.expand_dims(tf.expand_dims(v2, axis=-1), axis=-1),[1,2*N,2*M])   
        
        RR=tf.compat.v1.matrix_inverse(tf.multiply(v2M, tf.matmul(A_,A_,adjoint_b = True))+sigma2*IM/2)   #  2M*2M
#        print(RR)
        W=tf.multiply(v2N, tf.matmul(A_, RR, adjoint_a = True))  #LMMSE Matrix  2N*2M
#        W=tf.matmul(v2*A_,tf.matrix_inverse(tf.matmul(v2*A_,A_,adjoint_b = True)+sigma2*tf.eye(N)),adjoint_a = True)  #LMMSE Matrix
        tr=tf.compat.v1.trace(tf.matmul(W,A_))
  #      print(tr)
        tr=tf.tile(tf.expand_dims(tf.expand_dims(tr, axis=-1), axis=-1),[1,2*N,2*M])
        
        W=2*N/tr*W 
        
        z = tf.expand_dims(Y_,-1)-tf.matmul(A_,s)
        
        r = s + gamma[t]*tf.matmul(W,z)
        
        v2 = tf.maximum(tf.compat.v1.div(tf.norm(z, axis=[-2,-1])**2-M*sigma2,tf.compat.v1.trace(tf.matmul(A_,A_,adjoint_a=True))), eps)
#        tau2 = v2/(2*N)*(2*N+(gamma[t]**2-2*gamma[t])*2*M)+gamma[t]**2*sigma2/(2*N)*tf.trace(tf.matmul(W,W,adjoint_b=True))
        B=IN-gamma1[t]*tf.matmul(W, A_)
        
        tau2 = v2/2/N*tf.compat.v1.trace(tf.matmul(B,B,adjoint_b=True))+gamma1[t]*gamma1[t]*sigma2/4/N*tf.compat.v1.trace(tf.matmul(W,W,adjoint_b=True))
        
        s = shrink_bg_QPSK(r ,tau2)
        
        s=corre[t]*(s- expec[t]*r) 
        
    s=s[:,:,0]

    cost  = tf.nn.l2_loss(s-X_)  # l2 loss function

    err_temp = tf.compat.v1.to_float(tf.not_equal(tf.sign(s),tf.sign(X_)))

    err=tf.reduce_sum(err_temp)

    ser = tf.reduce_mean(err_temp)
         
    #learning_rate=0.001
    with tf.compat.v1.variable_scope('opt', reuse=reuse):
        optimizer = tf.compat.v1.train.AdamOptimizer(0.0001).minimize(cost)

    #Training DetNet
    if args.train: Train()
    update=Network_ini(sc.loadmat(weight_mat+'oamp.mat'))
#   loss_test, ser_test, s_, gamma_, gamma1_, corre_, expec_ =Test(snrdb_test)
    loss_test, ber_test, s_, gamma_, gamma1_, corre_, expec_ = Test(snrdb_test)
