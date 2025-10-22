# プリコンパイル専用スクリプト
# stan_pre_compilation.R

# ライブラリ
library(rstan)

# auto_write設定
rstan_options(auto_write = TRUE)

cat("Stanモデルのプリコンパイルを開始...\n")

# 全モデルをコンパイル（*.rdsファイルを生成）
cat("1/5: MCMC4f.stan\n")
stan_model(file = "MCMC4f.stan")

cat("2/5: MCMC4m.stan\n")
stan_model(file = "MCMC4m.stan")

cat("3/5: stan_MGS4e.stan\n")
stan_model(file = "stan_MGS4e.stan")

cat("4/5: stan_MHGD1a.stan\n")
stan_model(file = "stan_MHGD1a.stan")

cat("5/5: stan_Lang1a.stan\n")
stan_model(file = "stan_Lang1a.stan")

cat("完了！ *.rdsファイルが生成されました。\n")
