# 低流强 S3 逐位置 q 候选

2026-11-02。完整包覆盖173轮；新改动仅低流强40轮ring15–24，5,197个run–ring–sector位置更新。参考0616固定，其他位置继承0244；继承不自动表示已修好。尚未切换pars正式入口或生产ROOT。

`coefficients.json`使用`sm150_s3_pixel_energy_v2`。低流强采用24×32的`q_pixel`，按0起始的ring/sector索引；它是总倍率，取代`q_position`，不是再乘一次。两面共用q，先于25–100MeV和正背面能差2.5MeV门。高流强继续用原`q_position`，仅ring1–8。

E0=g(c0+c1*ADC)，Er=q(r,s)*(Ar+Br*E0r)，Es=q(r,s)*(As+Bs*E0s)。源刻度与A/B和0244相同，必须搭配`source_calibration/`三个文件。Si array、Ge、时间、门和alpha保持。

`coefficients.csv`为完整A/B表；`q_0244_ring_baseline`只是原逐ring基线，实际新q见完整JSON及`s3_energy_attenuation/parameters/low_pixel_20261102/pixel_factors.csv`。

本轮候选在旧分段上拟合并冻结后，通过E盘新增1,400万C2条目的代表轮谱检验；原生build在合成事件及3轮各前200万真实条目上核验。并非全量生产重建，也不是全部位置的物理误差认证。详细局限与未解决谱见`jupyter/codex_0245_150sm_low_beam_s3_pixel_energy_correction.ipynb`。

维护程序：`s3_energy_attenuation/conditioned_calibrate.cpp`、`pixel_fit.h`、`refit_pixels.py`、`pixel_results.py`；事件应用维护于`build_all/build.cpp`。记录在`codex/analysis/s3_low_pixel_energy_20261102/`。
