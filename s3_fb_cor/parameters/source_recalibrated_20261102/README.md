# 新源刻度配套 S3 A/B/q 候选（2026-11-02）

这是 0244 的完整候选，尚未切换 `pars/active.json`。全部173轮；高流强 ring1–8，低流强 ring1–24。高流强外环系数仅用于完整结构，禁止作为生产范围启用。

`coefficients.json` 与 `source_calibration/` 内三个基础刻度文件必须配套；只更新S3源系数，Si array保持现用系数。文本最后一列保留历史源拟合χ值，不是本次两点一次刻度的χ值。

E0=hardware_gain*(c0+c1*ADC)，Er=q[matched_ring]*(Ar+Br*E0r)，Es=q[matched_ring]*(As+Bs*E0s)。不再叠加旧正背面和衰减链。q作用于匹配ring的两面，不强迫不同ring等能。

全量系数和采纳状态见CSV；重拟未通过者保留等价换算。新基础刻度不独立决定高能绝对能标，参考轮0472/0616/0680的现用物理能标固定。数据为E盘保存的全173轮分段ADC配对样本，本次未扫描全量C2。

程序：`s3_fb_cor/recalibrate_source.cpp`、`source_refit.py`、`s3_energy_attenuation/refit_position.py`；完整日志/方案/核验/拟合在 `codex/analysis/s3_source_calibration_refit_20261102/`；阅读版 `jupyter/codex_0244_150sm_s3_new_source_calibration_abq.ipynb`。
