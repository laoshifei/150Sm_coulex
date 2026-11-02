# 低流强逐轮、逐 ring–sector 修正

2026-11-02。`pixel_factors.csv`完整列出40轮ring15–24的12,800项：5,197项本轮更新、320项固定参考、7,283项继承。使用q_selected，不能把q_fitted提案直接当已验证结果。保留所有拒绝原因；不需要改变的近参考位置和未解决的坏位置均可能在继承项中。

两面共用同一个q。完整A/B/q及配套源刻度位于`s3_fb_cor/parameters/low_pixel_20261102/`；`complete_package.json`保存路径和SHA。`run_metrics.csv`分别列历史检查和新E盘分段结果，D为固定总净计数累计谱距离，不是能量分辨率。

第一批20–90MeV内条件归一化存在移峰出窗后错配尾部的风险，保留为诊断；本包仅采用固定原始总净计数归一化的结果。q是相对0616同一像素的倍率，没有跨ring等能约束。未启用正式参数或生产ROOT。

阅读：`jupyter/codex_0245_150sm_low_beam_s3_pixel_energy_correction.ipynb`。计算源码在本目录上级，配套来源/核验在`codex/analysis/s3_low_pixel_energy_20261102/`。
