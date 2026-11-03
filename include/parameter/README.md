# 现用能量、时间、符合窗口与Doppler参数

现用数值文件直接放在本目录，不再套日期版本目录。184份 Ge 逐轮增益表集中在唯一的一层子目录 `run_k/`。

本目录提交的是已经采用的数值与校验清单；发布源码不启用新的Ge能标或效率候选。项目的pars/运行环境、原始数据和分析笔记仍保留本地，新环境运行时需配套提供pars/active.json、清单及实验配置。根目录的忽略规则保持，仅在本目录开放已采用的文本参数。

| 文件 | 用途 |
| --- | --- |
| [cali_ge.txt](cali_ge.txt) | Ge 基础能量刻度 |
| [run_k/](run_k/) | Ge 逐轮增益；保留184轮表，生产使用严格173轮名单 |
| [cali_si_r460.txt](cali_si_r460.txt) | 高流强Ⅰ的 Si array / S3 基础能量刻度 |
| [cali_si_r457_low_current.txt](cali_si_r457_low_current.txt) | 低流强的 Si array / S3 基础能量刻度 |
| [cali_si_r457_high_current.txt](cali_si_r457_high_current.txt) | 高流强Ⅱ的 Si array / S3 基础能量刻度 |
| [particle_source_parameters.csv](particle_source_parameters.csv) | 逐时期、通道的采用/保留明细 |
| [s3_coefficients.json](s3_coefficients.json) | 与上述新 S3 刻度配套的完整 A/B/q；程序实际读取 |
| [s3_ab_parameters.csv](s3_ab_parameters.csv)、[s3_q_parameters.csv](s3_q_parameters.csv) | A/B 和实际 q 的可读表 |
| [time_alignment.json](time_alignment.json) | `joint_zero_20261028` 的完整 Si/S3/Ge offset 与能量相关时间函数 |
| [time_alignment_metadata.json](time_alignment_metadata.json) | 时间参数正式采用记录 |
| [coincidence_windows.json](coincidence_windows.json) | 正式窄/宽prompt及左右random窗；0251第一次方案，原C++ build直接读取并一次应用 |
| [coincidence_windows.csv](coincidence_windows.csv) | 三时期逐Ge、Si/S3的240行窗口、宽度和alpha，供阅读核对 |
| [random_subtraction.json](random_subtraction.json) | 所有组附加c=1；与逐对窗宽alpha分开，供后续谱处理读取 |
| [doppler.json](doppler.json) | 已预先折算的最终β（v/c），按时期、Si/S3及ring读取；运行时不再乘倍率，原角度保持 |
| [doppler_beta.csv](doppler_beta.csv) | 全部三期96项实际β、反冲角和粒子参考能量，供阅读核对 |
| [doppler_beta_provenance.json](doppler_beta_provenance.json) | β导出的原设置、来源校验和转换说明；不参与事件修正 |
| [manifest.json](manifest.json) | 数值版本、207份文件 SHA256、来源路径及字节数 |

`time_alignment_*.txt/json` 是同一时间包的模式、格式、校验和来源附件；完整时间函数只对原生时间应用一次。`time_alignment_inherited_metadata.json` 中候选阶段的文字保留原文，实际采用以全局入口为准。

`cali_si_r457_wave_mapping.txt` 和 `particle_source_ge_snapshot.txt` 是粒子源包保留的接线/Ge快照；现用 Ge 基础表由 `cali_ge.txt` 选择。不要把这些快照当作新增修正。

统一选择入口仍为 [pars/active.json](../../pars/active.json)，现用选择版本为`doppler_tabulated_beta_20261103`，能量版本仍为 `energy_parameters_20261102`，时间版本仍为 `joint_zero_20261028`，窗口仍为`narrow_review_20261103`。`pars/` 中既有路径和旧 `include/parameters/` 数值路径均通过相对链接访问这里。窗口JSON、CSV及附加c也已原样迁入本目录；旧pars窗口路径保留兼容链接，不另维护一套现用数值。

2026-11-03按用户决定，将已采用的倍率预先折算为各时期、各环的最终β。`doppler.json`的`beta_by_period`直接保存β与配套反冲方向，程序只读取，不再乘倍率；原倍率仅保存在provenance来源记录中。冻结导出由`doppler/export_beta.cpp`生成，存于`doppler/parameters/tabulated_beta_20261103.json`。物理修正与上次采用结果完全相同，全部Ge角度保持，这次表示方式调整无需重新build或整批重做Doppler。

正常build直接从C2读取原始事件，完成刻度、时间修正和配对后，分别一次应用窄/宽窗；不再提供刚才增加的reselect入口。恢复前的代码和已生成ROOT的来源记录保留在[三轮直接build验证任务](../../codex/analysis/build_direct_window_validation_20261103/)。本次窗口目录迁移保持数值不变。

源刻度与 A/B/q 的配套关系、历史采用范围见[粒子采用记录](../parameters/particle_energy_20261102/README.md)；Ge 来源见[原迁移记录](../parameters/ge_energy_20261102/README.md)。此前目录整理保持数值；其后的Doppler采用不改变build刻度、时间与窗口。

维护源和冻结候选仍在 `cali_ge/`、`cali_si/`、`s3_fb_cor/`、`s3_energy_attenuation/`、`ge_time_walk/`、`doppler/`。后续采用新参数时，先保存冻结版本，再成套更新这里的文件、`manifest.json`、`pars/manifest.json` 和全局清单校验值。时间工具的 `install` 是候选安装步骤，不自动启用或覆盖现用文件。

从项目根目录核验：

```bash
python3 pars/manage.py verify
build_all/build 0617 --inspect-parameters
```

[本次迁移与读取核验](../../codex/maintenance/parameter_flat_layout_20261103/report.md)
