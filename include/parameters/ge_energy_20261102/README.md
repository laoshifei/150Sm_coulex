# Ge 基础刻度与逐轮增益

本目录保存正式使用的 Ge 能量参数。2026-11-02 按用户要求从 `pars` 原版本导入，**全部文件逐字节保持；没有重新拟合或采用新候选**。统一选择仍为 [pars/active.json](../../../pars/active.json)。

| 内容 | 文件 | 范围 |
| --- | --- | --- |
| 基础刻度 | [cali_ge.txt](cali_ge.txt) | 原表64行，包含37路物理Ge及未使用电子学地址 |
| 逐轮逐Ge增益 | [run_k/](run_k/) | 原有184轮完整表保留；build严格名单仍为173轮 |
| 文件来源与SHA256 | [manifest.json](manifest.json) | 1份基础表及184份逐轮表 |
| 迁移范围和采用记录 | [migration.json](migration.json) | 原路径、公式、严格轮次和已知局限 |

实际计算公式为 `E_keV = k_run,Ge × (p0 + p1×ADC + p2×ADC²)`。p0、p1、p2的单位分别为keV、keV/ADC、keV/ADC²；k无量纲，应用一次。原逐轮k由511keV峰位确定，基础表和逐轮k必须配套维护。

兼容入口 `pars/cali_ge.txt` 和 `pars/run_k/` 指向本目录。`pars/calibration/ge_energy_20261102` 是相对目录链接。C++仍在运行时读取文本，不把逐轮数值硬编码进头文件。

```bash
python3 pars/manage.py path ge_source_calibration
python3 pars/manage.py path ge_gain
python3 pars/manage.py verify
build_all/build 0617 --inspect-parameters
```

从项目根目录运行这些命令；只检查参数，不读取实验ROOT或重新拟合。

`pars/calibration/source_baseline_20261101/cali_ge.txt`、`pars/calibration/gain_baseline_20261101/`以及粒子参数包中的Ge副本保留为冻结来源。日常更新以本版本对应的正式选择为准；将来生成新版本后同步统一选择、清单与兼容链接，不手动改动两套数值。

本次只是改变正式文件位置，不改变Si/S3、时间修正、符合窗、效率或Doppler，也不需要为这次迁移重建已有ROOT。此前采用新Si/S3参数是否需要重建，是另一个任务。

既有Ge全能区残差、逐轮跳变及停用通道问题仍见[0247核查笔记](../../../jupyter/codex_0247_150sm_ge_calibration_and_run_gain_audit.ipynb)；迁移不表示这些科学问题已解决。[迁移核验](../../../codex/maintenance/ge_energy_include_20261102/completion.json)记录逐文件一致性和三个时期build的实际读取比较。
