# 共同时间零点正式采用版本

2026-10-28经五轮整轮C++ build对照后正式采用，随后完成全部173轮的窄/宽窗重建。完整函数与0219交付及生产运行副本逐字节相同：Ge与Si/S3配套，直接加到原生时间一次，不再叠加旧offset或新增残差。

- `calibration.json`：完整科学数值源，覆盖高流强I 101轮、低流强40轮、高流强II 32轮。
- `ge_parameters.csv`：173轮×20个Ge×normal/failed，共6920条函数在指定能量的可读摘要。
- `particle_parameters.csv`：173轮×152个Si/S3通道，共26296组完整修正的可读摘要。
- `package.json`：六个运行/来源文件的SHA256；`metadata.json`记录正式采用，`inherited_metadata.json`保留候选阶段的原始状态及限制，不把历史`production_active=false`解释为当前选择。

完整JSON的SHA256：`6ddcb304f1b861b5e73468effa1a9e6c9f4943e17eb7a759956d0448dd1b7618`。

在项目根目录运行`python3 ge_time_walk/parameters/manage.py verify`校验默认包；`install`将它复制到`pars/time_alignment/versions/joint_zero_20261028/`，不修改`active.txt`。运行此包不需要本地两个`candidate_*`目录；来源路径只用于追溯。CSV是查询摘要，时间求值使用完整JSON。

本地实际采用与全173轮重建记录位于`codex/analysis/build_joint_zero_adoption_20261028/`，窗口阅读版为0220/0221；这些资料及能量刻度`pars/`不随Git提交。宽左尾和局部残差仍保留，不把采用结论扩大为所有能区或CFD分支均已归零。原候选和0200保留其历史状态。
