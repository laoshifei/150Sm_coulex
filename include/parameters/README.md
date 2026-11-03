# 原参数路径的兼容入口

2026-11-03起，现用数值已直接放入 [include/parameter](../parameter/README.md)，其中也包含完整时间参数。此处保留旧路径链接及原采用记录，下面是原目录说明；日常查看请使用新入口。

按用户要求，Ge/Si/S3基础能量刻度、Ge逐轮增益及配套S3 A/B/q的数值文件集中在此目录。

- [particle_energy_20261102](particle_energy_20261102/README.md)：0242通过检查的基础刻度，配套0244→0245→0246的S3修正；含完整参数表、来源和保留旧值的通道。
- [ge_energy_20261102](ge_energy_20261102/README.md)：原Ge基础刻度与184份逐轮k原样迁入，数值未变；严格生产名单仍为173轮。
- 实际采用版本由[pars/active.json](../../pars/active.json)统一选择；`pars`中的路径通过相对链接指向这里，避免手工维护两套现用数值。
- 生成程序和冻结候选仍保留在`cali_ge/`、`cali_si/`、`s3_fb_cor/`、`s3_energy_attenuation/`。时间、符合窗和Doppler等其他参数继续按`pars`现有布局保存。

运行时加载参数文件，无需把逐轮数值硬编码进C++头文件。完整SHA和此次启用记录见版本目录及`pars/manifest.json`。
