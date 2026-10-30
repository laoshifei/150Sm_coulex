# S3 衰减与位置参数：joint_energy_20261030

本目录是同一套S3完整能量候选的衰减分解，涵盖全部173轮；不是独立选定的另一套参数。整理日期2026-10-30，数值继承已验证的0230候选，尚未切换正式参数。

- [attenuation.csv](attenuation.csv)：173×56=9,688行逐轮、逐通道`k_strip`。
- [ring_position_factors.csv](ring_position_factors.csv)：173×24=4,152行`q`，包括非生产高流强外环，明确标注`production`。
- [manifest.json](manifest.json)：与正背面目录共享的版本、主文件SHA256、来源及适用条件。
- [完整参数主文件](../../../s3_fb_cor/parameters/joint_energy_20261030/coefficients.json)：build实际只读取这个JSON，CSV仅供查看。
- [完整公式及使用说明](../../../s3_fb_cor/parameters/joint_energy_20261030/README.md)。

`B = p1_reference × k_strip`，`E = q[matched_ring] × (p0 + B × E0)`；同一q作用于ring与sector。E0只含基础刻度及硬件增益，不能重复使用旧衰减。p1/k分解依赖参考规范，不把k单独解释成已经测定的辐射损伤。

相对参考轮：低流强0616，高流强Ⅰ0472，高流强Ⅱ0680。高流强生产ring1–8、低流强ring1–24。44轮有非单位位置q；sector能量因此需要匹配ring，无法完整压成旧的单通道衰减文本表。

参数只有一份维护源：完整JSON。若将来重新求参，应发布新版本并同步再导出这两张表，不独立编辑CSV或混用历史正背面系数。现有逐轮质量与0634的折中情况见正背面目录的`run_status.csv`及完整说明。
