# S3 完整能量参数：joint_energy_20261030

这是一套已整理、可由现有 `build_all/build` 显式读取的完整替换候选，**尚未改变默认正式参数或重建现存 ROOT**。主文件为本目录的 [coefficients.json](coefficients.json)。它与0230定向复核候选、0231四轮全量试建使用的文件逐字节一致；2026-10-30只是本包的整理日期，没有重新求参。

覆盖173轮：低流强40、高流强Ⅰ101、高流强Ⅱ32；每轮24个ring和32个sector。高流强生产只使用ring1–8，外环9–24保留为诊断参数；低流强使用ring1–24。逐轮相对参考为0616、0472、0680，分别对应低、高Ⅰ、高Ⅱ。

## 一次完整应用

所有能量、截距均以keV计，倍率无量纲：

```text
E0 = hardware_gain * (c0 + c1*ADC + c2*ADC^2)
A = p0
B = p1_reference * k_strip
Er = q[run, matched_ring] * (A_ring   + B_ring   * E0_ring)
Es = q[run, matched_ring] * (A_sector + B_sector * E0_sector)
```

E0只含基础刻度和硬件增益，不含旧衰减或旧正背面项。完整JSON代替旧的衰减/正背面链，不能再叠加旧参数。q来自这次配对的ring，正背面各施加一次；sector不能脱离匹配ring使用一个固定q。

`p1_reference`和`k_strip`是在固定参考尺度下对总斜率B的分解，不是两种被独立测出的物理效应。p0按轮次变化；p1固定参考规范，k按轮次及通道给出。44轮含非单位q；整套修正130轮与旧正式值不同，43轮保持原值。这些是数值/覆盖统计，不表示每个更新通道均独立验证改善。

## 文件

- [coefficients.json](coefficients.json)：唯一供build读取的完整参数主文件，保留适用时期、参考轮和逐项来源。
- [channel_parameters.csv](channel_parameters.csv)：全部9,688行A/B及p0/p1/k分解，ring行列出q，sector标明按匹配ring查询；用于人工检查。
- [normalization.csv](normalization.csv)：全部9,688行正背面p0/p1分解。
- [run_status.csv](run_status.csv)：全部173轮来源、与旧值的差异、全量试建覆盖及局限。
- [衰减配套目录](../../../s3_energy_attenuation/parameters/joint_energy_20261030/)：k_strip和4,152项q表。这两张表是同一主文件的精确导出，不单独作为另一套参数维护。
- [manifest.json](manifest.json)：版本、SHA256、核验依据、适用条件。CSV采用完整浮点精度；build不读取CSV。

## 使用

在项目根目录执行，输出到新的目录；以下命令只是使用示例，本次整理没有执行事件build：

```bash
./build_all/build 0634 \
  --s3-energy-file s3_fb_cor/parameters/joint_energy_20261030/coefficients.json \
  --outdir rootfile/s3_joint_energy_20261030
```

只核对参数，不读取事件或写ROOT：

```bash
./build_all/build 0634 \
  --s3-energy-file s3_fb_cor/parameters/joint_energy_20261030/coefficients.json \
  --inspect-parameters
```

当前只支持直接读取原始C2的`windows`模式。默认时间参数和窄/宽符合窗仍由原build读取；本包只改变S3能量链。省略`--s3-energy-file`即使用现存正式旧能量参数；已有输出会拒绝覆盖。

## 验证及限制

0230完成173轮参数选择，其中31轮做了定向追加确认；其余142轮继承0229，并非都在0230重新测定。0231对0630、0634、0553、0718做完整原始事件试建，验证了能量应用、时间不变、Ge/Si不变及ROOT读回。其余169轮未在0231做全量试建。

0634是联合折中候选：共同事件统一旧尺度后的总体半宽略降，但ring11、sector30的半宽约增14%。保留原联合解，不能随意拼接几个旧系数。高流强谱形变化、相对参考尺度、缺少独立绝对能标和完整参数置信区间等局限仍在。

参数SHA256：`2ec47907e9281fabfe39468f644be3f4111f5ad7a5e4d72c60e3043db3e25f0b`。

原研究与本批核验见[整理记录](../../../codex/analysis/s3_replacement_parameter_package_20261030/report.md)、[0230来源](../../../codex/analysis/s3_calibration_targeted_review_20261029/report.md)、[0231全量试建](../../../codex/analysis/s3_build_candidate_pilot_20261030/report.md)。
