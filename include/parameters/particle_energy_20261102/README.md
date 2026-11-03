# Si/S3新基础刻度与配套A/B/q

这套运行参数按用户要求整理于2026-11-02。**S3 A/B/q已建立在0242的新S3刻度上，不能再配原S3基础刻度使用。**当前选择以[pars/active.json](../../../pars/active.json)为准；完整来源、采用范围和未更新项见[adoption.json](adoption.json)。已有生产ROOT不会因参数文件改变而自动更新。

## 从哪里查看参数

|内容|文件|含义|
|---|---|---|
|高流强Ⅰ基础刻度|[cali_si_r460.txt](source_calibration/cali_si_r460.txt)|0460源测量，应用于0472–0605时期的严格轮次|
|低流强基础刻度|[cali_si_r457_low_current.txt](source_calibration/cali_si_r457_low_current.txt)|0457源测量，已按低流强接线交换映射，应用于0616–0674|
|高流强Ⅱ基础刻度|[cali_si_r457_high_current.txt](source_calibration/cali_si_r457_high_current.txt)|0457源测量，应用于0680–0718|
|完整Si/S3刻度表及取舍|[source_parameters.csv](source_parameters.csv)|456个“时期–通道”实例；新旧系数、源通道、实际通道、采用状态和未采用原因|
|实际读取的A/B/q|[coefficients.json](s3/coefficients.json)|完整173轮、配套刻度SHA、物理通道映射、逐环/逐像素q|
|A/B数值表|[ab_parameters.csv](s3/ab_parameters.csv)|173×56=9,688行；A单位keV，B无量纲|
|实际应用的q数值表|[q_parameters.csv](s3/q_parameters.csv)|31,784行；sector=0表示该ring全部sector共用，正整数表示具体像素|

三个时期的运行表来自两次源测量，不是三次独立源刻度。表中基础刻度最后一列保留历史源拟合χ信息，**不是本次峰位一次刻度的χ²**。

## 基础刻度采用范围

|时期|源轮|Si array采用新值|Si array保留原值|S3采用新值|
|---|---|---:|---:|---:|
|高流强Ⅰ|0460|69/96|27/96|56/56|
|低流强|0457|74/96|22/96|56/56|
|高流强Ⅱ|0457|74/96|22/96|56/56|

这里采用0242中`window_checked=true`的主候选。图色变绿只说明图示源点残差改善，并不自动替代窗口/峰位审核；复核或仅诊断的Si通道保留原系数，完整列表在`source_parameters.csv`，未用邻道或任意替代值补齐。

物理S3为24ring+32sector，三期都使用配套新基础刻度。生产build仍只保留高流强ring1–8、低流强ring1–24，全部sector。高流强A/B表中的外环项用于完整结构和诊断，不改变生产筛选范围。Ge基础刻度、接线倍率及逐轮Ge增益保持。

## A/B/q与新刻度的关系

|阶段|工作|与本版本关系|
|---|---|---|
|0242|根据已拟合α峰中心重新求一次刻度|本次Si/S3基础系数来源|
|0244|把S3配对事件重新表达在新刻度上，处理A/B与逐ring q|完整173轮的基础A/B/q；部分重新拟合，支持不足处保留旧修正的等价换算|
|0245|固定新源刻度和A/B，检查低流强后部ring的逐像素q|40轮ring15–24的逐位置候选；未通过项继承原值|
|0246|按可重复的高能下降沿复核0641/0659/0674的ring20–24|只替换新增数据确认的166项q，其余继承0245|

本版本的全部A/B/q数值与0246完整包逐项相同，所有S3基础系数也逐项相同。现在同时启用Si新系数，因此三份完整源文件SHA发生变化；这是Si行更新导致，**没有重新修改S3的A/B/q或S3基础系数**。build核对配套完整源文件SHA，混用版本会在读取事件前报错。

部分A/B继承项经过等价变换，不是把旧A/B数字直接套在新基础刻度上。若旧、新一次刻度的斜率比为k=c1_old/c1_new，则保持同一ADC的旧修正能量可写为：

\[
B_{eq}=B_{old}k,\qquad
A_{eq}=A_{old}+B_{old}g(c0_{old}-k c0_{new}).
\]

其中g是该通道原接线倍率。重新拟合或采用等价项的逐通道状态保留在完整JSON的`channels`中。参考轮0472/0616/0680的能标锚定沿用0244条件；这次不是新的绝对重离子能量标定。

## 实际计算顺序

基础能量为

\[
E_0=g(c_0+c_1ADC+c_2ADC^2),
\]

单位keV；新一次刻度c2=0，保留旧值的通道仍按原系数。Si array使用自己的基础刻度及既有接线倍率，不应用S3的A/B/q。

对已经配对的S3 ring–sector：

\[
E_r=q_{rs}(A_r+B_rE_{0r}),\qquad
E_s=q_{rs}(A_s+B_sE_{0s}).
\]

- A/B替代旧S3正背面/衰减处理链；旧链不再次叠加。
- 低流强包保存`q_pixel`，两面使用同一个匹配像素值。高流强使用匹配ring的`q_position`，对应表中sector=0。
- `q_pixel`是完整剩余倍率，替代`q_position`，不把二者相乘。低流强内环和未确认位置的矩阵值可以来自原逐ring值，存成矩阵不意味着每格都独立测定。
- 时间修正仍从原ADC按已冻结的时间标定关系求值；此次能量参数更新不改变已采用的时间函数。
- `build_s3`、单道谱等前置诊断阶段保留其能量阶段。没有ring–sector配对时不能直接应用像素q；完整修正结果以`build_all/build`输出为准。

## 使用、追溯与局限

`pars/active.json`选择`particle_energy_20261102`后，正常执行`build_all/build`即使用本版本。`pars/calibration/particle_energy_20261102`、`pars/s3_energy/particle_energy_20261102`和原三个`pars/cali_si_*.txt`为相对兼容链接，实体仅在本目录。显式旧候选应配自己的源刻度和参数选择，不能只换A/B/q文件。

```bash
python3 pars/manage.py path source_calibration
python3 pars/manage.py path s3_energy
python3 pars/manage.py verify
build_all/build 0617 --inspect-parameters
```

这些命令只检查参数，不读E盘事件。此次更新没有重建或覆盖生产ROOT；已有ROOT仍包含创建时的旧参数版本，新旧数据不得只按文件名混用。

源刻度约5–6MeV的检查不独立认证15–100MeV外推；0246只确认三轮的部分像素，低能拖尾、部分下降沿及其他未解决位置仍保留。这些限制没有因正式采用而消失。

详见[本次核验报告](../../../codex/maintenance/particle_energy_install_20261102/report.md)。科学来源：[0242](../../../jupyter/codex_0242_150sm_si_s3_peak_position_recalibration.ipynb)、[0244](../../../jupyter/codex_0244_150sm_s3_new_source_calibration_abq.ipynb)、[0245](../../../jupyter/codex_0245_150sm_low_beam_s3_pixel_energy_correction.ipynb)、[0246](../../../jupyter/codex_0246_150sm_s3_energy_landmarks_and_q.ipynb)。这些历史笔记的“未启用”描述对应其整理时刻，现用版本由统一入口及本次采用记录确定。
