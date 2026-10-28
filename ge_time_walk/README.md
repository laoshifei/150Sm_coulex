# 当前采用：共同零点版本 joint_zero_20261028

2026-10-28已根据用户授权，在五轮整轮C++ build对照后切换。维护源为[parameters/joint_zero_20261028](parameters/joint_zero_20261028)，生产active及本目录默认选择均指向此版本。calibration.json与0219交付一致；完整Si/S3/Ge函数只加一次。

全173轮/346份窄宽ROOT已由原`build_all/build`完成重建、核验和替换；原346份旧ROOT及同一批数据的残留硬链接已按用户要求清理，旧回执、参数和核验记录保留。当前批次状态见`codex/analysis/build_joint_zero_adoption_20261028/`；中文窗口笔记为0220。窗口数值仍沿用0191，两端高度、随机归一和新谱检查单独列明；不声称所有CFD尾部或能区残差归零。

当前版本从项目根目录校验、复制：

```bash
python3 ge_time_walk/parameters/manage.py verify
python3 ge_time_walk/parameters/manage.py install
```

当前默认复制目标是`pars/time_alignment/versions/joint_zero_20261028/`；复制工具不修改`active.txt`。完整包含173轮的6920条Ge/CFD函数和26296组粒子通道函数，可读摘要分别为包内`ge_parameters.csv`和`particle_parameters.csv`。`inherited_metadata.json`保留候选阶段的原始状态，正式采用信息见同目录`metadata.json`。本地Notebook、研究记录、原始数据和能量刻度`pars/`不随本次Git提交；两个`candidate_*`目录仅为本地研究来源，运行当前包不需要它们。

以下保留0200历史说明，旧状态和默认选择只对应原记录日期；当前选择以`parameters/selected_profile.json`及`pars/time_alignment/active.txt`为准。

---

# Ge time-walk：完整参数、诊断谱与核心峰拟合

本目录同时维护C++程序和数值参数。已核验的0200完整包在[parameters/0200](parameters/0200)，覆盖173轮、20个Ge、normal/failed两类，共6920条完整函数；配套Si/S3函数保存在同一包，保持共同参考一致。`build`生成逐Ge、逐CFD、逐参考的诊断谱；生产符合事件由`build_all/build`生成。原多级ROOT宏已退役，旧实现保存在Git历史中。

## 参数维护与供build使用

- [calibration.json](parameters/0200/calibration.json)是完整科学数值源；[package.json](parameters/0200/package.json)保存六个原始文件的SHA和来源。此次归位未改变任何物理参数。
- [ge_parameters.csv](parameters/0200/ge_parameters.csv)逐轮列出Ge的常数C=F(334keV)、门外回退B、CFD和节点数，是可读导出；运行求值使用完整JSON。
- [residual_models.json](parameters/0200/residual_models.json)记录0200新增19分支及其模型/系数/支持域。这些项**已经合入完整F，不能再加一次**；未更新分支仍保留继承的完整函数。
- [selected_profile.json](parameters/selected_profile.json)选择本目录默认版本。`ge_time_walk/build`默认读取它；`--profile active`才读取生产全局active，`--profile NAME`读取已安装版本，`--profile-dir DIR`直接读取完整包。

复制到生产参数目录，从项目根目录执行：

```bash
python3 ge_time_walk/parameters/manage.py verify
python3 ge_time_walk/parameters/manage.py install
```

默认安装到`pars/time_alignment/versions/s3_common_reference_0200/`。复制逐文件核对SHA；已存在且完全相同可复用，不同内容拒绝覆盖。工具**不修改`active.txt`**。后续版本在本目录创建新版本、保留来源，更新选择清单，再复制到新runtime版本。不要只改runtime副本而遗漏维护源。

2026-10-26本次发布时，生产active仍为`ge_partial_improvements_20261026`，用户已选0200作为接下来全173轮重建的参数，整批替换尚未执行。历史`metadata.json`中的候选状态按原样保留，包的来源与采用阶段在`package.json`单独说明。Git包含源代码和参数，不包含编译后二进制、原始数据或本地能量刻度`pars/`。

## 阅读笔记与参数定义

本地中文阅读笔记：`jupyter/codex_0204_150sm_ge_time_walk_parameter_package.ipynb`。包含完整173轮参数、20个Ge的normal/failed三时期曲线、拟合项及可运行查询代码。按项目规则该Notebook保留本地，不随Git提交。已有0198记录C++方法；0200记录共同参考与选择依据；0202列offset；0203记录0617/0472/0680三整轮build对照及normal/failed/合计时间谱。

完整修正直接加到原生解码时间一次：`t_Ge_corrected = t_Ge_native + F(run,Ge,CFD,E)`。normal原生时间含CFD细分；failed沿用粗时间。定义`C=F(334keV)`、`W(E)=F(E)-C`，所以`F=C+W`。完整F由ln(E)上的节点插值与既有failed解析项组成；所有继承项已经包含，旧offset与新增残差不可再加。

0200共同参考为每时期S3环1–8等轮/等环平均`G(40MeV)=0`，Si相位已合入。334keV是Ge函数的拆分能量，不是参考探测器。三整轮对照支持S3内部对齐改善，但宽左尾仍存在，500–2000keV净谱没有普遍变干净；高流强Ge残差拟合主要使用逐轮约5%缓存，不能当作全部173轮Ge的独立全统计标定。

## 运行

从项目根目录执行：

```bash
make -C ge_time_walk
ge_time_walk/build 617 --project . --inspect-parameters
ge_time_walk/build 617 --project . \
  --input /mnt/e/Sm150/rootfile/data_C2_0617_0617.root \
  --output /tmp/walk_0617.root \
  --profile local --write-pairs
```

编译依赖C++17、ROOT（含Minuit2和nlohmann/json头）与OpenSSL开发库；原始树读取器位于本目录`raw.cpp/raw.h`，不依赖尚未提交的`get_ts_offset`改动。参数复制工具只需Python3.8+标准库；运行原始诊断还需要本项目本地`pars/`中的刻度、映射与衰减参数。

输出必须不存在。默认读取整轮；`--begin/--entries`只用于明确的局部核验。ROOT输出以`.partial`写入，关闭并检查后才发布正式文件，同时保存`.json`回执。失败部分保留，不能因有文件或旧进度文件就认定完成。默认不保存配对树，`--write-pairs`用于固定事件集合的重投影核验。

将输出ROOT的**绝对路径**写入JSON数组，例如`["/tmp/walk_0617.root"]`，然后：

```bash
ge_time_walk/fit --inputs inputs.json --output new_fit_directory --group low
```

`--group`取`low`、`high1`或`high2`。不同期、重复区间、不同参数版本或分箱不能混合。可用`--ge 304`检查单个Ge；默认全部20个、normal/failed都计算。多个同期期次共同拟合的是**现用逐轮完整函数之后的共同残差**，不自动证明各轮形状相同；逐轮检验时每个清单只放一轮。

批处理：`ge_time_walk/batch.sh RUN_LIST INPUT_DIR NEW_OUTPUT_DIR 2 PROFILE`。PROFILE默认`local`，整批固定同一参数目录；也可显式指定`active`或已安装版本名。名单每行一个轮号，可加`#`注释；并行数显式限制为1–8，默认2，持续补位，最终等待全部任务并输出退出码表。程序不会自动继续历史名单。

## 处理口径

- Ge为20个当前使用通道，40≤Eγ<2000keV，使用未做Doppler修正的实验室能量；Si为15<E<100MeV，S3为25<E<100MeV。pileup/outofr事件排除并计数，CFD normal/failed都保留。
- 复用`include/sm150_event_calibration.h`中的现行刻度，保留S3双面归一、逐轮衰减与环映射。时间函数复用`include/sm150_time_alignment.h`，粒子和Ge冻结为同一完整版本。
- 固定方向`dt = particle − Ge`。原生Ge时间保留当前normal CFD解码/整数化，failed使用粗时间。诊断中不加随机抖动，完整修正后保留double；生产绝对时间的整数化/抖动不由此文件替代。
- 配对使用原生粗时间差严格小于2000ns的全部唯一Ge–粒子raw-entry组合，没有最近时间选择。对同一组配对同时保存`native_ge`和`current`；前者**粒子已修正，Ge尚未加任何offset/time-walk**，不是两端均未修正。配对树另存`dt_native_both_ns`。
- 参考分别保存`si`、`s3_ring_all`、`s3_sector`、`s3_matched`。后者要求ring/sector修正时间差<80ns且能量差<2.5MeV、双向唯一匹配；低流强环1–24，高流强环1–8。双面匹配不参考Ge时间；两面不当作两个独立粒子相加。这些匹配门属于本批诊断设置。
- 时间谱默认±1000ns、4ns/bin；γ能量10keV/bin。两种CFD、两种时间坐标、各参考和前后半段均保留，包括空谱与越界计数。两半按Ge原始entry划分，半段边界配对不丢弃。
- 共享Ge或粒子的配对存在相关性；计数不是独立粒子数。宽的原生配对门独立于本轮Ge候选，但仍需在正式分析中核对门外损失。

## 拟合与参数

按保存的能量区间拟合局部`Gaussian核心+非负线性本底`，采用Poisson计数目标与bin积分，不用ProfileY全时间平均值。核心范围默认峰位±80ns，再放宽25%复核。失败、触界、弱峰、次峰歧义和窗口敏感性均保存原因；空谱不会写成零修正。共享命中使形式Poisson误差仅是条件诊断，不能直接作为最终标定不确定度。

曲线分别比较常数、1/√E、1/E、log(E)、log²(E)和四结点log分段线性。令`u=ln(E/334)`，各模型的系数基分别为`[1]`、`[1,√(334/E)−1]`、`[1,334/E−1]`、`[1,u]`、`[1,u,u²]`或各结点的插值权重；系数单位ns。

前后时间半段互相检验；两个方向的中心RMS均须改善至少`max(1ns,原RMS的10%)`，CV结果相差1ns内优先简单曲线。最终由合并谱求出的实际导出曲线也必须分别在两半达到同一改善要求：合并后可用能区可能改变，不能只验证两半各自的曲线。它是本批筛选规则与内部复核，不是独立盲验证。S3匹配参考有足够核心支持时优先使用，否则明确退到Si；另一参考的能量形状与零点偏差分别报告。模型只在测得支持域内拟合，域外新增残差夹到端点，不能当作域外也已测量；支持域内部缺点区间仍需用真实事件回放检查。

`fit_results.json`保存全部能区峰心、参数、失败原因、窗口检查、六类模型、双向检验和Si/S3对照；`branches.csv`便于阅读；`aggregate.root`保存实际合并谱。

若整轮数据有通过检查的分支，`candidate_profile/`输出与现有读取器兼容的完整JSON：

`F_new(run,Ge,CFD,E) = F_old(run,Ge,CFD,E) + r(E)`，故`dt_new = dt_old − r(E)`。

只更新本清单完整覆盖的轮次；其余Ge/轮次、粒子函数和既有failed模型继承。旧failed项在完整函数内部只计算一次，不重复相加。导出采用自适应log插值并由实际C++读取器密集网格复算，容差0.002ns。无可靠候选或只有局部输入时不生成完整候选目录。**候选不自动写入pars、不更改active、不重建生产ROOT**；仍需真实事件重投影、逐轮稳定性、参考、峰宽/左右尾和符合接受率核验。

本地详细核验证据位于`codex/maintenance/ge_time_walk_refactor_20261024/report.md`和`codex/maintenance/ge_time_walk_package_20261026/report.md`，研究记录不随Git提交。能量曲线只能平移相同能量的时间分布，不能把同能量的随机宽尾拟合消失。
