# I1-COLD-PHYSICS-V3.1 参数合同、电功归因与持久恢复

日期：2026-10-19。工程实施及限定数值恢复验收，不是器件物理标定。

## 1. 分项结论

| 分项 | 实际结果 | 不包含的结论 |
|---|---|---|
| WP1 参数合同、入口拒绝 | PASS；版本 I1-ION-SURFACE-1.0，显式来源/单位/域/不确定性与模式检查 | 没有新标定的离子或表面参数 |
| WP2 已有电功账本归因 | 完成 50 个共同接受区间离线审计；唯一根因 UNRESOLVED | 未修复或升级旧 19.24% 弱电功精度 |
| WP3 持久接受态恢复 | PASS，独立子进程、同构建串行 B0 短窗 | 非跨版本、分布式、任意源或完整电路恢复认证 |
| legacy 兼容 | 两个独立库的小型两步结果逐字节一致 | 未重跑旧科学算例 |
| 真实器件、全局端口、RF、自动冷热 | 保持未知/未验收/缺闭合 | 不由恢复成功升级 |

当前分支 `architecture-v2.1-ultrafast-rf`，HEAD `8ec25920adff63a0b1a875c620b9a32105458ad7`。
HEAD 不是执行身份。开始时 GENERIC_DEV diff、2965 个保护文件哈希、执行版源码、实际新库/驱动和结束时授权增量分别归档。
见 [来源索引](source_index.json)、[执行记录](build_test_and_execution_record.json)、[保护前](protection_before.json)、[保护后](protection_after.json)。

## 2. 参数合同与物理资格

| 参数 | 当前显式值 | 分类 | 来源资格 |
|---|---:|---|---|
| mu_p | 2e-4 m2/(V*s) | C | V3 合成输运测试配置，不是空气高场实测 |
| D_p | 5e-6 m2/s | C | 同一有效正离子测试配置，未认证 Einstein/高场关系 |
| gamma_i | 0.02，量纲 1 | C | 离子轰击有效产额测试参数，不是 gamma_critical |
| 实际离子组成、表面材料/状态、入射能量、器件产额及误差 | 未知 | D | 当前本地证据不足 |

本轮没有从本地文献索引找到能升级为 A 类器件适用参数的证据，也没有把 B 类部分匹配文献填成器件值。
没有新增外部经验数值；原 MorrowLowke、ReactionModel、ColdPhysics 物理公式及 SP3 吸收边界不变。
声明的 0--1e5 Td、290--310 K、1e5--1.1e5 Pa 仅为原合成配置的运行检查域，绝不是这么宽的物理标定有效域。

[参数包](ion_surface_parameter_contract.json)、[schema](ion_surface_parameter_schema.json)及[来源表](ion_surface_parameter_source_register.csv)
记录物种、气体、E/N/T/p、单位、常数近似、来源文件 SHA、来源等级、未知不确定性、表面状态和域外 REJECT。
`parameters.py` 做 schema 之外的语义校验：mu/D 物种/气体/有效域匹配，表面物种/气体匹配，值/单位/来源哈希正确，缺项和非有限值拒绝。
首版仅支持显式常参数；插值表并未实现，不自动外推。未知表面或入射能量没有填为零。

`REAL_DEVICE` 要求 A 类、显式器件证据及非 UNKNOWN 表面/入射条件；本包必然拒绝。
`CONDITIONAL_LITERATURE` 不接受 C 类参数；B 类必须给出部分匹配限制，仍非标定。
旧 V3 入口也拒绝外层 REAL_DEVICE 标签与内部合成合同不一致；LEGACY_ELECTRON_ONLY 默认物理行为不变。
成功读取用户证据也只表示 `USER_EVIDENCE_PROVIDED_NOT_INDEPENDENTLY_VALIDATED`，软件不能替代科学审查来源真实性。
当前环境没有 jsonschema，未安装依赖；执行了显式字段/常量/枚举及完整语义拒绝测试，**不声称获得通用 Draft 元模式验证**。

## 3. WP2：19.2% 残差的可归因部分

唯一输入是旧 V3 已归档的 `local_charge_and_work_budget.csv`、物种表和离子/壁面表；逐行核对 pre/post 身份。
只使用 0--1 ps 的 50 个完整共同区间；不使用后续不完整墙面片段，不运行旧轨迹或静态补算。
[逐区间审计](port_power_residual_attribution.csv)保留 H/Return 各自离子收集、二次发射、传导电流与电荷、反应和 floor/clip 数量。
粒子数修正不等于修正电功，缺少势加权配对时不换算成 J。

最差区间 0.06--0.08 ps：

| 已归档量 | 带符号数值 |
|---|---:|
| W_port | +6.5219135023835506e-19 J |
| Delta U_E | +3.2036083626184623e-26 J |
| W_cond | +5.4697532567591718e-19 J |
| 原生总势电功残差 | +1.0522327233286668e-19 J |
| 残差 / abs(W_cond) | 0.19237297807325857 |
| 500 V * Delta(q_B_active_V1 - q_B_stable_V2) | +1.0522238223686019e-19 J |

最后一项来自相邻端点两个已存电荷诊断的差，而非从电功残差反推。它与残差高度关联，支持优先审查有限总势求解/弱量消去误差，但不构成独立 r_phi/r_psi 收缩验证，不能宣布唯一根因已证明。

| 候选原因 | 本轮能确定什么 | 仍需什么 |
|---|---|---|
| 时间求积 | 固定 500 V 下宏/子步电压加权差为零；面功采用 pre 通量配端点平均势 | 不能据此证明全部分裂/输运时间误差为零 |
| 参考面不匹配 | 所用 Q 和 U 来自同 ActiveNetwork；未将 q_A 混入 | 不能由此认证全局电极面 |
| 空间算子/求解误差 | V1/V2 差分关联可直接计算；来源身份可核对 | 各区间独立 r_phi/r_psi、空间修正与势配对缺失，UNRESOLVED |
| 数值修正 | 实际 floor/clip、源数量与电极交换可读 | 旧表未存 W_outer/W_correction/W_reaction 及完整内部面功数组；留空而非置零 |
| 材料/发射能量 | 明确 UNAVAILABLE | 未证明是上述纯电网络恒等式残差的原因；禁止残差当热或发射能 |

本轮新恢复短窗保存更多功字段，只用于恢复一致性，不替旧轨迹补数据，也不用于宣称修复旧弱功误差。

## 4. 原生持久检查点实现

新增 [native_persistent_checkpoint.hpp](../../../adapters/i1_fp3_port_circuit/native_persistent_checkpoint.hpp)，
在原 [ColdFullProvider](../../../cpp/src/streamer/ColdFullProvider.cpp) 增加 `save_accepted_payload` / `load_accepted_payload`。
复用已有 `NativeFullSnapshot` 和事务，不复制 StreamerSolver，也不改旧 density/time checkpoint 的语义。

保存九场 ne/np/nn/rho/phi/Er/Ez/|E|/Sph、time、dt、上个接受时间、版本/步号、电子/离子面记录、H/Return 电荷/头部历史、两种电容缓存、权重场及 key、六个 SP3 warm 数组、求解诊断、累计 map、输出 offset。
参数、几何、网格、边界、源、源码/库/执行文件哈希由外层完整配置身份绑定；不同配置或构建不允许静默恢复。
累计 map 保存调用者提供的全部键，本测试实际含数值修正、电荷、电功、photo、离子收集和发射累计量。接受 trace 同时保存各步库存/电子收集/边界诊断。
这不是对旧 V3 丢失机制账本的恢复器；未提供的逐机制累计反应及初始化之前修正不能从检查点补造，首次接受之前的初始化修正也不重复累计。

PETSc 对象、指针、预条件器不序列化；同配置重建 Poisson/SP3 工作空间，保留必要 warm 和物理历史，后续实际两步结果检验其影响。
观察器指针不持久化，非空观察器必须明确重绑定，不能恢复悬空地址。
当前方案是固定电压 B0 同构建串行协议，不是任意电路状态、分布式提交或跨 ABI 格式。

外层 [checkpoints.py](../../../adapters/i1_cold_physics_v3_1/checkpoints.py) 使用 `.incomplete` 目录、fsync、SHA256 manifest、COMPLETE 校验及原子 rename。
只发布 ACCEPTED_CHECKPOINT；TRIAL_STATE/DIAGNOSTIC_SNAPSHOT、缺完成标志、版本/配置/电容归属不符、损坏 payload 或输出前缀均拒绝。
恢复到新输出目录，只复制已接受 offset 之前的内容，不重复追加旧试探尾部。
内部二进制另有校验、大小/形状/有限性和接受时间版本检查。低层格式仅同构建使用，不能脱离外层身份校验当通用 checkpoint。

## 5. 实际恢复及兼容验收

[事前计划](execution_plan.json)在运行前登记字段尺度、1e-10 场比较门、1e-8 固定小交换尺度账本门及严格时间/步号规则，没有事后放宽。
同一 B0 +500 V、SP3 ON、2.5 um、合成参数配置：

1. 第一个进程原生推进 1 步到 20 fs，发布接受 payload，继续 2 步到 60 fs。
2. 第一个进程退出后，第二个进程读取同一接受态，执行相同 2 步到 60 fs。
3. 合计 B0 `step()` 5 次。两个进程正常退出，计数和最终时刻均已保存。

[170 项比较](recovery_comparison.csv)全部 BITWISE_IDENTICAL，归一差最大 0；包括九场、面通量、SP3 warm、权重/历史、时间/步号、累计账本及功字段。
恢复输出没有重复接受行。PETSc 工作空间虽重建，其下游物理量在此测试中逐位一致，不将内部指针/耗时相同作为要求。
见 [恢复报告](recovery_reproducibility_report.md)及[机器摘要](recovery_summary.json)。

另用原小型 legacy probe 链接旧/新库，各执行 2 步，共 4 次，输出逐字节相同；未运行历史 CTest/长轨迹。
实际损坏 native payload 的新子进程在 load 阶段拒绝，0 次 step、无接受 trace。
最终隔离测试还覆盖参数/单位/域/身份拒绝、TOTAL/SPLIT 双计、版本/试探/损坏/不完整文件，以及真实子进程被终止时已发布 envelope 仍可校验。

**执行版本分层：** 原生恢复对照使用 `executed_sources` 中的 Python 驱动，当时在原生参考进程结束后封装 envelope。
随后将外层强化为子进程运行期间发布并处理超时，另用真实子进程复制既有 payload 后被终止的隔离测试验证；没有为 Python 强化再跑 B0。
原生保存/读取实现与已测试二进制未变。`final_sources` 与执行版差异显式保留，不把存储测试当第二次 PDE 验证。

## 6. 实际修改与复用

| 文件/模块 | 授权增量 |
|---|---|
| ColdFullProvider.hpp/.cpp | 新接受态持久保存/读取和 ledger 类型；旧 trial/commit 行为不变 |
| native_persistent_checkpoint.hpp | 同构建版本化完整 payload 序列化；PETSc 重建合同 |
| adapters/i1_cold_physics_v3_1/ | 参数语义校验、envelope、RecoverableI1ScenarioRunner、短窗驱动、测试及离线归档 |
| 原 adapters/i1_cold_physics_v3/scenario.py | 可选新合同接入；真实模式无合同/标签冲突拒绝；旧 legacy 保持 |
| CMakeLists.txt | 独立恢复驱动 target，新 build_i1_cold_physics_v3_1，不覆盖旧库 |
| 两份长期文档 | 追加本次证据，不改历史 Gate/Issue |

未修改 StreamerSolver、ColdPhysics、MorrowLowke、ReactionModel、Poisson、SP3、FP3 stable_V2 的物理或电功算法。
TOTAL/SPLIT 模式随 checkpoint 严格绑定且有拒绝测试；本次原生恢复实际测试 TOTAL，不冒称新跑了 SPLIT 电路等价测试。

## 7. 资源、入口和停止

可复现命令、环境、编译失败/修复日志均在 [运行记录](build_test_and_execution_record.json)。
本轮 B0 原生子进程约 50.94 s，整体原生恢复编排约 51.97 s；编译单列。
总原生调用 9 次：B0 5 次 + 小型 legacy 4 次；无新 1/10/250 ps 长轨迹，无参数扫描、热/TENG/RF/COMSOL。
额外隔离校验、离线归档为秒级，已登记执行费用远低于 1200 s；少量交互式小测试没有精确合计计时，未伪报精确总秒数。

`RecoverableI1ScenarioRunner.preflight/run_segment` 是统一入口的持久恢复扩展，当前驱动只支持上述有界恢复 fixture。
`run.py` 拒绝已有结果目录，不会自动复跑；新任务需要独立输出与授权预算，不能将此脚本当任意长时间运行器。

下一项科学依赖是有来源且覆盖实际高场/物种/材料/入射能量的离子与表面参数包。
弱电功精度另需独立线性残差和逐单元修正功配对记录；本次只登记缺项，不为补记录自动开启新轨迹。
全局电极端口、发射材料能、自动冷热、自持导电通道和 RF 均未升级。本轮完成并停止。
