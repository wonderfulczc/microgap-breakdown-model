# I1 当前状态与阶段路线图

建立日期：2026-10-13；最近审查：2026-10-15。本文维护任务边界和证据依赖，不授予新的计算权限。

## 来源与缺失输入

本次用户提及的附件《I1 当前状态与阶段路线图》未在工作区找到，已请求提供但尚未收到：`ROADMAP_ATTACHMENT=NOT_FOUND_REQUESTED_NOT_RECEIVED`。本文依据本轮用户任务、下列仓库报告和既有台账编写，不冒充附件内容。附件到达后的差异审查仍待完成。

- [通用几何、边界、端口合同](i1_geometry_boundary_port_contract.md)。状态定义与物理端口边界继续有效，不被本文替换。
- [长期问题台账](i1_open_issues_and_gates.md)。继续复用八个永久 Issue，不新增重复问题。
- [初始载流子位置不确定性与论文表述归档](i1_initial_carrier_position_uncertainty.md)。按2026-10-15用户附件保留理论备忘与论文草稿，P1–P8已对应实际报告；不提供新的实验概率或计算授权。
- [FP-2 单次加密恢复报告](../adapters/comsol_fp2/single_refinement_run/recovery_8p8/FP2_SINGLE_REFINEMENT_REPORT.md)。
- [V5 本轮状态充分性报告](../validation/i1_part_v/state_sufficiency/V5_STATE_SUFFICIENCY_REPORT.md)及[机器判定](../validation/i1_part_v/state_sufficiency/v5_state_screen_gate_summary.json)。
- [热阶段 Part IV-R](../thermal/iv_dual_state_refinement/IV_REFINEMENT_FINAL_REPORT.md)。

## 冷阶段当前证据

Part V-A 的解析、有限体积和原生预算证据仅支持其已测试范围内的 K1 数学/离散结构，不等于 K2 跨工况预测、状态充分性或物理端口验证。

FP-2 的 1.25/0.625 um 局部 10 ps 相邻网格比较已完成。COMMON_GAP 的 G、增长与总机制指标为 1.169067%、3.044452%、3.379323%；FULL_GAS 对应为 1.169248%、3.043824%、3.383173%。预注册相邻网格目标通过，但整体参考仍为 `CONDITIONAL`，不是连续误差认证。最细网格专属时间证据有限，W_rho 尚无专属动态准入。8.8 ps 恢复处的端子历史重置不作为本轮导数证据。

## Part V5 首轮筛选（历史保留）

本轮正式进入有界的状态充分性筛选，不再追加 FP-2 全窗口加密，也不训练 ROM。

主状态定义域为 FULL_GAS：C1=[G_R]，C2-Q=[G_R,W_rho]。COMMON_GAP 仅辅助诊断。几何、500 V、气体、零 Neumann 和 SP3 OFF 固定；初始宽度为 10.8/12/13.2 um，仅幅值匹配同一 G，电中性保证 W=0。完整原生反应、面输运和线性化 Poisson 迁移率项均计入 dotG；数值修正分开。

| 项目 | 本轮结果及限制 |
|---|---|
| 计算执行 | PASS；5 个形状/网格配置、9 个接受步，11 项只读测试通过 |
| 同网格匹配态 | G 匹配满足 1e-8，W=0；不同原生 dotG 可测 |
| 细档窄/宽 dotG 差 | 0.00584174785 S/s；主要为边界项与反应项部分抵消 |
| 保守数值敏感性 | 0.03730060408 S/s；预注册逐状态空间差之和占主导 |
| C1 / C2-Q 精确闭合筛选 | NUMERICALLY_UNRESOLVED；不是充分性 PASS，也不是整个 I1 失败 |
| 5% 瞬时近似目标 | NOT_FALSIFIED_IN_TESTED_PAIRS；原始半差/S_pair=0.28143%，不能证明通用近似可用 |
| 参数拟合/自由演化 | 本轮未执行，不以短步筛选替代 |

总电子数随形状及幅值改变，因此没有同时保持电子总数，也未隔离“仅位置”因素。若未来模型将形状显式列为持续输入，本次配对不再拥有相同全部输入。现有结果不证明 W_rho 无端口价值，不证明 C2-S/C3 足够，不自动增加第三状态。

## Part V5 零拟合局部预测（2026-10-14）

用户另行授权后完成一次 C1-L0 近似评价，见[报告](../validation/i1_part_v/local_prediction/V5_LOCAL_PREDICTION_REPORT.md)和[Gate JSON](../validation/i1_part_v/local_prediction/v5_local_prediction_summary.json)。不再针对旧微小配对差进一步加密，旧误差包络与结论原封不动保留。

在0.625 um标准参考首个floor后FULL_GAS pre态，完整A/T/B/M给出唯一固定率 `k_star=9.13463628041107e9 s^-1`。`dotG=A+T+B+M`是投影恒等式；`dotGhat=k_star*Ghat`是新闭合假设。后者有一个预计算有效率，但没有根据后续/测试曲线拟合的自由参数，不是无需参考信息的第一性定律。

基准时间外推达到门限后，完成两条粗网格窄/宽初态及一条宽初态细网格的完整0-10ps轨迹。细网格对象按事前规定的粗档增长预测误差最大规则选择，非无偏随机泛化评估。全部网格沿用相同连续初态幅值，不重新匹配G；预测只使用各条实际pre G0，不回灌未来值。

| 新轨迹 | G RMS误差 | 增长RMS误差 | G相对误差P95 |
|---|---:|---:|---:|
| 窄种子，1.25 um | 0.0259082% | 0.495502% | 0.0372963% |
| 宽种子，1.25 um | 0.0334317% | 0.640991% | 0.0462089% |
| 宽种子，0.625 um | 0.0858242% | 1.61356% | 0.166285% |

三条均达到本轮5%/10%/10%探索目标。宽种子原生G/增长的粗细敏感性为1.17186%/3.09045%，大于对应模型误差；整体判定为 **LOCAL_PREDICTION_SUPPORTED_CONDITIONALLY**。这只支持限定几何、500V、空气、SP3 OFF及所测初态下的局部输出预测，不能认证连续系统5%精度。

9.98ps的完整瞬时率仍比k_star*G高约3.91%，因此积分电导预测好不意味着恒定增长率精确成立，也不授权更长时间外推。W_rho只作诊断。该模型没有电子数、净电荷、能量或物理端口的低阶守恒方程。当前没有证据显示必须增加拟合参数才能满足这个限定窗口的目标，不自动训练ROM或增加状态。

## V5-G归档与V6单点电压迁移（2026-10-15）

[V5模型选择终审](../validation/i1_part_v/V5_MODEL_SELECTION_FINAL.md)冻结C1-L0为 `ENGINEERING_MODEL_ACCEPTED`，只适用于已有证据支持的限定500 V居中初态点。当前工程冻结值为 `k_star=9134636280.4 s^-1`；C1精确闭合 `NOT_PROVEN`、C2-Q `NOT_SELECTED`、C2-S `OBSERVATIONAL_CANDIDATE_ONLY`。V5-C1的±10 μm偏移增长误差13.072578%/14.077755%失败原样保留，记录为“已知适用范围限制，当前不再投入计算”，不是声称实验起始必然居中或偏移罕见。

用户另行前瞻授权后完成唯一一条居中600 V、1.25 μm、0–10 ps迁移测试，见[V6报告](../validation/i1_part_v6/centered_600v/V6_CENTERED_VOLTAGE_TRANSFER_REPORT.md)及[Gate](../validation/i1_part_v6/centered_600v/v6_gate_summary.json)。这项授权不受旧V5-C1“七组全过再自动运行Level2”的调度条件阻断，也不回改旧失败或未运行判定。使用相同连续中性Gaussian初态、解析cell-average和原生算法，不匹配G或电子数；预测只读取600 V首次floor后的G0，系数不变。

|本次FULL_GAS观测量|结果|原探索目标|
|---|---:|---:|
|G RMS预测误差|4.281991%|≤5%|
|增长RMS预测误差|45.914630%|≤10%|
|逐时相对G误差P95|6.730002%|≤10%|
|实际电导响应|17.919016%|独立报告|
|末值有符号误差|-7.084670%|独立报告|

唯一轨迹完成515个accepted steps；新增科学进程1348.822 s（动态1344.818 s、初始只读RHS重放4.005 s），未重跑500 V、未编译或改求解器。判定 **FIXED_K_VOLTAGE_TRANSFER_LIMITED**：总G和P95达标不能掩盖增长误差失败；保持500 V工程模型，不把600 V纳入固定k已接受适用域。初始完整原生k0_test约1.620898248e10 s^-1，仅作解释，未用于预测或重估系数。

数值资格仍为 `CONDITIONAL_SINGLE_GRID_NO_600V_ERROR_BOUND`。旧500 V相邻网格误差不自动成为600 V误差界；本轮只有一个新电压/初态/网格，不能宣称连续系统、实验或500–600 V整个区间精度。既有C1精确闭合/C2-Q/C2-S、COMSOL独立验证、全局端口、冷热、RLC/RF与K2状态均不升级。固定系数迁移已获得本次离散目标下的明确限制，不自动追加电压、空间状态、拟合、V7或Part VI。

## 位置不确定性归档与当前行动（2026-10-15）

[位置不确定性主文档](i1_initial_carrier_position_uncertainty.md)保留原附件及P1–P8来源映射。±10 μm改变的是中性载流子云的位置，非电极；实际初始质心约±9.958335 μm，W0均为0。旧V5-C1 `C1_INITIAL_STATE_LIMITED`及增长误差失败不变；当前行动分列为 **KNOWN_SCOPE_LIMITATION / NO_FURTHER_COMPUTE_AUTHORIZED**。这不表示实际偏移罕见、实验必然居中，也不以停止计算替代物理闭合。

当前没有可用本装置实验p0（`EXPERIMENTAL_P0=MISSING`）。条件确定性动力学加初态概率律是一种理论描述框架；附件概率公式、质心矩恒等式及随机模型示例均未作为新模型实施。V5-F没有证明ze不能作为状态或必须增加第三状态；C2-S仍为 `OBSERVATIONAL_CANDIDATE_ONLY`，不恢复开发。论文段落只作可引用草稿，外部文献只是方法背景，不替代本装置概率、击穿或实验验证。

V5本轮模型选择已结束，V6唯一600 V任务已结束且迁移受限；固定k_star仍为9134636280.4 s^-1，既非时间平均率也非普适气体常数。保持原500 V工程适用范围，不接受600 V固定k增长目标；保留10 ps、单网格和未提交GENERIC_DEV身份限制。其余V6几何/极性、V7、Part VI/FP-3、冷热/RF均未因此通过。原Gate不变，不因未通过的Gate自动创建新验证任务或执行400/800 V、几何、位置、随机初态测试。

## Part V7 共同适用域与既有冷输出对照（2026-10-15）

用户明确给出新V7协议后，仅对既有500 V居中五点与唯一600 V的FULL_GAS标量完成离线比较，见[V7报告](../validation/i1_part_v/v7_model_comparison/V7_COLD_BASELINE_COMPARISON_REPORT.md)、[Gate](../validation/i1_part_v/v7_model_comparison/v7_gate_summary.json)及[适用域矩阵](../validation/i1_part_v/v7_model_comparison/model_applicability.csv)。三种预测为G0、G0(1+k_star t)、G0 exp(k_star t)，统一使用实际首次floor后pre G0、51个0.2 ps节点，冻结k_star=9134636280.4 s^-1。未调用PDE/Poisson/COMSOL/StreamerSolver、编译、拟合或新增轨迹；偏移组只引用原报告，未重新处理。

500 V五点：常量保持增长误差均为100%；切线增长误差约3.227%–3.469%，C1约0.493%–0.641%，两者均达到原探索目标。C1相对常量增长误差减少约99.359–99.507个百分点，相对切线约2.587–2.974个百分点。常量只用G0，切线和指数都额外使用同一个全物理初始率信息，故相对常量的改善不能全部归因于指数形式。短10 ps中切线也有效，不能据此认定指数独有的物理优势或通用状态闭合。六条重算C1指标与旧归档一致，旧表和Gate不改。

600 V C1增长误差45.914630%，切线47.783675%，两者均未达到原10%目标；`FIXED_K_VOLTAGE_TRANSFER_LIMITED`和500 V既有工程范围不变。单网格与连续误差限制不升级，旧500 V敏感性不充当600 V误差界。没有以六个案例或相关时间点构造独立随机实验/显著性结论。

`V7_COLD_OUTPUT_BASELINES=COMPLETED_ON_AVAILABLE_DISCRETE_DATA`；物理常数R、R+V_f、RW、Toepler、Mayr、Cassie全部 `DEFERRED_NO_COMPATIBLE_PORT_REFERENCE`，原规划没有取消。G_R未验证为物理通道/端口电导，不使用V*G_R或V²*G_R生成电流/热输入，也不把RW/Toepler一概限定为LTE热弧。`V7_COMPLETE_PHYSICAL_BENCHMARK=DEFERRED_NOT_VALIDATED`，不能称完整V7科学验收通过。

[Part VI移交清单](../validation/i1_part_v/v7_model_comparison/PART_VI_HANDOFF.md)只将变量/参考面、局部外域、q_A/q_B/物理电荷、电容/位移电流防双计、未来反馈及后期通道依赖交网页端定型。`PART_VI_THEORY_READY=READY_FOR_WEB_THEORY_REVIEW_ONLY_NOT_IMPLEMENTATION`：无需为结束本任务而补跑全部经典模型，网页端可审查理论；Codex未开始Part VI/FP-3、电路、冷热或RF。本次按预定范围结束，不按模型排名追加候选或计算。

## Part VI 理论合同与源码映射归档（2026-10-15）

在用户本轮单独授权下，完成[统一端口、电路与能量理论合同 I1-PORT-1.0](i1_part_vi_port_circuit_energy_contract.md)与[源码映射/最小缺口](i1_part_vi_interface_mapping.md)，机器状态见[Part VI理论状态](../validation/i1_part_vi/part_vi_theory_status.json)。用户消息提供完整合同；独立同名附件文件未找到，分列 `FULL_CONTRACT_ATTACHMENT_MISSING` 与 `FULL_CONTRACT_TEXT_AVAILABLE_IN_USER_MESSAGE`，不冒称读取不存在的文件。

当前 `THEORY_DEFINED / MAPPING_DOCUMENTED`；实施 `IMPLEMENTATION_NOT_AUTHORIZED`，物理验收 `PHYSICAL_VALIDATION_PENDING`。只读核对C/G/F和FP-2观察器，未运行科学求解器、编译、拟合或任何测试。旧V7 Gate保持原归档值，不将新理论合同追溯写成V7已验证。

核心区别已映射：q_A是历史中心场诊断，q_B是局部同算子候选，物理Q_H仍待全局边界/参考面验收；原生i_disp_hv已包含空间电荷响应，不仅Cp*dV/dt。TOTAL_CHARGE与SPLIT_GEOMETRIC_CAPACITANCE在理论上互斥，运行时模式尚未实施。现有G2/G3的Cext并联拓扑不同于本合同用户串联Cext实例；保留旧实现和结果，不改名解释。旧C_A强制OpenCharge，C_B用当前配置边界，不能自动当成同一问题。

现有电压驱动FULL、原生吸收/面通量和逐物种修正可复用，但缺完整trial/commit/rollback、统一区间端口电荷/能量及能力拒绝；现有checkpoint加载会重置电极历史。C1只输出受限G_R，不能直接接成物理电阻、提供气体热或从一个G自动lifting到FULL。局部→全局、TENG源、冷热能量和Stage F/H源映射保留独立依赖。下一行动仅提交文档审查；**没有FP-3或数值实施许可**。

## 路线与准入

### 当前行动更新：FP-3 合成端口工程（2026-10-15）

用户新授权 I1-FP3-ENGINEERING-1.0 后，已在独立目录完成[合成端口工程报告](../adapters/i1_fp3_port_circuit/results_20261015_final/FP3_PORT_CIRCUIT_ENGINEERING_REPORT.md)、[工程Gate](../adapters/i1_fp3_port_circuit/results_20261015_final/fp3_engineering_gate_summary.json)及[FULL移交](../adapters/i1_fp3_port_circuit/FULL_PROVIDER_HANDOFF.md)。上轮 `IMPLEMENTATION_NOT_AUTHORIZED` 是历史状态，不追溯修改；本轮授权仅覆盖合成工程，不覆盖真实FULL连接。

独立实现不可变accepted状态、trial/commit/reject、TOTAL/SPLIT互斥归属及串联Cext隐式中点残差。每个耦合步实际经PortProvider试探，SPLIT由电路显式装配Cp电荷；矩阵指数仅作为独立连续线性ODE参考。最终25项隔离测试及20组固定合成矩阵通过；总接受步41464。最大无量纲代数残差2.36535e-14，累计能量残差1.29887e-14，TOTAL/SPLIT状态及累计量差0，最细主状态RMS相对误差0.0114583%，观测阶1.99479–2.00019。初次执行的浮点终点伪短步故障和失败测试记录保留；仅修复时钟调度，未放宽阈值或改物理参数。

结论 `FP3_SYNTHETIC_PORT_ENGINEERING_ACCEPTED`：合成骨架结案，不重复同一套R/C验收。C1能力/适用域拒绝已有独立wrapper，但不提供物理电流、热、反馈或lifting；FULL返回明确未实现，TENG源拒绝。`FULL_TRANSACTION_INTEGRITY=NOT_TESTED`，`REAL_PLASMA_CIRCUIT_COUPLING=NOT_IMPLEMENTED`。全局物理端口、B0、低场、冷热和RF不升级；旧G2/G3并联拓扑、Part VI理论状态JSON、V5/V6/V7历史Gate原样保留。下一步仅建议科学审查后另行授权真实FULL协议接入、同参考面状态回滚和区间预算；当前停止，不自动执行。

### 当前行动更新：原生FULL局部工程（2026-10-15）

用户另行授权I1-FP3-FULL-1.0后，完成[原生FULL报告](../adapters/i1_fp3_port_circuit/native_full_run_20261015/FP3_NATIVE_FULL_PROVIDER_REPORT.md)与[Gate](../adapters/i1_fp3_port_circuit/native_full_run_20261015/fp3_full_provider_gate_summary.json)。真实调用原有StreamerSolver::step，不复制求解器；只增friend状态访问及Poisson诊断快照接口，原离散/物性/时间推进不改。47次原生step含所有拒绝、受控失败及reference；两种归属分别从相同原始native快照执行3个20 fs宏步，实际根试探32次，场/历史和总账本归属差0。没有重跑旧合成矩阵或新10 ps轨迹。

保存/拒绝恢复精确；重复试探及observer ON/OFF一致。一接受步直接重放逐值相同；两步独立重算最大相对差1.187e-15，history相同，不伪称全部bitwise重放。原生串行事务Gate PASS，但不证明SP3 ON、进程持久事务或全局端口。

独立定义q_B_active_operator，消去原生Neumann替换行；q_A/q_B_native_links保留。C_local=1.725237679e-14 F，仅本局部数值网络。500 V恒压带权电荷残差/无抵消输运最大0.0165%；耦合接受步则31.13%，微弱载流子账本未解析。恒压电能残差/交换功0.5525%；耦合中残差/导电功最高10.43倍，不能用电容充电或预存场能掩盖。SPLIT微小剩余能量约-3.285e-20 J，原值保留，非负物理能量声明不成立；不截零或重定义C消除。

六宏步q根均达预定门，最大|residual|/(2e-16 C)=5.331e-6，Vp到500.034778 V；但根残差相对载流子输运及强L/dt²放大限制微小电流/功解释。因此决定 **FULL_PROVIDER_ACCEPTED_COUPLING_PENDING**：联通/事务确已实现，微弱电荷/电功精度准入待审，非“未执行”或模型失败。运行后7项只读/协议测试通过；封存时补全candidate累计字段后，最终8项API/只读检查通过，未增加PDE调用或改变已有accepted账本。已记录科学前端/内核合计约57.25 s，编译约9.73 s单列。初次stdout协议警告故障两进程均0次step，不重复已完成推进。

合成FP-3接受、V5/V6/V7和C1状态原样保留。全局物理端口、气体热、TENG/RF不升级；下一步最小待验仅是同网格端点Q/场能精度与外层残差放大的匹配，交科学审查，不自动修改求解容差、追加推进或扫描。

### 当前行动更新：FP-3 弱量稳定求值（2026-10-15）

用户另行授权 I1-FP3-SMALL-SIGNAL-1.0；前述“待审/不自动调容差”是历史停止记录，未回改。
本轮[报告](../adapters/i1_fp3_port_circuit/weak_port_run_20261015/FP3_WEAK_PORT_ACCURACY_REPORT.md)
及[Gate](../adapters/i1_fp3_port_circuit/weak_port_run_20261015/fp3_weak_port_gate_summary.json)
决定 **FULL_LOCAL_NUMERICAL_COUPLING_ACCEPTED**，仅该固定局部有效网络。
同源电荷/电功残差及旧外层 D*rq 放大得到独立收缩/真实端点支持；旧历史Gate不追溯升级。
新增默认OFF的 q_B_active_operator_stable_V2：源电荷投影与真实零电压响应链接交叉检查，
保留 q_A/native-links/V1；有限权重能量交叉项及负的remaining偏置不删除、不解释成物理热。
仍用原矩阵/原物性/ISG0/推进，只有运行时线性容差收紧、同网络稳定诊断与外层增量提交变化。
两恒压接受区间后，TOTAL/SPLIT各3个20fs宏步真实反馈，场/历史和总账本完全一致。
弱电荷余量/Qtr最大约1.59e-15，电能余量/导电功约3.50e-7；
外层保留增量|RV|约2.84e-14V，实际浮点提交端点重算约3.17e-12V；|Rq|/Qtr约6.08e-9、电路能量余量/导电功约2.04e-7，达到本轮预注册门。
24次真实step、14次显式静态RHS，科学执行约45.83s；16项隔离只读检查通过。
首次静态比较误用不同Cp的remaining，0次step即停；固定同一Cp后复用原参考和权重续行，
事件原记录保留，未重复参考求解或动态过程。实际执行源码与封存后mock默认值修订分别留存。
本限定数值接口工程结案，停止重复验收。全局物理端口、局部外域、气体热、低场、TENG/RF、
冷热及C1/V5/V6/V7状态均不升级；下一行动仅交科学审查，不自动计算。

### FP-4 冷—过渡—热交接接口审查（2026-10-15，仅文档）

[一次审查报告](../adapters/i1_fp4_handoff/FP4_COLD_THERMAL_HANDOFF_REVIEW.md)、
[字段能力矩阵](../adapters/i1_fp4_handoff/cold_to_thermal_field_map.csv)、
[schema](../adapters/i1_fp4_handoff/cold_thermal_handoff_schema.json)和
[分项就绪状态](../adapters/i1_fp4_handoff/fp4_handoff_readiness.json)已归档；未实施运行时交接。
FULL局部电荷、真实面电功、场能及数值修正可提供；W_cond不是W_gas_heat，rho_q不是rho_m，
i_port不是热柱I_col，U_remaining不是气体材料内能。最新FP3局部数值接口接受保留，
不重复合成矩阵、原生事务或弱量测试，也不重开C1/迁移/网格计算。

`AUTO_COLD_TO_THERMAL=MISSING_PHYSICAL_CLOSURE`：缺电功到电子/重粒子/激发化学储能及边界交换，
缺中性质量/动量与同EOS能量零点、二维到径向剖面/长度及净电荷/场能去向；G0连通性与QJ不足以生成热态。
不得从60fs记录、ne阈值或任意热化效率制造高温通道。
`GIVEN_THERMAL_INITIAL_STATE=READY_FOR_BOUNDED_IMPLEMENTATION`仅指独立来源的LTE参考、
明确T/p与u=0、长度/网格/EOS/边界/I_col及质量能量基准这一现有G1子集。
非零速度/守恒输入合同已定义，但接口未实现、数据待给定；真实B0热初态未提供。
独立给定热参考不证明自动击穿/交接，已有1mm人工参考不能改成70um并称标定。

当前G1只用T/p构造U=[rho,0,rho*e]，重置时间和累计量；公共U赋值不是已验初始化入口。
G2/G3真实存在，旧“未开始”文档按历史解释；G2给热端的是导电支路RMS，
其并联Cext拓扑不能替代FP3串联Cext。后续需全局时间锚与独立账本分段，不能丢弃速度/历史或重复计化学能。
下一次最小候选为独立GivenThermalInitialStateAdapter：输入有来源给定态和支路/边界/能量约定，
输出同EOS守恒态及库存审查或结构化拒绝；仅建议，需另行授权，不自动实施或运行。
热(H_c,M2)闭合、全局物理端口/外域、LTE、气体热、TENG、冷热/RF和实验Gate均不升级。

### FP-4 给定热态适配器实施与等价回放（2026-10-15）

用户另行授权I1-FP4-GIVEN-THERMAL-1.0后，完成
[给定热态报告](../adapters/i1_fp4_handoff/given_thermal_run_20261015/FP4_GIVEN_THERMAL_ADAPTER_REPORT.md)
及[分项Gate](../adapters/i1_fp4_handoff/given_thermal_run_20261015/fp4_given_thermal_adapter_gate_summary.json)，
决定 **GIVEN_THERMAL_ADAPTER_ENGINEERING_ACCEPTED**。上一轮只读报告/schema/readiness不追溯改写。
只新增独立GivenThermalInitialStateAdapter、runner和隔离测试，真实调用原G1构造与推进；
不复制PDE、不改物性/原G1/G2，不写U/time/累计量绕过构造器。
输入与EOS/源码/驱动来源绑定，严格原生网格中心T/p和显式u_r=0；
未知/非零初速、冷标量造热态、错误电流支路、历史续接、虚假B0或旧G2冒称串联均拒绝。

原150单元、0.5mm半径、1mm柱长、既有给定I_col、10ns主参考，两路径各103接受步、0拒绝，
合计206步；每步T/p/U/实际u/time/累计账本逐值一致，末值也复现匹配源码的旧主参考。
初态零速度不约束后续原生流体：最终径向速度169.93798m/s。
初始独立EOS质量/能量相对差约1.17e-16/2.60e-16；质量残差/M0最大1.17e-16，
能量残差/库存交换尺度3.14e-16、/真实QJ6.27e-10，Joule关系5.46e-16，均过本轮工程门。
回放函数计时约3.37s；30项隔离检查通过，额外非零人工时钟只做输入/函数测试，无第三条轨迹。
执行源码与封存后纯非法输入拒绝补强分别留存，未重跑PDE或冒称后修订版重新推进。

时间以tau和显式t_anchor分开；R_map=null+NOT_APPLICABLE_INDEPENDENT_REFERENCE，
新段累计归零不是继承冷历史。1mm人工热参考仍非70um真实B0、非自动热火花。
AUTO=MISSING_PHYSICAL_CLOSURE、非零速度/守恒初始化NOT_IMPLEMENTED、B0热态WAITING_DATA、
LTE实验NOT_VALIDATED、Hc/M2 NOT_RESOLVED及全局端口/TENG/RF/C1/V5/V6/V7状态均不升级。
仅交付[后续串联热接口清单](../adapters/i1_fp4_handoff/given_thermal_run_20261015/THERMAL_SERIES_COUPLING_HANDOFF.md)，
I_col/V_col/R_col、电容归属和完整热trial状态仍须另行审查；没有热RLC或其他后续运行许可。

### FP-4 给定LTE热柱单串联耦合（2026-10-15，本次另行授权）

[本轮报告](../adapters/i1_fp4_handoff/thermal_series_run_20261015/execution/FP4_THERMAL_SERIES_COUPLING_REPORT.md)
及[分项Gate](../adapters/i1_fp4_handoff/thermal_series_run_20261015/execution/thermal_series_gate_summary.json)
记录 `GIVEN_THERMAL_SERIES_COUPLING_ENGINEERING=PASS`，只限独立给定1mm LTE数值参考。
旧给定态适配器接受和FP3局部数值端口接受保持，不追溯修改旧交接清单或历史Gate。
新增独立热串联wrapper：原G1保持不变，复用G2无拓扑依赖的冻结R/柱电流RMS思路、
FP3串联参数与电容归属guard；101pF在串联回路，人工10V源不是TENG/800V实验。

A1/A2/B1/B2均完整10ns，200/400/200/200正式接受步；另5次被丢弃原生试探，合计1005次step。
TOTAL/SPLIT的每层完整热态/历史语义SHA、公共电路状态及账本逐值一致。
50/25ps对照峰值i、累计QJ、末R、轴温相对差分别0.14218%、0.08327%、0.02586%、0.01604%，
均过5%工程门；不声称连续误差界或独立物理验证。
独立区间能量残差按真实交换功归一最大1.44e-7（门1e-3），质量/M0最大3.50e-16。
源功、Rs耗散、柱体Joule、热边界通量及电感/两类电容储能分列；位移电流不加热。
热事务真实测试两种试探、原生step后受控失败、拒绝/重试/共同commit和stale/foreign保护。
本轮没有原生缩步，缩步重试实现不冒称已触发；串行内存事务不冒称持久恢复。
科学执行18.67s，14项隔离检查通过，旧源码/历史结果SHA及GENERIC_DEV diff不变。
首次关键字API调用失败在0个原生步时停止，原记录保留，修正后隔离执行，未增加第五条轨迹。

`GIVEN_THERMAL_ADAPTER_ENGINEERING_ACCEPTED`保持；
`AUTO_COLD_TO_THERMAL=MISSING_PHYSICAL_CLOSURE`保持。
真实B0热初态、LTE科学适用、近电极/轴向损失、全局端口、TENG、RF及热ROM均不升级。
本次结束，交科学审查；不自动追加热计算、经典模型拟合或用户实际电路预测。

### FP-4 轴向长度与阻尼敏感性（2026-10-15，单次另行授权）

[长度/阻尼报告](../adapters/i1_fp4_handoff/channel_length_run_20261015/FP4_CHANNEL_LENGTH_DAMPING_REVIEW.md)
及[真实指标](../adapters/i1_fp4_handoff/channel_length_run_20261015/thermal_damping_metrics.json)
独立用原150单元T/p和同D'Angola EOS复现G_length=4.14789481307e-6 S*m、1mm柱阻241.086152ohm。
相同径向物性下，70um假设柱阻16.876031ohm；100nH/串联101pF/Rs3ohm的
冻结二阶模型zeta由3.87859变0.315835；临界长度248.590294um，f0=50.079438MHz、
70um冻结fd=47.516077MHz。并联人工0.2pF另按三阶矩阵求极点，仅解析、不增热轨迹。
这些是长度敏感性，70um来自名义净气隙而非实测有效热柱长。

仅新增B70-P/H两条70um假设LTE参考，50/25ps cap均到100ns，2000/4000步，无拒绝；
材料库存与电阻均由原生length参数处理，新manifest明确HYPOTHETICAL_AXIAL_LENGTH_SCALING。
动态各只有一次电流反号（约10.911/10.914ns），没有完整重复周期；H负峰约-0.01793A可分辨。
热阻升高，Rtot约17.85ns穿越冻结临界值，100ns柱阻约785ohm；
不能以瞬时冻结判据声称时变系统已完成振铃。
50/25ps对照峰值、首零交、末R、QJ相对差0.04340%/0.03436%/0.06549%/0.009144%。
100ns FFT独立间隔10MHz，Hann/去均值/不补零；实际周期数不足，INSUFFICIENT_CYCLES。
不能把最大10MHz谱bin或初始47.5MHz复极点称实验频率；旧10ns窗口也不能认证50MHz。

解析长度阻尼敏感性PASS、限定热数值参考完成且一次反向波瓣已分辨；
EXPERIMENTAL_70UM_50MHZ_RF_VALIDATION=NOT_VALIDATED。
既有适配器、1mm四组和FP3 Gate保持，未重复这些测试；原生G1/G2/G3、物性及旧结果不改。
6000次原生接受步、科学进程约117.52s，11项隔离后处理/数学检查通过。
轴向均匀/50um径向尺度与70um长度同量级、端面热损失、近电极压降、LTE、TENG和寄生均待物理支持。
AUTO仍MISSING_PHYSICAL_CLOSURE，真实B0热态/全局端口/RF/热ROM均不升级。

### 源侧冻结位置V-Q-x接口（2026-10-15，单次另行授权）

[源侧报告](../adapters/i1_teng_source/results_20261015/TENG_SERIES_SOURCE_INTEGRATION_REPORT.md)
与[机器Gate](../adapters/i1_teng_source/results_20261015/teng_source_gate_summary.json)
记录SOURCE_INTERFACE_ENGINEERING=PASS；SOURCE_PHYSICAL_CALIBRATION与EXPERIMENTAL_50MHZ_VALIDATION
仍NOT_VALIDATED。真实净间隙70um不变，1mm热柱仍独立参考；旧70um热轨迹不重跑。
新增TENG_VQX_FROZEN_POSITION：VT=Voc-qT/CT、dqT/dt=i，复用FP3串联trial/commit协议，
区分源CT、串联Cext101pF、端口Cp和互斥TOTAL/SPLIT；未修改旧G2/G3。
缺源电荷/参数/单位拒绝，动态机械返回未实现，断开源的另一拓扑仅登记未选用。

四个纯电路SYNTHETIC_RATIO例（CT/Cext=0.1、1、10、理想源极限）用同一数学R_test，
无阻尼f0分别166.094705、70.823020、52.523757、50.079438MHz；非实测TENG频率。
完整源电荷矩阵的约束零根与真正动态模态分开；并联Cp事件按三独立模态解释。
唯一1V/Cp0.2pF/CT=Cext合成开关事件保持i、Vc、qT、CpVp连续，源/电感/物理电容能量不删除；
最初使用0.168760ps步长解析约3.375ps快过程。TOTAL/SPLIT公共状态与总账本差0。
矩阵指数参考的最细模态状态误差最大0.05214%、事件0.01438%，观测阶最小1.99926；
独立电荷/KVL/能量最大归一残差3.17e-14。22项新隔离检查通过；失败收集/逐位浮点断言记录保留，
未重跑科学曲线、未放宽事前门限。科学进程1.93s，所有热/PDE/COMSOL/RF调用为0。
固定源坐标能量qT²/(2CT)及Voc积分功不是完整TENG机械/静电能验证。
下一步真实计算需先给出实际闭合拓扑、同一机械位置和电荷基准下有来源的CT/Voc/qT0输入包，
不能用源比值或70um间隙猜测CT。真实热态、自动冷热、全局端口与接收链仍独立待验。

### FP-4 既有热阻与实验观测解释（2026-10-15，仅离线）

[离线归因报告](../adapters/i1_fp4_handoff/resistance_interpretation_run_20261015/FP4_THERMAL_RESISTANCE_AND_RINGDOWN_INTERPRETATION.md)
及[解释状态](../adapters/i1_fp4_handoff/resistance_interpretation_run_20261015/fp4_thermal_resistance_interpretation_status.json)
不是新模型Gate。读取6000个旧接受节点及初末150单元剖面；同源EOS只离线评价3次，无任何演化调用。
25ps组R从16.87603到785.00211ohm（46.5158倍）；全截面G_length为初值2.1498%，
初末EOS积分与ell/R一致。轴sigma为初值1.4853%，不能代替整个截面电导。
R_used/pre、R_end/post及有符号电流、区间RMS分开，原接受时序/历史误差不改写。

总Joule4.557763nJ，材料加宏观动能净增4.557763nJ，宏观动能增6.140176nJ，
EOS材料内能净减1.582413nJ；因此轴心降温不是全域能量流失。
r>50um初库存占96.193%（含热尾和大体积近环境），全域QJ/初库存0.031265%不能单独断言核心受热不足。
没有中间径向剖面或核心源/面通量预算，不能给膨胀/导热对R变化的贡献百分比，也不能唯一判实验根因。

Stage I实测合同仍MEASUREMENT_TEMPLATE_NO_DATA，已存示波器示例明确合成；
在审查的I1/H/I资料中未识别匹配70um装置的真实测点/原始波形/CT/Voc/qT0和高频开关连接。
直接回路多周期与仅接收端振铃必须分开；只有后者不能反演热阻，参考无多周期不排除实验接收振铃。
给定LTE+10V不符直接实测回路时也只能条件性不符，不能据不匹配激励/初态判G1 PDE数值错误。
最小下一输入先为原始波形测点/参考面/时钟负载及高频接线图；再依据真实源状态审查储能，
确有直接回路差异后才由科学审查选最小非LTE/近电极/有效通道改进，不自动拟合或扫描。
SOURCE_INTERFACE_ENGINEERING=PASS、所有FP3/FP4旧Gate及AUTO缺物理闭合保持；6项隔离离线检查通过。

### FP-4 早期驱动与自由松弛对照（2026-10-16，单次另行授权）

[对照报告](../adapters/i1_fp4_handoff/early_relaxation_run_20261016/FP4_EARLY_THERMAL_RELAXATION_CAUSAL_REVIEW.md)
及[决定](../adapters/i1_fp4_handoff/early_relaxation_run_20261016/fp4_early_relaxation_decision.json)
记录 `ARTIFICIAL_INITIAL_RELAXATION_DOMINANT_IN_THIS_REFERENCE`，不是新的实验或模型准入Gate。
两条同初态70um给定LTE参考各完成0–25ns、1000原生接受步，无拒绝；旧100ns结果不动。
10V串联驱动逐节点复现旧25ps路径；独立质量/动量/能量RHS逐单元差0。
固定r<=50um核心真实面通量预算合格，未改G1/EOS、提供器或时间积分。

5/10/20/25ns的驱动R为29.6590/42.2880/64.6677/76.0376ohm，
零Joule组为30.8305/43.7958/66.3230/77.8853ohm；相对驱动效应为
-3.7997%/-3.4427%/-2.4958%/-2.3723%，全窗最大幅度4.2589%。
自由松弛本身已产生主要上升；电驱动有可分辨的抑制作用，不是零反馈。
25ns驱动核心Joule4.55580nJ、导热外流50.4550nJ、流体数值能通量外流78.7468nJ，
核心总能下降124.646nJ；自由组核心无Joule而下降127.263nJ。
该面流体能量包含原生压力/平流与Rusanov数值通量，不将其直接命名为纯膨胀贡献。
核心逐步能量/质量交换归一残差最大约2.79e-13/1.12e-10；
全域弱边界质量及自由能量交换低于库存ulp，保留绝对与原始相对残差并标分辨率不足，不伪称通过。

用户本轮补充的“**发射端电压和电流也振荡，接收端周期更多**”保留为USER_REPORTED；
不能继续把实际陈述简化为仅接收端振铃。原始同参考面波形/拓扑/带宽及源CT/Voc/qT0仍未提供。
同EOS初中心质量密度约0.034625kg/m3，为环境的2.95855%；此人工高温/近2atm状态不是冷材料能映射。
下一步优先来源可信的预击穿/源储能、冷到热初态与测点数据，而非追加无来源的热轨迹。
FP3、FP4和TENG源接口原工程Gate保持；AUTO、真实B0/LTE、全局端口、TENG/RF与热ROM不升级。

### 动态火花电阻与RF输出优先级（2026-10-16，单次另行授权）

[文献与条件性敏感性报告](../validation/i1_part_vi/spark_rf_priority_20261016/I1_SPARK_RESISTANCE_LITERATURE_AND_RF_PRIORITY_REPORT.md)
及[分层状态](../validation/i1_part_vi/spark_rf_priority_20261016/rf_resistance_priority_gate.json)
结束人工10V/先验高温低密度LTE柱的进一步阻值扫描；旧G1/FP3/FP4工程与松弛对照不撤销。
10V可以是塌陷后柱压，缺的是实际预击穿储电/开关/热初态，不把旧功乘(800/10)^2。
论文主线优先可信回路波形与RF源量；完整气体热化学不是无条件前置门，但能量效率主张仍需独立材料/光/电磁账本。

六篇原出版来源按摘要/可访问正文段落分层核对；2025正文核对Toepler/RW公式与单位。
2010摘要与2018对既往研究的回顾有概括差异，未宣布普适优劣；无文献常数移植。
Mayr/Cassie保留适用阶段/参数缺口，经典模型物理比较仍延期，未借GR伪造电流。
新计算只有SYNTHETIC_EVENT：500/800V初始开路端电压、CT/Cext10、Cp0.2pF，串联101pF保持。
四条规定数学R曲线共八次事件尝试；S500/R20在旧1V内层尺度下求根中止且未重跑，
其余七条100ns在事前事件尺度下完成，不覆盖首失败或补造缺失曲线。
800V中R20与R200的首峰电流比4.20268、10–80MHz代理RMS比6.33897；
预定20→200ohm包络无反号，R20有四个同号完整周期，其他案例周期不足。
说明导通/存活/恢复时间尺度值得约束，不代表真实火花参数被标定。
恒阻主要状态expm误差最大1.1323e-5；变阻仅做局部一步/两半步抽检，未称全窗收敛认证。
TOTAL/SPLIT为实际同端点端口trial/commit分账重放差0，不新增第二套完整电路轨迹。

唯一RF代理为dt*abs(FFT(di_loop/dt))，不是辐射功率/接收电压。100ns独立bin间隔10MHz，
无补零；单调/单波瓣曲线即使最大bin落50/60/80MHz，也不叫识别出对应振荡。
现有Stage F空间J/电流矩与H传播/接收接口不被标量环流替代，旧频带信任不扩展。
源CT/Voc/qT0/Cp、HF开关、真实热态及发射/接收测点原始数据仍未知；
冷GR模型、AUTO、全局端口、TENG/RF实验全部原状态，不开展参数拟合或下一条轨迹。

### 受限形成阶段动态端口（2026-10-16，单次另行授权）

[Toepler/RW形成端口报告](../validation/i1_part_vi/spark_formation_20261016/SPARK_FORMATION_MODELS_ENGINEERING_REPORT.md)
及[四项正交判断](../validation/i1_part_vi/spark_formation_20261016/spark_formation_gate_summary.json)
记录形成律数学PASS、两个合成短事件耦合工程PASS，物理范围仅FORMATIVE_UNIPOLAR_ONLY，
DEVICE_PARAMETER_IDENTIFIABILITY仍NOT_CALIBRATED。不是本装置火花电阻或RF验收。
原七例100ns和S500/R20中止/PARTIAL不动、不补算；旧FP3/FP4/TENG及V5–V7 Gate保持。

新增独立端口复用旧FP3/TENG事务与串联电路；每个电压trial联立I_col和Q/S，
形成状态以C或A^2*s保存，不借x_V、不把总端口电流用于形成积分。
两个1V/5ns合成种子事件各1320接受步：R0均为100ohm，Toepler末67.485ohm、RW末83.345ohm。
K、名义70um、初始Q0/S0均为显式软件测试/重参数化，不是文献常数或实测历史。
独立连续ODE/DOP853与主隐式中点的最大相对差2.555e-5，独立端点账本最大归一残差1.56e-13；
73项隔离检查通过。实际拒绝/重复/篡改/stale/foreign及有历史共同commit已检查。
两种归属另以同端点端口trial/commit比较，无第三条电路事件；原生求解器/物性未改。

Toepler数学负流仍带符号，Q<=0拒绝；物理端口反向或形成时间过期即拒绝。
RW零流时S/R保持，固定K/ell无自然恢复；这两个形成基线不能替代完整双极通道/10–80MHz动态模型。
2026文章另找到作者稿v3，核对毫米级间隙及不同源/负载/探头拓扑，不称其为70um参数标定。
2010/2018概括差异和原文访问层级保留。无恢复项、下限、热化效率、Mayr/Cassie、拟合或RF谱图。

1. 当前停止点：形成律数学/合成工程交付完成，等待科学审查。真实K、R0/历史、形成时刻/支路、源/拓扑及恢复证据缺失时不开展器件预报或添加恢复律；不自动追加任何事件、热或RF计算。
2. 后续 V5 参数训练：需另行明确适用初态族、所需可分辨机制和数值参考。单恒压、单几何及本轮配对不建立可辨识性或自主闭合。当前 `ROM_TRAINING_READINESS=NOT_ESTABLISHED`。
3. V6/V7：V6居中600 V固定k迁移受限，其余电压/几何/极性未开展；V7冷输出对照完成不等于完整物理基准通过。经典模型仍保留规划并等待匹配参考，无自动推进或经验拟合。
4. Part VI/FP-3：合成工程及本轮限定局部弱量数值接口结案；全局物理端口仍未验证。不能将局部q_B、G_R或W_rho改名为全局端口量，不能重复计电容/吸收。原生事务及弱量证据仅在本限定SP3 OFF范围成立，不由旧G2或合成测试继承；后续物理阶段接口实施需另行授权。
5. 热阶段：Part IV-R 的设计矩阵可辨识性 PASS 和局部输出 CONDITIONAL_PASS 不抵消数值参考、覆盖、动态闭合及泛化的 NOT_RESOLVED；不冻结 H_c/M2，不继续热训练。
6. RLC/RF、冷热交接、COMSOL 最终独立验证和实验预测：本轮均不升级。NATIVE_RF_350MHZ 仍 NOT_RESOLVED。

### 实测CSV证据纳入（2026-10-16，单次离线授权）

[原始波形审查](../validation/stage_i/i1_waveform_review_20261016/I1_EXPERIMENTAL_WAVEFORM_REVIEW.md)
完成三份根目录CSV的全量质量/时间检查、原始尺度图、固定4/6/8sigma敏感性及分窗频谱。
资料可用性为 `EXPERIMENTAL_WAVEFORM_FEATURES_READY_ONLY`，不是物理模型验收。
CH1/CH3同设备/触发时间/逐点时间数组；CH2电流触发晚825.970939s，不能配成同步V/I。
CH1测点未确定；V/A为导出单位，探头校准、参考面、支路和deskew未齐。
CH1后100–500ns约22.5ns局部周期、Hann最大bin42.5MHz；电流后段最大bin约32.5MHz但周期离散。
初强后弱得到时域支持；不能据此确认接收端比同事件发射端多多少周期，也不认证精确50MHz。
CH3慢背景/起伏不自动等于气隙RF电压。不得直接V/I或形成电流积分。
真实火花参数保持 `NOT_CALIBRATED`；FP3/FP4/TENG/形成律工程Gate、旧七合成事件及中止记录全部保留。
下一步仅等待测点、探头和同步同事件资料；本轮没有科学求解器、参数辨识或新轨迹，完成即停止。

### K1/K2低阻物理需求与局部反应候选（2026-10-16，单次授权）

[条件物理桥接报告](../validation/i1_part_vi/conductive_channel_bridge_20261016/I1_CONDUCTIVE_CHANNEL_PHYSICAL_BRIDGE_REPORT.md)
完成低阻需求、当前原生物性/源项194项对照和恰好三条均匀物种积分。
三个计算Gate PASS不升级物理端口或RF。实际500V/1.25um Gaussian pre总电子数96.21741467；
在70um贯通均匀电子漂移柱假设下，20/60/200ohm需约3.7364e9–4.1939e10电子，
为种子的3.88e7–4.36e8倍。12/50um仅为已有剖面尺度转成的明示半径假设，非实测截面。
500/800V除净隙的恒场10ps数学检查分别增长12.0408%/52.9946%；零场1ns分支仅复合、
低场不连续限制保留。没有新空间FULL、热、电路或RF，也未读取实验CSV拟合。
下一步科学决定B：优先明确桥接空间态/自洽场到同面传导电流/电荷端口，并核对供能与高密度适用性。
本轮不实施后续连接；C1/V6/V7、FP3/FP4、TENG/形成律工程及AUTO缺口全部保持，完成即停止。

### K1/K2空间雪崩与冻结导电桥接（2026-10-16，单次授权）

[空间桥接报告](../validation/i1_part_vi/spatial_bridge_20261016/I1_SPATIAL_AVALANCHE_AND_CONDUCTION_BRIDGE_REPORT.md)
复用同一500V/1.25um/B0的pre0及post10ps，未初始化或推进FULL。
近轴70um冻结路径Lambda约2.64308、单次增殖14.0564、穿越247.615ps，不能解释为最终电子数。
恰好两个静态sigma网络及一个小均匀柱制造解完成；数学/数值检查PASS不升级物理端口。
G_R从1.11646e-10增至1.22356e-10S；条件G_frozen从1.26360e-12降至2.19487e-13S。
弱连续尾部有数值连接，但相对阈值改变会断开，未建立20–200ohm物理通道。
原G_R与半格接触电导网络不具严格算子上界同一性；仅同网络链接的psi试探能量变分界通过。
FP3 stable_V2的2.5um/状态/时段不匹配，保持其已接受局部能力但不拼接为本10ps物理端口。
唯一建议是有真实源/能量依据的空间通道生成/维持及同面电流-电荷连接，不接冻结G入动态电路。
历史Gate、全局端口、AUTO、实际TENG/RF全部保持；本次结束即停止，不授权后续计算。

### K1 B0光电离补给任务（2026-10-16，原生准入阻塞）

[B0/SP3报告](../validation/i1_part_vi/b0_photoionization_20261016/I1_B0_PHOTOIONIZATION_REPLENISHMENT_REPORT.md)
按本次授权核对当前源码，决定C：构造器在geometry_view+photoionization时明确抛出
`dual-geometry SP3 photon boundary not validated; use SP3 OFF`。SP3现有接口没有B0内部金属光子边界。
未删除guard、换针板或外部注入Sph；新初始化/step/SP3/Poisson/冻结电导/动态轨迹均0。
因此光电离活性、Return补给及250ps贯通均NOT_EVALUATED，不能写成光电离弱或无持续雪崩。
匹配声明配置的2.5um旧OFF10ps资料存在，但无ON配对，未生成差值；四个动态CSV仅表头，不填物理0。
原生电子输运/吸收存在、离子无空间输运、无二次发射；photo源对ne/np同加Sph的电荷结构存在。
仅对上轮1.25um OFF冻结路径做辅助B_ion积分，约13.4839、理想gamma_critical约0.0741624，
实际gamma和离子阴极通量仍未知，不能认证汤森自持或输入FULL。
唯一后续建议为科学审查B0内部光子边界合同及另行实施授权；不自行改核心、加阴极模型或扩大工况。
历史Gate及真实端口/AUTO/TENG/RF状态保持，本轮停止。

### I1-COLD-PHYSICS-V3（2026-10-16实施，2026-10-19中断恢复封存）

[实施报告](../validation/i1_cold_physics_v3/implementation_20261016/I1_COLD_PHYSICS_V3_IMPLEMENTATION_REPORT.md)
和[能力表](../validation/i1_cold_physics_v3/implementation_20261016/i1_cold_physics_v3_capabilities.json)
记录本轮显式授权的原生扩展：正离子保守漂移扩散、真实离子壁面驱动的阴极二次发射、
气体域六分量SP3内部金属吸收Robin边界、全启用物种的局部电荷/电功与串行事务、统一COLD_FULL_B0入口。
LEGACY_ELECTRON_ONLY默认行为保留；独立旧/新库小测试逐字节一致。36项原生、23项入口测试通过；
SP3六分量独立直接解最大相对差6.91e-12。离子/表面参数为SYNTHETIC_OR_CONDITIONAL_PHYSICS，未标定。

唯一新B0运行中断，完整共同账本及快照到1ps/50接受步；后续逐面片段独立保留，未充作完整结果，
没有从初态重启。B0_INTEGRATED_EXECUTION=CONDITIONAL_PASS，未完成10ps及正常退出记录。
最大物种源交换相对残差1.13e-13、稳定局部电荷/输运残差3.04e-15、SP3预算2.08e-9。
原总势电能残差/导电功最差0.1924，V3弱电功精度未认证；发射能量与气体材料能UNAVAILABLE。
上轮SP3 guard阻塞报告保留为历史事实，新版本仅在显式吸收边界能力下接通。
下一项架构依赖为有来源/有效域的离子输运与材料表面参数包及端口能力约束；
全局端口、自持通道、自动冷热、TENG/RF和全部历史Gate不升级。本轮封存停止，不追加轨迹。

### I1-COLD-PHYSICS-V3.1（2026-10-19，参数合同与持久恢复）

[V3.1报告](../validation/i1_cold_physics_v3_1/implementation_20261019/I1_COLD_PHYSICS_V3_1_PARAMETER_AND_RECOVERY_REPORT.md)
和[能力状态](../validation/i1_cold_physics_v3_1/implementation_20261019/i1_scenario_capability_status.json)
记录版本化离子/表面参数合同、统一入口拒绝与接受态持久恢复。当前mu_p/D_p/gamma_i仍全部为C类合成参数；
真实物种/高场范围/表面及入射能量为待补证据，REAL_DEVICE不自动采用合成值。
同构建串行B0在20fs接受态封存，连续两步与真实子进程重启后两步均到60fs，170项比较逐位一致；
共5次B0原生调用，另有旧/新库小型legacy各2步逐字节一致。只验证限定恢复，不延长旧1ps轨迹。
原V3的50个共同接受区间只读归因：最差0.06--0.08ps残差+1.05223e-19J，为导电功19.2373%。
与500V乘V1/V2电荷差分高度关联，但独立线性残差/逐单元修正功未归档，唯一根因仍UNRESOLVED，弱电功Gate不升级。
原生物理公式、SP3边界、历史库及Gate不改；原始执行版与后续Python发布/拒绝强化分别归档。
下一科学依赖是覆盖实际离子/表面及高场的来源参数包；全局端口、发射材料能、自动冷热和RF仍未认证。
本轮停止，不自动追加250ps或参数扫描。

### I1-COLD-PHYSICS-V3.2（2026-10-19，来源表与全物种空间源）

[V3.2报告](../validation/i1_cold_physics_v3_2/implementation_20261019/I1_COLD_PHYSICS_V3_2_IMPLEMENTATION_REPORT.md)
和[分项Gate](../validation/i1_cold_physics_v3_2/implementation_20261019/v3_2_gate_summary.json)
记录原生PositiveIonTransportProvider常数/EN表模式、来源/域拒绝、表面参数合同及Stage F实际消费。
未取得匹配高场原始物性表，本轮使用显式C类SYNTHETIC_TEST_TABLE；聚合np仍是EFFECTIVE_SPECIES_APPROXIMATION。
两个扩散方向可存储，但不等时原标量FV明确拒绝，不冒称实现完整各向异性扩散。
5次B0原生step完成0--60fs短窗及20fs接受态子进程恢复；场、面流、历史、累计量逐位一致。
真实电子净面流与正离子面流进入新Stage F读取器及既有轴对称电流矩诊断，不重复加入发射或位移电流。
新窗同网络稳定电功残差/有符号导电功最大约3.00e-11，原总势约0.1131；独立ubar*Delta r_phi解释其误差至弱功尺度5.58e-9。
该结果不追改旧V3 0--1ps的19.2%记录，不升级全局端口、材料发射能、RF传播或接收。
最小空间接口依赖仍是同参考面导体/引线回流、表面电荷和外域响应与气体J的一致源/端口包；
真实离子/材料来源缺口独立保留。本轮停止，不延长轨迹、不扫描参数。

### I1-COLD-PHYSICS-V3.3（2026-10-19，守恒全局源接口）

[V3.3报告](../validation/i1_cold_physics_v3_3/implementation_20261019/I1_V3_3_GLOBAL_SOURCE_AND_RF_INTERFACE_REPORT.md)
记录GlobalSourcePacket、I1GlobalSourceAssembler及缺全局资料时的LOCAL_GAS_ONLY/结构化拒绝。
只读V3.2 40--60fs：内部物种配对面严格抵消，全域扣除原生修正后的电荷残差相对面交换约2.30e-15；
逐单元修正/反应缺陷未归档，连续性资格保留CONDITIONAL_PASS，不由残差补造数组。
两类有限体积合成源实际进入旧Stage F Jefimenko传播和Stage H场投影；电/磁偶极独立检查最差约0.803%。
闭合回路总电流矩近零，空间布局由磁矩及方向场体现。没有接收电压或实验RF认证。
无新原生步、Poisson、热或电路演化，历史Gate不变。真实B0完整源仍NOT_READY：
唯一优先依赖为同参考面的真实电极/引线回流、支路及表面电荷空间映射；不猜几何，不自动追加计算。

### I1-COLD-PHYSICS-V3.4（2026-10-19，单元账本与几何/外域接口）

[V3.4实施报告](../validation/i1_cold_physics_v3_4/implementation_20261019/I1_COLD_PHYSICS_V3_4_IMPLEMENTATION_REPORT.md)
记录新增原生接受步机制数组、I1_NATIVE_CELL_CHARGE_LEDGER_V1、B0导体几何适配、外域Schur响应及GlobalSourcePacket 1.1。
唯一新原生窗为0--60fs、3次step；每步49340个气体单元，最大单元残差/无抵消交换约4.38e-14，
全域L1约1.08e-14，公共面严格抵消。旧V3.3逐单元CONDITIONAL_PASS不追改。
H/Return各297个局部金属面与旧轮廓来源核对；非零外域源制造例同时保持电势、端子电荷和能量。
显式有限面片、导体电流及回流进入既有Stage F源消费者；没有重复传播或新增RLC/热计算。
真实引线/回流、全局表面分布和B0外域算子仍缺失，FULL_GLOBAL_PHYSICAL_SOURCE继续拒绝。
下一唯一依赖是同H/Return参考面的真实外部几何与外域响应/表面映射输入；不自动追加计算。

### I1-COLD-PHYSICS-V3.5（2026-10-20，冷热独立模式统一入口）

[V3.5实施报告](../validation/i1_cold_physics_v3_5/implementation_20261020/I1_V3_5_MULTIPHASE_ARCHITECTURE_REPORT.md)
记录I1PhysicalAssemblyPackage、I1MultiphaseScenarioRunner及ScenarioStateEnvelope。
同一入口实际复用冷FULL三步至60fs，以及原1mm给定LTE热柱串联四步至0.2ns；没有复制或修改原生求解器。
热源仍为10V合成阶跃，0.2pF为合成端口电容；101pF仍串联，RMS柱电流与总回路/位移电流分列。
热预算残差/实际交换量最大约5.91e-9；冷单元电荷比约4.38e-14。独立段使用人工时间锚，不继承冷材料能。
冷源实际进入既有Stage F读取器/电流矩；热柱只输出条件性Jz和sigma，缺自由电荷/回流时不补造完整RF源。
装置输入、两模式和状态时间合同限定工程PASS；AUTO仍MISSING_PHYSICAL_CLOSURE，真实全局几何WAITING_DATA。
[中文测绘指南](I1_EXPERIMENTAL_ASSEMBLY_MEASUREMENT_GUIDE.md)列明实际引线、回流、测点、源初态及外域/表面映射需求。
下一唯一依赖是按同H/Return参考面取得来源绑定的实际装置输入包；不追加人工热轨迹或传播计算。
历史1mm参考、70um假设热态及实际70um名义净隙保持独立身份，旧Gate不追改。

### I1-P1（2026-10-20，装置来源及双源电路）

[P1实施报告](../validation/i1_p1/implementation_20261020/I1_P1_DEVICE_PARAMETERIZATION_AND_SOURCE_REPORT.md)
记录I1DeviceCasePackage、TENGContactSeparationSource及统一SOURCE_NETWORK_REFERENCE。
旧B0保持200um/30deg；R1=0.1406、DAG1=25.3646仅形成mm/全角假设下的140.6um候选，
D2和原单位未核实，没有新几何/PDE。工作区缺参数评估书，用户转录与合成介电参数严格分离。
几何CT/Voc真正接入原TENGFrozenPositionSource和串联电路；放大器另有Cext并联泄放及显式输出阻抗身份。
两条10ns参考完成，33项隔离检查通过；独立矩阵指数最大固定尺度状态差1.93e-5，
离线电能残差/实际交换尺度小于9.84e-13。无冷/热/传播新轨迹。
合成xT=0.1/1mm的CT约100/13.6pF，Cp=0无阻尼参考约71.0/145.4MHz，不是器件频率认证。
真实QT0、vC0、xT、表面电荷密度、胶层与寄生/回流仍缺，TENG初态SOURCE_INITIAL_STATE_UNRESOLVED。
冷/热入口只做共享来源兼容检查，未接通未知TENG物理耦合；AUTO及全局/RF资格保持。
下一唯一依赖是同参考面、同事件的真实高频拓扑与预击穿源/电容初态包；不追加扫描。

### I1-E2E-V1（2026-10-20，条件性端到端工具与轴向ROM）

[实施报告](../validation/i1_e2e/v1_20261020/I1_E2E_V1_IMPLEMENTATION_REPORT.md)
记录I1PredictivePipeline实际串接P1 TENG、原串联电路、守恒有限源、Stage F/H与负载接收电压。
一条100ns合成端口参考完成6个完整同极性周期；初始1ps RC快模态被20ps步长欠解析，
早期vg尺度误差约0.818，不能把能量闭合或晚窗口场误差0.051%升级为启动精度通过。
新一维保守电子/正离子、反应和Poisson ROM通过真实电荷端口耦合，完成两条Rs=3/300ohm、
70um参考及一条100um几何参考，各200ps；不是给定R(t)回放。
R变化引起弱内部物种变化；该低电压种子仍高阻，未形成实际低阻火花。
ROM短窗不足以分析10–80MHz；等效有限源/天线及所有TENG参数仍为条件性参考。
旧FULL/G1/传播内核不变，无新冷/热长轨迹。下一理论重点是FULL→ROM同参考面的
物种/面通量/QH/Wcond保守投影，G1材料能映射保持独立未闭合；不自动追加计算。

## 责任与开放问题

| Issue | 负责阶段 | 保留状态/边界 |
|---|---|---|
| I1-GEO-01 | B/C/E，Part V 审查 | BLOCKER，面向未限定几何推广；不要求重建已核验局部轮廓 |
| I1-BC-02 | B/C、Part VI | BLOCKER，全局边界/电容/物理端口 |
| I1-PORT-03 | Part V/VI、C | BLOCKER；q_A/q_B/物理端口不混用 |
| I1-NUM-04 | C、Part V | BLOCKER，未限定的参数训练；本轮筛选授权和执行证据单列 |
| I1-MAT-05 | C、Part V | DOMAIN_RESTRICTED；本轮非零电子支撑未触及低场奇异区不等于物性验证 |
| I1-THERM-06 | G1、Part VI | NOT_RESOLVED |
| I1-CIR-07 | Part VI、G/H | NOT_RESOLVED |
| I1-RF-08 | F/H、G、Part VI | DOMAIN_RESTRICTED |

上述状态保持各自作用域，不能使无关工作自动停滞，也不能凭执行成功关闭物理 Gate。旧报告、旧数值和未提交 GENERIC_DEV 保留；本轮不提交 Git。
