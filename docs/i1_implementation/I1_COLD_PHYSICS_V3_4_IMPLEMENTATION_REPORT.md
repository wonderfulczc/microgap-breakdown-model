# I1-COLD-PHYSICS-V3.4 实施报告

日期：2026-10-19。性质：原生接受记录扩展、几何/外域/空间源接口实施及有限工程检查。
实际分支为 architecture-v2.1-ultrafast-rf；HEAD仍为8ec25920adff63a0b1a875c620b9a32105458ad7。
HEAD不是执行身份，GENERIC_DEV基线、授权增量、执行源码、独立库/驱动及最终源码分别封存。

## 1. 分项结论

| 项目 | 结果 | 限定 |
|---|---|---|
| NATIVE_CELL_CHARGE_LEDGER_IMPLEMENTED | PASS | 真实RHS、物种端点及修正数组，接受后导出 |
| NATIVE_CELL_CHARGE_CONTINUITY | PASS | 本轮0--60 fs三步，不追改V3.3历史CONDITIONAL_PASS |
| B0_CONDUCTOR_GEOMETRY_ADAPTER | PASS | 既有参数化局部金属面、轮廓来源及固定H/Return |
| EXTERNAL_BOUNDARY_OPERATOR_INTERFACE | PASS | 来源/单位/分块/端子泛函及非零外域源；真实B0算子不可用 |
| SYNTHETIC_SCHUR_CHARGE_ENERGY_EQUIVALENCE | PASS | 小型离散制造网络，不是COMSOL结果 |
| GLOBAL_SOURCE_SURFACE_CHARGE_INTERFACE | PASS_ENGINEERING_LIMITED | 显式有限面片到FV节点、导体电流分布及接触交换 |
| I1_SCENARIO_ARCHITECTURE_INTEGRATION | PASS | 实际原生运行、局部源消费及合成完整源消费 |
| REAL_DEVICE_LEAD_RETURN_GEOMETRY | UNKNOWN_EXPERIMENTAL_GEOMETRY | 未提供实际引线、回流、接收结构 |
| REAL_B0_GLOBAL_PORT | NOT_READY | GLOBAL_BOUNDARY_OPERATOR_UNAVAILABLE，缺全局表面分布 |
| RF_RECEIVER_VALIDATION | NOT_VALIDATED | 本轮无传播或接收计算 |

历史物理参数、冷热交接、C1/V5/V6/V7、FP3/FP4、TENG和Toepler/RW状态不升级。

## 2. 真正修改与调用

- `StreamerSolver.hpp/.cpp`：扩展AcceptedStepRecord，只复制实际求值后的ReactionSources：电子/正离子/负离子RHS、离子复合及净电荷缺陷；保存pre-floor之前的原始三物种密度。不改反应、输运、Poisson、ISG0或SP3公式。
- `adapters/i1_cold_physics_v3_1/recovery_native.cpp`：新增默认关闭的I1_V34_HOOKS。新驱动包含原驱动，继续使用ColdFullProvider evaluate_trial/commit；只有commit之后才导出。
- `adapters/i1_cold_physics_v3_4/cell_ledger_export.hpp`：原生Geometry2DView金属面、mask、实际轮廓与每接受区间完整单元/面导出；没有第二套求解器。
- `python/streamer_rf/rf/source/cell_charge.py`：版本化读取、来源/时间/完整性检查，独立物种与净电荷账本，明确稀疏诊断选择。
- 同目录`conductor_geometry.py`、`external_boundary.py`、`surface_exchange.py`：实际局部几何适配、Schur响应及显式表面/接触/导线合同。
- `global_packet.py`：加法扩展1.1，保留1.0入口；局部源可要求完整单元账本，完整源检查新表面合同。Stage F/H传播内核未修改。
- `adapters/i1_cold_physics_v3_4/scenario.py`：继承既有SourceBoundI1ScenarioRunner，增加基于数据检查的五种源能力路径；FULL_GLOBAL_PHYSICAL_SOURCE仍拒绝。
- `CMakeLists.txt`：独立i1_cold_v34_smoke目标，使用build_i1_cold_physics_v3_4，不覆盖旧库。

实际链：原生三步 -> ColdFullProvider共同commit -> 单元/面记录 -> Stage F accepted_species读取器 ->
GlobalSourceAssembler完整单元检查 -> 原有轴对称旋转及current_moment。
合成有限面片/导线/回流另进入同一装配器与SourceSeries；未调用传播内核。

## 3. 原生逐单元账本

本轮新运行唯一一条+500 V、原2.5 um B0工程配置、SP3 ON轨迹。
离子E/N表、gamma=.02及激发身份严格沿用V3.2的C类合成配置，不是器件标定。
3次实际StreamerSolver::step，3次接受，节点0、20、40、60 fs，正常退出。
原生驱动在20 fs生成接受payload，但本轮没有另开恢复分支或重复恢复测试。

每步49,340个气体单元、197,332条有向面记录。全部气体单元覆盖，零轴面积未产生额外通量。
原始密度、体积、实际反应率、各阶段floor/clip及原始面流均落盘。
主CSV增加raw-pre/post物种库存N和净电荷C；原生CSV仍保留密度单位m^-3。
机制源单位m^-3 s^-1，面积m²，面粒子通量m^-2 s^-1，dt与体积换算到真实区间数量。

时间层：原始已接受pre -> 实际pre-floor -> 原生pre层通量/反应 -> 原生更新 -> clip/post-floor -> accepted post。
初始化已执行的floor不重复计入第一个区间；之后所有真实记录的修正按物种计入。
电子面流已经扣除二次发射；不把它再当体源。负离子本模型静止。

`DeltaQ + sum(face charge out) - actual reaction net charge - actual numerical correction = R`。
电离和SP3产生等量电子/正离子；附着和复合保留实际化学计量。
物理净电荷源为零，实际浮点反应缺陷单列并仅扣一次。
使用实际三物种RHS的带符号收缩；同时保存原生charge_balance字段，二者的算术舍入差单列，
不是从最终R反推源。没有热源或能量补偿。

| 指标（三步最大，可能不同步） | 实际值 |
|---|---:|
| 单元绝对残差 | 3.52736e-35 C |
| 单元残差/本单元无抵消交换尺度 | 4.37571e-14 |
| 全域L1残差/无抵消交换尺度 | 1.08433e-14 |
| 内部公共面配对误差 | 0 C |
| 反应机制重构与原生RHS最大舍入差 | 1.24386e-38 C |
| 原生charge_balance与扩展精度收缩差 | 4.85682e-42 C |

事前目标：逐单元1e-6，全域1e-8；分母为绝对物种面电荷、反应与修正交换之和，
不是预存大库存。三步交换尺度约7.90e-20 C，内部面按每单元两侧各计一次，定义不混用旧弱功尺度。
本轮不存在零交换单元；制造零交换检查返回NOT_EVALUABLE/绝对残差，不加epsilon。
后处理用long double累加已经舍入的原生double，不能恢复原生求解的额外精度。

[逐单元表](native_cell_charge_ledger.csv)、[逐单元残差](native_cell_residuals.csv)、
[逐步审计](native_cell_continuity_audit.csv)及native/cell_step_*封存索引给出完整证据。
可选selected_cells_diagnostic.csv只选8个单元，明确省略数量且不能作为完整全局源。

## 4. B0几何与参考面

从原生Geometry2DView导出H/Return各297个金属接触面，核对气体/金属mask、法向、面积和面中心。
法向为气体指向金属，面积已包含轴对称旋转，不再乘2*pi*r。
原生实例仍为LOCAL_HV_RETURN，H/Return身份不随极性改变。

与`adapters/comsol_fp2/b0_geometry_boundary_config.json`逐项核对：70 um间隙、200 um曲率、
15度半锥角、1 mm杆半径、30 mm尖端到背面长度。这些来自现有参数合同，不是由本Prompt新设。
本地场域只覆盖局部尖端，不代表整个30 mm电极都在此次网格内。
该轮廓资料为DERIVED、原可编辑CAD缺失、全局物理端子选择仍NOT_RESOLVED。
不同尺寸的旧FP1外域不能直接与本局部零Neumann网络拼接。

[映射](b0_electrode_geometry_mapping.json)保留这些来源和原生mask/面SHA。
真实引线、TENG连接、回流路径、材料和外域响应保持null及明确原因。
未均匀分配Q_H为sigma_s，也未将i_port复制到导线或气体。

## 5. 外域响应与能量

接口采用同一对称电势能网络：A=K_ff，B=K_fT，D=K_TT，q为自由节点电荷C，V为端子电势V。
约束自由块SPD，完整能量块非负；矩阵单位F，不能直接传入混合单位Poisson替换行。
自由方程 `A u = q-B V`；端子泛函 `Q_T=B^T u+D V`。
将自由节点分l/o后，同步得到A_eff、B_eff、D_eff、q_eff。
端子电荷保留 `B_o^T A_oo^-1 q_o` 偏置；不是仅用新的局部u套旧泛函。

全储能严格分为：
`0.5 [u_l,V]^T K_eff_full [u_l,V] + 0.5 q_o^T A_oo^-1 q_o`。
制造例q_o=2e-13 C，H电荷偏置=-5e-14 C，外域源能量5e-15 J；
局部约化能4.718e-12 J，总能4.723e-12 J。
全系统直接求解和消元后的最大电势/端子电荷/总能量差在本例计算中均为0。
该恰好舍入相等不证明任意网络都逐位一致，事前相对目标为1e-12。

TOTAL/SPLIT接口只把同一响应计入一次并拒绝双计，总电荷/总能量相同。
没有新做动态RLC两模式测试，也没有用2U/V²或猜测的外部C改变本地电容。
真实B0没有外域矩阵：GLOBAL_BOUNDARY_OPERATOR_UNAVAILABLE。制造矩阵绝不用于B0物理端口。

## 6. 表面与导线源接口

GlobalSourcePacket 1.1接收显式面片面积/位置/法向、sigma_s(t)、有限体积支撑映射、导体边电流分布、
接触面和支路连接。检查表面库存映射以及
`DeltaQ_metal=Deltaq_gas_into_metal+Deltaq_lead_into_metal`。
新增制造闭环为4个相邻有限体积单元、明确引线及回流、独立给定电流和表面储荷，
三个静态时间样本0/5/10 ps，无电路ODE、无新RF波形或传播。
每区间H/Return储荷分别-5/+5 pC，导线交换匹配，合成气体交换为定义上的0；
电极平衡残差和表面映射误差均为0。

当前生产结构的离散范围是显式面片到相邻笛卡尔FV支撑及刚体坐标变换；
不承诺任意曲线导线、CAD面元到传播网格的保守重映射。
独立surface_current字段缺失时仅在明确提供体导体电流分布的分支可用，缺项不伪填0。
来源绑定全局入口还要求匹配的数值分布payload、实际几何payload和外域算子及其端子电荷核对；
没有实际数据时只报告输入接口存在，不宣称SOURCE_BOUND_GLOBAL_CONDUCTOR已在真实器件上运行。
FULL_GLOBAL_PHYSICAL_SOURCE明确拒绝。

局部Stage F实际Mz：电子约-4.14657925075e-12、正离子约-2.12545514768e-14、总计约-4.16783380222e-12 A*m。
空间源不含第二份位移电流，也不重复发射。60 fs不能用于10--80 MHz谱或接收预测。
本轮没有重复V3.3传播；Stage H历史验证保持原状态。

## 7. 验证、版本与保护

最终26项隔离检查通过：Schur非零源、单位/来源/非SPD/重复电容拒绝、面交换、错误时间层、
参考面、接触、回流、分布SHA、位移/发射双计、能力冒用、稀疏完整性及零交换处理。
首次测试收集缺包标识失败，补齐新目录__init__后通过，未启动额外原生轨迹。
离线分析曾在NumPy整数JSON序列化处失败，失败产物保存在analysis_failed_1，修复后只重做离线分析。
另一次最终检查在计时变量遮蔽处立即失败，修复后通过，无原生或传播调用。

独立编译约13.51 s，唯一原生进程约44.76 s；成功离线审计约10.95 s，最终来源/库存检查约4.55 s。
第一次离线失败没有完整finally计时，明确保留未知实际耗时，并保守扣300 s预算，不伪称精确总时间。
合计计费预算及其余隔离测试/封存成本见build_and_test_record.json，低于1200 s；编译单列。
没有超过3次原生调用，没有新恢复分支、长窗、扫描、COMSOL、热、RLC或传播执行。

与旧V3.2末态只读比较不是逐位一致：ne相对旧最大尺度约3.76e-13、phi约2.78e-11、
场分量约3.5e-10；各字段绝对差均保留于read_only_V32_state_comparison.json。
不把跨构建差异冒称同构建持久恢复通过，也不据此追改旧Gate。
本轮守恒检查完全由本轮真实状态和机制记录完成。

执行时C++源码与最终C++源码相同；之后只完善Python来源/拒绝/面片检查和报告脚本。
执行源码与最终源码各自SHA封存，原生结果适用于executed_scenario所绑定版本；最终Python检查另有记录。
5982个开始时保护文件核对，授权改动仅列出的记录/驱动/源装配/CMake和两份长期文档；
旧库、结果、历史Gate及其他物理源码不得改变。未提交Git。

## 8. 唯一下一项依赖

取得并绑定**同一H/Return参考面下的真实电极接触、引线及回流几何，和由它们确定的外域边界响应/表面电荷空间映射**。
现有接口可接收这些明确数据；没有该输入，不值得再增加局部短窗来猜完整辐射源。
离子/表面标定、源/开关身份和接收链仍分别未认证。本次停止，不自动追加计算。
