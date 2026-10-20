# I1-COLD-PHYSICS-V3.2 实施报告

日期：2026-10-19。本轮完成真实代码扩展、原生短窗、子进程恢复和 Stage F 实际源消费；不是只读审查。
实际状态分别见 [Gate](v3_2_gate_summary.json)、[执行记录](native_build_and_test_record.json)及[来源索引](source_index.json)。

## 1. 分项结论

| 项目 | 结论 | 限定 |
|---|---|---|
| ION_FIELD_DEPENDENT_TRANSPORT_PROVIDER_IMPLEMENTED | PASS | 常参数兼容 + E/N线性插值 + 显式K0换算 + 域外拒绝 |
| PHYSICAL_ION_TABLE_COVERAGE | SOURCE_BACKED_HIGH_FIELD_TABLE_NOT_AVAILABLE | 未取得本次物种/气体/高场匹配的原始数值表 |
| SURFACE_YIELD_PROVIDER | PASS_ENGINEERING / NOT_CALIBRATED | 来源、物种、材料、表面、入射能量和未知能量状态检查 |
| SPECIES_RESOLVED_STAGE_F_SOURCE_EXPORT | PASS | 原生接受面通量到实际 Stage F SourceRecord 与电流矩 |
| V3_WEAK_ELECTRICAL_WORK | PASS_IN_TESTED_LOCAL_WINDOW | 仅新0--60fs同网络、合成输运参数的独立诊断；旧19.2%不追改 |
| B0_SHORT_INTEGRATION_AND_CHECKPOINT | PASS | 5次实际step、20fs保存、新进程恢复、60fs终点 |
| REAL_DEVICE_CALIBRATION | NOT_CALIBRATED | 实际离子混合物/表面/TENG源仍未知 |
| GLOBAL_PHYSICAL_PORT | NOT_VALIDATED | 当前仍为局部网络及气体源 |
| RF_PROPAGATION_AND_RECEIVER | NOT_RUN / NOT_VALIDATED | 没有调用传播、全波或接收计算 |

## 2. 源码与离子表

在原 `ColdPhysics.hpp/.cpp` 增加 `IonFieldTable` 和 `PositiveIonTransportProvider`。
`IonTransportModel` 保留原常参数字段和默认路径；新增表模式从实际单元 abs(E)/N 求系数。
气体面漂移速度取两侧 mu*E 的对称平均，扩散系数取对称平均，沿用原保守上风/Fick通量及正性dt控制。
金属面沿用真实半单元吸收及原法向，未修改电子ISG0、MorrowLowke、反应、Poisson或SP3公式。

表合同包含物种、背景气体、Td、mu或K0单位、D_parallel/D_perpendicular、T/p/EN范围、来源DOI/数据库ID/文件SHA、A/B/C/D、近似与拒绝规则。
线性插值不外推；K0=mu*N/Nref，Nref必须显式提供。缺扩散值时返回 `ION_DIFFUSION_UNAVAILABLE`，不补0或Einstein估值。
提供器可表示两个扩散方向，但当前FV仍为标量扩散；两方向不等时返回 `ANISOTROPIC_DIFFUSION_FV_NOT_IMPLEMENTED`，没有忽略张量交叉项后冒称实现。

原 `np` 是聚合正离子，不携带O2+/NO+/簇离子反应分布。
即使未来取得单离子资料，也必须标 `EFFECTIVE_SPECIES_APPROXIMATION`；不能把它当作已辨识混合物。
本轮统一入口只接受固定B0有效空气合成身份，拒绝把该fixture改标签为另一离子或气体。

[文献访问记录](literature_access_record.md)区分原始摘要、数据库平台说明与未取得的数值表。
Sinnott研究提供物种差异背景，Pitchford/LXCat说明数据类型，但本次没有取得适用的高场原始数组。
因此只运行 [SYNTHETIC_TEST_TABLE](synthetic_ion_table.csv)：0/100/100000 Td的mu为2e-4/2.1e-4/2.2e-4 m2/(V*s)，两个D均5e-6 m2/s。
这些节点为制造软件测试，不来自论文、不拟合B0，也未因结果调整。
运行中全气体单元实际 E/N 范围约 **2.46e-12--347.272 Td**；宽的合成表支持域不等于物理高场覆盖。
详见 [离子来源/域表](ion_field_transport_provider_source_and_domain.csv)。

## 3. 表面提供器

原 `SurfaceEmissionModel` 增加可选严格来源合同，V3.2入口必须启用。
绑定入射离子、气体、材料、表面状态、入射能量假设、来源SHA、证据等级和不确定性。
非C类缺材料/表面/有限入射能量域时拒绝；没有默认gamma或从理想gamma_critical替换产额。
实际测试gamma=.02仍为C类；材料、表面、入射能谱未知。

发射公式仍由真实壁面离子通量驱动。零入射、gamma=0、H/Return极性身份的独立测试通过。
电子面记录已经扣除二次发射，Stage F不再添加第二份发射电流。
发射材料能、电子初始能谱和近电极压降保持UNAVAILABLE，详见 [表面来源/域表](surface_yield_source_and_domain.csv)。

## 4. Stage F实际消费

新增 `python/streamer_rf/rf/source/accepted_species.py::load_accepted_species_packet`，不修改旧电子专用读取器的含义。
真实执行链为：ColdFullProvider原生接受记录 → 净电子/正离子面CSV → 分物种中心重构 → SourceRecord →
既有 `rotate_axisymmetric_record(8)` → 既有 `current_moment`。
**没有运行Jefimenko传播或全波求解**；旋转只用于轴对称电流矩积分，径向分量在完整方位角积分中抵消。

保留原始法向面流：Je=-e*Gamma_e_net，Jp=+e*Gamma_p，J=Je+Jp；负离子为IMMOBILE_BY_MODEL。
中心重构为Jr=(J_rplus_out-J_rminus_out)/2，Jz=(J_zplus_out-J_zminus_out)/2，逐物种独立计算。
这是中心诊断重构，不把中心差分重新称为原保守散度；原始面仍为连续性核对的权威量。
轴对称体积直接取原生环带体积，没有重复乘2*pi*r。

数据包保存geometry/mesh/operator/reference_plane/source SHA、单位与 `ACCEPTED_INTERVAL_PRE_FLUX`：
本次导出最终接受区间 **40--60 fs** 的pre层面流及同层rho，不伪称60fs瞬时通量。
终端总电流和位移电流不在空间J中，端口量另在FP3账本。位移响应不作为第二份外加J注入Stage F。

实际z向电流矩：

| 源 | Mz [A*m] |
|---|---:|
| 电子 | -4.146579250744071e-12 |
| 正离子 | -2.1254551476912927e-14 |
| 总传导 | -4.167833802220983e-12 |

分物种相加、原始面到中心重构、电极面电荷对原生收集/发射记录均独立核对。
见 [源包](stage_f_source_packet.json)、[消费审计](stage_f_source_consumer_audit.csv)、[原始面流](species_resolved_current_faces.csv)。
单个极短接受区间不足以评价10--80MHz谱；这些量不是接收电压、RF功率或完整全局辐射源。

## 5. 独立弱电功

不修改FP3旧算法。复用已支持电子及正离子的 `ActiveNetwork` / `StableNetwork`，由新只读接受记录钩子输出
`ACTIVE_NETWORK_SPECIES_STABLE_WORK_V3_2` 版本化诊断，不回写原生场/密度或旧q_A/q_B。
同一链接网络独立计算 r_phi=q+bV-Au、r_psi=b-Apsi及e_Q收缩，保留V1与stable_V2电荷。
对同一rho另解零导体电压的u_rho；比较电荷投影和电极链接响应，不用u-Vpsi作为唯一参考。
电功来自实际电子/正离子面电荷；发射已经在净电子面中。数值修正和反应电荷缺陷分别与势配对。
场能、交叉项、C_E-Cp均按原同网络定义保留，不换电容、不用残差反推传导功。

[事前计划](execution_plan.json)对弱电荷及有符号/无抵消面功均登记1e-3目标。
实际新窗结果如下，最大值不一定来自同一区间：

| 诊断 | 最大相对值 |
|---|---:|
| stable电荷预算 / 无抵消原生输运电荷 | 2.36245e-15 |
| stable电功残差 / abs(Wcond_signed) | 2.99921e-11 |
| stable电功残差 / Wcond_abs | 2.99592e-11 |
| 原总势电功残差 / abs(Wcond_native) | 0.113058 |
| abs(R_E_native-ubar*Delta r_phi) / abs(Wcond_native) | 5.58329e-9 |

本窗每区间无抵消电荷约3.66e-20 C，Wcond约5.47e-19 J；stable残差约-1.64e-29 J，未用预存大场能归一。
在40--60fs区间，原残差-6.1844205182e-20 J，独立ubar*Delta r_phi为-6.1844204766e-20 J。
这给出了**新窗同网络原总势残差的独立线性求解误差归因**，不是从最终预算倒造修正项。
零电压链接与源投影电荷端点差最多约1.21e-31 C；原总势与响应重构势最大差约1.61e-8 V。
完整带符号项、参考误差与非零分母见 [独立残差表](v3_weak_work_independent_residuals.csv)。

旧V3的0--1ps/19.2373%记录没有重算或追改。新窗不同参数/时间不能宣称旧轨迹已通过。
新出口是明确精度身份的同网络诊断，不是本装置的气体热、发射能量或全局端口物理认证。

## 6. 集成、恢复与保护

采用与V3.1相同原生短窗驱动，在默认关闭的V3.2钩子下接入表配置、接受记录诊断与源输出。
没有复制StreamerSolver、Poisson或RLC。`SourceBoundI1ScenarioRunner`复用V3.1身份/资源/原子接受checkpoint协议。
外层完整参数包和CSV哈希随checkpoint配置身份绑定；恢复重新验证输入后构造同一模型。
一次连续进程3步，20fs时保存后继续2步；新进程恢复该接受态后2步，最终60fs。
**合计5次step、3个唯一接受节点**。九场、历史、warm、权重、电子/离子面、累计账本及接受trace均逐位一致，无重复累计。

20项新C++数学/域/通量测试、11项运行前Python检查先通过，才启动原生烟测。
后续仅增加Python固定fixture身份与Stage F拒绝硬化及其隔离测试，未重跑原生场。
最终拒绝测试曾因项目与/tmp跨文件系统、夹具硬链接失败而未执行到拒绝断言；改用独立副本后复查，原失败日志保留，未修改科学参数或源数据。
运行版源码与最终Python层分别封存；当前source hash改变时旧执行输入应拒绝，不静默冒称同版本。
TOTAL/SPLIT互斥规则和持久模式身份保留；实际原生测试TOTAL，本轮不重跑旧双模式RLC矩阵。

科学执行约51秒，消费/隔离分析秒级，编译约13秒单列；精确已记录成本见执行记录。
只有原生5次step及同网络诊断的7次独立空间电荷RHS，无250ps、热/G1/G2、TENG演化或RF传播。
新建独立 `build_i1_cold_physics_v3_2`；旧库、报告、结果和历史Gate均按保护清单核对，不提交Git。

## 7. 最小后续依赖

当前已有可实际消费的**局部气体全运动物种J**及限定同网络弱电功诊断。
要进入可信RLC/RF，最小空间接口缺口是把同参考面的电极/引线回流、表面电荷和外域响应与气体J组成一致的全局电荷连续源/端口包。
不能以gas-only电流矩替代整根电极/回路源，也不能将i_port直接复制为J。
物理参数方面仍须来源匹配的高场离子输运、真实离子组成及材料表面/入射能量依据；这与软件源包验收分开。
本轮停止，不自动增加传播、参数扫描或形成/恢复模型。
