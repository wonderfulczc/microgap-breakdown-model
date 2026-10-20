# I1-COLD-PHYSICS-V3.3 全局守恒源与 F/H 接口实施报告

日期：2026-10-19。性质：真实接口实施、既有B0数据离线消费、两类合成源传播制造检查。
没有新增StreamerSolver/Poisson/SP3时间推进，没有重跑物理电路、热模型或实验波形拟合。
本次新场量仅为 `SYNTHETIC_ENGINEERING_REFERENCE`，不是70 um装置的RF预测。

## 1. 结论分层

| 项目 | 结论 | 范围 |
|---|---|---|
| GLOBAL_SOURCE_PACKET_CONTRACT_IMPLEMENTED | PASS | 版本化数据、语义拒绝、完整合成源及局部源两种入口 |
| SPECIES_FACE_TO_SOURCE_CONTINUITY | CONDITIONAL_PASS | 原生配对面/全域修正账本通过；逐单元修正及反应缺陷未归档 |
| ELECTRODE_LEAD_RETURN_ASSEMBLY_SYNTHETIC | PASS | 有限端点储荷与闭合回流两类源，独立节点和电极平衡 |
| STAGE_FH_ACTUAL_SOURCE_CONSUMPTION | PASS | Stage F实际时域传播、Stage H场投影；不是接收电压 |
| DISPLACEMENT_AND_EMISSION_NO_DOUBLE_COUNT | PASS | 净电子面流包含发射，空间J不含第二份位移电流 |
| REAL_B0_GLOBAL_SOURCE_READINESS | NOT_READY | 真实引线/回流、全局表面电荷及参考面映射缺失 |
| EXPERIMENTAL_RF_READINESS | NOT_VALIDATED | 没有实验源/传播/接收链标定 |

详见[分项Gate](v3_3_gate_summary.json)、[能力范围](global_vs_local_source_capability.json)。
历史V3/V3.1/V3.2、FP3/FP4、TENG、Toepler/RW、C1/V5/V6/V7结论均不追改。

## 2. 实际新增代码与调用链

新增 `python/streamer_rf/rf/source/global_packet.py`：

- `GlobalSourcePacket`：实例/网格/算子/边界/参考面/时间/源参数与构建身份。
- `ConservativeGraph`：原始有向FV面电流、区间电荷与独立节点储荷。
- `I1GlobalSourceAssembler`：局部诊断或完整合成自由源装配、同面电极/导线账本检查。
- `audit_native_faces`：逐单元原始余额、内部配对及外边界交换；未知修正保持null。
- `save_global_source_bundle` / `load_global_source_bundle`：新目录发布、SHA与完成标记、拒绝损坏/不完整源包。

源schema位于同模块目录；执行和拒绝检查位于 `adapters/i1_cold_physics_v3_3/`。
这是独立薄适配层，没有复制FULL、Poisson、RLC或电磁传播求解器。
没有改动旧 `accepted_species.py`、ColdPhysics、ColdFullProvider、FP3稳定电功、Stage F/H内核。
既有文件只修改两份长期文档。

实际合成链：
`ConservativeGraph -> GlobalSourcePacket -> I1GlobalSourceAssembler -> SourceSeries`
`-> jefimenko.solver.evaluate_waveform -> fullwave.native_receiver.project_field`。

Stage H本轮只投影E，输出仍为V/m。没有调用天线、H_rx传递函数、S参数或openEMS；
没有生成接收电压。旧H1/H2并联G3参考面及历史频带合同没有被替换为I1串联拓扑。

## 3. 守恒与空间表示

[完整合同](electrode_surface_lead_return_contract.md)定义各符号和缺项。
有向边a到b的流出incidence为a处正、b处负；节点电荷独立提供：

`Delta q_node + D * Delta q_edge = 0`。

对H和Return永久身份，分别检查：
`Delta Q_metal = Delta q_gas_into_metal + Delta q_lead_into_metal`，
且 `Delta q_Hg = -Delta q_gas_into_metal`，满足Part VI端口方向。
源装配不会因电流反号重命名电极。

本版制造源采用相邻、等截面笛卡尔有限体积单元的有向面网络。
最终入口显式拒绝非邻接/斜线图冒充该FV离散；不承诺任意导线网格的保守重映射。
每单元为1 mm立方体，端点电荷分布在有限体积端帽中，无零半径或未正则化点源。
`GRAPH_HALF_EDGE_VOLUME_V1`把每条边的I*位移向量各一半分配至相邻单元，
产生体积平均J，保留原始面流为守恒权威。这是有限体积中心求积，并非近导线连续解认证。

完整自由源使用rho与传导J，电荷积累的位移响应由场解体现，不能另加partial_tD。
TOTAL/SPLIT保持同一份物理自由电荷和J；本轮只检查所属枚举/双计拒绝及同源一致性，
没有重做历史RLC两模式矩阵，也没有从制造几何反推真实Cp。

## 4. 两类制造源及独立传播检查

事前计划见[execution_plan](execution_plan.json)。
两类均为1 A幅值、50 MHz正弦导电电流，源时间0--100 ns，1001个源节点，dt=0.1 ns。
Q(t)=(1-cos(omega*t))/omega；区间面电荷使用解析Q差，独立给定节点库存。
梯形电流积分与解析区间积分还单独核对，不用改节点电荷来消除时间离散误差。
这里的50 MHz是制造输入，不是反演的实验频率，不做FFT或接收谱。

### A. 有限端点短导电段

21个节点、20条边，中心跨度20 mm，沿z方向；端点储荷为-Q和+Q。
它有明确的电容性端点回流/电荷释放，不伪造一条未给出的导电回路。
解析电偶极矩 p=0.02 Q，pdot=0.02 I，pddot=0.02 dI/dt。
独立赤道场公式：

`Ez = -1/(4*pi*epsilon0) * [p/r^3 + pdot/(c*r^2) + pddot/(c^2*r)]_(t-r/c)`。

与既有Jefimenko的空间迟滞积分对比，包含近场/感应项，不只比一条辐射公式。

### B. 明确导线与回流的闭合回路

60个有限单元构成x-z面10 mm乘20 mm矩形闭环，回流边明确标记RETURN。
各节点入出电荷相等，节点净储荷为零，所有边参与空间源。
独立磁偶极矩 m=-y*(0.01*0.02)*I；用磁偶极迟滞场检查赤道E。

需要纠正一个容易混用的指标：闭合无积累回路的 `integral J dV` 应接近0，
不能要求改变布局后它必然改变。该回路的磁偶极矩与方向场才反映空间布局。
本轮在同一传播批次中计算赤道和法向观测点；利用刚体旋转协变关系，
这对应在固定测点旋转整套回路的对照，不声称额外执行了第二个回路布局。
电流相同而方向场不同，不能把端口电流单值视为完整空间辐射源。

### 实际结果

| 指标 | 短段A | 闭环B |
|---|---:|---:|
| 最大单元区间电荷残差 [C] | 0 | 0 |
| 全局净电荷变化 [C] | 0 | 0 |
| 10 m独立解析E相对差 | 2.87780e-4 | 2.04009e-3 |
| 20 m独立解析E相对差 | 2.20131e-4 | 8.02628e-3 |
| 10 m / 20 m辐射项幅值比 | 1.99986738 | 2.00575134 |
| 10 m赤道峰值E [V/m] | 0.06309745 | 0.0006610852 |
| 10 m赤道最大B [T] | 2.10459e-10 | 2.18569e-12 |

所有场比较满足事前1%探索目标；B的20 m误差约0.803%，不是机器精度认证。
有限体积、紧致偶极近似及原内核时间插值误差共同存在；本轮不追加网格/时间收敛扫描。
A以已知输入频率检查辐射项相位，最大等效时延误差约6.11e-15 s，低于0.11 ns目标。
距离10/20 m的光行时分别33.3564/66.7128 ns。B未另作单独时延拟合，字段为null；
其独立迟滞公式对比包含相位误差，不虚构另一项延迟结果。

B峰值磁矩为-2.0000000000000025e-4 A*m2，与解析值-2e-4相符；
总电流矩最大约1.53e-18 A*m，法向E峰值2.77e-17 V/m，
法向/赤道比4.19e-14。瞬态B及近场不是无线接收电压。

完整数值见[制造源指标](manufactured_source_continuity_and_field_metrics.csv)、
[A场时序](segment_propagation_fields.csv)、[B场时序](loop_propagation_fields.csv)。
每类仅调用一次已有evaluate_waveform，三个观测点各601个采样，共两批3606个FieldSample。
输出观察时间70--130 ns包含传播延迟；源激励仍严格限于0--100 ns。

## 5. 真实V3.2局部数据接入

只读旧40--60 fs接受区间、200x400原网格、pre层净电子与正离子FV面流，
以及step2/step3九场二进制端点；文件SHA由旧数据包和本轮保护清单约束。
独立重构中心J与旧中心CSV逐值一致；轴对称旋转只分配环带体积，不再乘2*pi*r。

| 分量 | z电流矩 [A*m] |
|---|---:|
| 电子 | -4.146579250744071e-12 |
| 正离子 | -2.1254551476912927e-14 |
| 总传导 | -4.167833802220983e-12 |

分量相加最大绝对差5.87e-28 A*m。内部相向面电荷配对误差为0。
物种面常规电流恒等式最大绝对舍入差5.68e-14 A/m2；未把该中心重构当保守散度。
逐面二次发射已在Gamma_e_net中；Return发射数量1.090466531981215e-9与原生账本相符，
没有再次注入体电流。H/Return吸收、离子收集、发射见[电极审计](native_electrode_exchange_audit.csv)。

全域原始余额3.8028547950850696e-26 C；独立累计数值修正差3.802854811928518e-26 C；
反应电荷缺陷-4.72817e-39 C。扣除已知项后残差-1.68430e-34 C，
相对无抵消面电荷7.32409e-20 C为2.29967e-15。
该分母按每气体单元外向面求绝对值，内部面出现两次；与V3.2独立面尺度不同，未偷换旧门限。

**逐单元**floor/clip和实际反应电荷缺陷数组未随旧输出落盘，
因此[逐单元表](native_species_face_charge_mapping.csv)只列原始余额，
`physical_source_C / numerical_correction_C / qualified_residual_C`保持空值及明确状态。
光电离配对净电荷为零是化学计量关系，不代替实际浮点源缺陷的逐单元记录。
没有由残差反算一个“守恒修正”，也没有重跑FULL补数据。

同局部参考面的stable Deltaq_port=1.0940115004928722e-21 C，
DeltaQ_rho=1.0915126816203065e-21 C，Deltaq_Hg=2.4988188725658097e-24 C；
端口恒等式舍入差约-1.31e-37 C。电势恒500 V，所以本区间几何电容端点增量为零。
这些数字来自旧独立账本，不是新全局物理电流。
历史q_A=8.487244748367312e-12 C、q_B_native=8.655161854172055e-12 C、
active_V1=8.626188398774207e-12 C仍保留不同身份。
稳定V2采用Cp*V与空间电荷投影及稳定区间增量；不把局部q_B当全局表面电荷。

局部数据实际读取链与三类电流矩见[消费审计](stage_fh_actual_consumer_audit.csv)。
未从20 fs区间生成10--80 MHz谱、场传播、气隙辐射功率或接收电压。

## 6. 拒绝、来源与保护

运行前22项隔离检查通过。首次测试收集因测试装饰器语法错误停止，未执行传播；
修复后的重启仅限preflight，原失败日志保留，没有原生/PDE重试。
两次传播完成后只新增FV邻接、源包SHA/完成标记等守卫与隔离检查，不重做传播。
实际运行版与最终守卫版分别保存在executed_sources和final_sources，SHA不混称同版。

拒绝覆盖：错误参考面、缺回流/导线/表面电荷、错误单位、发射双计、位移双计、
TRIAL时间层、重复电容、错误节点电荷/面电流、错误电极端口交换、
真实装置资格冒用及不支持频谱请求。完整源包必须通过数组与空间映射检查，不接受只有标签的JSON。

执行约309.77 s，另有秒级离线逐面核对和最终守卫/保护检查；总成本见
[执行记录](build_and_test_record.json)与finalization_record.json。
两次传播各约152/155 s；峰值进程RSS约373300 KiB。未编译，原生step次数0。
900 s预算没有扩展。没有全波、大规模FDTD、热/TENG/Toepler/RW或历史测试重跑。

[来源索引](source_index.json)绑定实际GENERIC_DEV、输入、旧二进制及新执行代码；
[保护前](protection_before.json)/[保护后](protection_after.json)检查旧文件字节。
除两份授权长期文档外，原有源码、库、结果和历史Gate保持不变；不提交Git。

## 7. 下一项唯一依赖

优先取得并定义**与局部H/Return参考面一致的真实电极、引线、回流及表面电荷空间映射**，
同时明确高频源/开关拓扑与支路电流身份。现有真实资料不足时只输出LOCAL_GAS_ONLY，
不能猜导线形状、复制i_port为J或把本次合成端帽替代装置电极。
逐单元物理源/数值修正数组的持久导出是后续真实全局源连续性资格的一项明确接口缺口，
本轮不为补齐它新增轨迹。
离子/表面参数标定、自动冷热、全局端口及接收传递函数仍为独立物理依赖。
本次完成即停止，不自动开启下一轮计算。
