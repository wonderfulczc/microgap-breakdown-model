# I1-COLD-PHYSICS-V3.5 冷热独立模式统一架构

日期：2026-10-20。HEAD只是历史提交身份，实际使用GENERIC_DEV源码与既有V3.4库；
完整执行/最终哈希及增量见source_index和保护记录。没有重建或覆盖旧二进制。

## 1. 结论

PHYSICAL_ASSEMBLY_INPUT_CONTRACT、COLD_FULL_UNIFIED_MODE、GIVEN_THERMAL_SERIES_UNIFIED_MODE、
UNIFIED_STATE_AND_TIME_CONTRACT均为限定工程PASS。
AUTO=MISSING_PHYSICAL_CLOSURE；真实全局几何WAITING_DATA；全局物理端口NOT_READY；RF接收NOT_VALIDATED。
历史V3.3单元CONDITIONAL_PASS和V3.4/FP3/FP4/TENG/Toepler/RW等Gate均不追改。

## 2. 真正新增代码与复用

`adapters/i1_cold_physics_v3_5/assembly.py`提供I1PhysicalAssemblyPackage，引用旧I1ConductorGeometryAdapter。
外部几何、拓扑、元件位置、第三导体/地、测点负载、静电响应、动态表面映射和器件标定分开保存。
已知值绑定文件SHA和内容；未知保持null+reason。实际B0局部包引用旧原生金属面与轮廓合同，
不复制源文件、不假设实验引线。合成完整连接只验证输入，不冒充外域算子或表面电荷解。

`scenario.py::I1MultiphaseScenarioRunner`继承V3.4入口。冷模式调用既有run_segment及其
ColdFullProvider原生驱动；`cold_output.py`消费新接受记录、单元账本和现有Stage F读取器、
GlobalSourceAssembler、轴对称旋转及current_moment。没有新传播。

热模式调用原GivenThermalInitialStateAdapter构造，然后原ThermalSeriesProvider真实
evaluate_trial/commit；原RadialThermalSolver执行热步。未直接覆写U、time或累计账本。
未复制Poisson、SP3、化学、G1/EOS、RLC或传播RHS。
首版热入口仅接受原提供器已有的10V理想阶跃、100nH/3ohm/串联101pF和显式零回路初值；
其他源/初值，包括请求TENG热耦合，返回UNSUPPORTED_SERIES_FIXTURE，不被悄悄替换。
TENGFrozenPositionSource的旧V-Q-x接口原样保留，真实CT/Voc/qT0仍未知。

`state.py`提供ScenarioStateEnvelope和切换拒绝。新段统一时间坐标而非统一物理能量库存。
输入/状态schema是声明结构，实际Python语义检查负责来源、单位、支路和域校验；
未安装依赖，也不声称通用Draft元模式认证。

## 3. 实际有界运行

| 项目 | 冷模式 | 给定热模式 |
|---|---|---|
| 实际接受/调用 | 3/3原生步 | 4/4 G1试探步，无拒绝 |
| 本地终点 | 6e-14 s | 2e-10 s |
| 初态与几何 | 原B0+500V、SP3 ON、C类离子/表面输入 | 原150单元、0.5mm径向、**1mm**人工LTE柱 |
| 电路 | 恒压局部场，无新RLC | 原10V合成单串联回路，Cp=0.2pF，TOTAL |
| 场/材料 | 同局部网络场能；材料能不可用 | EOS内能和径向动能；无解析静电空间场 |

两段分别标V35_COLD_INDEPENDENT和V35_THERMAL_INDEPENDENT，时间锚均为人工1us，
不是物理冷热交接时刻。global_time=anchor+local_time逐值核对。
源包保留原生PRE_FLUX时间层，不能把40--60fs平均流当60fs瞬时流。

冷步最大逐单元无抵消电荷残差比4.37571e-14；稳定弱电功比3.11697e-11。
原总势与稳定出口同时保留。源消费者实际读取电子和正离子，发射已在净电子流内，未重复位移流。
不把局部q_B当全局Q，也不由G_R生成热源。

热柱R由241.086152变为248.148222 ohm；
末回路电流0.0170709643 A，区间末Icol=0.0134125629 A，
区间末Cp支路=0.00365840133 A；更新热态的Icol=0.013316639 A。
R_used、R_end、区间RMS、区间末与更新后量全部分列，不能交叉使用时间层来核对KCL。
径向速度最大7.196697 m/s，未被清零。
累计QJ=2.54051166e-12 J。电路/热/跨域/总预算最大残差相对实际交换量
5.91311e-09，低于事前1e-3；预算代数恒等式误差
7.14714e-17。没有用初始大热能掩盖交换误差。

只新增这一条四步热轨迹；不重跑两模式等价矩阵。TOTAL/SPLIT拒绝语义复用历史合同。
独立热参考R_map仍NOT_APPLICABLE，不能写成冷到热能量误差为零。

## 4. 热柱空间源及持久状态边界

输出THERMAL_COLUMN_REFERENCE：实际sigma(r)、G_length、ell、Icol和Jz=sigma*Icol/G_length。
面积积分电流相对误差2.60535e-16。
单位S*m与A/m2明确；环带面积不再乘2*pi*r。
它是给定轴向均匀LTE局部柱的条件性电流分布，不是完整自由电荷/导体回流源。
G1没有自由电荷状态，因此rho_free=null，未填0以强行生成传播SourceRecord。
本轮实际Stage F消费限冷模式局部全物种源；热模式输出条件性柱分布，未做传播/接收。

冷Envelope引用实际九场、机制/面流、累计账本及接受payload。原V3.1 seal/verify恢复入口仍可用；
本轮未重测进程重启，也未把未封装的payload或Envelope本身当成可直接恢复的发布检查点。
热Envelope引用U/T/p/ur/sigma及原生账本，kind=DIAGNOSTIC_SNAPSHOT；完整热持久恢复未实施。
这不影响既有进程内热trial/commit隔离。

## 5. 装置输入和唯一后续依赖

实际名义净隙70um、既有参数化B0局部面和参考面可读取。1mm人工热柱与它**无物理映射**。
真实材料、引线、回流、屏蔽/接地、元件位置、探头负载、CT/Voc/qT0、外域响应和表面分布仍缺失。
完整合成装置测试只证明有限图连接与参考面校验。完整真实源申请被准确拒绝，
不是因为资料不足而拒绝整个工具运行。

详见中文测绘指南。下一项唯一优先工作：按同一H/Return参考面收集并绑定实际装置几何、
高频连接及测点输入包，再交由现有外域/源装配合同检验；不追加人工热态扫描。
自动冷热仍缺冷材料能/质量动量、空间投影、LTE依据、共同能量零点与源/电容/场库存匹配。

## 6. 复现、保护和停止

运行入口：`env PYTHONPATH=python:. .venv/bin/python -m adapters.i1_cold_physics_v3_5.scenario <scenario.json> <新目录>`。
当前执行包绑定确切源码/库；源码变化时拒绝，不能用旧包冒充新版本运行。
具体命令及成本见build_test_record.json，原始日志、输入、执行源码、最终SHA均保留。
本轮两模式都由同一入口启动；23项隔离前置检查通过，后续仅离线检查真实输出与schema。
没有新编译、长轨迹、参数扫描、RF传播或依赖安装。MPI沙箱socket告警保留于stderr，原生正常退出。
历史文件仅允许两份长期文档追加证据，其他保护文件必须字节不变。
任务至此停止，不提交Git、不自动开启后续计算。
