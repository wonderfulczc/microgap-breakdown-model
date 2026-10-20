# I1-P1 装置参数化与双源电路实施

日期：2026-10-20。身份：条件性工程参考，不是本装置TENG、气隙电阻或RF认证。
本轮完成实际代码连接与两条有界源—电路计算；历史Gate不追改。

## 1. 来源与几何

实际分支为architecture-v2.1-ultrafast-rf，HEAD为8ec25920adff63a0b1a875c620b9a32105458ad7；
运行身份还包含未提交GENERIC_DEV和适配器。执行源码副本、哈希、授权增量及保护记录独立保存。
工作区未找到`I1_parameter_assessment_20261020.md`。本轮没有伪称读取评估书，
以本次用户消息形成`input/user_transcription.json`，不是原始量测报告。

`I1DeviceCasePackage`逐项绑定原始JSON内容、单位、来源分类及SHA。
名义70um、3ohm、100nH、101pF与R1/DAG1/D2转录值分开；没有MEASURED_VERIFIED器件参数。
未知材料牌号/表面、第二电极对称性、支座完整形状、引线/回流、探头及外域/表面映射保持null+reason。
包引用V3.5 I1PhysicalAssemblyPackage与原FP2/FP1 B0来源，不复制几何内核。

| 分支 | 尖端半径 | 圆锥全角 | 身份 |
|---|---:|---:|---|
| B0_NOMINAL | 200um | 30deg | 保留原来源模型 |
| B0_MEASUREMENT_CANDIDATE | 条件140.6um | 条件25.3646deg | 仅当R1单位为mm、DAG1为全角 |

候选半角12.6823deg；R1原值0.1406与旧0.1263分别保存。D2=1.0114没有尺寸映射。
候选只生成参数比较，不生成CAD/mask或改用原生FULL几何；未经明确假设调用会拒绝。
真实70um名义净隙不是TENG机械分离距离，也不能认证有效热柱长度。

## 2. 真正实现与复用

新增`adapters/i1_p1/device.py`、`sources.py`、`amplifier.py`、`network.py`及输入/测试/归档脚本。
旧`adapters/i1_cold_physics_v3_5/scenario.py`只增加SOURCE_NETWORK_REFERENCE分派和可选共享装置/源匹配守卫。
没有修改原TENGFrozenPositionSource、TENGSeriesCircuit、FP3端口/求根/时间调度、G1或原生C++。

`TENGContactSeparationSource.freeze`从几何计算CT/Voc，真正构造原TENGFrozenPositionSource；
原TENGSeriesCircuit联合推进QT、电感/外接电容和数学气隙端口，调用evaluate_trial/commit。
不是预先给电压轨迹后回放。每条10,000接受步、30,000次残差/端口评价。

放大器`AmplifierSeriesCircuit`继承原SeriesPortCircuit的Newton、trial和调度，
仅增加Cext并联泄放支路与可选串联Rout+Lout的残差/账本，不复制RLC积分器。
放大器波形支持来源绑定的分段线性时间表、禁止外推；跨波形节点的区间须先拆分，否则明确拒绝。
动态输出阻抗首版限Z=Rout+sLout，有效域显式声明；未知时只能请求IDEAL_AMPLIFIER_LIMIT，
不宣称支持任意传递函数或实测S参数。R+sL仅作单次试探检查，无第三条完整轨迹。

## 3. TENG制造合同及静态参考

重叠面积1.6e-3m2，上/下铜片4x4cm和4x4.5cm，铜厚10um；
额外暴露铜片不计有效摩擦面积。PET/FEP各50um。胶10um的真实位置未知。
测试明示PET/FEP在介质路径、胶不在路径；另作胶在路径的公式测试，不认证实际层叠。

本地缺评估书介电先验，故运行前登记epsPET=3、epsFEP=2、epsadh=3、
sigma_eff=1e-6C/m2，全部SYNTHETIC_ENGINEERING，**不是文献先验或实测值**。
没有安装物性数据、没有拟合这些数值。常介电假设不提供真实频率色散资格。

采用用户给定理想平行板关系：
`d_eff=sum(t/epsr); CT=eps0*A/(x+d_eff); Voc=sigma*x/eps0; Qsc=CT*Voc; VT=Voc-QT/CT`。
CT与Voc不是相互独立的拟合输入；QT为外电路转移电荷坐标，实际由电路积分。
首版冻结x与sigma；泄漏显式关闭，只属参考。动态机械、接触充电历史或真实漏电模型请求拒绝。

| 合成xT | CT | Voc | Ceff(串联101pF) | 无阻尼f0，Cp=0 |
|---|---:|---:|---:|---:|
| 0.1mm | 100.000239pF | 11.294091V | 50.248817pF | 70.999814MHz |
| 1mm | 13.600032pF | 112.940907V | 11.986064pF | 145.372360MHz |

理想电压源极限f0=50.079438MHz。上述f0不适用于含0.2pF的三状态系统的精确模态，
也不是两条10ns波形的频谱识别，更不是实际40–50MHz标定。
`0.5 CT Voc^2`分别6.377839nJ和86.738616nJ，
仅为固定源坐标从QT=0到Qsc的理想净交付尺度/对应坐标库存，不能称实验可用能量或初始机械能。
完整两位置表保留Qsc、Woc及能量坐标；真实QT0和vC0仍为null。

## 4. 两套实际运行

共同数学端口：200ohm并联0.2pF，不是本装置电阻/电容；L=100nH、Rs=3ohm、串联Cext=101pF。
两条均0--10ns，步长1ps；总计20,000个电路接受步。没有新冷/热步。

| 项目 | TENG目标拓扑参考 | 放大器替代拓扑参考 |
|---|---|---|
| 源 | 上述xT=0.1mm冻结几何源 | 1V恒值IDEAL_AMPLIFIER_LIMIT |
| 连接 | TENG+→GAP→L→Rs→Cext→TENG- | AMP+→GAP→L→Rs→(Cext并联Rbleed)→AMP- |
| Rbleed | 不存在 | 1Mohm，显式合成，绝非Rs |
| Cp归属 | TOTAL_CHARGE | SPLIT_GEOMETRIC_CAPACITANCE |
| 初态 | i=Vc=QT=0，Vp=11.294091V | i=Vc=0，Vp=1V |
| 峰值绝对环路电流 | 0.050005835A | 0.004612051A |
| 累计数学柱支路耗散 | 2.694200898nJ | 0.030033505nJ |
| Voc/放大器源功 | 3.989616533nJ | 0.038102906nJ |

这是显式合成的开路平衡→接入数学电阻事件。事件不改变电感电流、Cext和Cp电压或QT；
Icol可跳变，Cp电流即时补偿以保持i=0。源/气隙电压由代数源关系生成，
不将气隙初始500–800V等同Cext预充，也没有默认为真实QT0=0。
初始Cp板电荷和节点偏置显式记录；QT并不是TENG全部束缚表面电荷。
无充电历史时不假设所有导体节点中性，不声称合成初态来自实验。

每个接受步保存i、Icol、Cp位移流、Cext电流、泄放流、Vp、Vc、VT、QT及电功。
原放大器记录中VT为区间中点；汇总审计另列VT端点，两时间层不混用。
源坐标能QT²/(2CT)、L/Cext/Cp库存、柱功、Rs损耗及泄放损耗各自列出。
TOTAL组件Cp能与SPLIT外侧Cp能只计一次；同输入单次试探的两归属对照通过，
未额外运行TOTAL/SPLIT完整轨迹。放大器的TENG源坐标字段为null/NOT_APPLICABLE。

## 5. 独立核对与门限

33项隔离检查通过，包含拒绝/重复/外来提交、胶层、极性、来源哈希、单位、
未知初态、未匹配源、AUTO拒绝及冷/热兼容检查。没有调用旧历史测试入口。
独立连续系统矩阵指数采用同初态、不同方程实现，核对早期快模态及全窗有限节点；
矩阵指数是两条轨迹的解析参考评价，不是第三/第四条ODE轨迹。

| 指标 | TENG | 放大器 |
|---|---:|---:|
| 最大固定尺度状态误差 | 1.929984e-5 | 1.929962e-5 |
| 最大KVL绝对残差/V | 3.499423e-13 | 4.352074e-14 |
| 最大端口电荷残差/C | 1.956375e-28 | 1.262177e-29 |
| 最大单步能量残差/J（原账本） | 5.177957e-25 | 4.827040e-27 |
| 累计能量残差/J（原账本） | -3.308722e-24 | -2.121619e-24 |
| 离线独立能量残差/实际无抵消交换尺度 | 9.834836e-13 | 7.645861e-13 |

预注册状态目标1%，KVL/电荷/能量固定尺度1e-8，互斥归属1e-10；均满足。
尺度明确为V0、Cext*V0、Cext*V0²、V0/sqrt(L/Cext)，另保存绝对残差及实际交换尺度，
不用epsilon处理零分母。独立账本重算中点Icol²R、源功与各储能，不从残差构造热源。
能量只属于数学电路，W_column不能升级为气体热、RF、光或材料内能。

## 6. 统一入口及物理限制

真实运行入口仍是`I1MultiphaseScenarioRunner.run_mode`。
SOURCE_NETWORK_REFERENCE可执行上述两种源。冷FULL和给定热串联模式可读取同一装置来源包，
本轮只作preflight，没有重跑V3.5轨迹。候选几何不能冒充旧B0原生网格。
旧冷模式保持规定电压，旧热模式保持来源明确的1mm人工LTE+10V独立参考；
请求不具备连接的TENG/放大器时明确拒绝，不能静默替换源。
没有实现真实TENG-FULL或TENG-G1双向物理耦合；这是保留的能力缺口，不是软件标签可授予的资格。
AUTO=MISSING_PHYSICAL_CLOSURE，真实源标定NOT_VALIDATED，全局端口NOT_READY，RF未验证。

**下一项唯一优先工作：**形成同一H/Return参考面、同一预击穿事件的来源绑定源初态包：
核实实际高频连接与探头负载，给出xT及CT/Voc的测量关系、QT0或可换算充电历史、vCext0，
并核实介质/胶层路径和R1/DAG1单位定义。先让已有源/初态合同可被物理约束，
不为匹配频率调sigma、引线寄生或人工热阻。

## 7. 成本、复现与保护

完整动态与矩阵参考执行约4.52s；33项pytest进程约0.95s，离线账本处理不足1s；
精确时间、内存和后续静态封存检查记录见执行文件。总新增数值/检查成本远低于900s。
首次启动在归档Python可执行文件符号链接时失败，发生于任何动态步之前，
保留未执行源码快照与错误记录；修正为绝对路径+SHA后仅执行上述两条轨迹，无重复参考。
无编译、依赖安装、FULL、G1、传播或实验CSV调用。

复现单条模式（输出必须是新目录；本轮已停止，不自动重跑）：
`PYTHONPATH=python:. PYTHONDONTWRITEBYTECODE=1 .venv/bin/python -m adapters.i1_cold_physics_v3_5.scenario <input/teng_scenario.json或amplifier_scenario.json> <新输出目录>`。
完整一次性编排为`python -m adapters.i1_p1.execute`，已存在轨迹则拒绝自动重放。
输入文件绑定当前来源SHA；修改来源后须显式重新绑定，不能沿用旧身份。
执行版/最终版、冷库/驱动只读身份与保护结果见source_index及protection_after。
仅授权入口分派与两份长期台账修改；所有历史结果和物理内核保持。无Git提交、reset、clean或stash。
