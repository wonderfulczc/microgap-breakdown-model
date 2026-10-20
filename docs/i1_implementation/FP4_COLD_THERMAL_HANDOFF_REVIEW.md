# FP-4 冷—过渡—热通道交接的最小状态与能量接口审查

版本 I1-FP4-HANDOFF-1.0；2026-10-15。来源：本轮用户提供的网页版合同全文。
本文件是只读源码审查与接口设计，不是冷热转换实现、热求解或物理验收。

## 1. 分项决定与历史保护

|对象|当前决定|准确范围|
|---|---|---|
|交接合同/字段映射|THEORY_DEFINED / MAPPING_DOCUMENTED|本文件、CSV与声明式schema；无运行时提供器|
|AUTO_COLD_TO_THERMAL|MISSING_PHYSICAL_CLOSURE|材料能量分配、质量/动量/组分、LTE及剖面尚缺|
|GIVEN_THERMAL_INITIAL_STATE|READY_FOR_BOUNDED_IMPLEMENTATION|仅独立、来源明确的LTE参考；当前G1的T/p且显式u_r=0子集|
|一般非零速度/守恒热初态|INPUT_CONTRACT_DEFINED_WAITING_DATA|需合法数据及显式初始化接口；现有构造器不支持|
|真实B0热初态|INPUT_CONTRACT_DEFINED_WAITING_DATA|本轮没有给定真实热剖面或转换证据|
|FP-3局部数值接口|FULL_LOCAL_NUMERICAL_COUPLING_ACCEPTED|继承原报告范围，不重复测试，不推广全局|
|热降阶(H_c,M2)|NOT_RESOLVED|保持Part IV-R，未冻结、未拟合|

这里的READY只说明可以另行授权实施一个有明确输入的受限薄适配器，**不是本轮实施许可或物理PASS**。
独立热参考不依赖自动交接先通过，不以一个总NOT_RESOLVED阻断它；也不能把独立参考写成自动击穿预测。

实际分支 architecture-v2.1-ultrafast-rf，HEAD
8ec25920adff63a0b1a875c620b9a32105458ad7；已有9个tracked文件修改及未跟踪GENERIC_DEV/适配器。
完整当前文件SHA、tracked diff身份与原工作树见[保护清单](protection_before.json)；
旧原生库SHA为d09e9d8f80d4c805d8b9509d159f57d1fc801c9af846edbb59faea1565a3438f，
最新已运行驱动SHA为9eb23a121618fd1db7ce7acb5a14edd52a884db440badd1af48e1efd3c755ff9。
二者只是本次只读核实的已有文件，未编译或执行。
网页端给出的G1文档/radial Git blob SHA只是基线来源；当前文件另记SHA-256和实际Git blob身份，不能混同。

## 2. 真正读到的代码与结果

仓库路径均相对项目根；`路径:函数`表示源码定位，非已调用。

|证据|实际入口及用途|
|---|---|
|FP3|[弱量报告](../i1_fp3_port_circuit/weak_port_run_20261015/FP3_WEAK_PORT_ACCURACY_REPORT.md)、[Gate](../i1_fp3_port_circuit/weak_port_run_20261015/fp3_weak_port_gate_summary.json)、[来源](../i1_fp3_port_circuit/weak_port_run_20261015/source_index.json)；读取接受CSV字段，未重算指标|
|FULL|[StreamerSolver.hpp](../../cpp/include/streamer_rf/streamer/StreamerSolver.hpp):StreamerState/AcceptedStepRecord；[StreamerSolver.cpp](../../cpp/src/streamer/StreamerSolver.cpp):fields/step/electron_transport_current_source；[ReactionModel.cpp](../../cpp/src/streamer/ReactionModel.cpp):evaluate_reactions/electron_temperature|
|局部端口|[native_full.py](../i1_fp3_port_circuit/native_full.py):NativeFullProvider；[native_full_stable.hpp](../i1_fp3_port_circuit/native_full_stable.hpp):StableNetwork；[native_full_snapshot.hpp](../i1_fp3_port_circuit/native_full_snapshot.hpp):StreamerTransactionAccess；[native_coupling.py](../i1_fp3_port_circuit/native_coupling.py):IncrementalLocalFullCircuit|
|C2/C3/C-R2/C-R4|[C2](../../docs/stage_c2_dynamic_electrode_streamer.md)、[C3](../../docs/stage_c3_electrode_current_resistance.md)、[C-R2](../../docs/stage_cr2_reaction_source_decomposition.md)、[C-R4](../../docs/stage_cr4_joule_handoff.md)；[JouleHandoff.cpp](../../cpp/src/streamer/JouleHandoff.cpp):compute_joule_handoff_diagnostics/JouleHandoffAccumulator::sample|
|G0|[说明](../../docs/stage_g0_cold_thermal_handoff.md)、[handoff.py](../../python/streamer_rf/circuit/handoff.py):conductive_percolation/state_from_cr4_row/g1_initial_state_contract；真实结果在[circuit/g0_handoff](../../circuit/g0_handoff/g0_handoff_summary.json)，不是不存在的thermal/g0|
|G1|[文档](../../docs/stage_g1_thermal_channel.md)、[radial.py](../../python/streamer_rf/thermal/radial.py):RadialThermalSolver/thermal_initialization_from_g0/validate_calibrated_radial_state/required_thermal_energy；[properties.py](../../python/streamer_rf/thermal/properties.py):EquilibriumAirProperties|
|G2|[说明](../../docs/stage_g2_thermal_circuit_coupling.md)、[coupling.py](../../python/streamer_rf/thermal/coupling.py):ThermalCircuitCoupling/circuit_interval；[实际摘要](../../thermal/g2/g2_reference_summary.json)|
|G3|[说明](../../docs/stage_g3_transient_port.md)、[port.py](../../python/streamer_rf/thermal/port.py):port_frame/PortTransientContract；[实际摘要](../../thermal/g3_port/g3_port_summary.json)|
|热降阶|[Part IV-R](../../thermal/iv_dual_state_refinement/IV_REFINEMENT_FINAL_REPORT.md)：数值参考、联合覆盖、自主推广仍未解决|
|旧理论映射|[I1-PORT-1.0](../../docs/i1_part_vi_port_circuit_energy_contract.md)、[归档映射](../../docs/i1_part_vi_interface_mapping.md)：是此前时间点记录，后续FP3能力按最新报告叠加，不回改旧文件|

G1文档中的“G2尚未开始”、G2文档中的“没有G3”是历史时点说明：当前G2/G3源码和结果确实存在。
保留旧文档，不据它们把已完成工程判成未实现。
本轮所需核心文件已找到；缺失的是后文明确列出的物理状态/闭合数据，而非通过新计算补找文件。

## 3. 变量、单位与控制面合同

|本合同名|单位|不能混同|
|---|---|---|
|rho_q|C/m³|StreamerState.rho是净电荷，不是气体质量|
|rho_m|kg/m³|G1 props.rho_kg_m3，不等于ne+np+nn乘任意质量|
|ne/np/nn|m^-3|聚合带电物种计数，不提供完整中性组分|
|e_charge|C|元电荷，不是比内能|
|e_m / h_specific|J/kg|G1 e_J_kg / h_J_kg；不是总热能H_c或网格间距h|
|u_r|m/s|重粒子宏观径向速度，不是电子漂移速度|
|U_field / U_material / K_gas|J|分别为场储能、定义过能量零点的材料内部库存、宏观动能|
|W_cond / W_gas_heat|J|带电输运电功不自动等于重粒子热输入|
|i_port / I_col|A|总端口响应与热柱传导支路，不自动相等|
|C_p / C_local / C_global / C_ext|F|端口归属、局部诊断、全局结构与外部元件分别声明|
|ell_ch|m|所建模柱段长度，不自动取净气隙或计算域长度|
|G_length|S·m|截面积分，不是单位S的G_R或柱电导G_length/ell_ch|

每个量绑定geometry_id、boundary_id、operator_id、reference_plane_id、control_volume_id、
时间层和能量归属；轴对称dV=2pi*r*dr*dz，边界为真实旋转面积。
H/Return身份固定；Vp=VH−VR，i_port正向为外电路进入H，I_Hg为H向气体的常规电流。
几何、长度、气体、初态都作为参数实例，不将B0尺寸写入通用接口方程。

## 4. 完整交接理论合同

以下为网页版本轮理论结构；标签为MATHEMATICAL_IDENTITY/INTERFACE_THEORY_DEFINED，
**不代表仓库已实现材料能方程或已验证热转换**。

### 4.1 电路与场域电功账本

保留FP3实际原生面量、有效网络和修正项的符号：

```text
W_port + W_outer + W_correction + W_reaction_defect
    = Delta U_field + W_cond + R_E
```

correction/reaction_defect来自实际更新的数值项，不是气体热源；R_E是残差，不作为附加输入。
改变控制面正向时必须显式给出变换。物理连续极限、不计数值修正时：

```text
W_port + W_outer,in = Delta U_field + W_cond
```

此式没有给出电功在材料自由度之间的分配。
同一端口只选TOTAL_CHARGE或SPLIT_GEOMETRIC_CAPACITANCE之一；SPLIT场能扣除同一Cp项、
电路侧记一次，不能在R_map中两边各加一次。有限权重下remaining也保留交叉项/偏置。

### 4.2 材料能量账本

固定气体控制体、相同能量零点下，U_material由互不重叠的电子/重粒子平动、
激发、解离/电离等**实际建模**储能组成；宏观流体动能K_gas另记：

```text
Delta(U_material + K_gas)
    = W_cond + E_material,in - E_material,out
```

in/out包括同边界上的物质携能、表面交换和辐射等，按同一符号分类，不能重复计损失。
电子—重粒子内部交换在总材料账本中抵消，但在重粒子热预算中必须出现。
现冷模型无材料能/输运方程，缺电子能库存、激发/化学能库存与组分、重粒子能、动量、
电极热交换及辐射能账本，逐项为UNAVAILABLE，不能当0。
不能设任意热化效率eta来声称映射已验证，不能从一个总量恢复唯一T(r)、p(r)、u_r(r)或组分。
G1 LTE EOS的e_m=h_specific−p/rho_m使用其平衡热化学约定；
若已含平衡激发/解离/电离贡献，不再把同一化学能另外加进G1总能量。
冷阶段没有可直接匹配该EOS零点的材料能输出；不能仅凭ne/np/nn反应计数补出它。

### 4.3 热柱空间状态与导电支路

一般守恒初态：

```text
U_h(r,t_h) = [rho_m, rho_m*u_r, rho_m*(e_m+u_r^2/2)]
```

当前G1构造器由T,p经EOS得到rho_m/e_m，并明确构造U=[rho_m,0,rho_m*e_m]。
没有显式非零速度或守恒状态初始化入口。本轮不改构造器或将未知速度填零。
轴向均匀、柱状、所选LTE电导模型中：

```text
G_length = integral_A sigma(T,p) dA            [S*m]
R_col = ell_ch/G_length                       [ohm]
E_z = I_col/G_length                          [V/m]
P_J = ell_ch*integral_A sigma*E_z^2 dA
    = I_col^2*R_col                           [W]
```

这是G1柱模型关系，不是任意70um放电的严格定律。原参数current是I_col，
只有位移/感应/其他支路已分离并有相应条件时，才能近似取i_port=I_col。
如需近电极压降，先定义Vp=V_col+V_near_electrode，再分配柱体/电极功；
当前无该物理关系，不凭经验填V_f，不把全部Vp*I_col沉积进柱体。

### 4.4 零时长换模型的必要条件

相同物理控制体、参考面、真实电路元件与能量零点下，纯数值表示切换不得创造电荷、质量、能量。
无规定脉冲注入/交换时，真实电感电流连续、同一物理电容电荷连续；
被表示总电荷、质量、能量匹配；重要动量须匹配或明确其去向：

```text
R_map = [W_circuit + U_field + U_material + K_gas]_after
      - [W_circuit + U_field + U_material + K_gas]_before
```

纯表示变更的R_map应在明确投影/数值误差内；电极、源、外域库存若在范围内，两侧同范围补入。
没有冷材料能时R_map=UNAVAILABLE，不是0；不得反推“未知加热”或拟合损失。
电容表示改变时不能强制同时保留不相容的Q、V和能量。
局部到全局须同时映射端口、边界、储能，禁止C_global−C_local补差。
热G1本身不演化rho_q；冷端非零净电荷/场能若未由明确组件继续承载，必须拒绝自动切换，不能丢弃。
标量守恒是必要条件，不能证明径向剖面唯一、轴向均匀或LTE成立。

## 5. FULL、G0、G1/G2的逐字段结论

完整[字段能力CSV](cold_to_thermal_field_map.csv)使用六种规定状态，
把“原生已有”“已有数据可推但本轮未算”“物理闭合”“外部初态”“未实现”“不适用”分开。

**FULL可提供的是真实数值电荷/电功，而非气体热。**
ne/np/nn/rho_q、原生电子面通量、物种floor/clip、q_A和两种q_B保留；
stable_V2给同局部活跃网络的感应电荷、Deltaq_port、Deltaq_Hg及其预算。
平均局部i_port可由Deltaq_port/dt得到，不等于I_col。
U_field、稳定remaining/DeltaU、原生/外层端口功、导电面功及无抵消功、外域和数值修正分别可追溯。
remaining含有限权重偏置，负值不是“可用负热能”。这些只在匹配局部控制面成立，非全局Q_H。
本轮只引用最新FP3 Gate中的已记录数字：电荷残差/Qtr约1.59e-15，
电能残差/有符号导电功约3.50e-7，实际提交端点KVL约3.17e-12V，未重新评价这些指标。

**G0不是热化器。** 现有g1_initial_state_contract携带候选时间、
channel_length/effective_radius/volume、sigma_eff、PJ/QJ、E/ne统计、Gb/Rb；
thermal_profile_status始终NOT_AVAILABLE。state_valid=true只表示冷标量对象可用，不表示热初态有效。
Q_required/Pi_H/阈值无生产标定，缺校准时HANDOFF_CALIBRATION_PENDING。
conductive_percolation需要显式gas/HV/ground邻接mask；ground在这个旧接口表示Return材料，不证明物理接地。
当前G0历史结果是旧针板、SP3 ON算例，不与最新双针SP3 OFF的60fs记录拼成同一事件。
摘要实际为PRODUCTION_HANDOFF_NOT_REACHED、候选state_valid=false；latest冷对象true但仍无热剖面。
sigma相对阈值/轴线ne阈值只是数值连通诊断，不是LTE或热火花判据。

**C-R4的Joule名字不证明加热。** 它用原电子面通量的对面平均中心J与E做J·E，
采样功率再梯形累计QJ；可有负贡献，ROI按sigma阈值随状态变动。
因此它既非固定物质库存，也不等于FP3同网络逐面能量账本；两者不可无控制面/求积审查直接互换。
C3的p_cond/Gb=P_drift/V²又是另一历史诊断，不能替代G_R或完整传导功。
ReactionModel的T_e=e_charge*D/(mu*k_B)只是重组率使用的代数有效温度，未求电子能方程，
不能把它当已演化电子能量库存或证明电子—气体温度相等。

**G1/G2真的需要的输入。** G1至少需RadialGrid、ell_ch、T(r)、p(r)、properties、
ambient/closed边界、ambient_T/p及current(t)的物理支路/支持范围。
LTE拟合输入域为50–60000K、0.01–100atm，求解器另要求初始e_m>0；这是物性定义域，不是LTE已成立证据。
rho_m与e_m由同EOS导出，e_m=h−p/rho_m；质量和能量须绑定同径向体积及长度。
边界有原生导热及流体能/质量交换；ambient ghost不是通用无反射边界。
辐射关闭，电极/轴向损失、黏性、磁力、非平衡与双温没有实现；诊断中radiation_loss=0
是“未启用”，不能作真实辐射为0的测量。

**current的实际身份。** G1 current(t)进入I²R，排除位移电流。
G2原拓扑是Vs→Rs→L→并联(Cext,Cgap,Gsp)，不是FP3的串联Cext；
其I_A=V_channel/Rsp，热步使用sqrt(Qsp/(Rsp_n*dt))这个导电支路RMS，非带符号平均电流、I_L或总端口。
G3的I_port=I_L−I_Cext=I_sp+I_Cgap保持历史身份，不改成原生FULL i_port。
今后接用户串联拓扑需薄适配与新范围验收，不能重解释原G2结果。

## 6. 构造器、能量零点和历史丢失风险

RadialThermalSolver.__init__没有u_r/U_h/time/history参数；不支持的关键字会拒绝，
若未来适配器只调用T/p入口而忽略输入速度，就会丢失动量和动能，因此合同要求提前拒绝。
公共self.U虽可被Python赋值，**不等于受验证的守恒初始化入口**：必须同步T/p反演、
initial_mass/energy、时间及预算，不能由外部直接改数组冒充完整初始化。
validate_calibrated_radial_state只检查T/p维数、来源及EOS有效性，不检查u/电荷/历史/映射守恒或LTE实证。
required_thermal_energy要求调用方给定质量、cold_e、target_e，不自行选择温度，不含未知损失或膨胀功。

构造器把time、steps/rejected_steps及QJ/Qconduction/Qhydro/mass_out归零，以给定热态建立新初值库存。
它不自动保留FULL的t_h、端口Q/电压历史、回路i/Vc、累计源功/修正或材料能零点对应关系。
G2构造器还要求thermal.time=0且steps=0，重建热/电路初始预算；已有浅copy共同接受依赖
G1替换数组的实现细节，不是跨冷/热模型的全状态事务。
未来必须保留全局物理时间t=t_h+tau，tau只是热求解局部时钟；旧累计账本与新段增量分开，
不能用归零累计量掩盖已有电功，也不能假定G1的EOS零点与冷能量已匹配。

## 7. 两入口的最小依赖与下一次最小模块

### A. AUTO_COLD_TO_THERMAL：MISSING_PHYSICAL_CLOSURE

明确缺四条关系，不是笼统“需要更多验证”：
1. 传导电功到电子/重粒子/激发化学储能、边界和电极交换的空间时间分配及守恒方程；
   当前无dU_e/dt、电子—重粒子交换或完整化学/中性组分能量。
2. 冷非平衡材料库存到同EOS零点下LTE质量/内能/组分的映射与适用性依据；不允许任意eta。
3. 二维轴向非均匀冷场到一维柱体rho_m/u_r/e_m剖面及ell_ch的有来源投影，含动量、非零净电荷/场能去向。
4. 同参考面、同外域与相同物理电容/回路状态的同时映射，包括近电极压降和损耗归属。

现有60fs电荷/电功记录不足以补出上述自由度；不构造高温通道，不以ne或固定时刻触发切换。
G0连接信息可以保留，但单独不能证明热化或确定T/p/u。

### B. GIVEN_THERMAL_INITIAL_STATE：受限独立参考

已有G1/G2合成参考清楚给定Rmax=0.5mm、ell_ch=1mm、
T=300+9700*exp(-(r/50um)^2)、p=101325*(1+exp(-(r/50um)^2))、u=0。
它保持NUMERICAL_REFERENCE_INITIAL_CONDITION及LTE_APPLICABILITY_PENDING_CALIBRATION，
**不把1mm改成70um后声称微间隙已标定**。
G1原参考有给定导电电流；G2历史参考有自己的并联电路与明确初值，均可作为独立受限参考证据。
新给定参考需提供来源、T/p/u与空间采样约定、ell_ch/网格、质量与能量基准、
EOS/组分适用性、边界、导电支路电流支持范围；不能让适配器自动给定这些物理量。
任意非零u或守恒输入需显式初始化接口；不兼容则返回NOT_IMPLEMENTED，未知数据保持null+原因。
该入口不是从冷端来，R_map冷→热不适用；仍应独立检查其自身质量/热能初值与EOS一致性。

### 下一次可实施的唯一最小模块（本轮不编码、不运行）

**GivenThermalInitialStateAdapter，先只承接当前T/p+显式u=0的独立参考。**
输入：本目录schema的GIVEN_THERMAL_INITIAL_STATE包、来源/能量零点、
径向单元/采样约定、T/p/u、ell_ch、EOS、边界、明确I_col驱动及独立起始时间；
若请求电路连接，还需确切拓扑、支路映射、状态与电容归属，不默认复用FP3串联到G2并联。
输出：带身份的初始化请求、同EOS导出的U_h/质量/内能/动能与检查报告、
I_col身份、全局时间锚与分段预算标识；输入不满足时结构化拒绝，不输出伪造热态。
通过后才允许由另行授权的运行器调用原G1；本模块不负责推断热化、自动触发或拟合。

最小触点为新增独立薄适配器，复用radial的网格/EOS校验；原冷FULL提供器/stable_V2和G0均不用改。
如后续确需非零速度/守恒初态，再审查一个显式G1初始化入口，携带U_h并同步T/p与初始账本；
不得本轮直接写self.U。热阶段电路trial/commit与冷热瞬时映射不是同一个验收，不能相互继承。
本轮没有授权下一次实施或科学计算预算，清单仅交科学审查。

## 8. 归档与停止

[readiness](fp4_handoff_readiness.json)逐项保存状态；[schema](cold_thermal_handoff_schema.json)
只是声明式输入合同，不是运行时验证器，跨数组长度、EOS反演、积分守恒及LTE证据须语义检查。
schema可表达未知量，未知不填0；schema可表达未来保守态不代表现有G1支持它。
本地未安装jsonschema元模式验证包，未安装新依赖；仅检查JSON语法、声明结构与本地$ref解析，
不把这称为完整Draft2020-12元模式认证或任何运行时/科学测试。
路线图与I1-THERM-06/CIR-07/PORT-03等既有Issue仅补证据，不增加重复Issue、不改历史Gate。

仅执行Markdown/JSON/schema/CSV结构、链接、文件身份及git diff静态检查；
本轮PDE、Poisson、StreamerSolver、G1/G2、RLC、TENG、RF调用、编译、拟合及科学测试均为0。
保护记录与源码/已有结果身份见source_index.json及protection_after.json；
历史事故中原smoke JSON未知字节仍UNKNOWN，不恢复猜测，不执行旧CTest。
不提交、不reset/clean/stash；完成此审查后停止，不追加诊断轮次或科学计算。
