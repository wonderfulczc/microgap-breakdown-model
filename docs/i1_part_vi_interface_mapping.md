# I1 / Part VI：现有接口只读映射与最小缺口

日期：2026-10-15；对应[合同 I1-PORT-1.0](i1_part_vi_port_circuit_energy_contract.md)。
`MAPPING_DOCUMENTED`仅代表本轮源码/文档清单，不代表接口实施、数值或物理验收。
未调用求解器、编译、拟合或历史测试；下列代码位置来自真实未提交 GENERIC_DEV。
完整文件身份与保护检查见[状态JSON](../validation/i1_part_vi/part_vi_theory_status.json)。

## 1. 身份与证据入口

分支 `architecture-v2.1-ultrafast-rf`，HEAD `8ec25920adff63a0b1a875c620b9a32105458ad7`。
开始时9个tracked文件已有修改，且有大量未跟踪研究文件；不是干净HEAD运行身份。
tracked binary diff SHA256为 `2fa08169f4b46eb4954079dfb97a46367fa2fc5beffe212ed3dc6f0cea26ac04`。
只读哈希核对到现有库 `build/libstreamer_rf_numerics.a`：
`d09e9d8f80d4c805d8b9509d159f57d1fc801c9af846edbb59faea1565a3438f`；
既有 `transfer_driver`：`70aa4933769f85b8b2a6082405944602fbd2bf55844219d808c5f10d5cfe8b34`。
这些只是当前文件/既有报告身份，不是本轮编译证明或新运行记录。

|证据代号|真实文件入口|
|---|---|
|C-H|[StreamerSolver.hpp](../cpp/include/streamer_rf/streamer/StreamerSolver.hpp)，StreamerConfig、AcceptedStepRecord、NativeFaceFlux、StreamerSolver|
|C-S|[StreamerSolver.cpp](../cpp/src/streamer/StreamerSolver.cpp)，fields、step、electrode_surface_diagnostics、conjugate_hv_charge_for|
|GEO|[geometry2d.hpp](../cpp/include/streamer_rf/geometry2d.hpp)，Geometry2DView、FaceRole2D；[electrode_profile.hpp](../cpp/include/streamer_rf/electrode_profile.hpp)|
|P|[poisson.hpp](../cpp/include/streamer_rf/poisson.hpp)、[poisson.cpp](../cpp/src/poisson.cpp)，assemble、electrode_operator_identity、FixedElectrodePoisson|
|V|[voltage_waveform.hpp](../cpp/include/streamer_rf/voltage_waveform.hpp)，VoltageWaveform::value、ConstantVoltage、SampledVoltage|
|OBS|[原生accepted-step合同](../adapters/comsol_fp2/FP2B_V5_EXPORT_CONTRACT.md)、[local_dynamic.cpp](../adapters/comsol_fp2/local_dynamic.cpp)|
|C3|[历史Stage C3说明](stage_c3_electrode_current_resistance.md)、[旧端口合同](i1_geometry_boundary_port_contract.md)|
|G-CIR|[circuit/model.py](../python/streamer_rf/circuit/model.py)，SeriesRLCGapCircuit；[旧电路内核文档](stage_g1_circuit_kernel.md)|
|G2|[thermal/coupling.py](../python/streamer_rf/thermal/coupling.py)，circuit_interval、ThermalCircuitCoupling；[G2说明](stage_g2_thermal_circuit_coupling.md)|
|G3|[thermal/port.py](../python/streamer_rf/thermal/port.py)，port_frame、PortTransientContract；[G3说明](stage_g3_transient_port.md)|
|G0|[JouleHandoff.cpp](../cpp/src/streamer/JouleHandoff.cpp)、[circuit/handoff.py](../python/streamer_rf/circuit/handoff.py)、[G0说明](stage_g0_cold_thermal_handoff.md)|
|G1|[thermal/radial.py](../python/streamer_rf/thermal/radial.py)，joule_heating、RadialThermalSolver、thermal_initialization_from_g0；[热G1说明](stage_g1_thermal_channel.md)|
|F|[rf/source/schema.py](../python/streamer_rf/rf/source/schema.py)、[adapters.py](../python/streamer_rf/rf/source/adapters.py)、[integrals.py](../python/streamer_rf/rf/source/integrals.py)|
|V5–V7|[V5-G](../validation/i1_part_v/V5_MODEL_SELECTION_FINAL.md)、[V6](../validation/i1_part_v6/centered_600v/V6_CENTERED_VOLTAGE_TRANSFER_REPORT.md)、[V7](../validation/i1_part_v/v7_model_comparison/V7_COLD_BASELINE_COMPARISON_REPORT.md)、[旧移交](../validation/i1_part_v/v7_model_comparison/PART_VI_HANDOFF.md)|

状态词分开：`NUMERICALLY_IMPLEMENTED_PHYSICAL_PENDING`为已有数值字段但物理解释待验；
`DERIVABLE_NOT_EXPORTED`为已有数据可推导但没有该直接字段（本轮不计算）；
`FIELD_MISSING`为缺字段；`FORMULA_NOT_IMPLEMENTED`为缺所需公式/调用链。
代码存在不等于其所有支路验证通过。

## 2. 场、端口与接受步字段映射

|字段/函数|理论物理量|单位|积分域/参考面|时间层|当前实现状态|证据|最小缺口|负责阶段|
|---|---|---|---|---|---|---|---|---|
|Geometry2DView::geometry_id/revision/reference_plane_id、classify/face|几何/材料/边界面身份|标识、m、m²|轴对称规则网格；gas/HV/return，Axis/Outer分开|固定|IMPLEMENTED_SUPPORTED_SUBSET|GEO|命名参考面非物理验收；需绑定全局端口/边界实例|B/C/VI|
|AxisymmetricGrid、AcceptedStepRecord::mesh_id|体积/网格|m、m³|旋转体积，导体不参与气体积分|固定|IMPLEMENTED|C-H、C-S|持续携带mesh_hash与坐标/单位声明|C/VI|
|StreamerConfig::poisson_boundary；electrode_operator_identity|边界/算子|边界势V或梯度V/m；标识|Neumann与OpenCharge不同；Return=0|固定配置，每次场求解|IMPLEMENTED_LOCAL|P|observer的operator_id是固定族名，缺直接携带完整boundary_id/算子指纹字段|C/VI|
|FixedElectrodePoisson::solve、StreamerSolver::fields|总势、总场、rho|V、V/m、C/m³|当前局部域、同geometry_view|每步更新；不是冻结静态场|IMPLEMENTED_VOLTAGE_DRIVEN|P、C-S|不是物理全局端口；固定epsilon0，非通用epsilon分布后端|C|
|record.weight_psi/er/ez；WeightCache|psi、Ew|1、m^-1|同配置1V零电荷问题|固定参考，接受记录附带|IMPLEMENTED|C-S observer组装|缓存已绑定operator key；需统一元数据出口，不能加COMSOL总场|C/VI|
|electrode_surface_diagnostics；diag.q_hv、record.q_native_hv_C|q_A，非已验Q_H|C|气体中心场求和，历史阶梯HV界面|当前采样；record为post|NUMERICALLY_IMPLEMENTED_PHYSICAL_PENDING|C-S、C3|非同算子共轭；不可唯一物理端口|C/VI|
|conjugate_hv_charge_for；record.q_conjugate_hv_C|q_B候选|C|epsilon0*A/h*(VH−phi_gas)，完整中心链接距离|当前；record为post|NUMERICALLY_IMPLEMENTED_PHYSICAL_PENDING|C-S|独立物理参考面/外域和误差仍缺；不改名Q_H|C/VI|
|vacuum_gap_capacitance；native_capacitance_F|历史C_A诊断|F|该函数**强制OpenCharge**的1V问题|固定缓存|IMPLEMENTED_DIFFERENT_REFERENCE_POSSIBLE|C-S|不能默认与本次Neumann C_B同一参考；保持旧字段不替换|C/VI|
|conjugate_vacuum_capacitance；conjugate_capacitance_F|C_B局部同算子候选|F|使用c_.poisson_boundary的1V问题|固定缓存|NUMERICALLY_IMPLEMENTED_PHYSICAL_PENDING|C-S|不是全局C_p；需capacitance_scope/owner|C/VI|
|NativeFaceFlux::gamma_e_m2s、area_m2、role|Gamma_e,out|m^-2 s^-1、m²|每气体单元向外，内面两侧反号；轴面零面积|实际pre源/输运更新|IMPLEMENTED_NATIVE|C-S capture_face、C-H|永久保留时间/邻接/法向；不是中心J|C|
|hv_absorbed_electrons、return_absorbed_electrons|Delta N吸收|粒子数|HV/Return气体侧实际吸收|接受区间dt积分|IMPLEMENTED_NATIVE|C-S step、OBS|Deltaq_Hg=e*HV计数可得；实际当前仅电子输运|C/VI|
|diag.i_cond_hv=e*hv_abs|I_Hg电子分量|A|H金属向气体常规电流为正|本步pre通量形成的区间率|NUMERICALLY_IMPLEMENTED_PHYSICAL_PENDING|C-S:370、C3|统一候选端子控制面；不补虚构离子运动项|C/VI|
|outer_lost_electrons、diag.outer_boundary_current|外边界电子损失及e*损失率|粒子数、A|整个外气体边界|接受区间|IMPLEMENTED_NOT_FO|C-S、OBS|Fo须psi加权且J·n=−e Gamma_out；现有正e损失率**不是**Fo|C/VI|
|Qo、dotQo、Po_in|外边界静电修正/电功|C、A、W|场域外边界，非粒子通量|端点/区间|FIELD_MISSING|C-H outer_electrostatic_correction_status|状态字符串不是数值积分；零Neumann截断理论Qo=0不认证真实外域|B/C/VI|
|diag.i_disp_hv、i_total_hv|历史Delta q_A/dt及其与I_Hg之和|A|历史中心场参考面|前次历史到当前post|NUMERICALLY_IMPLEMENTED_PHYSICAL_PENDING|C-S:252–258、370–374|已包含电压及空间电荷场变化，非仅Cp dotV；不再叠加完整Iw|C/VI|
|record.G_R_pre/post_S、W_rho_pre/post_C|FULL_GAS投影、权重电荷|S、C|全部局部气体；COMMON_GAP另名|真实floor后pre、post|IMPLEMENTED_INTERNAL_ONLY|C-H、C-S、OBS|Wrho非端子Q；未获独立W动态/全局端口准入|V/VI|
|delta_G_reaction/transport/boundary/mobility/correction_S|Delta G各机制|S，非直接S/s|原生面通量和权重差；局部气体|整个accepted步|IMPLEMENTED_DISCRETE_BUDGET|C-S observer、OBS|有限乘积M非全时序严格瞬时M；极短对齐步不能假精度|C/V|
|impact/photo/attachment2/3/recombination_m3s|原生电子源分量|m^-3 s^-1|气体单元|实际pre RHS|IMPLEMENTED|C-S step；ReactionModel由step调用|正/负离子由原生反应更新，非离子输运；完整反应净电荷审计不可漏物种|C|
|pre_floor/post_floor/clip的ne/np/nn delta数组|数值密度及电荷改变量|m^-3；积分后C|gas/mask更新|初始化floor与pre→post修正分开|IMPLEMENTED_NATIVE_ARRAYS|C-H；OBS local_dynamic|deltaQ=e*sumV(delta_np−delta_ne−delta_nn)；不从残差定义修正|C/VI|
|连续两个q_B端点、pre.phi/post.phi、HV计数|候选DeltaQ、integral I_Hg dt|C|同一局部离散端子|accepted区间|DERIVABLE_NOT_EXPORTED_AS_UNIFIED_PORT|C-H、C-S|record仅直接给post q_B；需显式配对pre端点/首步/拒绝事务，不能继承重启丢失历史|C/VI|
|Deltaq_port、port_energy|物理区间电荷/端口功|C、J|经验证的同一H端口|accepted区间|FIELD_MISSING|C-H对照新合同|候选数值电荷≠物理Q；缺统一port求积和物理准入|VI|

### 符号与不重复记账

当前电子进入H，Gamma_out>0，J·n_gas=−e Gamma_out；合同I_Hg=−integral J·n_gas为正，
与原生 `i_cond_hv=e*hv_abs`一致。等式是符号映射，不是物理参考面已一致的证明。
反应不自动含离子运动；原生dnp/dnn是源项，没有离子面输运。
VI-4用完整Q_H加一次吸收；VI-9用Cp dotV+Iw+Fo+dotQo，**不再加吸收**。
外边界电子计数、加权Fo、静电Qo、外边界功Po_in是四个不同对象。

## 3. 能量、电容、事务与能力映射

|字段/函数|理论物理量|单位|积分域/参考面|时间层|当前实现状态|证据|最小缺口|负责阶段|
|---|---|---|---|---|---|---|---|---|
|conductance_diagnostics::p_cond、gb、rb|历史electron drift功率及P/V²、倒数|W、S、ohm|气体sigma_e*E²|当前状态|IMPLEMENTED_NOT_FULL_PCOND|C-S、C3|不含全部扩散功；gb非GR，非物理I/V|C/VI|
|compute_joule_handoff_diagnostics::PJ_gas/channel|重构电子J·E功率|W|gas及sigma阈值定义的channel ROI|采样状态|IMPLEMENTED_DIAGNOSTIC|G0|保留正负功；非原生面量的精确离散端口能量，非自动热源|G0/VI|
|JouleHandoffAccumulator::QJ_gas/channel|上述采样功率梯形积分|J|同ROI，非完整端口|相邻采样区间|IMPLEMENTED_DIAGNOSTIC|G0|未闭合粒子/化学/热/场能分配；不能全赋H_c|G0/G1/VI|
|U_E、E_port、E_outer、DeltaU_sc|统一FULL场能/端口/外域预算|J|同端口场域及电容分区|端点/accepted区间|FORMULA_NOT_IMPLEMENTED_AS_UNIFIED_CONTRACT|C-H、G0对照|现有phi/E可供后续适配，不称已有经验证预算|C/VI|
|CircuitParameters::Cext_F/Cgap_F/cgap_source；SeriesRLCGapCircuit.rhs|历史外部并联Cext、局部gap支路|F|Vs→Rs→L→并联Cext/Cgap/Gb节点到ground|连续ODE/输出|IMPLEMENTED_DIFFERENT_TOPOLOGY|G-CIR|非新合同串联Cext；必须独立拓扑配置/审查，不能改名复用原方程|VI/G|
|evaluate_solution::cext_energy_J/cgap_energy_J；port_frame|旧电路储能及gap端口|J、V、A|I_port=I_L−I_Cext=I_sp+I_Cgap|旧电路输出|IMPLEMENTED_REFERENCE_ONLY|G-CIR、G3|不是原生C端口；H_GAP_CAPACITANCE_PARTITION_REQUIRED保留|G/H/VI|
|capacitance_owner、capacitance_scope|TOTAL/SPLIT互斥归属|枚举|local/global/source/external|固定实例|FIELD_MISSING|G-CIR、C-H|现有C数值/来源字段不足以阻止重复装配；先明确owner再实施|VI|
|VoltageWaveform::value、StreamerConfig::voltage_waveform|给定V_p|V|HV、Return0|原生pre/post时刻|IMPLEMENTED_PRESCRIBED_INPUT|V、C-S fields|不是未知Vp试探反馈；不能复制500V物性到新输入|C/VI|
|snapshot/evaluate_trial/commit/rollback、R_I|双向电路事务|状态token、A|PDE+电路同accepted状态|trial→accepted|FORMULA_NOT_IMPLEMENTED_FOR_STAGE_C|C-H/C-S|缺事务、候选Vp入口和残差wrapper；step会改变内部状态/历史|C/VI|
|save_checkpoint/load_checkpoint|已有重启数据|s、m^-3|time、ne/np/nn与grid维数；加载重求场|文件重启|IMPLEMENTED_INCOMPLETE_FOR_TRANSACTION|C-S:40–41|不保留完整源/场/电极电荷历史、head历史或外部累积；load会reset_electrode_history，非事务快照|C/VI|
|ThermalCircuitCoupling::step|热G2已有拒绝/共同commit|s、J|热流体+历史并联电路|trial copy→accepted|IMPLEMENTED_G2_ONLY|G2|可借鉴而非声称Stage C已有；不接FULL，不替代所需端口定义|G/VI|
|C1-L0固定率预测脚本、V5-G域说明|内部GR|S、s^-1|FULL_GAS、500V居中已测试点0–10ps|初值→输出|IMPLEMENTED_ANALYSIS_BASELINE|V5–V7|非可连接的物理电阻组件|V/VI|
|capability_flags、OUT_OF_VALIDATED_SCOPE、UNAVAILABLE_CAPABILITY、STATE_LIFTING_UNAVAILABLE|能力拒绝/有效域/lifting|枚举|未来统一ROM/FULL入口|调用前/切换|NOT_IMPLEMENTED_AS_RUNTIME_WRAPPER|在python/cpp目标入口检索未发现这些合同标识|仅文档限制；需最小wrapper，不能伪装已部署|VI|
|TENG F_T、q_T、mechanical input|源V-Q-x及能量|V、C、m、J|用户源与回流|连续/accepted|FORMULA_NOT_IMPLEMENTED_IN_INSPECTED_CIRCUIT_KERNEL|G-CIR|仅source_voltage给定函数，无本装置TENG参数/源能量账本|VI/源模型|
|joule_heating、RadialThermalSolver|热LTE导电、qJ、内能/动能/损耗|ohm、W/m³、J|柱对称热通道，指定长度|热accepted步|IMPLEMENTED_THERMAL_REFERENCE_ONLY|G1|i²Rsp身份不提供冷G_R→热源映射；LTE适用域与网格限制保留|G1|
|thermal_initialization_from_g0；g1_initial_state_contract|冷热初态/能量接口|状态、J、K、Pa|匹配物质/控制体|候选handoff|UNRESOLVED_PHYSICAL_HANDOFF|G0、G1|明确返回CALIBRATED_RADIAL_PROFILE_REQUIRED等；不能将冷标量初始化为可信热剖面|G0/G1/VI|
|H_c/M2低维闭合|热降阶状态|J、J·m²|固定核心，非全部冷能量|热动态|NOT_RESOLVED|[Part IV-R报告](../thermal/iv_dual_state_refinement/IV_REFINEMENT_FINAL_REPORT.md)|既有热FULL与G2能力不证明双状态自主闭合|G1/VI|

## 4. Stage F 的真实源与无线限制

|字段/函数|理论物理量|单位|积分域/参考面|时间层|当前实现状态|证据|最小缺口|负责阶段|
|---|---|---|---|---|---|---|---|---|
|electron_transport_current_source::jr/jz|电子传导J_RF=−e Gamma_e的对向面平均|A/m²|gas，金属为0；中心重构|调用时state，非直接接受步face数组|IMPLEMENTED_ELECTRON_SOURCE|C-S:157–217|与原生输运公式相容不等于该中心J可精确替代所有面账本|C/F|
|load_petsc_stage4_field_source_csv|rho和Jr_RF/Jz_RF空间源|C/m³、A/m²|2D轴对称来源，旋转体积元数据|field输出时间|IMPLEMENTED_TAGGED_SOURCE|F adapters|不是端口电流，不含虚构离子；独立端口到空间源映射缺失|F/VI|
|Afivo源两种source_definition|旧漂移J或连续性面通量重构J|A/m²|其3D源cells|导出时间|IMPLEMENTED_DISTINCT_LABELS|F adapters常量及load_afivo_stage_e_source_csv|漂移旧列未含扩散，不可与新面源混写|E/F|
|SourceMetadata/SourceRecord；current_moment|源身份、M=integral J dV|标识、A·m|声明的空间支撑域|series采样|IMPLEMENTED_SOURCE_SCHEMA|F schema、integrals|已有geometry/time/units，不具统一physical_port/capacitance_owner合同|F/H/VI|
|terminal_current→空间J/接收电压|端口到辐射/接收映射|A、A/m²、V|回路/电极/天线/接收链|未来耦合|NOT_VALIDATED|G3 FULL_WAVE_LOADING_FEEDBACK_NOT_COUPLED；旧移交|需参考面/外域/Maxwell源分区与加载反馈，防重复位移响应|F/H/VI|

## 5. 最小实施清单（建议位置，不授权编码）

1. **端口数据薄适配层**：复用C-H的AcceptedStepRecord和C-S既有q_A/q_B/吸收；补pre/post电荷配对、Deltaq_Hg、明确候选/物理资格、boundary/operator/build身份以及TOTAL/SPLIT互斥枚举。先不改ISG0/反应/Poisson。
2. **Stage C完整事务快照和试探Vp**：在C-H/C-S现有state/step周围增加所有可变历史、SP3 warm state和外部observer/预算的可恢复边界；试探不向正式observer累计，统一commit一次。当前文件checkpoint不足；拒绝步已有局部行为不证明完整事务。外层以VI-17残差联立，而非重写输运为电流边界。
3. **电路/储能归属**：复用G-CIR的参数和预算工具，但新增合同串联实例的独立装配/状态Vc，不覆盖旧并联benchmark；选TOTAL或SPLIT一次，区分C_local/global/Cext/source。已有G2 copy/commit可借鉴，不能直接声称可用于冷FULL。
4. **能量与外域出口**：复用phi/E/面通量与G0有符号功诊断，补同参考面的U_E、端口/外边界区间功与求积；未验证的粒子能/气体热保持不可用。局部→全局需要经审查边界响应及同步电荷/能量映射，不能用补偿电容省略。
5. **C1薄能力门**：仅封装已接受的内部预测及严格域标签，返回合同拒绝码；不新增模型、不变k、不把越界送入未经验证lifting。无需新通用框架。

这些是未来最少的工程触点，不是完成列表或新科学任务预算。物理端口未解决时，合成端口事务可另行授权审查，但不能宣称真实B0电路接入。

## 6. 保持的Gate与停止点

V7冷输出已完成，经典物理模型仍延期；V5-G受限工程接受、600V迁移受限、位置限制均原样保留。
本轮仅 `THEORY_DEFINED / MAPPING_DOCUMENTED`；`IMPLEMENTATION_NOT_AUTHORIZED / PHYSICAL_VALIDATION_PENDING`。
I1-PORT-03与BC-02仍阻塞全局物理端口；CIR-07仍NOT_RESOLVED；THERM-06与RF-08不升级。
未新增Issue，不补造缺字段，不重算历史误差，不把同算子q_B互易性当作真实端口证据。
完成链接/JSON/文档一致性与文件保护核查后停止，交网页版审查；不自动执行FP-3。
