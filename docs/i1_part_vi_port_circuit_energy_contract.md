# I1 / Part VI：ROM–FULL 统一端口、电路与能量接口理论合同

版本：**I1-PORT-1.0**；日期：2026-10-15。状态：`INTERFACE_THEORY_DEFINED`。
来源：用户本轮消息提供的完整合同正文。工作区未找到独立附件
`I1_PartVI_统一端口电路能量接口理论合同_V1.0.md`，因此文件来源标记为
`FULL_CONTRACT_ATTACHMENT_MISSING`；**消息中的完整正文可用**，本文件按其全部章节整理归档，非声称读取了不存在的附件。
本稿为网页版理论定义、连续恒等式及接口设计，**不是仓库已实施或已验证的结果**。
本轮只读源码并维护文档；不实施 FP-3，不计算新科学数值。
配套：[实际接口映射](i1_part_vi_interface_mapping.md)、[机器状态](../validation/i1_part_vi/part_vi_theory_status.json)、
[旧几何/边界/端口合同](i1_geometry_boundary_port_contract.md)、[问题台账](i1_open_issues_and_gates.md)。

## 0. 依据与本轮决定

[V7 原报告](../validation/i1_part_v/v7_model_comparison/V7_COLD_BASELINE_COMPARISON_REPORT.md)记录：离线执行 PASS；
冷输出对照 `COMPLETED_ON_AVAILABLE_DISCRETE_DATA`；经典模型
`DEFERRED_NO_COMPATIBLE_PORT_REFERENCE`；完整物理对照未完成。
500 V 五个居中单因素点的增长误差为：常量保持 100%，切线 3.227%–3.469%，C1-L0 0.493%–0.641%。
切线和指数共享初始有效率，不能据此证明指数独有的物理优势。
600 V C1 增长误差仍为 45.914630%；这些是**旧报告抄录**，未重新计算。

保留 [V5-G](../validation/i1_part_v/V5_MODEL_SELECTION_FINAL.md) 的 C1-L0
`ENGINEERING_MODEL_ACCEPTED`，仅内部 G_R、B0 局部域、500 V、已支持居中点、10 ps。
冻结 `k_star=9134636280.4 s^-1`，非时间平均率或普适气体常数。
C1 精确闭合 `NOT_PROVEN`；C2-Q `NOT_SELECTED`；C2-S `OBSERVATIONAL_CANDIDATE_ONLY`。
[V6](../validation/i1_part_v6/centered_600v/V6_CENTERED_VOLTAGE_TRANSFER_REPORT.md)保持
`FIXED_K_VOLTAGE_TRANSFER_LIMITED`。
[位置限制](i1_initial_carrier_position_uncertainty.md)保持
`KNOWN_SCOPE_LIMITATION / NO_FURTHER_COMPUTE_AUTHORIZED`，不声称偏移罕见或实验必然居中。
物理 R、R+V_f、RW、Toepler、Mayr、Cassie延期，非失败。
不重开位置/电压扫描、精确闭合搜索、加密或 B0 来源修复。

## 1. 物理对象、坐标与条件

\(\Omega\)为所选端口结构的场域，包括气体及需要的线性绝缘介质；
\(\Omega_g\)为移动载流子气体域。导体内部不参与气体积分。
边界为高压接触面 \(S_H\)、回流面 \(S_R\)、截断/其他外部面 \(S_o\)。
\(\mathbf n\)从场域向外；在电极处由气体/介质指向金属。
轴对称 \(dV=2\pi r\,dr\,dz\)，边界用真实旋转面积。
H/R 是永久导体身份，不随电压反号交换。
每次运行几何与标量介电率固定；移动电极、时变/各向异性介质、表面自由电荷须扩展合同。
合同允许线性介质**不代表当前 Stage C 支持介质区域**。

\[
V_p=V_H-V_R. \tag{VI-1}
\]

\(I_p\)为外电路进入 H 的常规电流，\(V_pI_p\)为流入组件的电功率。
仅在没有漏记第三端/环境位移/外边界通量时，单串联回路电流才可识别为 \(I_p\)。
\(V_R=0\)可以只是电势规范，不等于真实接地线；对地寄生、第三导体须单独节点/电容矩阵。
准静电局部端口假设不自动适用于 RF 引线感应或传播路径电压。

## 2. 主端口关系

令 \(\mathbf D=\epsilon\mathbf E\)，\(\mathbf J_{cond}\)只含实际实现的运动物种：

\[
Q_H=-\int_{S_H}\mathbf D\cdot\mathbf n\,dA,\tag{VI-2}
\]
\[
I_{Hg}=-\int_{S_H}\mathbf J_{cond}\cdot\mathbf n\,dA,\tag{VI-3}
\]
\[
I_p=\dot Q_H+I_{Hg}.\tag{VI-4}
\]

\(I_{Hg}\)以 H 金属向气体为正；电子进入 H 时贡献正的常规电流。
电荷/电流必须覆盖同一导体控制体与参考面，额外导线接触面须纳入或单列。
历史中心场 \(q_A\)、同算子链接候选 \(q_B\)、\(Q_{H,physical}\)保持不同身份，禁止静默替换。

## 3. 权重互易式与收集电流抵消

固定线性 \(\epsilon\)，外边界明确且与实际问题相容：

\[
\nabla\cdot(\epsilon\nabla\psi)=0,\quad\psi_H=1,\quad\psi_R=0,
\quad\mathbf E_w=-\nabla\psi.\tag{VI-5}
\]

\(\psi\)无量纲，\(\mathbf E_w\)单位 m^-1。
\(C_p=\int_{S_H}\epsilon\partial_n\psi\,dA\)单位 F；
\(W_\rho=\int_\Omega\psi\rho_q\,dV\)单位 C。仅气体带体电荷时可限于 \(\Omega_g\)；自由面电荷须补面项。
对 \(\nabla\cdot(\epsilon\nabla\phi)=-\rho_q\)，Green 恒等式给出
\(\int\psi\rho_qdV=\int_{\partial\Omega}\epsilon(\phi\partial_n\psi-\psi\partial_n\phi)dA\)，故

\[
Q_H=C_pV_p-W_\rho+Q_o,\tag{VI-6}
\]
\[
Q_o=\int_{S_o}\epsilon(\phi\partial_n\psi-\psi\partial_n\phi)\,dA.\tag{VI-7}
\]

\(Q_o\)是独立边界积分，不是拟合残差。相容齐次零 Neumann 可使它为零，
只对该截断模型成立，不认证全局电极。
在守恒反应系统 \(\partial_t\rho_q+\nabla\cdot\mathbf J_{cond}=0\)中，令
\(I_w=\int_{\Omega_g}\mathbf J_{cond}\cdot\mathbf E_wdV\)，
\(F_o=\int_{S_o}\psi\mathbf J_{cond}\cdot\mathbf n\,dA\)，则固定 \(\psi\) 下

\[
\dot W_\rho=I_{Hg}-F_o-I_w,\tag{VI-8}
\]
\[
I_p=C_p\dot V_p+I_w+F_o+\dot Q_o.\tag{VI-9}
\]

式 VI-9 **不能再加电极吸收电流**，它已在推导中抵消。
floor/clip 数值电荷及其电场响应须来自原生真实记录，独立评价，不能由残差倒推。
若物理反应净电荷不平衡，先核对物种是否漏记，不能把残差当真实源。

## 4. G_R 不等于完整端口电流

\(\boldsymbol\Gamma_e=-\mu_en_e\mathbf E-D_e\nabla n_e\)，
\(\mathbf J_e=-e\boldsymbol\Gamma_e\)。同边界分解 \(\mathbf E=V_p\mathbf E_w+\mathbf E_{sc}\)且无额外驱动时：

\[
G_R=\int_{\Omega_g}e\mu_en_e|\mathbf E_w|^2dV,
\qquad I_w=V_pG_R+I_{e,sc}+I_{e,diff}+I_{ion}.\tag{VI-10}
\]

\(I_{e,sc}=\int e\mu_en_e\mathbf E_{sc}\cdot\mathbf E_wdV\)，
\(I_{e,diff}=\int eD_e\nabla n_e\cdot\mathbf E_wdV\)。离子项仅对已实现输运求和。
独立外场须另加投影，不能藏在残差。各电流可有符号；\(G_R\ge0\)不保证 \(I_w\)与 \(V_p\)同号。
C1 只预测内部 G_R，缺少非欧姆、电荷与能量闭合；不能接成 \(R=1/G_R\)。
600 V 已知失败也不支持将冻结率用于未验证反馈电压轨迹。

## 5. 电容互斥归属

**TOTAL_CHARGE**：组件输出完整 Q_H 与 I_Hg，以 VI-4 装配；电路侧不能再并联同一 C_p。
**SPLIT_GEOMETRIC_CAPACITANCE**：固定 C_p，令 \(Q_{pl}=Q_H-C_pV_p\)，则

\[
I_p=C_p\dot V_p+\dot Q_{pl}+I_{Hg}.\tag{VI-11}
\]

也可采用 VI-9，但同一贡献只记一次。场能分账须与电流模式一致。
C_ext 是按网表连接的外部元件，仍有 \(C_{ext}V_C^2/2\)储能；
C_p 属于所选端口；C_local 是截断域诊断，非自动等于 C_p；源电容、对地/引线寄生归各自源或环境。
未经域分割和参考面匹配，禁止以 C_global−C_local 补偿。

## 6. 局部 FULL 到全局端口

对同一全局离散系统分块为
\(K_{ll}\phi_l+K_{lo}\phi_o=b_l(\rho)+B_lV\)，
\(K_{ol}\phi_l+K_{oo}\phi_o=b_o+B_oV\)。外域固定、线性、K_oo 可逆时：

\[
(K_{ll}-K_{lo}K_{oo}^{-1}K_{ol})\phi_l
=b_l-K_{lo}K_{oo}^{-1}b_o+(B_l-K_{lo}K_{oo}^{-1}B_o)V.\tag{VI-12}
\]

这是既有线性消元结构，非新物理。端子电荷和能量须同步消元，不能只改势场仍用局部 q_B。
未来可用经验证的等效边界算子、匹配域或同物理全域求解；本次不实现 Schur 补，不要求 COMSOL 全矩阵。
局部零 Neumann 保留原参考身份，不升级全局端口。

## 7. 电路双向耦合与事务

用户参考实例为单串联 TENG–R_s–L–C_ext–气隙，无额外地线；
R_s=3 ohm、L=100 nH、C_ext=101 pF、约800 V为用户参考参数，非本轮测量认证。

\[
L\dot i+R_si+V_C+V_p=V_T,\tag{VI-13}
\]
\[
C_{ext}\dot V_C=i,\tag{VI-14}
\]
\[
i=\dot Q_H(X,V_p)+I_{Hg}(X,V_p),\tag{VI-15}
\]
\[
\dot X=F_{FULL}(X,V_p;\theta),\quad\text{并解 Poisson 约束}.\tag{VI-16}
\]

TENG 保留源关系 \(F_T(V_T,q_T,\dot q_T,x_T,\dot x_T;\theta_T)=0\)，
无源侧分流时 \(\dot q_T=i\)。不能因约800 V就使用无限供电源；参数和机械输入另定。
初始 i、V_C、场和源电荷须符合预充历史，不默认101 pF预充800 V，也不由“调频”推出 V_C(0)=0。
现有电压驱动内核可保留，外层给试探 V_p，使

\[
R_I=i_{circuit}-[\dot Q_H+I_{Hg}]=0.\tag{VI-17}
\]

未来事务为 checkpoint/snapshot → evaluate_trial → report_residual → commit/rollback。
每次试探须起于同一 accepted 状态，不能重复永久推进。名称是设计，不代表代码已有。
通用网表可使用同一残差；上式串联实例不承诺自动验证任意拓扑。

## 8. 场能、电功与热

固定几何/epsilon，守恒连续方程下：
\(U_E=\tfrac12\int_\Omega\mathbf E\cdot\mathbf DdV\)，
\(P_{cond}=\int_{\Omega_g}\mathbf J_{cond}\cdot\mathbf EdV\)，
\(P_{o,in}=-\int_{S_o}\phi(\mathbf J_{cond}+\partial_t\mathbf D)\cdot\mathbf ndA\)。

\[
V_pI_p+P_{o,in}=\dot U_E+P_{cond}.\tag{VI-18}
\]

P_o,in 可有符号，开放域不能未经审计删去。粒子能、激发/电离化学能、气体内/动能、电极和辐射另记。
当 \(\phi=V_p\psi+\phi_{sc}\)、\(\phi_{sc}\)电极上为0、外边界使交叉项消失时，
\(U_E=C_pV_p^2/2+U_{sc}\)。进一步 \(F_o=Q_o=P_{o,in}=0\)时：

\[
V_pI_w=P_{cond}+\dot U_{sc}.\tag{VI-19}
\]

\(P_{cond}=P_{e,drift}+P_{diff}+P_{ion}\)，
\(P_{e,drift}=\int e\mu_en_e|\mathbf E|^2dV\)。
V_p I_w 不恒等于 V_p²G_R，更不恒等于气体加热。
外电路 \(W_{circ}=Li^2/2+C_{ext}V_C^2/2\)：

\[
\dot W_{circ}=V_Ti-R_si^2-V_pi,\tag{VI-20}
\]
\[
V_Ti+P_{o,in}=R_si^2+\frac{d}{dt}(W_{circ}+U_E)+P_{cond}.\tag{VI-21}
\]

SPLIT 若将 C_pV_p²/2 放电路侧，场侧只记剩余 U_sc。TENG机械功、内部储能和损耗不在 W_circ 中，另有源预算。

## 9. Accepted 区间合同

每个 [t_n,t_n1] 需交换下表，标明 pre/trial/accepted、求积、源码/库/网格/算子身份。

|字段|单位|要求|
|---|---|---|
|t_n,t_n1,dt|s|实际接受步，极短浮点对齐步不代表导数精度|
|V_p端点及区间表示|V|明确试探/实际输入|
|DeltaQ_H|C|同参考面真实端点差|
|Deltaq_Hg|C|原生带电面通量的区间积分|
|Deltaq_port|C|DeltaQ_H+Deltaq_Hg|
|E_port|J|integral V_p I_p dt，注明求积；一般不能端点V乘平均I冒充精确能量|
|DeltaU_E,E_cond,E_outer及可获得热/粒子预算|J|互斥储能分区，不补造缺失项|
|floor/clip及电荷/投影影响|m^-3、C、S|来自实际更新，非反推残差|

平均 I_p=Deltaq_port/dt。PDE/电路失败均回滚，不重复累计电荷或功。

## 10. ROM/FULL 能力与越界

当前 C1-L0 的**合同能力**为：

```json
{"predict_internal_GR":true,"provide_physical_terminal_current":false,
 "provide_gas_heating":false,"accept_validated_circuit_feedback":false,
 "provide_full_state_lifting":false}
```

合同要求越界返回 OUT_OF_VALIDATED_SCOPE（600 V、未知波形/反向过零、不支持初态、超10 ps），
请求电流/热返回 UNAVAILABLE_CAPABILITY，而非伪装0；这些拒绝是否已编码见映射，不能仅凭文档声称实现。
未来 ROM 需同时提供 \(\dot z=f(z,V_p;\theta)\)、Q_H、I_Hg、U_model和能量/适用域；
FULL 的 X 包含实际动态场变量，phi为约束，数值 q_B 不等于物理全局端口。
单个 G_R 无法唯一恢复 ne/np/nn。无全状态检查点、可信lifting或授权重放时，
ROM→FULL 要求 STATE_LIFTING_UNAVAILABLE；不能宣称免费自动切换，也不能每10 ps重置时钟续用。

## 11. 冷热、切换与 Stage F/H

冷段 Stage C/FULL负责输运/反应/空间电荷；现有工程窗口10 ps。
G0负责连接和粒子/电荷/能量映射，未验证条件不补经验阈值。
G1/G2负责各自热流体/电导/电路域；H_c/M2降阶仍未解决。
10–80 MHz对应12.5–100 ns，后期通道和更长衰减需独立源/回路证据，不能由10 ps指数外推。
切换保持回路电感电流、物理电容电荷及可追踪总能量；电容模型改变时不强迫Q/V同时不变。
H_c仅是特定控制体相对内能，非冷等离子体全部能量；不得把累计端口功全部设为H_c。

Stage F/H须标源类型 electron_conduction / total_conduction / terminal_current、参考面、支撑域和时间单位。
原生电子面电流、中心重构J、总传导J、端口I不能混用。
完整Maxwell已有介质位移响应时，未经审计不能再次注入 J_cond+partial_tD。
无线接收需天线/传播/接收链，端口预测不等于无线预测；无空间源映射的ROM仅可驱动明确的回路/辐射模型。

## 12. 最小接口身份

|类别|合同字段|
|---|---|
|实例|geometry_id,boundary_id,reference_plane_id,operator_id|
|维度|coordinate_system,dimension,units|
|能力|capability_flags,validity_domain|
|电容|capacitance_owner,capacitance_scope，TOTAL/SPLIT互斥、local/global分开|
|事务|accepted_state_token,trial_input,commit_token|
|端口|Q_H,I_Hg,Deltaq_port,port_energy|
|外域|Q_o,F_o,P_o_in或等价经验证算子|
|内部|G_R,W_rho,P_cond,U_E|
|预算|charge_budget,energy_budget,numerical_corrections|
|溯源|source_hash,mesh_hash,build_identity,time_layer|

## 13. 未来实施验收（本轮不执行）

后续另行授权的FP-3可依次审查合成理想R/C端口、TOTAL/SPLIT等价和浮置两节点、
试探回滚及零/反向有限性、C1能力拒绝，之后才是有物理前提的FULL反馈与经典模型比较。
合成电阻不得由当前G_R冒充。本稿不授予新科学预算或编译/执行许可。

## 14. 标签与相对旧合同的变更记录

|内容|身份/状态|
|---|---|
|连续电荷、Green及能量关系|EXISTING_PHYSICS / MATHEMATICAL_IDENTITY，须满足列明假设|
|q_B、接受步接口、外域消元适配|NUMERICAL_METHOD_ADAPTATION，工程/适用域尚需验证|
|冻结初始有效率|NEW_CLOSURE_HYPOTHESIS，保留旧受限工程接受|
|本轮统一接口|INTERFACE_THEORY_DEFINED，非IMPLEMENTED或实验验证|

2026-10-15：旧合同保留不改。新增完整符号/外边界/能量及事务定义；
旧G3并联Cext模型与本稿用户串联Cext**不是同一拓扑**，旧结果不重解释。
q_A/q_B/物理Q、电容互斥、局部/全局及不可lifting进一步细化，不追溯关闭Gate。
实施和数值验证 `IMPLEMENTATION_NOT_AUTHORIZED / PHYSICAL_VALIDATION_PENDING`。

## 15. 来源与核验层级

- R1：[V7原报告](../validation/i1_part_v/v7_model_comparison/V7_COLD_BASELINE_COMPARISON_REPORT.md)，本轮本地读取；数值仅抄录。
- R2：用户本轮提供的串联TENG–R–L–C_ext–气隙、浮置回流及初始储能说明，非新增测量。
- R3：[COMSOL电路连接说明](https://doc.comsol.com/6.3/doc/com.comsol.help.acdc/acdc_ug_modeling.06.39.html)。
- R4：[COMSOL External I-Terminal](https://doc.comsol.com/6.3/doc/com.comsol.help.acdc/acdc_ug_electrical_circuit.09.25.html)。
- R5：[COMSOL External I vs. U](https://doc.comsol.com/6.4/doc/com.comsol.help.semicond/semicond_ug_acdc.7.52.html)。
- R6：[Riegler、Windischhofer (2020)](https://arxiv.org/abs/2001.10592)，用户合同说明只核对摘要适用层级。
- R7：[Niu等 (2013)，DOI 10.1002/adma.201302808](https://advanced.onlinelibrary.wiley.com/doi/abs/10.1002/adma.201302808)，用户合同说明只核对摘要V-Q-x结构。

R3–R7为**合同作者提供的外部来源登记**，本轮未联网/未独立重审原文，未借用本装置参数。
VI-1–VI-21按本稿法向和假设归档，不称其为外文逐字公式，也不以文献替代项目离散验证。
本轮完成后停止，不自动启动FP-3或任何后续计算。
