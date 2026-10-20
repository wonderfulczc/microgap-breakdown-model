# FP-4 给定热初态薄适配器与原生 G1 等价验证

2026-10-15；I1-FP4-GIVEN-THERMAL-1.0。
最终决定：**GIVEN_THERMAL_ADAPTER_ENGINEERING_ACCEPTED**。
范围仅独立、来源明确的LTE数值参考、原生径向单元中心T/p、显式初始u_r=0。
不是自动冷→热、真实B0、LTE实验、全局端口或热降阶验收。
本次新授权覆盖实施与两条原生回放；旧FP4只读合同及其NOT_IMPLEMENTED历史记录保持原样。

## 1. 实际输入与版本

分支architecture-v2.1-ultrafast-rf，HEAD
8ec25920adff63a0b1a875c620b9a32105458ad7；实际仍为未提交GENERIC_DEV及适配器。
原库、冷FULL/stable_V2、G0/G1/G2/物性、历史测试与所有旧数据只读。
真正调用的是安装在当前工作区的streamer_rf.thermal.radial.RadialThermalSolver，
没有复制PDE/EOS或覆写U、时间、初始库存和累计账本。
保护、模块实际路径、源码/系数/输入/结果SHA见source_index.json和protection_after.json。

原始来源为[原G1生成脚本](../../../thermal/g1/generate_g1_reference.py)及
[真实归档摘要](../../../thermal/g1/g1_reference_summary.json)，没有执行该脚本的写回入口。
原生源文件SHA与该历史摘要中的source_sha256一致后才继续。
从脚本AST只提取网格、初态、电流表达式；直接路径采用其原lambda，
适配器采用等价的显式参数信号，不以测试轨迹拟合或猜测波形。

|量|实际主参考|
|---|---|
|径向网格|150单元；Rmax=0.5mm；原生单元中心采样|
|柱长|1mm，不是70um，也未作长度缩放|
|T(r,0)|300+9700 exp(-(r/50um)^2) K|
|p(r,0)|101325*(1+exp(-(r/50um)^2)) Pa|
|u_r(r,0)|显式全零数组，仅为初态|
|I_col(tau)|0.05*(1+0.5*sin(pi*tau/10ns)^2) A|
|原生选项|ambient；300K/101325Pa；CFL=0.3；conduction=True|
|时间|tau从0到10ns；按原脚本在5ns及10ns对齐|
|来源标签|NUMERICAL_REFERENCE_INITIAL_CONDITION / LTE_ASSUMED_FOR_REFERENCE / NOT_B0_THERMAL_CALIBRATION|

冻结输入见[input/handoff_package.json](input/handoff_package.json)、
[profile](input/profile.json)、[current](input/current.json)、[runtime config](input/runtime_config.json)。
通用适配器没有写死这些尺寸、温度或电流幅值；它们只属于这一个有来源的参考包。

## 2. 实际实现及语义边界

[given_thermal.py](../given_thermal.py):GivenThermalInitialStateAdapter沿用上一轮schema字段；
运行参数和信号分别放入其既有source_manifest_uri/current.source_uri指向的带哈希sidecar，
未改历史schema，也未另造不兼容主字段。
construct调用原构造器，返回原生solver、独立EOS库存、语义审计、全局时间锚和segment_id。
可运行入口是[run_given_thermal.py](../run_given_thermal.py)，不会调用G2/G3或电路求解。

实际检查来源文件及当前SHA、EOS/系数/材料/能量参考、原生网格面、数组维度与有限性、
原生EOS范围及e_m>0、边界、柱长、信号支路/单位/支持范围、时间基准和独立新段身份。
只支持原生CELL_CENTER，不插值、平滑、外推或作未经定义的守恒重映射。
缺u_r为UNKNOWN_INITIAL_VELOCITY；任意可表示非零值（测试含1e-300m/s）明确拒绝。
字符串/布尔值不强制转成速度0；JSON非零数下溢也拒绝，不把1e-400解析成合法显式零。

结构化拒绝见[rejection_audit.json](rejection_audit.json)、isolated_tests.xml与input_edge_tests.xml：

|请求|拒绝身份|
|---|---|
|仅G_R/W_cond/W_port/ne/np/nn或AUTO|MISSING_PHYSICAL_CLOSURE|
|非零初速、守恒态入口|UNSUPPORTED_INITIAL_STATE_MODE|
|缺初速|UNKNOWN_INITIAL_VELOCITY|
|i_port/I_L/总电流当柱体电流|CURRENT_BRANCH_MISMATCH|
|缺来源/哈希不符|SOURCE_UNAVAILABLE / SOURCE_HASH_MISMATCH|
|缺时间基准/支持不足|REQUIRED_FIELD_MISSING（明确缺项）/ SIGNAL_SUPPORT_INSUFFICIENT|
|非法网格/维数/有限性/长度|INVALID_GRID / GRID_MAPPING_UNSUPPORTED / INVALID_PROFILE / INVALID_NUMERIC_INPUT|
|原生EOS域外或e_m<=0|INITIAL_PROPERTY_INVALID / INITIAL_ENERGY_NONPOSITIVE|
|冷库存继承/检查点恢复|COLD_HISTORY_MAPPING_UNAVAILABLE / RESUME_NOT_IMPLEMENTED|
|旧G2并联冒称新串联/真实B0科学声明|THERMAL_SERIES_COUPLING_NOT_IMPLEMENTED / UNSUPPORTED_PHYSICAL_CLAIM|

没有jsonschema依赖，不安装新包、不宣称Draft元模式认证；上述支持子集的语义检查真实执行，
不是只验证JSON可解析。测试用新的隔离输入，不改历史数据。
非正EOS返回的拒绝单测是显式测试替身，不称为新的物理初态或原生EOS失效案例。

## 3. 独立初始库存与单位

从输入及同一EOS独立计算，不读取构造器库存再回填预期值：

```text
A_i=pi*(r_right^2-r_left^2) [m²] = native RadialGrid.volumes
V_i=ell_ch*A_i [m³]         (不再乘2*pi*r)
rho_m=rho(T,p); e_m=h_specific-p/rho_m
mass0=sum rho_m*V_i
material_energy0=sum rho_m*e_m*V_i
kinetic_energy0=0，仅因显式初始u_r=0
```

|库存|独立EOS求和|原生构造器|相对差|
|---|---:|---:|---:|
|质量/kg|8.877255034826210e-10|8.877255034826212e-10|1.165e-16|
|材料内能/J|2.0825416216199168e-4|2.0825416216199174e-4|2.604e-16|
|气体动能/J|0|0|EXACT_ZERO，不定义0/0|

差别来自sum与原生dot的运算次序，低于1e-12门限；原值见initial_inventory_audit.csv。
另有initial_column_electrical_audit.csv独立核对G_length/G_col/R_col及体积分PJ。
初始G_length约4.14789481307e-6 S*m，R_col约241.0861521484ohm，PJ约0.602715380371W。
公式沿用G_length=sum(sigma*A)、R_col=ell/G_length、Ez=I_col/G_length、
PJ=sum(sigma*Ez²*V)=I_col²R_col。没有引入epsilon电导，原生零电导/零电流行为不修改。

能量ID为NATIVE_DANGOLA_FIT_H_MINUS_P_OVER_RHO加实际系数SHA，
表示这个版本的原生EOS参考约定，不是已匹配冷阶段的绝对能量零点。
没有另加已包含在平衡EOS中的化学/激发能；材料库存不等于冷W_cond、W_port或H_c。
R_map=null，状态NOT_APPLICABLE_INDEPENDENT_REFERENCE，不生成冷→热误差0。

## 4. 时间、当前支路及新段预算

tau是原G1局部时间，t_global=t_anchor+tau。主参考t_anchor显式为0，
未从冷FULL的60fs末端推断任何交接时刻。
prior_cold_history=null+NOT_APPLICABLE。原生steps/rejected_steps、QJ/Qconduction/Qhydro/mass_out
从新段开始，绝不声称继承冷历史。

信号支持THERMAL_LOCAL；GLOBAL明确在回调中取I_global(t_anchor+tau)，
两者均检查完整计划时间域。不自动推断基准或外推。
人工2us非零锚测试仅作函数取值/时间输出检查，标SYNTHETIC_CLOCK_TEST，
未新增第三条热轨迹或任何物理交接含义。
current始终是COLUMN_CONDUCTION_CURRENT；无i_port/I_L反演、无C1或RLC调用。

## 5. 两条真正原生回放

初始化先核对T/p/U、时间和库存、G/R/PJ；输入内容不变，两个实例不共享可变T/p/U数组。
再按相同原生dt限制、相同5ns/10ns对齐逐步推进，无插值、真值回灌或重初始化。
两条各完成103接受步、0拒绝，共206接受步，终点10ns。
每个接受时间、T/p/U各分量、实际u_r、G/R/PJ、累计边界/Joule/质量账本均逐值一致。
见thermal_replay_equivalence.csv和thermal_mass_energy_budget.csv。
最终原生最大向外速度169.93798040774922m/s，**没有每步清零速度**。
原生场快照保存在5ns/10ns文件，最终剖面另存；没有每步保存全场。

在源码SHA和本参考配置核对相符的前提下，归档主参考末值也逐值复现：
T_axis=7870.204792238168K、R_col=620.1303486974772ohm、
原生total_energy=2.0827145299770567e-4J、QJ=1.7290835713937548e-8J。
这是封装等价与该历史数值回归，不是独立PDE物理认证或新网格收敛证明。

## 6. 质量、材料能和Joule账本

独立由当前U和原生实际累计项计算：

```text
R_mass=M(t)-M(0)+mass_out
R_energy=E_total(t)-E_total(0)-QJ+Qconduction_out+Qhydro_out
E_total=sum U[:,2]*V_i
```

共同库存/交换尺度固定为max(abs(E0),abs(E),abs(QJ)+abs(Qconduction)+abs(Qhydro))。
有符号边界量保持原号；绝对值只用于诊断尺度，不用于改写守恒式。

|最大诊断|实际|预注册目标|
|---|---:|---:|
|初始独立库存相对差|2.604e-16|1e-12|
|两路径同字段差|0，逐值一致|1e-12|
|abs(R_mass)/M0|1.167e-16|1e-10|
|abs(R_energy)/库存交换尺度|3.137e-16|1e-10|
|abs(R_energy)/abs(QJ)，QJ非零时|6.267e-10|1e-6|
|PJ体积分与I_col²R_col相对差|5.465e-16|1e-10|

独立求和的最大绝对R_mass约1.036e-25kg，R_energy约6.532e-20J；
与原生dot运算的最后几位差别保留，不以巨大背景内能掩盖新增Joule尺度。
QJ=0的初始行该比值留空NOT_EVALUABLE，不加epsilon。
radiation_loss=0同时标MODEL_DISABLED，非实测辐射为0；
轴向/电极损失、近电极压降、磁力、非平衡等保持UNAVAILABLE/NOT_MODELED。
没有从残差反推热源或损失；本LTE柱内的I²R不是从冷FULL获得的全部气体热输入。

## 7. 资源、源身份与保护

本轮只执行新目录测试与两条参考；原生回放函数计时3.3713s，首次24项隔离测试报告0.88s。
封存时补强输入拒绝（布尔/字符串、JSON非零下溢、畸形容器），再执行6项纯输入检查，
没有再次构造/推进完整轨迹。30项均通过；该最后检查的合法参考只做validate，不调用step。
执行时源码在executed_source，封存后源码另有SHA；后修订只收紧非法输入，
接受的实值数组、原生构造及推进路径不变，不冒称修订后重跑过PDE。
解释器导入/prepare未单独计时，不把函数计时称为完整OS进程记账；成本见cost_record.json。
原生执行命令受timeout580s限制，总预算600s；实际206接受步远低于2000上限。

源/输入/输出身份、可重现命令、增量diff与保护记录均在本目录。
旧原生/热源、G0/G1/G2/G3结果、FP3及FP4历史Gate未变；tracked GENERIC_DEV diff保持。
未安装新依赖、未编译、未跑全仓/历史测试、未提交或执行reset/clean/stash。

## 8. 分层Gate与停止

本轮八项Gate均PASS，见fp4_given_thermal_adapter_gate_summary.json。
AUTO仍MISSING_PHYSICAL_CLOSURE；非零速度/守恒初始化仍NOT_IMPLEMENTED；
真实B0热态WAITING_DATA、LTE实验NOT_VALIDATED、H_c/M2闭合NOT_RESOLVED；
串联热电路NOT_TESTED_THIS_TASK，全局端口/TENG/RF/实验及C1/V5/V6/V7均不升级。

下一步只有[串联连接移交清单](THERMAL_SERIES_COUPLING_HANDOFF.md)供科学审查，
不自动实施热RLC、冷转热或任何后续计算。给定热态适配器本限定工程验收到此结案。
