# Toepler / Rompe–Weizel 形成阶段受限动态端口

2026-10-16；PHYSICAL_SPARK_RESISTANCE_FORMATION-1.0。数学和受限工程，不是本装置物理电阻/RF验收。

## 1. 四项正交判断

- FORMATIVE_LAW_MATHEMATICS：PASS，单位、历史积分、负流数学行为及域拒绝已检验。
- COUPLED_PORT_ENGINEERING：PASS，仅下述两个合成1V/5ns事件。
- FORMATION_ONLY_PHYSICAL_SCOPE：FORMATIVE_UNIPOLAR_ONLY；不含自然恢复、熄弧或再击穿。
- DEVICE_PARAMETER_IDENTIFIABILITY：NOT_CALIBRATED。

PHYSICAL_SPARK_RESISTANCE_VALIDATED=false、RF_VALIDATED=false。数学/事务/能量通过不升级实际器件。
前轮七条100ns、S500/R20中止及执行PARTIAL、800V20/200Ω首峰和代理比6.338973全部保持，未重跑补齐。
不把旧源接口、FP3/FP4工程或冷C1的内部GR改名为已标定物理电阻。

## 2. 原文可访问边界

[文献合同](executed_literature_contract.md)逐项记访问层级。
[2018](https://doi.org/10.1541/ieejfms.138.316)只重新读取出版摘要，
500/4000V平行盘结果不移植到双针TENG；[2010摘要](https://doi.org/10.1587/transcom.E93.B.1801)
与2018对既往人体放电的概括差异保留，不宣布500/800V普适排序。
[2025正文§3.3–3.4式2–3](https://doi.org/10.1002/tee.24152)可核对积分式和单位，
不是原始历史推导全文。[2026出版DOI](https://doi.org/10.1016/j.elstat.2026.104253)直接打开仍失败；
本轮找到[作者稿v3](https://arxiv.org/html/2509.05791)，阅读方法2/理论3等可访问段落，
核对毫米级间隙及电容充电/串联负载/探头支路不同于本拓扑，未声称已逐字核对最终出版全文。
这些来源只支持候选结构和限制，不提供本装置常数、能量比例或RF认证。
源码实际在adapters/i1_teng_source，不是预期标题中的validation路径；通过源索引定位，未猜API。

## 3. 数学合同和离散式

Q=Q0+int Icol dt[C]；R_T=K_T ell/Q，K_T[V*s/m]。
S=S0+int Icol²dt[A²*s]；R_RW=K_RW ell/sqrt(S)，K_RW[V*sqrt(s)/m]。
ell[m]，两式结果均ohm。Q/S是形成历史，不是源电荷qT、端口CpVp或热库存。
数学单测保留带符号负电流：Toepler Q下降而非取abs(I)，Q<=0拒绝；RW平方历史不下降。
物理端口的pre/post电压均须非负，从而线性区间电压/正R不出现逆向传导；反向请求返回
FORMATIVE_UNIPOLAR_ONLY，超过显式5ns合同返回FORMATION_TIME_LIMIT_REACHED。
独立连续参考也设置首次反向事件停止；本次实际窗口内未触发反向。

给定同一接受态和试探端点V1，以Vbar=(V0+V1)/2、a=K ell同时解：

```text
Toepler: j = Vbar*Qn/(a-dt*Vbar/2)
         Q1=Qn+dt*j
RW:      j = Vbar*sqrt(Sn)/sqrt(a²-dt*Vbar²/2)
         S1=Sn+dt*j²
d(j)/d(V1) analytic -> existing outer circuit Newton Jacobian
Deltaq_port=Cp*(V1-V0)+dt*j
W_column=dt*Vbar*j=dt*j²*R(history_bar)
```

该闭式仅为原形成律的**隐式中点离散联立求解**，不是另一个物理电阻律。
非正离散分母拒绝TRIAL_INTERVAL_OUT_OF_MODEL_DOMAIN；它是本试探积分区间限制，
不称连续正流Toepler模型的有限时物理奇点。没有分母epsilon、R floor或新增恢复项。

## 4. 独立实现、源身份及事务

新增spark_formation_models.py只继承原FP3 IdealCapacitorPort，复用TrialInput/TrialResult、
PortDescriptor和原TENGSeriesCircuit的残差、源电荷及共同commit。旧源码未改。
FormationAccepted扩展不可变接受态，显式formation_value/unit/law；原x_V保持0，不跨单位借用。
每次trial从同一Qn/Sn出发，联立j和历史；并非先积分给定Icol再回灌R。
参数包必须单位/身份齐全，仅SYNTHETIC_NOT_CALIBRATED被当前运行入口授权；冒称实测会拒绝。

预注册两个1V数学事件，CT=Cext=101pF、Cp0.2pF、L100nH、Rs3Ω，i0=Vc0=qT0=0、Vp0=1V。
R0=100Ω、K_T=1e-4、K_RW=1（单位如上）纯软件参数，不来自论文或测量。
ell70um只是名义几何；Q0=7e-11C、S0约4.9e-13A²s由R0/K重参数化。
这不是从绝缘气体自发生成形成历史；外部事件把开路支路接为已给正历史种子的候选支路，
未解释种子的真实来源，亦未为其伪造储能。源/电感/物理电容状态连续；101pF未预充800V。
源模型VT=Voc-qT/CT，没有TENG机械变化。端口i_loop=j+Cp*dVp/dt，形成积分只用j。

在t=0和已有799步历史的接受态，执行两种不同电压试探、重复第一个、reject以及完整
电路trial后受控失败；Q/S、源电荷、Cp、环路态、时钟和累计功/输出计数逐值不变。
唯一commit后拒绝duplicate/stale/foreign及篡改形成候选；不靠只检查无状态电阻证明事务。
串行内存事务，不称进程崩溃持久恢复或原生FULL回滚验证。
Toepler主事件TOTAL，RW主事件SPLIT；每个区间另以相反归属端口实际trial/commit重放同端点，
形成状态、总电荷、总Cp能/传导功差0。没有第三条独立电路事件或新全套归属矩阵。

## 5. 两个实际短程及独立参考

|模型|实际窗口/ns|接受步|末R/Ω|通道功/pJ|最大参考相对误差|
|---|---:|---:|---:|---:|---:|
|TOEPLER|5|1320|67.484994297|19.744274716|2.55454e-05|
|RW|5|1320|83.345217630|19.583231159|2.43332e-05|

每例前0.4ns最大1ps、随后最大5ps，原有界Newton≤8次；没有参数/步长扫描或完整20ns追加。
独立参考是另行手写的连续9变量ODE：i/Vc/qT/Vp/归一Q或S及独立源功、Rs功、通道功、传导电荷。
它不用端口隐式更新式，不调用FormationPort.resistance；DOP853的误差控制1e-11/1e-12，max_step5ps。
参考与主积分属于同一两个合成物理实例，不是另选第三参数事件。形成增量以增量RMS归一，
未用大Q0/S0掩盖小历史变化。源电荷、各状态、形成增量与三个独立功积分均核对。
最大参考误差2.55454e-05，低于预登记1e-3；不是连续物理误差或模型适用性认证。
73项隔离检查通过；另有6段解析/独立积分量纲审计，含正、零、负流和初值/时钟/单位拒绝。

独立从保存端点重算KVL、串联电容/源电荷、Cp加传导电荷、形成更新与导电功，
最大固定尺度残差1.55986e-13，目标1e-8。原始绝对V/C/J误差见独立审计CSV。
W_column仅通道电功/模型耗散，W_source/储能分别保存，不全称气体热；
Q/S无假想dR/dt储能，gas_heat=UNAVAILABLE，无RF/光/电离分配或谱图。

## 6. 可证伪边界与下一决议

RW固定K/ell且无恢复时dR/dt<=0，Icol=0则S/R保持；任何真实无驱动回升要求额外恢复闭合，
不是本模型数值积分能产生的行为。Toepler有符号反向可能令Q减小并达奇点，端口在反向即拒绝，
不解释双极性环振。模型边界明确不等于已知真实形成终点；5ns只是本测试授权窗口。
这两个形成基线不能取代完整10–80MHz存活/熄弧/恢复/再击穿模型。
不新增Mayr/Cassie、经验指数、S/tau、退电离时间常数、效率或物理状态。

[来源与初态要求](SOURCE_AND_INITIALIZATION_REQUIREMENTS.md)列出真实K、有效长度、R0/历史、
形成时刻/阶段、传导支路与源/拓扑/波形均未标定。
[最简路线决议](PHYSICAL_RESISTANCE_NEXT_DECISION.md)：已有敏感性足够，不再任意R(t)包络；
先有真实源/电路和实验端口数据，再审查恢复/热化/非LTE/分段状态的可辨识性。
发射端电压/电流振荡与接收端周期更多的用户陈述保留，不能当相同电流或反演R。
冷GR、空间辐射源、天线传播、接收振铃保持分层，真实RF/TENG/全局端口/AUTO Gate不升级。

## 7. 身份、资源与保护

分支architecture-v2.1-ultrafast-rf；HEAD8ec25920adff63a0b1a875c620b9a32105458ad7，
实际仍含未提交GENERIC_DEV，未冒称Git包含全部接口。
原生科学/合成进程1.030053s，两个主事件共2640接受区间；
参考积分、事务拒绝均包含在此计时，隔离测试与静态后处理另记，远低于600s。
没有PDE、Poisson、G1/G2、冷FULL、EOS、机械TENG、COMSOL、全波RF、编译或历史测试。
保护清单、环境库/执行副本哈希与授权diff单列；仅两份长期文档更新，不提交Git。
旧七案例与中止记录、FP3/FP4、V5/V6/V7等历史Gate保持，本任务完成后停止。

文件：[模型四项状态](spark_formation_gate_summary.json)、[量纲/解析审计](law_dimension_and_domain_audit.csv)、
[冒烟指标](synthetic_formation_coupling_metrics.csv)、[接受步](synthetic_formation_accepted.csv)、
[独立参考](independent_formation_reference_accuracy.csv)、[独立端点账本](independent_accepted_charge_energy_audit.csv)、
[事务](formation_transaction_audit.json)、[执行](execution_record.json)、[来源](source_index.json)。
