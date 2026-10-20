# I1PredictivePipeline：入口、模块与替换边界

## 运行入口

```bash
PYTHONPATH=python:. PYTHONDONTWRITEBYTECODE=1 .venv/bin/python \
  -m adapters.i1_e2e.pipeline scenario.json NEW_OUTPUT_DIRECTORY
```

输出目录必须不存在；现有四案例已经完成，不自动重跑。
示例输入位于本目录`synthetic_scenario.json`、`rom_R3_scenario.json`、
`rom_R300_scenario.json`、`rom_gap100um_scenario.json`。
改变参数需新文件、新输出、明确授权；不是静默改写已有结果。

## 实际调用图

| 模块 | 代码/复用 | 实际能力 |
|---|---|---|
| GeometryProvider | i1_e2e/pipeline.py | 正值几何检查、有效面积；一维长度真实进入K/FV |
| TENGSourceProvider | i1_p1/sources.py + i1_teng_source/source.py | 制造几何→CT/Voc→动态转移电荷源 |
| ExternalCircuitProvider | 原CircuitParameters/TENGSeriesCircuit | 可配置正值Rs/L/Cext，仅既有串联拓扑 |
| DischargeChannelProvider | 原ReferenceParallelRCPort；新rom.py | 数学测试端口与自主物种状态ROM严格分开 |
| CircuitPlasmaCoupler | 原FP3/TENG Newton与事务；ROMSeriesCircuit | 实际trial/commit，数值电荷Jacobian；只扩展场能记账 |
| RadiationSourceProvider | i1_e2e/radiation.py + 原GlobalSourcePacket | 有限尺寸偶极/闭合回路，rho与J及回流同时保存 |
| PropagationProvider | 原Stage F Jefimenko | 四案例实际场计算，源/时钟身份保持 |
| ReceiverProvider | 原Stage H project_field + ReceiverVoltageModel | 有效高度、负载及可选因果FIR |
| AnalysisAndComparison | review.py | 已接受CSV对照、独立矩阵/偶极参考、窗口合格性 |
| 共用状态 | interfaces.py | 电压、不同时间层电流、电荷、电功及缺失库存 |

没有复制冷FULL、G1、旧RLC或传播内核。新一维Poisson属于明确授权的ROM，
并不与二维原生算子互换。旧参数包以内容SHA读取。

## 三类提供器的共同合同

ROM和数学端口直接复用FP3 PortDescriptor/TrialInput/TrialResult及accepted token。
`I1_CHANNEL_PROVIDER_RESULT_1.0`保存physics_mode、time、voltage、currents、
charge、energy、state_source与scope。
端点环路i、区间I_Hg、几何位移流和空间电荷感应不能混同。
材料能不可用时保留null+reason，不由电功残差补造。

外层请求`mode=COLD_FULL/GIVEN_THERMAL_SERIES`时，必须附真实旧入口配置
`legacy_scenario`且模式匹配；直接调用原I1MultiphaseScenarioRunner。
本轮只用已归档结果检查共同Envelope，没有新原生/热步。
这些模式保留旧源与物理限制，不谎称已接通任意TENG物理端口。
AUTO_COLD_TO_THERMAL不需要构造求解器就返回MISSING_PHYSICAL_CLOSURE。

## 时间、初态及数值限制

所有E2E案例源、电路、场包使用同一秒制时间；接收时间包含空间传播延迟。
独立参考初态i=Vc=QT=0、Vp=Voc，明确SYNTHETIC_EVENT，
不从500/800V猜测Cext预充。其他充电历史首版拒绝，而非伪造兼容初态。
机械状态冻结；动态QT实际积分。网表只有首版串联结构，其他结构明确拒绝。

ROM粒子显式稳定子步+外电路中点，非全隐式二阶PDE。
数学端口有dt/(R Cp)快模态诊断。已执行合成20ps相对1ps快模态欠解析，
最新入口明确输出FAST_INITIAL_MODE_UNDERRESOLVED；守恒不能取代波形精度。
不自动缩步重算已有案例。晚时间场另有独立参考，但不能替代启动段认证。

## 可复现身份

每例保留scenario_inputs、provenance、源graph.npz/manifest、接受CSV和场/接收输出。
执行版数值模块副本见executed_sources与executed_versions.json。
最终仅补充共同Envelope、时间精度警告、频谱资格/奇数PSD归一及源适配器哈希；
这些不改变已经接受的轨迹/传播场。原数值文件不覆盖，reviewed文件另存。
首个fixture构建脚本曾提取公共ROM配置函数，旧哈希已登记但该脚本旧字节未另存；
合成案例的完整执行JSON及实际数值内核版本均已保存，不声称所有生成器逐字节封存。
