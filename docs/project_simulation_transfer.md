# Project Simulation Transfer

Stage 5 transfer purpose: identify which parts of the verified streamer workflow can be reused for later doctoral micro-gap breakdown simulations, and which parts must be replaced before the model is physically appropriate for electrode-driven gaps.

## 1. 已跑通的通用链路

- 气体输运和反应：Stage 3/4/5 使用 Morrow-Lowke analytic baseline，能驱动三粒子 streamer 求解。
- 三粒子连续性方程：电子、正离子、负离子已耦合显式推进。
- Poisson：复用 Stage 2 OpenCharge Poisson。
- SP3：复用 Stage 2 三组 SP3 与 Robin 边界耦合。
- ISG-0：电子通量推进中已使用 Stage 2 ISG-0。
- 双流注传播与碰撞：Stage 4 paper-like 20 um 工况已完成碰撞流程；Stage 5 高场工况复用了该流程。
- 电流矩：C++ 每个 accepted step 直接积分 drift/electron current moment。
- FFT/ESD：Stage 1 频谱/ESD 约定在 Stage 4/5 后处理中复用。
- 频带能量：VHF/UHF/SHF 输出带 trusted/partial/untrusted 标记。

## 2. 可直接复用到课题的模块

- 网格与场数据结构：`AxisymmetricGrid`、`ScalarField2D`。
- PETSc/MPI 框架：Poisson/SP3 和短窗口 rank consistency 已验证。
- 自适应时间步：Stage 3–5 均使用同一控制逻辑。
- checkpoint：长运行和中断恢复流程已建立。
- 电流矩积分：`stage4_run` 直接输出 `current_moment.csv`。
- 辐射后处理：Delta current、导数、FFT、ESD、频带能量。
- 证据登记和验证器：run registry、artifact manifest、validation matrix、closure validator。

## 3. 课题中必须替换的内容

- 自由空间 Gaussian 种子：微间隙应由电极附近初始电离、背景电离或表面发射模型替代。
- 均匀背景场：必须替换为实际电极几何产生的非均匀场。
- Morrow-Lowke 输运模型：当前只适合 workflow baseline；目标气体、压力、湿度和微间隙高场需要重新校准或替换为 Boltzmann 表。
- 无电极边界：必须加入电极 Dirichlet/材料边界。
- 固定空气气氛：课题气体组分、湿度和压力需参数化。
- 纯流注碰撞场景：微间隙击穿更关注电极间发展、电流闭合和外电路耦合。

## 4. 微间隙课题下一步模型

建议顺序：

1. 加入实际针/箔/球电极几何。
2. 加入电极 Dirichlet 边界。
3. 建立电极附近种子电子或背景电离。
4. 引入实际电压波形。
5. 校准微间隙气体输运参数。
6. 输出时变放电电流。
7. 计算电流矩和原生宽频辐射。
8. 再与外部 RLC 网络耦合。
9. 区分原生谱与 RLC 调制谱。

## 5. 当前模型不能直接回答的问题

- 电极表面发射；
- 阴极鞘层；
- 热化与火花通道；
- 电极侵蚀；
- 等离子体温度；
- 外部 RLC 反作用；
- 击穿后的完整热等离子体阶段。

## 6. 推荐的课题仿真最小闭环

电极几何
→ 电场增强
→ 种子形成
→ 流注/击穿发展
→ 时变电流
→ 电流矩
→ 原生辐射
→ RLC 调制
→ 接收信号

当前代码没有完成上述全部课题模型；它完成的是 solver、碰撞、辐射链路、趋势逻辑和证据审计框架，可作为微间隙课题的工程起点。
