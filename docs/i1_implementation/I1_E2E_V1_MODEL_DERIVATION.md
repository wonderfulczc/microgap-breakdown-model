# I1-E2E-V1 模型方程、单位与适用域

这是现有物理方程的一维有限体积适配候选，不宣称新物理定律或已经验证的二维FULL严格降阶。
旧C++、电子系数、G1、FP3、TENG、Stage F/H内核未修改。

## 1. 来源和状态

电子系数实际调用`adapters/i1_conductive_channel/kinetics.py`，其逐项原生等价已有
`validation/i1_part_vi/conductive_channel_bridge_20261016/native_reaction_rate_parity.csv`证据。
本轮读取当前`cpp/src/streamer/MorrowLowke.cpp`、`ReactionModel.cpp`，没有重新编译或更换系数。
反应化学计量为原有电离、二体/三体附着、电子-正离子复合、正负离子复合。
光源本ROM明确关闭；没有把二维SP3结果复制为均匀光源。

状态是每个轴向单元的ne、np、nn [m^-3]。负离子有真实反应和空间电荷，但本版不输运。
气体温度/压力固定，电子能量、重粒子热、辐射化学、激发态和表面能谱均未建模。
电子局部场近似及其低场原生分支不是新的双极性/高密度适用性认证。
离子mu=2e-4m2/(V s)、D=5e-6m2/s与gamma=0均是显式C类合成输入，不是器件物性。

## 2. 轴向几何和保守输运

单元体积Vj=Aj*dz；有效半径及可选正值area_factors是明确模型输入，
不是从针尖曲率或实验电导反演。内部面面积取相邻有效面积的正对称调和平均；
端面取邻近单元面积。首版只建立轴向平面等效映射。
尖端半径/锥角可保存，但不在该映射中制造场增强；本轮改变的是实际进入模型的间隙长度。

连续形式：

    d(A ne)/dt + d(A Gamma_e)/dz = A Se
    d(A np)/dt + d(A Gamma_p)/dz = A Sp
    d(nn)/dt = Sn
    Gamma_e = -mu_e ne E - De d(ne)/dz
    Gamma_p = +mu_p np E - Dp d(np)/dz

代码中面量F=Af*Gamma，单位粒子/s；不是电流密度。
内部面漂移用一阶迎风，扩散用中心差分，公共面只构造一次：
`dn_j/dt=-(F_{j+1/2}-F_{j-1/2})/Vj+S_j`。
固定金属吸收边界仅允许向外漂移及到零浓度边界的扩散；距离是半单元。
阴极由端子电势判定，H/Return永不改名。发射为gamma乘真实正离子向外入射面流，
只扣入电子净边界面流，不再加体源。gamma=0/无入射的独立测试严格为零。

显式Euler按实际移出率和反应损失率限制dt，安全系数0.15；不用密度floor。
负值、非有限值或无法解析的步长直接拒绝试探。资源计数不因拒绝回退。
粒子使用每个子步pre层电场；trial的电压在宏步端点间线性连接。
这是一阶粒子推进与外电路中点的混合离散，不能继承合成RC的整体二阶资格。
本次2ps宏步各接受一个粒子子步，所以粒子显式使用该步起始电压，
已接受电压的变化在下一步进入输运/反应；不存在伪称的全隐式粒子求解。

## 3. 同网络一维Poisson与电荷

epsilon=epsilon0。节点电势为[Vp,phi_cells,0]，链接电容
`kf=epsilon*Af/df [F]`，df内部为dz、端面为dz/2。
自由块K对称正定，`K phi_cells=q+b Vp`，
`qj=e_charge*(np-ne-nn)*Vj [C]`。
这是新ROM的一维小矩阵，不是复制/替换原二维PETSc Poisson。

权重`K psi=b`，几何电容`Cp=1/sum(1/kf)`。
同参考面电极电荷有两个等价表达：

    QH = k_left*(Vp-phi_first)
       = Cp*Vp - psi^T q
    Ufield = 0.5 sum_f kf*(phi_left-phi_right)^2

稳定端点增量使用实际各物种更新差，不用大总电荷相减。
正z从H指向Return。常规电流`Jface=e*(Fp-Fe)`单位A；
H金属进入气体的I_Hg即左端正z常规电流。

    Deltaq_port = DeltaQH + integral I_Hg dt

它包含几何电容、空间电荷感应和电极带电粒子交换。
几何Cp只在组件记一次，ROM首版TOTAL；未实现ROM SPLIT时不能宣称已验证。
原FP3和TENG的旧TOTAL/SPLIT机制保持原样。
单元净电荷预算独立核对物种更新、公共面电荷和实际反应净电荷舍入缺陷，
不从最终残差构造补偿电流。

## 4. 传导诊断不替代端口

诊断sigma=e*(mu_e ne+mu_p np) [S/m]，
在有效轴向欧姆近似下`Rdiag=sum dz/(sigma*Aj)`、Gdiag=1/Rdiag。
若出现零电导单元则报告断路，不加电导地板。
这个冻结诊断忽略扩散、电极吸收及位移响应，绝不用于替代真实trial端口电荷。
它不是内部G_R，也不是器件火花电阻。
本轮均匀1e16m^-3种子约5498电子（70um）；它不是历史Gaussian约96电子的同一空间初态。
100um案例保持相同密度/截面，库存随长度增加；没有暗中保持错误体积。

## 5. TENG-RLC双向连接与事务

源调用P1几何源和原TENGFrozenPositionSource，冻结x/sigma，
`VT=Voc-QT/CT; dQT/dt=i`。R/L/C均来自配置，不硬编码为旧算例值。
原TENGSeriesCircuit/SeriesPortCircuit负责残差、Newton、时间调度和联合commit：

    L di/dt = VT - Rs*i - Vc - Vp
    Cext dVc/dt = i
    i*dt = Deltaq_port

ROM提供同一accepted态的纯trial。对试探Vp做小扰动重新推进以求端口电荷导数，
不用合成RC解析Jacobian冒充粒子导数。拒绝不改变物种、时间、场能及累计账本。
来源电荷、端口token和电路状态一次共同提交。
ROMSeriesCircuit只替换场能库存记账，复用原残差/求根/commit，不复制RLC RHS。
冷FULL和G1通过旧统一入口独立运行；共同结果合同不授予自动冷热转换或新的物理TENG耦合。

## 6. 电功与缺失能量

每个子步用真实常规面电荷与前后平均节点势差：
`Wcond=sum_f Deltaq_face*(phibar_left-phibar_right)`，有符号。
`Wport_native=sum Vpbar_sub*Deltaq_port_sub`。
反应净电荷舍入缺陷的功由phibar与实际源电荷配对。
独立检查`Wport_native+Wreaction_defect-DeltaUfield-Wcond`。
外层功为宏步平均电压乘真实端口电荷，时间求积差单列。
原电路兼容字段`W_column/resistor_dissipation`在ROM出口明确改标为有符号面电功，
不解释为正定电阻热或气体热。没有伪造dR/dt储能项。

`QT²/(2CT)`仍是指定源坐标库存，Voc*i是该坐标驱动功，
不是完整机械能/摩擦充电能。没有计算RF转换效率。

## 7. 等效空间源与接收器

COMPACT_ELECTRIC_DIPOLE：有限正体积相邻单元，端点电荷为±积分i_loop dt，
有显式电容性回流端点；闭合回路模式有明确回流边和零节点积累。
使用真实接受的Cext电荷增量交叉约束转移量；不从守恒残差补电荷。
这是外部导体的条件性空间映射，不把i_port复制为气体J。
没有再次外加位移流或二次发射。rho/J同时保存，FV图连续性精确；
中心重构仍是GRAPH_HALF_EDGE_VOLUME_V1的有限体积求积近似，不冒称逐点连续性。
本版支持坐标轴方向有限源，任意CAD/导线及任意旋转映射仍可替换，不能由标签获得。

实际调用旧GlobalSourceAssembler、SourceSeries、Jefimenko.evaluate_waveform、
Stage H.project_field；无FDTD/传播内核重写。
已知零前史下，旧中心时间导数在启动处有提前插值泄漏；
保留原始场，额外因果支撑门在最短传播距离/c之前输出严格零，
且单列泄漏幅值。该修正不认证启动边缘的时间精度。
只在存在完整延迟源支持的窗口解释场，不延长源历史造尾波。

`Voc_rx=h_eff dot E`，`Vrx=Voc_rx*Rload/(Rant+Rload)`，
最小模型显式Rant/Rload；可选有给定采样间隔的因果FIR，对应频率相关Hrx。
未标定天线不是实测单位增益。接收负载功Vrx²/Rload仅是条件电路量，
没有反向加载发射源，也不是RF辐射效率。
FFT采用去均值Hann窗、真实采样数、不补零；频谱与同极性周期数量分别报告。

## 8. 最值得继续验证的理论关系

应建立二维FULL→一维ROM的保守投影：轴向截面、物种库存、面通量、
同参考面QH及Wcond同时约束，而不是拟合单一R(t)。
G1只在有来源的材料能、质量/动量及LTE适用性可映射时成为热端状态闭合。
现有结果没有给出这些闭合，亦不阻止独立ROM条件性E2E运行。
