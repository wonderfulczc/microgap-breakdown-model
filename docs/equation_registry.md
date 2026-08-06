# Equation Registry

No PDE equation is implemented in Stage 1.1.

## EQ-0001 electron continuity

- 模块：fluid_model
- 数学形式：$∂n_e/∂t+∇·Γ_e=S_i+S_ph-S_att$
- 变量定义：符号遵循来源；离散下标为网格节点，波浪号为 Fourier 变换。
- 每个变量的 SI 单位：长度 m；时间 s；密度 m^-3；电场 V/m；电势 V；电流 A；电流矩 A m；频率 Hz 或角频率 rad/s；能量 J。
- 方程左右量纲：已按 SI 核对；SP3 使用 Table 3 A_j 单位族，ESD 按角频率定义。
- 来源：SRC-0007
- 来源位置：model equations/current definitions
- provenance_type：direct
- 在 Shi 2019 主复刻中的作用：Stage 2 以后候选关系；Stage 1.1 仅登记。
- 是否直接实现：否
- 计划数值离散：待 A2 决策，本阶段不实现。
- 所需边界条件：按模块来源；未公开细节列入 unknowns。
- 所需初始条件：见 parameter_registry.csv。
- 验证算例：Liu 2007 SP3、Kulikovsky transport 或 Shi/Luque radiation as applicable。
- 当前状态：registered_not_implemented
- 已知歧义：None beyond cited sign convention.

## EQ-0002 positive-ion continuity

- 模块：fluid_model
- 数学形式：$∂n_p/∂t=S_i+S_ph$
- 变量定义：符号遵循来源；离散下标为网格节点，波浪号为 Fourier 变换。
- 每个变量的 SI 单位：长度 m；时间 s；密度 m^-3；电场 V/m；电势 V；电流 A；电流矩 A m；频率 Hz 或角频率 rad/s；能量 J。
- 方程左右量纲：已按 SI 核对；SP3 使用 Table 3 A_j 单位族，ESD 按角频率定义。
- 来源：SRC-0007
- 来源位置：model equations/current definitions
- provenance_type：direct
- 在 Shi 2019 主复刻中的作用：Stage 2 以后候选关系；Stage 1.1 仅登记。
- 是否直接实现：否
- 计划数值离散：待 A2 决策，本阶段不实现。
- 所需边界条件：按模块来源；未公开细节列入 unknowns。
- 所需初始条件：见 parameter_registry.csv。
- 验证算例：Liu 2007 SP3、Kulikovsky transport 或 Shi/Luque radiation as applicable。
- 当前状态：registered_not_implemented
- 已知歧义：None beyond cited sign convention.

## EQ-0003 negative-ion continuity

- 模块：fluid_model
- 数学形式：$∂n_n/∂t=S_att$
- 变量定义：符号遵循来源；离散下标为网格节点，波浪号为 Fourier 变换。
- 每个变量的 SI 单位：长度 m；时间 s；密度 m^-3；电场 V/m；电势 V；电流 A；电流矩 A m；频率 Hz 或角频率 rad/s；能量 J。
- 方程左右量纲：已按 SI 核对；SP3 使用 Table 3 A_j 单位族，ESD 按角频率定义。
- 来源：SRC-0007
- 来源位置：model equations/current definitions
- provenance_type：direct
- 在 Shi 2019 主复刻中的作用：Stage 2 以后候选关系；Stage 1.1 仅登记。
- 是否直接实现：否
- 计划数值离散：待 A2 决策，本阶段不实现。
- 所需边界条件：按模块来源；未公开细节列入 unknowns。
- 所需初始条件：见 parameter_registry.csv。
- 验证算例：Liu 2007 SP3、Kulikovsky transport 或 Shi/Luque radiation as applicable。
- 当前状态：registered_not_implemented
- 已知歧义：None beyond cited sign convention.

## EQ-0004 electron drift-diffusion flux

- 模块：fluid_model
- 数学形式：$Γ_e=-μ_e n_e E-D_e∇n_e$
- 变量定义：符号遵循来源；离散下标为网格节点，波浪号为 Fourier 变换。
- 每个变量的 SI 单位：长度 m；时间 s；密度 m^-3；电场 V/m；电势 V；电流 A；电流矩 A m；频率 Hz 或角频率 rad/s；能量 J。
- 方程左右量纲：已按 SI 核对；SP3 使用 Table 3 A_j 单位族，ESD 按角频率定义。
- 来源：SRC-0007
- 来源位置：model equations/current definitions
- provenance_type：direct
- 在 Shi 2019 主复刻中的作用：Stage 2 以后候选关系；Stage 1.1 仅登记。
- 是否直接实现：否
- 计划数值离散：待 A2 决策，本阶段不实现。
- 所需边界条件：按模块来源；未公开细节列入 unknowns。
- 所需初始条件：见 parameter_registry.csv。
- 验证算例：Liu 2007 SP3、Kulikovsky transport 或 Shi/Luque radiation as applicable。
- 当前状态：registered_not_implemented
- 已知歧义：None beyond cited sign convention.

## EQ-0005 Poisson equation

- 模块：fluid_model
- 数学形式：$∇²φ=-ρ/ε0$
- 变量定义：符号遵循来源；离散下标为网格节点，波浪号为 Fourier 变换。
- 每个变量的 SI 单位：长度 m；时间 s；密度 m^-3；电场 V/m；电势 V；电流 A；电流矩 A m；频率 Hz 或角频率 rad/s；能量 J。
- 方程左右量纲：已按 SI 核对；SP3 使用 Table 3 A_j 单位族，ESD 按角频率定义。
- 来源：SRC-0007
- 来源位置：model equations/current definitions
- provenance_type：direct
- 在 Shi 2019 主复刻中的作用：Stage 2 以后候选关系；Stage 1.1 仅登记。
- 是否直接实现：否
- 计划数值离散：待 A2 决策，本阶段不实现。
- 所需边界条件：按模块来源；未公开细节列入 unknowns。
- 所需初始条件：见 parameter_registry.csv。
- 验证算例：Liu 2007 SP3、Kulikovsky transport 或 Shi/Luque radiation as applicable。
- 当前状态：registered_not_implemented
- 已知歧义：None beyond cited sign convention.

## EQ-0006 electric field definition

- 模块：fluid_model
- 数学形式：$E=-∇φ$
- 变量定义：符号遵循来源；离散下标为网格节点，波浪号为 Fourier 变换。
- 每个变量的 SI 单位：长度 m；时间 s；密度 m^-3；电场 V/m；电势 V；电流 A；电流矩 A m；频率 Hz 或角频率 rad/s；能量 J。
- 方程左右量纲：已按 SI 核对；SP3 使用 Table 3 A_j 单位族，ESD 按角频率定义。
- 来源：SRC-0007
- 来源位置：model equations/current definitions
- provenance_type：direct
- 在 Shi 2019 主复刻中的作用：Stage 2 以后候选关系；Stage 1.1 仅登记。
- 是否直接实现：否
- 计划数值离散：待 A2 决策，本阶段不实现。
- 所需边界条件：按模块来源；未公开细节列入 unknowns。
- 所需初始条件：见 parameter_registry.csv。
- 验证算例：Liu 2007 SP3、Kulikovsky transport 或 Shi/Luque radiation as applicable。
- 当前状态：registered_not_implemented
- 已知歧义：None beyond cited sign convention.

## EQ-0007 Gaussian neutral plasma seed

- 模块：fluid_model
- 数学形式：$n=n0 exp(-distance²/σ²)$
- 变量定义：符号遵循来源；离散下标为网格节点，波浪号为 Fourier 变换。
- 每个变量的 SI 单位：长度 m；时间 s；密度 m^-3；电场 V/m；电势 V；电流 A；电流矩 A m；频率 Hz 或角频率 rad/s；能量 J。
- 方程左右量纲：已按 SI 核对；SP3 使用 Table 3 A_j 单位族，ESD 按角频率定义。
- 来源：SRC-0007
- 来源位置：model equations/current definitions
- provenance_type：direct
- 在 Shi 2019 主复刻中的作用：Stage 2 以后候选关系；Stage 1.1 仅登记。
- 是否直接实现：否
- 计划数值离散：待 A2 决策，本阶段不实现。
- 所需边界条件：按模块来源；未公开细节列入 unknowns。
- 所需初始条件：见 parameter_registry.csv。
- 验证算例：Liu 2007 SP3、Kulikovsky transport 或 Shi/Luque radiation as applicable。
- 当前状态：registered_not_implemented
- 已知歧义：None beyond cited sign convention.

## EQ-0008 Zheleznyak integral

- 模块：photoionization
- 数学形式：$S_ph=∫I(r′)g(R)/(4πR²)dV′$
- 变量定义：符号遵循来源；离散下标为网格节点，波浪号为 Fourier 变换。
- 每个变量的 SI 单位：长度 m；时间 s；密度 m^-3；电场 V/m；电势 V；电流 A；电流矩 A m；频率 Hz 或角频率 rad/s；能量 J。
- 方程左右量纲：已按 SI 核对；SP3 使用 Table 3 A_j 单位族，ESD 按角频率定义。
- 来源：SRC-0003
- 来源位置：Bourdon PDF pp.3-9, Eqs.(1)-(27)
- provenance_type：direct
- 在 Shi 2019 主复刻中的作用：Stage 2 以后候选关系；Stage 1.1 仅登记。
- 是否直接实现：否
- 计划数值离散：待 A2 决策，本阶段不实现。
- 所需边界条件：按模块来源；未公开细节列入 unknowns。
- 所需初始条件：见 parameter_registry.csv。
- 验证算例：Liu 2007 SP3、Kulikovsky transport 或 Shi/Luque radiation as applicable。
- 当前状态：registered_not_implemented
- 已知歧义：None beyond cited sign convention.

## EQ-0009 photon emission source

- 模块：photoionization
- 数学形式：$I=ξ[pq/(p+pq)]νu ne$
- 变量定义：符号遵循来源；离散下标为网格节点，波浪号为 Fourier 变换。
- 每个变量的 SI 单位：长度 m；时间 s；密度 m^-3；电场 V/m；电势 V；电流 A；电流矩 A m；频率 Hz 或角频率 rad/s；能量 J。
- 方程左右量纲：已按 SI 核对；SP3 使用 Table 3 A_j 单位族，ESD 按角频率定义。
- 来源：SRC-0003
- 来源位置：Bourdon PDF pp.3-9, Eqs.(1)-(27)
- provenance_type：direct
- 在 Shi 2019 主复刻中的作用：Stage 2 以后候选关系；Stage 1.1 仅登记。
- 是否直接实现：否
- 计划数值离散：待 A2 决策，本阶段不实现。
- 所需边界条件：按模块来源；未公开细节列入 unknowns。
- 所需初始条件：见 parameter_registry.csv。
- 验证算例：Liu 2007 SP3、Kulikovsky transport 或 Shi/Luque radiation as applicable。
- 当前状态：registered_not_implemented
- 已知歧义：None beyond cited sign convention.

## EQ-0010 absorption function

- 模块：photoionization
- 数学形式：$g/pO2=[exp(-χmin pO2R)-exp(-χmax pO2R)]/[pO2R ln(χmax/χmin)]$
- 变量定义：符号遵循来源；离散下标为网格节点，波浪号为 Fourier 变换。
- 每个变量的 SI 单位：长度 m；时间 s；密度 m^-3；电场 V/m；电势 V；电流 A；电流矩 A m；频率 Hz 或角频率 rad/s；能量 J。
- 方程左右量纲：已按 SI 核对；SP3 使用 Table 3 A_j 单位族，ESD 按角频率定义。
- 来源：SRC-0003
- 来源位置：Bourdon PDF pp.3-9, Eqs.(1)-(27)
- provenance_type：direct
- 在 Shi 2019 主复刻中的作用：Stage 2 以后候选关系；Stage 1.1 仅登记。
- 是否直接实现：否
- 计划数值离散：待 A2 决策，本阶段不实现。
- 所需边界条件：按模块来源；未公开细节列入 unknowns。
- 所需初始条件：见 parameter_registry.csv。
- 验证算例：Liu 2007 SP3、Kulikovsky transport 或 Shi/Luque radiation as applicable。
- 当前状态：registered_not_implemented
- 已知歧义：None beyond cited sign convention.

## EQ-0011 SP3 equation 1

- 模块：photoionization
- 数学形式：$∇²φ1-(λpO2)²φ1/κ1²=-(λpO2/κ1²)νu/(cτu)$
- 变量定义：符号遵循来源；离散下标为网格节点，波浪号为 Fourier 变换。
- 每个变量的 SI 单位：长度 m；时间 s；密度 m^-3；电场 V/m；电势 V；电流 A；电流矩 A m；频率 Hz 或角频率 rad/s；能量 J。
- 方程左右量纲：已按 SI 核对；SP3 使用 Table 3 A_j 单位族，ESD 按角频率定义。
- 来源：SRC-0003
- 来源位置：Bourdon PDF pp.3-9, Eqs.(1)-(27)
- provenance_type：direct
- 在 Shi 2019 主复刻中的作用：Stage 2 以后候选关系；Stage 1.1 仅登记。
- 是否直接实现：否
- 计划数值离散：待 A2 决策，本阶段不实现。
- 所需边界条件：按模块来源；未公开细节列入 unknowns。
- 所需初始条件：见 parameter_registry.csv。
- 验证算例：Liu 2007 SP3、Kulikovsky transport 或 Shi/Luque radiation as applicable。
- 当前状态：registered_not_implemented
- 已知歧义：None beyond cited sign convention.

## EQ-0012 SP3 equation 2

- 模块：photoionization
- 数学形式：$∇²φ2-(λpO2)²φ2/κ2²=-(λpO2/κ2²)νu/(cτu)$
- 变量定义：符号遵循来源；离散下标为网格节点，波浪号为 Fourier 变换。
- 每个变量的 SI 单位：长度 m；时间 s；密度 m^-3；电场 V/m；电势 V；电流 A；电流矩 A m；频率 Hz 或角频率 rad/s；能量 J。
- 方程左右量纲：已按 SI 核对；SP3 使用 Table 3 A_j 单位族，ESD 按角频率定义。
- 来源：SRC-0003
- 来源位置：Bourdon PDF pp.3-9, Eqs.(1)-(27)
- provenance_type：direct
- 在 Shi 2019 主复刻中的作用：Stage 2 以后候选关系；Stage 1.1 仅登记。
- 是否直接实现：否
- 计划数值离散：待 A2 决策，本阶段不实现。
- 所需边界条件：按模块来源；未公开细节列入 unknowns。
- 所需初始条件：见 parameter_registry.csv。
- 验证算例：Liu 2007 SP3、Kulikovsky transport 或 Shi/Luque radiation as applicable。
- 当前状态：registered_not_implemented
- 已知歧义：None beyond cited sign convention.

## EQ-0013 SP3 isotropic combination

- 模块：photoionization
- 数学形式：$Ψ=(γ2φ1-γ1φ2)/(γ2-γ1)$
- 变量定义：符号遵循来源；离散下标为网格节点，波浪号为 Fourier 变换。
- 每个变量的 SI 单位：长度 m；时间 s；密度 m^-3；电场 V/m；电势 V；电流 A；电流矩 A m；频率 Hz 或角频率 rad/s；能量 J。
- 方程左右量纲：已按 SI 核对；SP3 使用 Table 3 A_j 单位族，ESD 按角频率定义。
- 来源：SRC-0003
- 来源位置：Bourdon PDF pp.3-9, Eqs.(1)-(27)
- provenance_type：direct
- 在 Shi 2019 主复刻中的作用：Stage 2 以后候选关系；Stage 1.1 仅登记。
- 是否直接实现：否
- 计划数值离散：待 A2 决策，本阶段不实现。
- 所需边界条件：按模块来源；未公开细节列入 unknowns。
- 所需初始条件：见 parameter_registry.csv。
- 验证算例：Liu 2007 SP3、Kulikovsky transport 或 Shi/Luque radiation as applicable。
- 当前状态：registered_not_implemented
- 已知歧义：None beyond cited sign convention.

## EQ-0014 total photoionization source

- 模块：photoionization
- 数学形式：$S_ph=Σj Aj pO2 Ψj$
- 变量定义：符号遵循来源；离散下标为网格节点，波浪号为 Fourier 变换。
- 每个变量的 SI 单位：长度 m；时间 s；密度 m^-3；电场 V/m；电势 V；电流 A；电流矩 A m；频率 Hz 或角频率 rad/s；能量 J。
- 方程左右量纲：已按 SI 核对；SP3 使用 Table 3 A_j 单位族，ESD 按角频率定义。
- 来源：SRC-0003
- 来源位置：Bourdon PDF pp.3-9, Eqs.(1)-(27)
- provenance_type：direct
- 在 Shi 2019 主复刻中的作用：Stage 2 以后候选关系；Stage 1.1 仅登记。
- 是否直接实现：否
- 计划数值离散：待 A2 决策，本阶段不实现。
- 所需边界条件：按模块来源；未公开细节列入 unknowns。
- 所需初始条件：见 parameter_registry.csv。
- 验证算例：Liu 2007 SP3、Kulikovsky transport 或 Shi/Luque radiation as applicable。
- 当前状态：registered_not_implemented
- 已知歧义：None beyond cited sign convention.

## EQ-0015 SP3 boundary conditions

- 模块：photoionization
- 数学形式：$coupled Robin conditions; Liu 2007 Eqs.(1)-(2)$
- 变量定义：符号遵循来源；离散下标为网格节点，波浪号为 Fourier 变换。
- 每个变量的 SI 单位：长度 m；时间 s；密度 m^-3；电场 V/m；电势 V；电流 A；电流矩 A m；频率 Hz 或角频率 rad/s；能量 J。
- 方程左右量纲：已按 SI 核对；SP3 使用 Table 3 A_j 单位族，ESD 按角频率定义。
- 来源：SRC-0003
- 来源位置：Bourdon PDF pp.3-9, Eqs.(1)-(27)
- provenance_type：direct
- 在 Shi 2019 主复刻中的作用：Stage 2 以后候选关系；Stage 1.1 仅登记。
- 是否直接实现：否
- 计划数值离散：待 A2 决策，本阶段不实现。
- 所需边界条件：按模块来源；未公开细节列入 unknowns。
- 所需初始条件：见 parameter_registry.csv。
- 验证算例：Liu 2007 SP3、Kulikovsky transport 或 Shi/Luque radiation as applicable。
- 当前状态：registered_not_implemented
- 已知歧义：None beyond cited sign convention.

## EQ-0016 ordinary SG flux

- 模块：numerical_flux
- 数学形式：$j=D/(hI0)[nk-exp(α)n(k+1)]$
- 变量定义：符号遵循来源；离散下标为网格节点，波浪号为 Fourier 变换。
- 每个变量的 SI 单位：长度 m；时间 s；密度 m^-3；电场 V/m；电势 V；电流 A；电流矩 A m；频率 Hz 或角频率 rad/s；能量 J。
- 方程左右量纲：已按 SI 核对；SP3 使用 Table 3 A_j 单位族，ESD 按角频率定义。
- 来源：SRC-0005
- 来源位置：Kulikovsky PDF pp.1-5, Eqs.(4),(20)-(31)
- provenance_type：direct
- 在 Shi 2019 主复刻中的作用：Stage 2 以后候选关系；Stage 1.1 仅登记。
- 是否直接实现：否
- 计划数值离散：待 A2 决策，本阶段不实现。
- 所需边界条件：按模块来源；未公开细节列入 unknowns。
- 所需初始条件：见 parameter_registry.csv。
- 验证算例：Liu 2007 SP3、Kulikovsky transport 或 Shi/Luque radiation as applicable。
- 当前状态：registered_not_implemented
- 已知歧义：None beyond cited sign convention.

## EQ-0017 virtual-node width

- 模块：numerical_flux
- 数学形式：$hv=sqrt(2εDh/|ΔW|)$
- 变量定义：符号遵循来源；离散下标为网格节点，波浪号为 Fourier 变换。
- 每个变量的 SI 单位：长度 m；时间 s；密度 m^-3；电场 V/m；电势 V；电流 A；电流矩 A m；频率 Hz 或角频率 rad/s；能量 J。
- 方程左右量纲：已按 SI 核对；SP3 使用 Table 3 A_j 单位族，ESD 按角频率定义。
- 来源：SRC-0005
- 来源位置：Kulikovsky PDF pp.1-5, Eqs.(4),(20)-(31)
- provenance_type：direct
- 在 Shi 2019 主复刻中的作用：Stage 2 以后候选关系；Stage 1.1 仅登记。
- 是否直接实现：否
- 计划数值离散：待 A2 决策，本阶段不实现。
- 所需边界条件：按模块来源；未公开细节列入 unknowns。
- 所需初始条件：见 parameter_registry.csv。
- 验证算例：Liu 2007 SP3、Kulikovsky transport 或 Shi/Luque radiation as applicable。
- 当前状态：registered_not_implemented
- 已知歧义：None beyond cited sign convention.

## EQ-0018 virtual-node locations

- 模块：numerical_flux
- 数学形式：$xLR=midpoint∓hv/2$
- 变量定义：符号遵循来源；离散下标为网格节点，波浪号为 Fourier 变换。
- 每个变量的 SI 单位：长度 m；时间 s；密度 m^-3；电场 V/m；电势 V；电流 A；电流矩 A m；频率 Hz 或角频率 rad/s；能量 J。
- 方程左右量纲：已按 SI 核对；SP3 使用 Table 3 A_j 单位族，ESD 按角频率定义。
- 来源：SRC-0005
- 来源位置：Kulikovsky PDF pp.1-5, Eqs.(4),(20)-(31)
- provenance_type：direct
- 在 Shi 2019 主复刻中的作用：Stage 2 以后候选关系；Stage 1.1 仅登记。
- 是否直接实现：否
- 计划数值离散：待 A2 决策，本阶段不实现。
- 所需边界条件：按模块来源；未公开细节列入 unknowns。
- 所需初始条件：见 parameter_registry.csv。
- 验证算例：Liu 2007 SP3、Kulikovsky transport 或 Shi/Luque radiation as applicable。
- 当前状态：registered_not_implemented
- 已知歧义：None beyond cited sign convention.

## EQ-0019 exponential density interpolation

- 模块：numerical_flux
- 数学形式：$n(x)+1=(nk+1)exp[a(x-xk)]$
- 变量定义：符号遵循来源；离散下标为网格节点，波浪号为 Fourier 变换。
- 每个变量的 SI 单位：长度 m；时间 s；密度 m^-3；电场 V/m；电势 V；电流 A；电流矩 A m；频率 Hz 或角频率 rad/s；能量 J。
- 方程左右量纲：已按 SI 核对；SP3 使用 Table 3 A_j 单位族，ESD 按角频率定义。
- 来源：SRC-0005
- 来源位置：Kulikovsky PDF pp.1-5, Eqs.(4),(20)-(31)
- provenance_type：direct
- 在 Shi 2019 主复刻中的作用：Stage 2 以后候选关系；Stage 1.1 仅登记。
- 是否直接实现：否
- 计划数值离散：待 A2 决策，本阶段不实现。
- 所需边界条件：按模块来源；未公开细节列入 unknowns。
- 所需初始条件：见 parameter_registry.csv。
- 验证算例：Liu 2007 SP3、Kulikovsky transport 或 Shi/Luque radiation as applicable。
- 当前状态：registered_not_implemented
- 已知歧义：None beyond cited sign convention.

## EQ-0020 ISG-0 flux

- 模块：numerical_flux
- 数学形式：$j=D/(hvI0)[nL-exp(αv)nR]$
- 变量定义：符号遵循来源；离散下标为网格节点，波浪号为 Fourier 变换。
- 每个变量的 SI 单位：长度 m；时间 s；密度 m^-3；电场 V/m；电势 V；电流 A；电流矩 A m；频率 Hz 或角频率 rad/s；能量 J。
- 方程左右量纲：已按 SI 核对；SP3 使用 Table 3 A_j 单位族，ESD 按角频率定义。
- 来源：SRC-0005
- 来源位置：Kulikovsky PDF pp.1-5, Eqs.(4),(20)-(31)
- provenance_type：direct
- 在 Shi 2019 主复刻中的作用：Stage 2 以后候选关系；Stage 1.1 仅登记。
- 是否直接实现：否
- 计划数值离散：待 A2 决策，本阶段不实现。
- 所需边界条件：按模块来源；未公开细节列入 unknowns。
- 所需初始条件：见 parameter_registry.csv。
- 验证算例：Liu 2007 SP3、Kulikovsky transport 或 Shi/Luque radiation as applicable。
- 当前状态：registered_not_implemented
- 已知歧义：None beyond cited sign convention.

## EQ-0021 ISG-0 zero-width limit

- 模块：numerical_flux
- 数学形式：$lim j=nmid[Wmid-(D/h)ln((n(k+1)+1)/(nk+1))]$
- 变量定义：符号遵循来源；离散下标为网格节点，波浪号为 Fourier 变换。
- 每个变量的 SI 单位：长度 m；时间 s；密度 m^-3；电场 V/m；电势 V；电流 A；电流矩 A m；频率 Hz 或角频率 rad/s；能量 J。
- 方程左右量纲：已按 SI 核对；SP3 使用 Table 3 A_j 单位族，ESD 按角频率定义。
- 来源：SRC-0005
- 来源位置：Kulikovsky PDF pp.1-5, Eqs.(4),(20)-(31)
- provenance_type：direct
- 在 Shi 2019 主复刻中的作用：Stage 2 以后候选关系；Stage 1.1 仅登记。
- 是否直接实现：否
- 计划数值离散：待 A2 决策，本阶段不实现。
- 所需边界条件：按模块来源；未公开细节列入 unknowns。
- 所需初始条件：见 parameter_registry.csv。
- 验证算例：Liu 2007 SP3、Kulikovsky transport 或 Shi/Luque radiation as applicable。
- 当前状态：registered_not_implemented
- 已知歧义：None beyond cited sign convention.

## EQ-0022 cylindrical divergence

- 模块：fluid_model
- 数学形式：$∇·Γ=(1/r)∂(rΓr)/∂r+∂Γz/∂z$
- 变量定义：符号遵循来源；离散下标为网格节点，波浪号为 Fourier 变换。
- 每个变量的 SI 单位：长度 m；时间 s；密度 m^-3；电场 V/m；电势 V；电流 A；电流矩 A m；频率 Hz 或角频率 rad/s；能量 J。
- 方程左右量纲：已按 SI 核对；SP3 使用 Table 3 A_j 单位族，ESD 按角频率定义。
- 来源：SRC-0007
- 来源位置：model equations/current definitions
- provenance_type：direct
- 在 Shi 2019 主复刻中的作用：Stage 2 以后候选关系；Stage 1.1 仅登记。
- 是否直接实现：否
- 计划数值离散：待 A2 决策，本阶段不实现。
- 所需边界条件：按模块来源；未公开细节列入 unknowns。
- 所需初始条件：见 parameter_registry.csv。
- 验证算例：Liu 2007 SP3、Kulikovsky transport 或 Shi/Luque radiation as applicable。
- 当前状态：registered_not_implemented
- 已知歧义：None beyond cited sign convention.

## EQ-0023 axial current density

- 模块：radiation/current_moment
- 数学形式：$Jz=Σs qs Γs,z (included terms unresolved)$
- 变量定义：符号遵循来源；离散下标为网格节点，波浪号为 Fourier 变换。
- 每个变量的 SI 单位：长度 m；时间 s；密度 m^-3；电场 V/m；电势 V；电流 A；电流矩 A m；频率 Hz 或角频率 rad/s；能量 J。
- 方程左右量纲：已按 SI 核对；SP3 使用 Table 3 A_j 单位族，ESD 按角频率定义。
- 来源：SRC-0007
- 来源位置：model equations/current definitions
- provenance_type：direct
- 在 Shi 2019 主复刻中的作用：Stage 2 以后候选关系；Stage 1.1 仅登记。
- 是否直接实现：否
- 计划数值离散：待 A2 决策，本阶段不实现。
- 所需边界条件：按模块来源；未公开细节列入 unknowns。
- 所需初始条件：见 parameter_registry.csv。
- 验证算例：Liu 2007 SP3、Kulikovsky transport 或 Shi/Luque radiation as applicable。
- 当前状态：registered_not_implemented
- 已知歧义：Whether diffusion, ion, or displacement current is included remains unresolved.

## EQ-0024 cross-sectional current

- 模块：radiation/current_moment
- 数学形式：$I(z,t)=∫Jz 2πr dr$
- 变量定义：符号遵循来源；离散下标为网格节点，波浪号为 Fourier 变换。
- 每个变量的 SI 单位：长度 m；时间 s；密度 m^-3；电场 V/m；电势 V；电流 A；电流矩 A m；频率 Hz 或角频率 rad/s；能量 J。
- 方程左右量纲：已按 SI 核对；SP3 使用 Table 3 A_j 单位族，ESD 按角频率定义。
- 来源：SRC-0007
- 来源位置：model equations/current definitions
- provenance_type：direct
- 在 Shi 2019 主复刻中的作用：Stage 2 以后候选关系；Stage 1.1 仅登记。
- 是否直接实现：否
- 计划数值离散：待 A2 决策，本阶段不实现。
- 所需边界条件：按模块来源；未公开细节列入 unknowns。
- 所需初始条件：见 parameter_registry.csv。
- 验证算例：Liu 2007 SP3、Kulikovsky transport 或 Shi/Luque radiation as applicable。
- 当前状态：registered_not_implemented
- 已知歧义：Whether diffusion, ion, or displacement current is included remains unresolved.

## EQ-0025 current moment

- 模块：radiation/current_moment
- 数学形式：$ICM=∫I dz$
- 变量定义：符号遵循来源；离散下标为网格节点，波浪号为 Fourier 变换。
- 每个变量的 SI 单位：长度 m；时间 s；密度 m^-3；电场 V/m；电势 V；电流 A；电流矩 A m；频率 Hz 或角频率 rad/s；能量 J。
- 方程左右量纲：已按 SI 核对；SP3 使用 Table 3 A_j 单位族，ESD 按角频率定义。
- 来源：SRC-0007
- 来源位置：model equations/current definitions
- provenance_type：direct
- 在 Shi 2019 主复刻中的作用：Stage 2 以后候选关系；Stage 1.1 仅登记。
- 是否直接实现：否
- 计划数值离散：待 A2 决策，本阶段不实现。
- 所需边界条件：按模块来源；未公开细节列入 unknowns。
- 所需初始条件：见 parameter_registry.csv。
- 验证算例：Liu 2007 SP3、Kulikovsky transport 或 Shi/Luque radiation as applicable。
- 当前状态：registered_not_implemented
- 已知歧义：Whether diffusion, ion, or displacement current is included remains unresolved.

## EQ-0026 collision current-moment increment

- 模块：radiation/current_moment
- 数学形式：$ΔICM=ICM_collision-ICM_isolated$
- 变量定义：符号遵循来源；离散下标为网格节点，波浪号为 Fourier 变换。
- 每个变量的 SI 单位：长度 m；时间 s；密度 m^-3；电场 V/m；电势 V；电流 A；电流矩 A m；频率 Hz 或角频率 rad/s；能量 J。
- 方程左右量纲：已按 SI 核对；SP3 使用 Table 3 A_j 单位族，ESD 按角频率定义。
- 来源：SRC-0007
- 来源位置：model equations/current definitions
- provenance_type：direct
- 在 Shi 2019 主复刻中的作用：Stage 2 以后候选关系；Stage 1.1 仅登记。
- 是否直接实现：否
- 计划数值离散：待 A2 决策，本阶段不实现。
- 所需边界条件：按模块来源；未公开细节列入 unknowns。
- 所需初始条件：见 parameter_registry.csv。
- 验证算例：Liu 2007 SP3、Kulikovsky transport 或 Shi/Luque radiation as applicable。
- 当前状态：registered_not_implemented
- 已知歧义：Whether diffusion, ion, or displacement current is included remains unresolved.

## EQ-0027 finite-antenna radiation

- 模块：radiation/current_moment
- 数学形式：$Bφ=sinθ/(4πε0c³R)∫∂I/∂t dz$
- 变量定义：符号遵循来源；离散下标为网格节点，波浪号为 Fourier 变换。
- 每个变量的 SI 单位：长度 m；时间 s；密度 m^-3；电场 V/m；电势 V；电流 A；电流矩 A m；频率 Hz 或角频率 rad/s；能量 J。
- 方程左右量纲：已按 SI 核对；SP3 使用 Table 3 A_j 单位族，ESD 按角频率定义。
- 来源：SRC-0001
- 来源位置：Shi 2019 PDF pp.4-6, Eqs.(2)-(6)
- provenance_type：direct
- 在 Shi 2019 主复刻中的作用：Stage 2 以后候选关系；Stage 1.1 仅登记。
- 是否直接实现：否
- 计划数值离散：待 A2 决策，本阶段不实现。
- 所需边界条件：按模块来源；未公开细节列入 unknowns。
- 所需初始条件：见 parameter_registry.csv。
- 验证算例：Liu 2007 SP3、Kulikovsky transport 或 Shi/Luque radiation as applicable。
- 当前状态：registered_not_implemented
- 已知歧义：None beyond cited sign convention.

## EQ-0028 far-field current-moment relation

- 模块：radiation/current_moment
- 数学形式：$Bφ∝∂ICM(t-R/c)/∂t$
- 变量定义：符号遵循来源；离散下标为网格节点，波浪号为 Fourier 变换。
- 每个变量的 SI 单位：长度 m；时间 s；密度 m^-3；电场 V/m；电势 V；电流 A；电流矩 A m；频率 Hz 或角频率 rad/s；能量 J。
- 方程左右量纲：已按 SI 核对；SP3 使用 Table 3 A_j 单位族，ESD 按角频率定义。
- 来源：SRC-0001
- 来源位置：Shi 2019 PDF pp.4-6, Eqs.(2)-(6)
- provenance_type：direct
- 在 Shi 2019 主复刻中的作用：Stage 2 以后候选关系；Stage 1.1 仅登记。
- 是否直接实现：否
- 计划数值离散：待 A2 决策，本阶段不实现。
- 所需边界条件：按模块来源；未公开细节列入 unknowns。
- 所需初始条件：见 parameter_registry.csv。
- 验证算例：Liu 2007 SP3、Kulikovsky transport 或 Shi/Luque radiation as applicable。
- 当前状态：registered_not_implemented
- 已知歧义：None beyond cited sign convention.

## EQ-0029 lifecycle current-moment model

- 模块：radiation/current_moment
- 数学形式：$ICM=ICM0 exp[α(t-T0)]/[1+exp((α+β)(t-T0))]$
- 变量定义：符号遵循来源；离散下标为网格节点，波浪号为 Fourier 变换。
- 每个变量的 SI 单位：长度 m；时间 s；密度 m^-3；电场 V/m；电势 V；电流 A；电流矩 A m；频率 Hz 或角频率 rad/s；能量 J。
- 方程左右量纲：已按 SI 核对；SP3 使用 Table 3 A_j 单位族，ESD 按角频率定义。
- 来源：SRC-0001
- 来源位置：Shi 2019 PDF pp.4-6, Eqs.(2)-(6)
- provenance_type：direct
- 在 Shi 2019 主复刻中的作用：Stage 2 以后候选关系；Stage 1.1 仅登记。
- 是否直接实现：否
- 计划数值离散：待 A2 决策，本阶段不实现。
- 所需边界条件：按模块来源；未公开细节列入 unknowns。
- 所需初始条件：见 parameter_registry.csv。
- 验证算例：Liu 2007 SP3、Kulikovsky transport 或 Shi/Luque radiation as applicable。
- 当前状态：registered_not_implemented
- 已知歧义：None beyond cited sign convention.

## EQ-0030 Fourier derivative relation

- 模块：radiation/current_moment
- 数学形式：$|F[dICM/dt]|=ω|ĨCM|$
- 变量定义：符号遵循来源；离散下标为网格节点，波浪号为 Fourier 变换。
- 每个变量的 SI 单位：长度 m；时间 s；密度 m^-3；电场 V/m；电势 V；电流 A；电流矩 A m；频率 Hz 或角频率 rad/s；能量 J。
- 方程左右量纲：已按 SI 核对；SP3 使用 Table 3 A_j 单位族，ESD 按角频率定义。
- 来源：SRC-0001
- 来源位置：Shi 2019 PDF pp.4-6, Eqs.(2)-(6)
- provenance_type：direct
- 在 Shi 2019 主复刻中的作用：Stage 2 以后候选关系；Stage 1.1 仅登记。
- 是否直接实现：否
- 计划数值离散：待 A2 决策，本阶段不实现。
- 所需边界条件：按模块来源；未公开细节列入 unknowns。
- 所需初始条件：见 parameter_registry.csv。
- 验证算例：Liu 2007 SP3、Kulikovsky transport 或 Shi/Luque radiation as applicable。
- 当前状态：registered_not_implemented
- 已知歧义：None beyond cited sign convention.

## EQ-0031 electric-field spectrum

- 模块：radiation/current_moment
- 数学形式：$Ẽ=sinθ(ωĨCM)/(4πε0c²R)$
- 变量定义：符号遵循来源；离散下标为网格节点，波浪号为 Fourier 变换。
- 每个变量的 SI 单位：长度 m；时间 s；密度 m^-3；电场 V/m；电势 V；电流 A；电流矩 A m；频率 Hz 或角频率 rad/s；能量 J。
- 方程左右量纲：已按 SI 核对；SP3 使用 Table 3 A_j 单位族，ESD 按角频率定义。
- 来源：SRC-0001
- 来源位置：Shi 2019 PDF pp.4-6, Eqs.(2)-(6)
- provenance_type：direct
- 在 Shi 2019 主复刻中的作用：Stage 2 以后候选关系；Stage 1.1 仅登记。
- 是否直接实现：否
- 计划数值离散：待 A2 决策，本阶段不实现。
- 所需边界条件：按模块来源；未公开细节列入 unknowns。
- 所需初始条件：见 parameter_registry.csv。
- 验证算例：Liu 2007 SP3、Kulikovsky transport 或 Shi/Luque radiation as applicable。
- 当前状态：registered_not_implemented
- 已知歧义：None beyond cited sign convention.

## EQ-0032 Poynting relation

- 模块：radiation/current_moment
- 数学形式：$S=cε0|Ẽ|²$
- 变量定义：符号遵循来源；离散下标为网格节点，波浪号为 Fourier 变换。
- 每个变量的 SI 单位：长度 m；时间 s；密度 m^-3；电场 V/m；电势 V；电流 A；电流矩 A m；频率 Hz 或角频率 rad/s；能量 J。
- 方程左右量纲：已按 SI 核对；SP3 使用 Table 3 A_j 单位族，ESD 按角频率定义。
- 来源：SRC-0001
- 来源位置：Shi 2019 PDF pp.4-6, Eqs.(2)-(6)
- provenance_type：direct
- 在 Shi 2019 主复刻中的作用：Stage 2 以后候选关系；Stage 1.1 仅登记。
- 是否直接实现：否
- 计划数值离散：待 A2 决策，本阶段不实现。
- 所需边界条件：按模块来源；未公开细节列入 unknowns。
- 所需初始条件：见 parameter_registry.csv。
- 验证算例：Liu 2007 SP3、Kulikovsky transport 或 Shi/Luque radiation as applicable。
- 当前状态：registered_not_implemented
- 已知歧义：None beyond cited sign convention.

## EQ-0033 energy spectral density

- 模块：radiation/current_moment
- 数学形式：$P=(ωĨCM)²/(6πε0c³)$
- 变量定义：符号遵循来源；离散下标为网格节点，波浪号为 Fourier 变换。
- 每个变量的 SI 单位：长度 m；时间 s；密度 m^-3；电场 V/m；电势 V；电流 A；电流矩 A m；频率 Hz 或角频率 rad/s；能量 J。
- 方程左右量纲：已按 SI 核对；SP3 使用 Table 3 A_j 单位族，ESD 按角频率定义。
- 来源：SRC-0001
- 来源位置：Shi 2019 PDF pp.4-6, Eqs.(2)-(6)
- provenance_type：direct
- 在 Shi 2019 主复刻中的作用：Stage 2 以后候选关系；Stage 1.1 仅登记。
- 是否直接实现：否
- 计划数值离散：待 A2 决策，本阶段不实现。
- 所需边界条件：按模块来源；未公开细节列入 unknowns。
- 所需初始条件：见 parameter_registry.csv。
- 验证算例：Liu 2007 SP3、Kulikovsky transport 或 Shi/Luque radiation as applicable。
- 当前状态：registered_not_implemented
- 已知歧义：None beyond cited sign convention.

## EQ-0034 frequency-band energy

- 模块：radiation/current_moment
- 数学形式：$W=2∫P dω$
- 变量定义：符号遵循来源；离散下标为网格节点，波浪号为 Fourier 变换。
- 每个变量的 SI 单位：长度 m；时间 s；密度 m^-3；电场 V/m；电势 V；电流 A；电流矩 A m；频率 Hz 或角频率 rad/s；能量 J。
- 方程左右量纲：已按 SI 核对；SP3 使用 Table 3 A_j 单位族，ESD 按角频率定义。
- 来源：SRC-0001
- 来源位置：Shi 2019 PDF pp.4-6, Eqs.(2)-(6)
- provenance_type：direct
- 在 Shi 2019 主复刻中的作用：Stage 2 以后候选关系；Stage 1.1 仅登记。
- 是否直接实现：否
- 计划数值离散：待 A2 决策，本阶段不实现。
- 所需边界条件：按模块来源；未公开细节列入 unknowns。
- 所需初始条件：见 parameter_registry.csv。
- 验证算例：Liu 2007 SP3、Kulikovsky transport 或 Shi/Luque radiation as applicable。
- 当前状态：registered_not_implemented
- 已知歧义：None beyond cited sign convention.

## EQ-0035 time-step restrictions

- 模块：fluid_model
- 数学形式：$Δt=min(advection,diffusion,dielectric,reaction restrictions); exact controller unresolved$
- 变量定义：符号遵循来源；离散下标为网格节点，波浪号为 Fourier 变换。
- 每个变量的 SI 单位：长度 m；时间 s；密度 m^-3；电场 V/m；电势 V；电流 A；电流矩 A m；频率 Hz 或角频率 rad/s；能量 J。
- 方程左右量纲：已按 SI 核对；SP3 使用 Table 3 A_j 单位族，ESD 按角频率定义。
- 来源：SRC-0007
- 来源位置：model equations/current definitions
- provenance_type：direct
- 在 Shi 2019 主复刻中的作用：Stage 2 以后候选关系；Stage 1.1 仅登记。
- 是否直接实现：否
- 计划数值离散：待 A2 决策，本阶段不实现。
- 所需边界条件：按模块来源；未公开细节列入 unknowns。
- 所需初始条件：见 parameter_registry.csv。
- 验证算例：Liu 2007 SP3、Kulikovsky transport 或 Shi/Luque radiation as applicable。
- 当前状态：registered_not_implemented
- 已知歧义：Exact controller and safety factors unpublished.
