# 第三方依赖说明

本项目源码包声明但不 vendoring Python 运行依赖；其许可证与 NOTICE 由各自发行包提供。
本项目也不分发 Afivo-streamer、openEMS、CSXCAD、fparser、AppCSXCAD 或 COMSOL 源码/二进制。

已核验的主要许可证包括：PETSc `BSD-2-Clause`；Open MPI 的 BSD-style Open MPI
license；NumPy 的复合宽松许可证表达式；SciPy/pandas/h5py 的 BSD 类许可证；PyYAML、
pytest 和 PyPA build 的 `MIT`；openEMS/Afivo/AppCSXCAD 的 `GPL-3.0`；CSXCAD 的
`LGPL-3.0`。Matplotlib 采用其 PSF-compatible license，并包含字体/样式附加许可证。

fparser 在本地固定源码树中未发现顶层许可证文件，因此标记为
`UNKNOWN_REQUIRES_REVIEW`。完整逐项来源与分发边界见
`packaging/third_party_license_inventory.json`。当前状态：
`THIRD_PARTY_NOTICES_STATUS=INCOMPLETE_REVIEW_REQUIRED`。
