# 第三方依赖说明

状态：`THIRD_PARTY_LICENSE_REVIEW = PASS`。

本项目自身软件代码采用 Apache-2.0；这是项目许可证，不是下列第三方后端的许可证。

本项目 wheel 和 sdist 不包含第三方依赖源码或二进制。Python 依赖由包管理器独立安装；
PETSc 与 Open MPI 是系统依赖；Afivo-streamer、openEMS、CSXCAD、fparser、AppCSXCAD
与 COMSOL 均保持外部边界。若未来分发预装二进制环境或外部后端，必须重新审查许可证
正文、NOTICE、源码提供及动态链接义务。当前 wheel/sdist 不分发 Afivo、openEMS、
CSXCAD、fparser 或 AppCSXCAD 的源码或二进制；若分发模式改变，必须重新审计。

| 依赖 | 固定版本/提交 | 许可证 | 当前关系 |
|---|---|---|---|
| PETSc | 3.24.4 | BSD-2-Clause | 系统动态链接，不随包分发 |
| Open MPI | 5.0.10 | BSD-style Open MPI license | 系统动态链接/命令调用，不随包分发 |
| NumPy | 2.5.1 | BSD-3-Clause 等复合宽松许可证 | PyPI 运行依赖 |
| SciPy | 1.18.0 | BSD-3-Clause，含组件许可证 | PyPI 运行依赖 |
| pandas | 3.0.3 | BSD-3-Clause | PyPI 运行依赖 |
| Matplotlib | 3.11.1 | Matplotlib License，含数据许可证 | PyPI 运行依赖 |
| PyYAML | 6.0.3 | MIT | PyPI 运行依赖 |
| h5py | 3.16.0 | BSD-3-Clause | PyPI 运行依赖 |
| pytest | 9.1.1 | MIT | 开发依赖 |
| Afivo-streamer | `a50b550...` | GPL-3.0，only/or-later 未由顶层证据明确 | 外部进程与文件接口 |
| openEMS | `8f480d0...` | GPL-3.0-or-later | 外部 Python 运行时/动态链接/文件接口 |
| CSXCAD | `0458ee1...` | LGPL-3.0-or-later | 外部 Python 运行时与动态链接 |
| fparser | `4b9c845...` | LGPL-3.0-only | openEMS 的外部动态链接依赖 |
| AppCSXCAD | `9249ab7...` | GPL-3.0-or-later | 可选 GUI，当前未安装/未分发 |

fparser 的证据为固定源码提交中 `docs/fparser.html` 的 Usage license 声明及同一提交跟踪的
`docs/lgpl.txt`。逐项证据、分发模式和风险见
`packaging/third_party_license_inventory.json`。本文件不是法律意见。
