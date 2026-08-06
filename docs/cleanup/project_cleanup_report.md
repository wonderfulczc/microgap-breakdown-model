# Project cleanup report

## 1. 清理目标

本次清理目标是在不影响 Stage 1 至 Stage 5 最终结果、证据链、源码、配置、测试、验证器和后续复用能力的前提下，归档失效/中断运行并删除缓存、临时文件和无效展开目录。

## 2. 清理原则

执行原则为先审计、后删除；不确定文件默认保留；被 closure report、validation matrix、run registry、derived artifact manifest、最终图表或 HTML 总结引用的文件不得删除；清理先进入 quarantine，验证通过后才 purge。

## 3. 清理前文件数

安全审计记录的清理前文件数为 13676。

## 4. 清理后文件数

清理、二次测试缓存清理、最终 `.pyc` 缓存清理并生成总结文档后的最终文件数为 8675。

## 5. 清理前磁盘占用

安全审计记录的清理前磁盘占用为 5.5G。

## 6. 清理后磁盘占用

最终磁盘占用为 5.1G。

## 7. 节省空间

多轮清理 manifest 候选原始体积为 1,232,339,811 bytes；此外最终删除测试和验证过程中再次生成的 `.pyc` 缓存文件 354 个。审计归档压缩包保留体积为 228,150,815 bytes。按 `du -sh` 粗略统计，项目目录占用减少约 0.4G；差异来自归档保留、新增总结图文文件以及文件系统块大小和四舍五入。

## 8. 删除文件分类统计

| cleanup_pass | classification | final_action | count | size_bytes |
| --- | --- | --- | --- | --- |
| post_validation_cache | cache | quarantined | 979 | 28268506 |
| post_validation_cache | temporary | quarantined | 2 | 58239 |
| primary | audit_only | quarantine_directory | 4 | 1091553299 |
| primary | cache | quarantined | 4904 | 112208499 |
| primary | temporary | quarantined | 2 | 251268 |

manifest 登记的单文件缓存/临时候选 5887 个，最终补充删除 `.pyc` 缓存 354 个；audit-only 展开目录 4 个，归档 manifest 内部文件 143 个。

## 9. 归档文件统计

| archive | original_dirs | size_bytes | sha256 | reason | replacement | created |
| --- | --- | --- | --- | --- | --- | --- |
| archive/audit_only/stage2_invalidated_results.tar.zst | results/stage2/forensic/invalidated_results | 108703 | be9b54f5a80923a4e812eed4df86339fd976a524fb32f96af5aaf07141d8d4b0 | Stage 2 invalidated forensic results; not valid final evidence. | Stage 2 closure evidence and forensic summary | True |
| archive/audit_only/stage4_interrupted_runs.tar.zst | results/stage4/runs/right_isolated_20um_paused_20260723 | 9371664 | f3263c3c88591f548406085c46ee7933993d280af7d94753671ddee10e1849f4 | Interrupted/corrupt/paused run directories replaced by formal evidence. | S4-RIGHT-ISOLATED and other formal Stage 4 runs | True |
| archive/audit_only/stage5_interrupted_runs.tar.zst | results/stage5/runs/highfield_right_isolated_20um_interrupted_corrupt_20260728;results/stage5/runs/highfield_collision_20um_partial_hung_20260727 | 541036466 | 779a49e39d66eb57416cce5d031329183c2c51f5eae0d1f4ec2ddf75626b67d8 | Interrupted/corrupt/paused run directories replaced by formal evidence. | S5 formal reruns and resource-aware terminations | True |

归档目录：`archive/audit_only/`。每个归档均包含 `.sha256`、`.manifest.csv` 和 `.README.md`。

## 10. 保留的 uncertain 文件

`results/cleanup/retained_uncertain.csv` 中保留 uncertain 文件 116 个。它们未被删除，后续如需进一步清理必须人工确认。

## 11. 正式 run_id 完整性

Stage 1 至 Stage 5 的正式 run registry、validation matrix、derived artifact manifest 和 closure report 仍保留。`tools/validate_cleanup_integrity.py` 已检查登记结果路径存在且 SHA256 匹配。

## 12. Stage 1 至 Stage 5 回归结果

清理后验证记录见 `results/cleanup/final_validation.log`。结果：构建通过、ctest 通过、pytest 通过、Stage 1 closure 通过、Stage 2 closure 通过、Stage 3 resource-aware validator 通过、Stage 4 closure validator 通过、Stage 5 closure validator 通过、cleanup integrity validator 通过、simulation summary validator 通过。二次和最终缓存清理后又以 `PYTHONDONTWRITEBYTECODE=1` 运行 cleanup integrity 和 simulation summary validator，均通过。

## 13. SHA256 验证

run registry 和 derived artifact manifest 中登记的正式结果 SHA256 已由 cleanup integrity validator 检查通过。审计归档 SHA256 记录在 `archive/audit_only/*.sha256`。

## 14. 最终项目目录结构

最终目录结构说明见 `docs/final_project_structure.md`。

## 15. 恢复审计归档的方法

审计归档只用于历史追溯，不得用于最终结果。恢复时可执行：

```bash
tar --zstd -xf archive/audit_only/stage4_interrupted_runs.tar.zst -C /tmp/recovered_stage4_audit
```

若系统 tar 不支持 zstd，可使用 Python 3.14+ 的 `tarfile.open(path, "r:zst")` 解压。恢复后必须对照对应 `.manifest.csv` 和 `.sha256` 校验。

## 16. 清理日期

2026-07-29 03:53:54 UTC
