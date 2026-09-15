# 文献证据维护

本目录只登记可追溯文献证据，不保存无明确再分发许可的论文全文。
`literature_registry.yaml` 是文献索引；`evidence_cards/` 记录经专家审阅后
对模型有实际影响的证据；`change_requests/` 保存经决策批准前的变更请求。

固定流程为：新论文 -> evidence registration -> expert review -> impact analysis
-> decision record -> optional change request -> minimal implementation -> targeted
regression -> global non-regression。禁止从论文标题直接推导代码变更。

