# Environment authority

The validated core runtime was R 4.5.1. `renv.lock` is preserved exactly but covers only part of the direct dependency set used across the current workflow. Do not treat it as a complete one-file restoration contract.

Exact version evidence is available for the R dependencies not represented in the lock: MCPcounter 1.2.0, Seurat 5.5.1, SeuratObject 5.4.0, UCell 2.14.0, cowplot 1.2.0, future 1.75.0, patchwork 1.3.2, ragg 1.5.2, readxl 1.4.5, systemfonts 1.3.2, tiff 0.1-12, and zip 2.3.3. The final Figure 6 script also uses `png` and `readr`; their exact versions were not established by the retained lock/session evidence. `estimate` occurs only in an optional, nonexecuted branch.

Final supplementary-figure engineering uses Python packages including Pillow, matplotlib, NumPy, pandas, and python-docx. Exact package versions were not frozen. The released figures, source workbooks, generating scripts, and hashes are the final authority.

Status: `RESTORABLE_ENVIRONMENT = PARTIAL`; `ENVIRONMENT_DOCUMENTATION = PASS`.
