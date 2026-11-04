# getBio

RNA-seq（bulk + 单细胞）分析常用函数集合，由个人脚本 `get_mamin.r` 整理成 R 包。

## 安装

```r
# 依赖中有 Bioconductor 包，建议先装 BiocManager
install.packages(c("remotes", "BiocManager"))
remotes::install_github("YOUR_GITHUB_USERNAME/getBio",
                        repos = BiocManager::repositories(),
                        dependencies = TRUE)
```

## 使用

```r
library(getBio)
tpm <- get_tpm(fpkm_matrix)
design <- Get_design(tpm, group_list)
```

## 函数一览

| 文件 | 内容 | 函数 |
|---|---|---|
| `R/utils.R` | 表达矩阵预处理、TPM、基因查找等通用工具 | `get_tpm`, `Top`, `Get_Med`, `Get_logic`, `Chect_gene`, `get_Normalize`, `get_rankMatrix`, `get_reannotated`, `get_numMatrix`, `get_genes`, `get_correlated_genes` |
| `R/heatmap.R` | 分组/样本热图 | `get_groupheatmap`, `get_heatmap` |
| `R/differential.R` | limma 差异分析（含 limma-trend / voom） | `Get_design`, `get_diffecence`, `run_limmatrend`, `run_voomlimma` |
| `R/gsea_gsva.R` | GSEA、ssGSEA/GSVA、GMT 基因集 | `get_gsea`, `get_gsva`, `get_immune_gsva`, `get_broadinstituteORG`, `GeneSet`, `createGeneSetCollection` |
| `R/plots.R` | 相关性、PCA、火山图、小提琴图、散点图、Venn 图 | `get_correlation`, `get_pca`, `get_volcano`, `get_violinplot`, `Get_Scatterplot2`, `get_vennplot` |
| `R/enrichment.R` | GO / KEGG 富集与作图 | `myGO_KEGG`, `myGO_KEGG2`, `get_GoKegg`, `get_goseletCSV`, `get_gobarplot`, `get_gopathway` |
| `R/survival.R` | 生存相关基因 | `Get_associated_gene` |
| `R/seurat.R` | Seurat 单细胞流程（降维、marker、伪细胞、CellChat、monocle3 等） | `GatherData.Seurat`, `PseudoCell`, `get_Cellchat`, `get_Topmarkers`, `modify_vlnplot`, `StackedVlnPlot`, `get_newReduction`, `get_localReduction`, `myfindmarkers`, `get_Identreplace`, `Get_cytometryiod`, `get_pseudotime`, `get_harmony`, `get_seuratGO`, `get_HTOassay`, `get_percentage`, `myfindmarkers2`, `get_seuratgroupGO`, `get_predcelltype` |

另外导出：颜色向量 `color_used`、`color_crc`，S4 类 `GeneSet`。

## 说明

- 函数体基本保持原脚本不变，仅做了拆分和打包；原脚本里的 `library()` 保留在各函数内。
- 原脚本中 `get_gsea` 定义了两次，只保留了后一个（带 `datatype` 参数）。
- 非 CRAN 依赖（`idmap3`、`CellChat`、`monocle3`）放在 Suggests，用到对应函数时再安装。
