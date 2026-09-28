# 固定版本 MuSiC 运行环境

核查与配置日期：2026-09-17。**完整官方 MuSiC 包已安装且合成技术控制通过；后续固定协议的真实数据模拟比较也已完成**，见 `research/MUSIC_COMPARISON_REPORT.md`。以下环境配置与技术控制记录保留其执行时事实；不是自写替代函数。当前研究状态以 `STATUS.md` 为准。

## 调用

现有 R 位于 `C:/Program Files/R/R-4.5.2/bin/Rscript.exe`，版本为 R 4.5.2 (2025-10-31 ucrt)，x86_64 Windows。初始 PATH 未找到 R 不表示电脑未安装；本轮只读检查常见路径发现它，没有安装新的 R、改 PATH、修改注册表或写入用户既有 R 包库。

在项目根目录运行：

```powershell
$env:LC_ALL = 'C'
$env:LANG = 'C'
$env:R_LIBS_USER = (Join-Path (Get-Location) '.tools/R-library')
$env:R_LIBS_SITE = $env:R_LIBS_USER
& 'C:/Program Files/R/R-4.5.2/bin/Rscript.exe' --vanilla environment/music/smoke.R
```

这些环境变量只作用于当前 shell 及子进程，不保存为系统配置。分析 R 脚本开头用 `source("environment/music/runtime.R")`；它只把项目 `.tools/R-library` 和 R 自带基础库放入 `.libPaths()`，并加载真实 MuSiC 与 SingleCellExperiment 包。无需 `devtools`，无需 Rtools。

## 官方来源和精确版本

- MuSiC 1.0.0，官方仓库固定 commit `f21fe67f5670d5e9fca0ad7550abaae3423eb59c`，原始 GPL >= 3 许可；从固定 [GitHub commit archive](https://codeload.github.com/xuranw/MuSiC/zip/f21fe67f5670d5e9fca0ad7550abaae3423eb59c) 获取，不是浮动 master。
- 原始 archive 大小 68,200,648 字节，SHA256 `b092e1c17dc4ebe4c90af9729d40bdc1bcc08e74e1e7792d3d67c9b3bf4c0126`。
- 原始源码目录：`data/raw/MuSiC_runtime/20260917T090134491841Z/MuSiC-f21fe67f5670d5e9fca0ad7550abaae3423eb59c`。`R/utils.R` 定义 `music_prop`；`R/construct.R` 定义 `music_basis`；`R/analysis.R` 定义 `music.iter`、`music.basic`。`source_manifest.json` 和 `runtime_manifest.json` 保存 URL、时间和 SHA256。
- Bioconductor 固定 3.22，以适配既有 R 4.5；没有把 R 4.6 的当前 Bioconductor 3.23 包混入。
- 核心实际版本：nnls 1.6、TOAST 1.24.0、Biobase 2.70.0、SingleCellExperiment 1.32.0、Matrix 1.7-4、MCMCpack 1.7-1、ggplot2 4.0.3。全部安装版本/位置/Built 记录见 `music/installed_packages.tsv`；实际已加载版本见 `music/session_info.txt`。
- 74 个 CRAN/Bioconductor 依赖全部使用官方 Windows 二进制包安装于项目库，加官方源码安装的 MuSiC 共 75 个项目包。R 自带 Matrix 等基础推荐包保持现有安装不变。
- 原安装过程在 R 临时目录下载包；之后专门把相同版本的全部 74 个 binary archive 缓存至 `data/raw/MuSiC_runtime/20260917T090134491841Z/dependency_binaries`（87,382,423 字节）。`music/snapshot_runtime.py` 逐一验证缓存内 2,879 个文件与实际安装文件字节相同，并保存 archive URL/SHA256。缓存是**安装后获取并比对**，不是声称首次安装使用了此缓存路径。
- 部分包 Built 字段为 R 4.5.3，在 R 4.5.2 下发出补丁版本警告；原始警告保留在安装/smoke 日志，加载与以下技术控制实际通过。不要把环境描述为所有包都在 4.5.2 编译。

参考：官方 [MuSiC tutorial](https://xuranw.github.io/MuSiC/articles/MuSiC.html)、[Bioconductor installation guidance](https://bioconductor.org/install/)、[R Windows FAQ](https://cran.r-project.org/bin/windows/base/rw-FAQ.html)。版本以本地锁定记录为准，未来网页变化不能静默升级本项目环境。

## 输入接口与重要尺度

```r
result <- MuSiC::music_prop(
  bulk.mtx = bulk_counts,       # gene x multiple target columns
  sc.sce = reference_sce,      # full-gene original counts assay
  markers = NULL,             # final comparison: no preselected marker panel
  clusters = "cell_type",     # colData column name
  samples = "donor",          # colData column name
  select.ct = fixed_cell_types,
  cell_size = NULL, ct.cov = FALSE, verbose = FALSE,
  iter.max = 1000, nu = 0.0001, eps = 0.01,
  centered = FALSE, normalize = FALSE
)
```

`reference_sce` 是 `SingleCellExperiment::SingleCellExperiment(assays=list(counts=...))`，行名为统一的基因 ID，列名为唯一 cell ID，`colData` 与列严格对齐。参考 counts 保留完整基因，以便库大小按完整计数计算；最终比较用 `markers=NULL`。若另一个版本明确研究marker面板，应通过 `markers` 参数而不是截断SCE，且需要另版协议，不能改变本次运行。bulk 保持 gene x 多目标矩阵，官方代码对单一目标的部分索引没有 `drop=FALSE`，因此研究批次至少两列。

返回的 `Est.prop.weighted` 与 `Est.prop.allgene` 都是 target x cell-type 矩阵，分别为 MuSiC 加权估计及其同设计矩阵 NNLS。两者均使用 MuSiC 的供者汇总和细胞大小处理；`allgene` 这个历史输出名不表示无视传入 marker，也不能把它称为旧试验中的 cell-pooled NNLS。输出中的 `Var.prop` 不是本项目供者总体不确定性或抽样 SD。这里仅记录 API，不替代正式科学协议。

## 已执行技术检查

`music/smoke.R` 用 60 个基因、6 类人工线性独立表达谱、3 位人工供者、72 个人工参考细胞、3 个人工混合目标测试完整包。**全部是代数控制，不是真实生物学数据或研究发现。**

- 从固定官方 R 文件重新解析函数，`music_prop`、`music_basis`、`music.iter`、`music.basic` 的函数体与安装包逐一相同。
- weighted 与 NNLS 相对已知组成的最大绝对比例误差均为 `4.440892098500626e-16`。
- 两次执行 weighted 的最大绝对差为 0；所有比例有限、非负，最大行和误差为 0。
- 小型完整拟合用时约 0.16 秒；这不是实际 PBMC 全网格的运行时间预估。
- 结果见 `music/smoke.json`；完整控制代码、warning 和会话信息均保留。

## 重建与版本边界

`music/install_dependencies.R` 和安装日志记录本次最初联网安装步骤；该脚本选择当时仓库包，不是未来严格锁版本安装器。严格重建应先验证 `music/runtime_manifest.json` 中每个缓存 archive 的 SHA256，然后以 `install.packages(cached_zip_paths, repos=NULL, type="win.binary", lib=project_library)` 安装全部缓存，最后在相同 R 4.5.2 环境执行：

```powershell
& 'C:/Program Files/R/R-4.5.2/bin/R.exe' CMD INSTALL --library=.tools/R-library --no-multiarch data/raw/MuSiC_runtime/20260917T090134491841Z/MuSiC-f21fe67f5670d5e9fca0ad7550abaae3423eb59c
```

本地已有环境不必重装。若搬迁电脑或升级 R，请建立新版本环境与新 manifest，重跑技术控制，不能覆盖本次版本记录。原有 R 的执行文件与 x64 DLL SHA256 也记录于 `runtime_manifest.json`，没有复制或修改系统安装。
