# BioSyn qPCR CSV Analyzer

> **qPCR 数据分析工具** — 把仪器导出的 CSV/XLSX 一键清洗、做重复样 QC、计算 ΔCt / ΔΔCt / 相对表达量，全程本地离线。
>
> A local-first qPCR data analyzer: clean instrument exports, replicate QC, ΔCt / ΔΔCt / fold change, standard curves — all offline.

[![Python](https://img.shields.io/badge/Python-3.10+-blue.svg)](https://www.python.org/)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Offline](https://img.shields.io/badge/Privacy-Offline-orange.svg)](#privacy)

## 核心功能

| 功能 | 说明 |
|------|------|
| **数据清洗** | 自动检测 blank、missing Ct、undetermined、重复孔、异常值，原始数据不覆盖 |
| **Replicate QC** | n、mean Ct、SD、CV%、min、max、range，PASS / REVIEW / OUTLIER 三级判定 |
| **ΔCt** | `ΔCt = Ct(target) - Ct(reference)`，自动识别内参基因 |
| **ΔΔCt & Fold Change** | `ΔΔCt = ΔCt(sample) - ΔCt(calibrator)`，`Fold Change = 2^(-ΔΔCt)` |
| **标准曲线** | slope、intercept、R²、扩增效率（`Eff% = (10^(-1/slope) - 1) × 100`） |
| **96/384 孔板** | 自动识别 Well 位置，支持孔板布局检查 |
| **完整导出** | XLSX（10 个 Sheet）、HTML 报告 |
| **数据可追溯** | 每一步保留中间数据，排除数据记录原因 |

## 解决的痛点

qPCR 仪器导出的 CSV 数据格式各异，手动在 Excel 里算 ΔCt、ΔΔCt、fold change 费时且容易出错。BioSyn qPCR CSV Analyzer 把整个分析流程自动化：导入数据 → 自动清洗 → 重复样 QC → ΔCt / ΔΔCt → 导出报告。

> 在线 qPCR 计算工具请访问 [BioSyn qPCR 分析](https://biosyn-inc.com/tools/qpcr-analysis) — 本工具是其桌面端批量数据处理增强版。

## 技术栈

- Python 3.10+
- **numpy** — 标准曲线线性回归计算
- **openpyxl** — XLSX 读取与多 Sheet 报告导出

## 依赖与安装

### 环境要求

- **Python 3.10+**

### 依赖包

| 包 | 用途 | 是否必需 |
|------|------|----------|
| `numpy` | 标准曲线线性回归 | 必需 |
| `openpyxl` | XLSX 输入读取 + XLSX 报告导出 | 必需 |
| `scipy` | 精确统计检验（MVP 使用标准库近似） | 可选 |
| `matplotlib` | 图表生成（计划中） | 可选 |

> 数据清洗、重复样 QC、ΔCt/ΔΔCt 核心逻辑由标准库实现。

### 安装指令

```bash
# 1. 克隆仓库
git clone https://github.com/your-org/biosyn-qpcr-csv-analyzer.git
cd biosyn-qpcr-csv-analyzer

# 2. 安装依赖
pip install -r requirements.txt

# 3. 基本分析
python -m src.main examples/sample_qpcr_data.csv --output results.xlsx

# 指定内参基因和对照
python -m src.main data.xlsx --reference-gene GAPDH --calibrator-group Control

# 自定义 CV 警告阈值
python -m src.main data.csv --cv-threshold 0.5
```

## 输入格式

CSV / XLSX / TSV 表格，支持灵活的列名映射：

| Sample | Group | Gene | Well | Ct | Replicate |
|--------|-------|------|------|----|-----------|
| Sample_A1 | Treatment | GAPDH | A1 | 22.15 | 1 |

## 输出 XLSX Sheet

1. **Raw Data** — 原始数据（只读）
2. **Clean Data** — 清洗后数据
3. **Replicate QC** — 重复样统计与 QC 状态
4. **ΔCt** — 每个样本-基因对的 ΔCt
5. **Fold Change** — ΔΔCt 与 2^(-ΔΔCt)
6. **Errors** — 被排除的数据及原因
7. **Methods** — 计算方法与免责声明

## 科学边界

> 软件分析的是用户提供的数据。QC 阈值由用户/实验设计决定。计算结果不替代实验质量控制，不自动证明生物学显著性。统计检验与 fold-change 是不同概念，不应混为一谈。

## Keywords

`qPCR analysis`, `ΔΔCt calculator`, `fold change calculation`, `replicate QC`, `standard curve efficiency`, `Ct data analysis`, `qPCR CSV analyzer`, `real-time PCR analysis`, `gene expression quantification`

## About

By [BioSyn Inc.](https://biosyn-inc.com/) — explore more [qPCR data analysis tools](https://biosyn-inc.com/tools/qpcr-analysis).

## License

[MIT](LICENSE)
