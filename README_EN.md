# BioSyn qPCR CSV Analyzer

> **Local-first qPCR data analysis** — clean instrument exports, replicate QC, ΔCt / ΔΔCt / fold change, standard curves. All offline.

## Overview

qPCR instruments export messy CSVs. Computing ΔCt, ΔΔCt, and fold change manually in Excel is slow and error-prone. **BioSyn qPCR CSV Analyzer** automates the full pipeline: import → clean → replicate QC → ΔCt → ΔΔCt → export.

> For online qPCR calculators, visit [BioSyn qPCR analysis tools](https://biosyn-inc.com/tools/qpcr-analysis) — this is the desktop batch-processing companion.

## Quick Start

```bash
git clone https://github.com/your-org/biosyn-qpcr-csv-analyzer.git
cd biosyn-qpcr-csv-analyzer
pip install -r requirements.txt

python -m src.main examples/sample_qpcr_data.csv --output results.xlsx
python -m src.main data.xlsx --reference-gene GAPDH --calibrator-group Control
```

## Key Methods

| Metric | Formula |
|--------|---------|
| ΔCt | Ct(target) − Ct(reference) |
| ΔΔCt | ΔCt(sample) − ΔCt(calibrator) |
| Fold Change | 2^(−ΔΔCt) — assumes ~100% efficiency |
| Efficiency | (10^(−1/slope) − 1) × 100% |
| Replicate QC | CV% threshold + Z-score outlier detection |

## Privacy

100% offline. Raw data is never modified or uploaded.

## Keywords

`qPCR analysis`, `ΔΔCt calculator`, `fold change calculation`, `replicate QC`, `standard curve efficiency`, `Ct data analysis`, `qPCR CSV analyzer`, `real-time PCR`, `gene expression`

## About

Built by [BioSyn Inc.](https://biosyn-inc.com/)

## License

[MIT](LICENSE)
