# Changelog

## [1.0.0] - 2024-01-01
### Added
- CSV / XLSX / TSV import with flexible column mapping
- Data cleaning: blanks, missing Ct, undetermined, duplicate wells
- Replicate QC: n, mean, SD, CV%, min, max, range with PASS/REVIEW/OUTLIER
- ΔCt = Ct(target) − Ct(reference)
- ΔΔCt = ΔCt(sample) − ΔCt(calibrator)
- Fold Change = 2^(−ΔΔCt)
- Standard curve: slope, intercept, R², amplification efficiency
- 96/384 well plate support
- XLSX export (10 sheets) and HTML report
- Raw data immutable — all exclusions logged with reasons
- 100% offline operation
