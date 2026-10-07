"""
CLI entry point for BioSyn qPCR CSV Analyzer.

Usage:
    python -m src.main data.csv --output results.xlsx
    python -m src.main data.xlsx --reference-gene GAPDH --calibrator-group Control
"""

from __future__ import annotations

import argparse
import csv
import logging
import sys
from pathlib import Path
from typing import List, Dict, Any

from .analyzer import qPCRAnalyzer, qPCRRecord
from .exporter import export_qpcr_results

logger = logging.getLogger(__name__)


def setup_logging(verbose: bool = False) -> None:
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )


def read_qpcr_file(path: Path) -> List[qPCRRecord]:
    """Read qPCR data from CSV / XLSX / TSV."""
    records: List[qPCRRecord] = []

    if path.suffix.lower() in (".xlsx", ".xls"):
        try:
            import openpyxl
            wb = openpyxl.load_workbook(path, read_only=True)
            ws = wb.active
            rows = list(ws.iter_rows(values_only=True))
            if not rows:
                return records
            headers = [str(h).strip() if h else "" for h in rows[0]]
            for row in rows[1:]:
                rec = dict(zip(headers, row))
                records.append(_row_to_record(rec))
        except ImportError:
            raise RuntimeError("openpyxl required for XLSX. pip install openpyxl")
    else:
        delim = "\t" if path.suffix.lower() == ".tsv" else ","
        with open(path, "r", encoding="utf-8", errors="replace") as fh:
            reader = csv.DictReader(fh, delimiter=delim)
            for row in reader:
                records.append(_row_to_record({k: (v or "").strip() for k, v in row.items()}))

    logger.info("Read %d qPCR records from %s", len(records), path)
    return records


def _row_to_record(row: Dict[str, Any]) -> qPCRRecord:
    """Convert a raw row dict to a qPCRRecord, with flexible column mapping."""
    def get(*keys: str, default: str = "") -> str:
        for k in keys:
            for actual in row:
                if actual.lower() == k.lower():
                    return str(row[actual])
        return default

    ct_raw = get("Ct", "ct", "CT", "Threshold Cycle", "Cp", default="")
    invalid_ct_values = {"undetermined", "nan", "na", ""}
    try:
        ct = float(ct_raw) if ct_raw and ct_raw.lower() not in invalid_ct_values else None
    except ValueError:
        ct = None

    rep_raw = get("Replicate", "rep", "replicate", default="1")
    try:
        rep = int(rep_raw)
    except ValueError:
        rep = 1

    return qPCRRecord(
        sample=get("Sample", "sample", "Sample_Name", default="Unknown"),
        group=get("Group", "group", "Condition", "Treatment", default=""),
        gene=get("Gene", "gene", "Target", "Assay", default="Unknown"),
        well=get("Well", "well", "Position", default=""),
        ct=ct,
        replicate=rep,
        reference_gene=get("Reference_Gene", "Ref_Gene", default=""),
        notes=get("Notes", "notes", default=""),
        raw_row=row,
    )


def run_analysis(
    input_path: Path,
    output_path: Path,
    reference_gene: str = "",
    calibrator_group: str = "",
    cv_threshold: float = 1.0,
) -> Dict[str, Any]:
    """Run the full qPCR analysis pipeline."""
    records = read_qpcr_file(input_path)
    if not records:
        raise ValueError("No records found in input file.")

    analyzer = qPCRAnalyzer(
        reference_gene=reference_gene,
        calibrator_group=calibrator_group,
        cv_warn_threshold=cv_threshold,
    )

    # Step 1: Clean
    cleaned, exclusions = analyzer.clean_data(records)

    # Step 2: Replicate QC
    rep_stats = analyzer.replicate_qc(cleaned)

    # Step 3: ΔCt
    dct_results = analyzer.compute_dct(rep_stats)

    # Step 4: ΔΔCt
    ddct_results = analyzer.compute_ddct(dct_results)

    # Step 5: Export
    output_path.parent.mkdir(parents=True, exist_ok=True)
    export_qpcr_results(
        records, cleaned, rep_stats, dct_results, ddct_results,
        exclusions, output_path,
    )

    summary = {
        "total_wells": len(records),
        "excluded": len(exclusions),
        "replicate_groups": len(rep_stats),
        "dct_pairs": len(dct_results),
        "ddct_results": len(ddct_results),
    }
    logger.info("qPCR analysis complete: %s", summary)
    return summary


def main(argv: List[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="biosyn-qpcr",
        description="BioSyn qPCR CSV Analyzer — local qPCR data analysis.",
        epilog="Built by BioSyn Inc. — https://biosyn-inc.com/",
    )
    parser.add_argument("input", type=Path, help="Input CSV / XLSX / TSV file")
    parser.add_argument("--output", "-o", type=Path,
                        default=Path("qpcr_results.xlsx"))
    parser.add_argument("--reference-gene", "-r", default="",
                        help="Reference / housekeeping gene name")
    parser.add_argument("--calibrator-group", "-c", default="",
                        help="Calibrator / control group name")
    parser.add_argument("--cv-threshold", type=float, default=1.0,
                        help="CV%% warning threshold (default: 1.0)")
    parser.add_argument("--verbose", "-v", action="store_true")
    args = parser.parse_args(argv)

    setup_logging(args.verbose)

    try:
        summary = run_analysis(
            args.input, args.output,
            reference_gene=args.reference_gene,
            calibrator_group=args.calibrator_group,
            cv_threshold=args.cv_threshold,
        )
        print(f"\nqPCR analysis complete: {summary}")
        print(f"Results saved to: {args.output.resolve()}")
        return 0
    except Exception as exc:
        logger.error("Fatal error: %s", exc)
        print(f"Error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
