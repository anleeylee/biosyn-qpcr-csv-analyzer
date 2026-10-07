"""
qPCR analysis engine — data cleaning, replicate QC, ΔCt / ΔΔCt, standard curve.

Raw data is never modified. Every step produces new intermediate data.
"""

from __future__ import annotations

import math
import statistics
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np


@dataclass
class qPCRRecord:
    """A single qPCR well/observation."""
    sample: str
    group: str
    gene: str
    well: str
    ct: Optional[float]
    replicate: int = 1
    reference_gene: str = ""
    notes: str = ""
    raw_row: Dict = field(default_factory=dict)


@dataclass
class CleanRecord:
    """Record after cleaning — includes exclusion reason."""
    record: qPCRRecord
    excluded: bool = False
    exclusion_reason: str = ""


@dataclass
class ReplicateStats:
    """Statistics for a replicate group."""
    sample: str
    gene: str
    group: str
    n: int = 0
    mean_ct: float = 0.0
    sd: float = 0.0
    cv: float = 0.0
    min_ct: float = 0.0
    max_ct: float = 0.0
    range_ct: float = 0.0
    status: str = "PASS"  # PASS | REVIEW | OUTLIER
    outliers: List[float] = field(default_factory=list)


@dataclass
class DCtResult:
    """ΔCt result for a sample-gene pair."""
    sample: str
    group: str
    target_gene: str
    reference_gene: str
    mean_target_ct: float = 0.0
    mean_ref_ct: float = 0.0
    delta_ct: float = 0.0


@dataclass
class DDCtResult:
    """ΔΔCt and fold change result."""
    sample: str
    group: str
    gene: str
    delta_ct: float = 0.0
    delta_delta_ct: float = 0.0
    fold_change: float = 0.0
    is_calibrator: bool = False


@dataclass
class StandardCurve:
    """Standard curve fit results."""
    gene: str = ""
    concentrations: List[float] = field(default_factory=list)
    cts: List[float] = field(default_factory=list)
    slope: float = 0.0
    intercept: float = 0.0
    r_squared: float = 0.0
    efficiency_pct: float = 0.0


class qPCRAnalyzer:
    """Main qPCR analysis engine."""

    def __init__(
        self,
        cv_warn_threshold: float = 1.0,   # %
        outlier_z_threshold: float = 2.0,
        reference_gene: str = "",
        calibrator_group: str = "",
    ) -> None:
        self.cv_warn_threshold = cv_warn_threshold
        self.outlier_z_threshold = outlier_z_threshold
        self.reference_gene = reference_gene
        self.calibrator_group = calibrator_group

    def clean_data(self, records: List[qPCRRecord]) -> Tuple[List[CleanRecord], List[Dict]]:
        """
        Clean raw qPCR data. Returns cleaned records and exclusion log.

        Rules:
        - Blank / missing Ct → excluded
        - Undetermined / NaN Ct → excluded
        - Duplicate well → flagged
        - Negative Ct → flagged for review
        """
        cleaned: List[CleanRecord] = []
        exclusions: List[Dict] = []
        seen_wells: Dict[str, str] = {}

        for rec in records:
            cr = CleanRecord(record=rec)

            # Check Ct validity
            if rec.ct is None or (isinstance(rec.ct, float) and math.isnan(rec.ct)):
                cr.excluded = True
                cr.exclusion_reason = "Missing or undetermined Ct"
            elif rec.ct < 0:
                cr.excluded = True
                cr.exclusion_reason = f"Negative Ct value: {rec.ct}"

            # Check duplicate well
            well_key = f"{rec.sample}_{rec.well}"
            if well_key in seen_wells:
                cr.excluded = True
                cr.exclusion_reason = f"Duplicate well {rec.well}"

            if cr.excluded:
                exclusions.append({
                    "sample": rec.sample, "gene": rec.gene,
                    "well": rec.well, "reason": cr.exclusion_reason,
                })
            else:
                seen_wells[well_key] = rec.well

            cleaned.append(cr)

        return cleaned, exclusions

    def replicate_qc(
        self, cleaned: List[CleanRecord]
    ) -> List[ReplicateStats]:
        """Compute replicate statistics per sample-gene group."""
        # Group by (sample, gene)
        groups: Dict[Tuple[str, str], List[float]] = {}
        group_meta: Dict[Tuple[str, str], str] = {}

        for cr in cleaned:
            if cr.excluded:
                continue
            rec = cr.record
            key = (rec.sample, rec.gene)
            groups.setdefault(key, []).append(rec.ct)
            group_meta[key] = rec.group

        results: List[ReplicateStats] = []
        for (sample, gene), cts in sorted(groups.items()):
            n = len(cts)
            mean_ct = statistics.mean(cts)
            sd = statistics.stdev(cts) if n > 1 else 0.0
            cv = (sd / mean_ct * 100) if mean_ct != 0 else 0.0

            rs = ReplicateStats(
                sample=sample, gene=gene, group=group_meta[(sample, gene)],
                n=n, mean_ct=round(mean_ct, 4), sd=round(sd, 4),
                cv=round(cv, 2), min_ct=min(cts), max_ct=max(cts),
                range_ct=round(max(cts) - min(cts), 4),
            )

            # QC status
            if cv > self.cv_warn_threshold:
                rs.status = "REVIEW"
            if n < 2:
                rs.status = "REVIEW"

            # Simple outlier detection (Z-score)
            if n >= 3 and sd > 0:
                for ct in cts:
                    z = abs((ct - mean_ct) / sd)
                    if z > self.outlier_z_threshold:
                        rs.outliers.append(ct)
                        rs.status = "OUTLIER"

            results.append(rs)

        return results

    def compute_dct(
        self, replicate_stats: List[ReplicateStats]
    ) -> List[DCtResult]:
        """
        Compute ΔCt = Ct(target) - Ct(reference).
        """
        # Build lookup: (sample, gene) -> mean_ct
        lookup: Dict[Tuple[str, str], float] = {
            (rs.sample, rs.gene): rs.mean_ct for rs in replicate_stats
        }
        sample_groups: Dict[str, str] = {
            rs.sample: rs.group for rs in replicate_stats
        }

        results: List[DCtResult] = []
        # Determine reference gene
        ref_gene = self.reference_gene
        if not ref_gene:
            # Use the most common gene as reference
            gene_counts: Dict[str, int] = {}
            for rs in replicate_stats:
                gene_counts[rs.gene] = gene_counts.get(rs.gene, 0) + 1
            ref_gene = max(gene_counts, key=gene_counts.get)

        for rs in replicate_stats:
            if rs.gene == ref_gene:
                continue
            ref_key = (rs.sample, ref_gene)
            if ref_key not in lookup:
                continue
            results.append(DCtResult(
                sample=rs.sample,
                group=sample_groups.get(rs.sample, ""),
                target_gene=rs.gene,
                reference_gene=ref_gene,
                mean_target_ct=rs.mean_ct,
                mean_ref_ct=lookup[ref_key],
                delta_ct=round(rs.mean_ct - lookup[ref_key], 4),
            ))

        return results

    def compute_ddct(
        self, dct_results: List[DCtResult]
    ) -> List[DDCtResult]:
        """
        Compute ΔΔCt = ΔCt(sample) - ΔCt(calibrator).
        Fold Change = 2^(-ΔΔCt).
        """
        if not dct_results:
            return []

        # Determine calibrator: use calibrator_group if set, else first sample
        calibrator_delta: Dict[str, float] = {}
        for d in dct_results:
            if self.calibrator_group and d.group == self.calibrator_group:
                calibrator_delta[d.target_gene] = d.delta_ct
            elif not self.calibrator_group:
                # Use first sample per gene as calibrator
                if d.target_gene not in calibrator_delta:
                    calibrator_delta[d.target_gene] = d.delta_ct

        results: List[DDCtResult] = []
        for d in dct_results:
            cal = calibrator_delta.get(d.target_gene, d.delta_ct)
            ddct = round(d.delta_ct - cal, 4)
            fc = round(2 ** (-ddct), 4)
            results.append(DDCtResult(
                sample=d.sample, group=d.group, gene=d.target_gene,
                delta_ct=d.delta_ct, delta_delta_ct=ddct,
                fold_change=fc,
                is_calibrator=(d.delta_ct == cal),
            ))

        return results

    def standard_curve(
        self, concentrations: List[float], cts: List[float]
    ) -> StandardCurve:
        """
        Fit a standard curve: Ct = slope × log10(conc) + intercept.
        Efficiency = 10^(-1/slope) - 1.
        """
        log_conc = [math.log10(c) for c in concentrations]
        n = len(log_conc)

        if n < 2:
            return StandardCurve()

        # Linear regression
        x_mean = statistics.mean(log_conc)
        y_mean = statistics.mean(cts)
        ss_xy = sum((x - x_mean) * (y - y_mean) for x, y in zip(log_conc, cts))
        ss_xx = sum((x - x_mean) ** 2 for x in log_conc)
        ss_yy = sum((y - y_mean) ** 2 for y in cts)

        slope = ss_xy / ss_xx if ss_xx != 0 else 0.0
        intercept = y_mean - slope * x_mean
        r_squared = (ss_xy ** 2) / (ss_xx * ss_yy) if ss_xx * ss_yy > 0 else 0.0
        efficiency = (10 ** (-1.0 / slope) - 1) * 100 if slope != 0 else 0.0

        return StandardCurve(
            concentrations=concentrations, cts=cts,
            slope=round(slope, 4), intercept=round(intercept, 4),
            r_squared=round(r_squared, 6),
            efficiency_pct=round(efficiency, 2),
        )
