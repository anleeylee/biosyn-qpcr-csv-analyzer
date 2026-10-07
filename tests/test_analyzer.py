"""Tests for BioSyn qPCR CSV Analyzer."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from analyzer import qPCRAnalyzer, qPCRRecord


def test_clean_missing_ct():
    a = qPCRAnalyzer()
    records = [
        qPCRRecord("S1", "G1", "GAPDH", "A1", 22.0),
        qPCRRecord("S1", "G1", "GAPDH", "A2", None),
    ]
    cleaned, exclusions = a.clean_data(records)
    assert len(exclusions) == 1


def test_replicate_qc():
    a = qPCRAnalyzer(cv_warn_threshold=1.0)
    records = [
        qPCRRecord("S1", "G1", "GAPDH", "A1", 22.0),
        qPCRRecord("S1", "G1", "GAPDH", "A2", 22.2),
        qPCRRecord("S1", "G1", "GAPDH", "A3", 22.1),
    ]
    cleaned, _ = a.clean_data(records)
    stats = a.replicate_qc(cleaned)
    assert len(stats) == 1
    assert stats[0].n == 3
    assert stats[0].mean_ct == 22.1


def test_dct():
    a = qPCRAnalyzer(reference_gene="GAPDH")
    records = [
        qPCRRecord("S1", "T", "GAPDH", "A1", 20.0),
        qPCRRecord("S1", "T", "Target_X", "B1", 26.0),
    ]
    cleaned, _ = a.clean_data(records)
    stats = a.replicate_qc(cleaned)
    dct = a.compute_dct(stats)
    assert len(dct) == 1
    assert dct[0].delta_ct == 6.0


def test_standard_curve():
    a = qPCRAnalyzer()
    sc = a.standard_curve(
        concentrations=[1e6, 1e5, 1e4, 1e3, 1e2],
        cts=[15.0, 18.5, 22.0, 25.5, 29.0],
    )
    assert sc.r_squared > 0.99
    assert 80 < sc.efficiency_pct < 120


if __name__ == "__main__":
    test_clean_missing_ct(); test_replicate_qc(); test_dct(); test_standard_curve()
    print("All tests passed!")
