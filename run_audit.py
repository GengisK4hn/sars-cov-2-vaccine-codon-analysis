#!/usr/bin/env python3
"""
run_audit.py — One-stop forensic audit
=======================================
Runs the full SupraCode forensic pipeline on all data and prints
a live audit table as results come in.

Usage:
    python3 run_audit.py              # Full audit (scorer + vector + detectors)
    python3 run_audit.py --quick      # Quick audit (scorer + vector screen only)
    python3 run_audit.py --scorer     # Scorer only (Z-scores + RSCU + contamination)

Each phase prints its results immediately so you can audit as it runs.
"""

import subprocess
import sys
import time
from pathlib import Path

BASE = Path('/home/dad/supracode-tool')


def run_phase(name, cmd, timeout=300):
    """Run a phase, print live output, return success."""
    print(f"\n{'='*80}")
    print(f"  PHASE: {name}")
    print(f"  CMD:   {' '.join(cmd)}")
    print(f"{'='*80}\n")
    t0 = time.time()
    try:
        result = subprocess.run(
            cmd, cwd=str(BASE), timeout=timeout,
            capture_output=False, text=True)
        dt = time.time() - t0
        ok = result.returncode == 0
        status = 'OK' if ok else 'FAIL'
        print(f"\n  [{status}] {name} ({dt:.1f}s)")
        return ok
    except subprocess.TimeoutExpired:
        print(f"\n  [TIMEOUT] {name} (>{timeout}s)")
        return False
    except Exception as e:
        print(f"\n  [ERROR] {name}: {e}")
        return False


def print_summary(results):
    """Print final audit summary."""
    print(f"\n{'='*80}")
    print('  AUDIT SUMMARY')
    print(f"{'='*80}")
    for phase, (ok, dt) in results.items():
        status = 'PASS' if ok else 'FAIL'
        print(f"  [{status}] {phase:40s}  {dt:.1f}s")

    # List output files
    print(f"\n  Output files:")
    outputs = [
        'calibrated_composite_scores.csv',
        'calibrated_baseline_features.csv',
        'rscu_codon_usage_comparison.csv',
        'combined_contamination_by_lab.csv',
        'combined_per_brand_summary.csv',
        'vaccine_vector_screen.csv',
        'comprehensive_vector_analysis.json',
        'COMBINED_SYNTHETIC_CONTAMINATION_REPORT.txt',
    ]
    for f in outputs:
        p = BASE / f
        if p.exists():
            size = p.stat().st_size
            print(f"    OK   {f:50s}  {size:>8,} bytes")
        else:
            print(f"    ---  {f:50s}  (not generated)")

    print(f"\n  Read the master report:")
    print(f"    cat COMBINED_SYNTHETIC_CONTAMINATION_REPORT.txt")
    print(f"{'='*80}\n")


def main():
    import argparse
    parser = argparse.ArgumentParser(
        description='One-stop SupraCode forensic audit')
    parser.add_argument('--quick', action='store_true',
                        help='Quick mode: scorer + vector screen only')
    parser.add_argument('--scorer', action='store_true',
                        help='Scorer only (Z-scores + RSCU + contamination)')
    args = parser.parse_args()

    py = sys.executable
    results = {}

    # ---- Phase 1: Calibrated composite scorer (always runs) ----
    t0 = time.time()
    ok = run_phase(
        '1. Calibrated Z-score + RSCU + contamination join',
        [py, 'calibrated_composite_scorer.py'],
        timeout=300)
    results['Calibrated scorer'] = (ok, time.time() - t0)

    if args.scorer:
        print_summary(results)
        return

    # ---- Phase 2: Vaccine vector screen ----
    t0 = time.time()
    ok = run_phase(
        '2. Vector backbone screen (pCMV-BD, SV40, KanR, f1, ColE1)',
        [py, 'vaccine_vector_screen.py'],
        timeout=120)
    results['Vector screen'] = (ok, time.time() - t0)

    if args.quick:
        print_summary(results)
        return

    # ---- Phase 3: Comprehensive vector analysis ----
    t0 = time.time()
    ok = run_phase(
        '3. Deep vector feature annotation',
        [py, 'comprehensive_vector_analysis.py'],
        timeout=120)
    results['Vector analysis'] = (ok, time.time() - t0)

    # ---- Phase 4: Key detectors on vaccine sequences ----
    test_seqs = [
        ('Pfizer 2P spike', 'validation_data/2P_Vaccine.fasta'),
        ('Pfizer vector',   'data/all_sequences/pfizer_bnt162b2_vector_OR134577.fasta'),
        ('Moderna vector',  'data/all_sequences/moderna_mrna1273_vector_OR134578.fasta'),
        ('Wuhan Hu-1',      'data/sequences/NC_045512_Wuhan-Hu-1_reference.fasta'),
    ]

    detectors = [
        ('FCS detector',           'fcs_detector.py'),
        ('Cloning scar detector',  'cloning_scar_detector.py'),
        ('Codon optimization',     'codon_optimization_verifier.py'),
        ('Stabilizing mutations',  'stabilizing_mutations_detector.py'),
    ]

    for det_name, det_script in detectors:
        det_path = BASE / det_script
        if not det_path.exists():
            print(f"\n  [SKIP] {det_name} ({det_script} not found)")
            results[det_name] = (False, 0)
            continue

        t0 = time.time()
        # Run on first available test sequence
        ran = False
        for seq_label, seq_path in test_seqs:
            full_path = BASE / seq_path
            if full_path.exists():
                ok = run_phase(
                    f'4. {det_name} on {seq_label}',
                    [py, det_script, str(full_path)],
                    timeout=60)
                results[det_name] = (ok, time.time() - t0)
                ran = True
                break
        if not ran:
            results[det_name] = (False, 0)

    print_summary(results)


if __name__ == '__main__':
    main()
