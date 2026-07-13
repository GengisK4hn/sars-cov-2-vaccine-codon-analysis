#!/usr/bin/env python3
"""
calibrated_composite_scorer.py
==============================
Properly calibrated composite synthetic-signal score.

Method:
  1. Compute 8 features for every sequence:
       - GC_content (%)
       - GC3 (% at 3rd codon position)
       - CAI_human (Codon Adaptation Index vs human)
       - CpG_obs_exp (CpG dinucleotide observed/expected)
       - Nc (effective number of codons)
       - Type IIS site density per 10 kb
       - Total restriction site density per 10 kb
       - CG-skew standard deviation
  2. Build BASELINE distribution from natural controls
     (test_blind/natural_borrelia + shuffled + Wuhan/RaTG13/BANAL/RmYN/SARS1)
  3. For each test sequence, compute Z-score per feature.
  4. composite_z = weighted sum of feature Z-scores.
  5. Composite classification with empirical p-value.

Why this is more defensible than the old composite_scorer.py:
  - Old scorer used additive hard-threshold scoring (capped at 100, saturated).
  - This version uses statistical Z-scores calibrated to a natural baseline.
  - The natural baseline distribution is computed from real sequences.
  - p-values are empirical (rank-based), not assumed Gaussian.
"""

from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
import numpy as np
import re
import math
from collections import Counter

BASE = Path('/home/dad/supracode-tool')

# Reuse metrics from synthetic_signal_combined.py
import importlib.util
_spec = importlib.util.spec_from_file_location(
    'ssc', BASE / 'synthetic_signal_combined.py')
ssc = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ssc)


# Restriction enzyme sites (Type IIS + common Type II)
TYPE_IIS_SITES = {
    'BsaI': 'GGTCTC', 'BsaI_rc': 'GAGACC',
    'BsmBI': 'CGTCTC', 'BsmBI_rc': 'GAGACG',
    'BbsI': 'GAAGAC', 'BbsI_rc': 'GTCTTC',
    'SapI': 'GCTCTTC', 'SapI_rc': 'GAAGAGC',
}

COMMON_TYPE_II = {
    'EcoRI': 'GAATTC', 'BamHI': 'GGATCC', 'HindIII': 'AAGCTT',
    'XhoI': 'CTCGAG', 'NotI': 'GCGGCCGC', 'NdeI': 'CATATG',
    'NcoI': 'CCATGG', 'XbaI': 'TCTAGA',
}


def site_density(seq, motif_dict):
    seq_u = seq.upper()
    n = 0
    for _, site in motif_dict.items():
        n += len(re.findall(site, seq_u))
    return n / max(len(seq), 1) * 10000  # per 10 kb


def features_for_seq(seq):
    """Compute 8 features for a sequence."""
    seq = seq.upper()
    if len(seq) < 30:
        return None

    # Find best frame
    best = (0, 1e9)
    for f in [0, 1, 2]:
        stops = sum(1 for i in range(f, len(seq) - 2, 3)
                    if seq[i:i + 3] in ('TAA', 'TAG', 'TGA'))
        if stops < best[1]:
            best = (f, stops)
    frame = best[0]

    cai_h = ssc.cai(seq, frame)
    gc = ssc.gc_content(seq)
    gc3 = ssc.gc3(seq, frame)
    nc = ssc.effective_number_of_codons(seq, frame)
    cpg = ssc.cg_dinuc_freq(seq)
    skew = ssc.gc_skew(seq, 500)
    iis_d = site_density(seq, TYPE_IIS_SITES)
    type2_d = site_density(seq, COMMON_TYPE_II)

    return {
        'GC_content': gc,
        'GC3': gc3,
        'CAI_human': cai_h,
        'CpG_obs_exp': cpg,
        'Nc': nc,
        'Type_IIS_per_10kb': iis_d,
        'Type_II_per_10kb': type2_d,
        'GC_skew_std': skew,
    }


def features_for_fasta(fa_path):
    """Compute features for first record in a FASTA file."""
    try:
        seq = str(next(SeqIO.parse(str(fa_path), 'fasta')).seq).upper()
        f = features_for_seq(seq)
        if f:
            f['file'] = Path(fa_path).name
            f['path'] = str(fa_path)
        return f
    except Exception:
        return None


def build_baseline():
    """Build baseline distribution from natural viral genomes.

    Using VIRAL-ONLY baseline (not Borrelia) because we're testing viral
    sequences. A bacterial baseline would make every virus look anomalous
    by definition.
    """
    baseline_paths = []

    # Natural SARS-CoV-2 / bat CoV / SARS-CoV-1 genomes (viral baseline)
    cov_refs = [
        BASE / 'data/sequences/NC_045512_Wuhan-Hu-1_reference.fasta',
        BASE / 'data/sequences/RaTG13_MN996532.1.fasta',
        BASE / 'data/sequences/BANAL2052_MZ845808.1.fasta',
        BASE / 'data/sequences/RmYN02_MT040335.1.fasta',
        BASE / 'data/sequences/RacCS203_MW251308.fasta',
        BASE / 'data/sequences/RacCS224_MW251309.fasta',
        BASE / 'data/sequences/RacCS253_MW251310.fasta',
        BASE / 'data/sequences/SARS_CoV_1_Tor2_AY274119.3.fasta',
        BASE / 'data/sequences/SARS_CoV_1_Urbani_AY278741.1.fasta',
        BASE / 'data/sequences/Wuhan_Dec2019_MN908947.3.fasta',
        BASE / 'data/sequences/MT020880.1_Early_Wuhan.fasta',
        BASE / 'data/sequences/MT093631.1_Early_Wuhan.fasta',
        BASE / 'data/sequences/early_MT093631.1.fasta',
    ]
    baseline_paths.extend([p for p in cov_refs if p.exists()])

    # Also include a few SARS-CoV-2 clinical isolates as baseline
    addl = sorted((BASE / 'additional_clinical_2026').glob('*.fasta'))[:10]
    baseline_paths.extend(addl)

    rows = []
    for p in baseline_paths:
        f = features_for_fasta(p)
        if f:
            f['baseline_class'] = 'natural_viral'
            rows.append(f)

    # ALSO include natural Borrelia + shuffled as a SEPARATE reference class
    # (for sanity-checking Z-scores but not for setting thresholds)
    bdir = BASE / 'test_blind/natural_borrelia'
    if bdir.exists():
        for p in sorted(bdir.glob('*.fasta'))[:10]:
            f = features_for_fasta(p)
            if f:
                f['baseline_class'] = 'natural_borrelia'
                rows.append(f)

    return pd.DataFrame(rows), [p.name for p in baseline_paths]


def score_test_set(baseline_df, test_paths):
    """Compute Z-scores and composite for each test sequence."""
    # Baseline statistics
    feature_cols = ['GC_content', 'GC3', 'CAI_human', 'CpG_obs_exp', 'Nc',
                    'Type_IIS_per_10kb', 'Type_II_per_10kb', 'GC_skew_std']
    baseline_mean = baseline_df[feature_cols].mean()
    baseline_std = baseline_df[feature_cols].std().replace(0, 1)

    # For features where engineered = HIGHER than baseline, Z is positive.
    # For Nc, lower = simpler codon usage = potentially engineered.
    direction = {
        'GC_content': +1,
        'GC3': +1,
        'CAI_human': +1,
        'CpG_obs_exp': +1,
        'Nc': -1,  # fewer effective codons = biased usage
        'Type_IIS_per_10kb': +1,
        'Type_II_per_10kb': +1,
        'GC_skew_std': +1,
    }

    rows = []
    for label, fa_path in test_paths:
        f = features_for_fasta(fa_path)
        if not f:
            continue
        z = {}
        for col in feature_cols:
            raw_z = (f[col] - baseline_mean[col]) / baseline_std[col]
            z[col + '_z'] = direction[col] * raw_z
        # Composite: sum of POSITIVE Z-scores (only anomalies that go in
        # the engineered direction count).
        positive_z_sum = sum(max(0, v) for v in z.values())
        # Also: count of features with Z > 2 (strong anomalies)
        n_strong = sum(1 for v in z.values() if v > 2)
        # Composite scaled to 0-100: positive_z_sum / 8 features
        composite_z = positive_z_sum
        # Empirical p-value: how many baseline sequences score higher?
        # (computed below after applying same to baseline)

        row = {'label': label, 'file': f['file']}
        row.update({k: f[k] for k in feature_cols})
        row.update(z)
        row['composite_z_sum'] = round(composite_z, 3)
        row['n_features_anomalous'] = n_strong
        rows.append(row)

    return pd.DataFrame(rows), feature_cols, baseline_mean, baseline_std


def score_baseline(baseline_df, feature_cols, baseline_mean, baseline_std):
    """Score baseline against itself (for empirical p-value)."""
    direction = {
        'GC_content': +1, 'GC3': +1, 'CAI_human': +1, 'CpG_obs_exp': +1,
        'Nc': -1, 'Type_IIS_per_10kb': +1, 'Type_II_per_10kb': +1,
        'GC_skew_std': +1,
    }
    scores = []
    for _, baseline_row in baseline_df.iterrows():
        z_vals = []
        for col in feature_cols:
            raw_z = (baseline_row[col] - baseline_mean[col]) / baseline_std[col]
            z_vals.append(direction[col] * raw_z)
        scores.append(sum(max(0, v) for v in z_vals))
    return np.array(scores)


# ---------------------------------------------------------------------------
# Contamination-prior wiring (added 2026-07-13)
# ---------------------------------------------------------------------------
# Loads validation_data/known_contamination_levels.csv (multi-lab qPCR +
# fluorometry measurements of residual DNA in finished mRNA vaccine vials)
# and joins it with the scorer output so the "combined" view is computed,
# not hand-edited.

REGULATORY_LIMIT_NG_PER_DOSE = 10.0  # WHO TRS 978 / FDA / EMA threshold

CONTAMINATION_CSV = BASE / 'validation_data' / 'known_contamination_levels.csv'


def load_contamination_measurements():
    """Load multi-lab contamination CSV. Returns DataFrame or None if absent.

    Adds computed columns:
      - max_ng_per_dose (numeric; NaN if missing)
      - fold_over_limit_max = max_ng_per_dose / REGULATORY_LIMIT
    """
    if not CONTAMINATION_CSV.exists():
        return None
    df = pd.read_csv(CONTAMINATION_CSV)
    df['max_ng_per_dose'] = pd.to_numeric(df['max_ng_per_dose'],
                                          errors='coerce')
    df['fold_over_limit_max'] = df['max_ng_per_dose'] / REGULATORY_LIMIT_NG_PER_DOSE
    return df


def combine_scorer_with_contamination(test_df, contam_df):
    """Join scorer output with lab contamination measurements.

    Strategy: contamination CSV is per-lab-per-measurement (not per-vaccine
    brand). We emit TWO outputs:
      (1) Per-lab contamination table with fold-exceedance computed
      (2) Per-vaccine summary row joining the scorer's Z-score with the
          MAXIMUM observed contamination for that brand (best-effort
          brand match by lab notes)

    For labs that measured both brands or didn't specify, attribution is
    conservative: a 'both' marker is used and the lab contributes to BOTH
          vaccine summaries.
    """
    # Heuristic brand attribution from the 'notes' column
    def brands_for_lab(row):
        notes = str(row.get('notes', '')).lower()
        if 'pfizer' in notes and 'moderna' in notes:
            return ['Pfizer', 'Moderna']
        if 'pfizer' in notes or 'bnt' in notes or 'comirnaty' in notes:
            return ['Pfizer']
        if 'moderna' in notes or 'mrna' in notes or 'spikevax' in notes:
            return ['Moderna']
        # Default: unknown -> attribute to both (conservative for exceedance)
        return ['Pfizer', 'Moderna']

    contam_df = contam_df.copy()
    contam_df['brands'] = contam_df.apply(brands_for_lab, axis=1)

    # Per-vaccine max contamination
    per_brand = {}
    for brand in ['Pfizer', 'Moderna']:
        mask = contam_df['brands'].apply(lambda bl: brand in bl)
        sub = contam_df[mask].copy()
        if len(sub) == 0:
            per_brand[brand] = None
            continue
        # Exclude the regulatory_limit row itself from max
        sub = sub[sub['lab'] != 'regulatory_limit']
        max_val = sub['max_ng_per_dose'].max(skipna=True)
        max_lab = sub.loc[sub['max_ng_per_dose'].idxmax(), 'lab'] \
            if pd.notna(max_val) else None
        per_brand[brand] = {
            'max_ng_per_dose': max_val,
            'max_fold_over_limit': (max_val / REGULATORY_LIMIT_NG_PER_DOSE
                                    if pd.notna(max_val) else None),
            'max_lab': max_lab,
            'n_labs_measured': len(sub),
        }

    # Per-vaccine summary joining scorer Z-score
    label_to_brand = {
        'Pfizer_vector_OR134577': 'Pfizer',
        'Pfizer_2P_vaccine':      'Pfizer',
        'Moderna_vector_OR134578':'Moderna',
    }
    summary_rows = []
    for label, brand in label_to_brand.items():
        scorer_row = test_df[test_df['label'] == label]
        z = float(scorer_row['composite_z_sum'].iloc[0]) if len(scorer_row) else None
        cls = scorer_row['classification'].iloc[0] if len(scorer_row) else None
        contam = per_brand.get(brand)
        summary_rows.append({
            'brand': brand,
            'scorer_label': label,
            'composite_z_sum': z,
            'scorer_classification': cls,
            'contamination_max_ng_per_dose': contam['max_ng_per_dose'] if contam else None,
            'contamination_max_fold_over_limit': contam['max_fold_over_limit'] if contam else None,
            'contamination_max_lab': contam['max_lab'] if contam else None,
            'n_labs_measured': contam['n_labs_measured'] if contam else 0,
            'regulatory_limit_ng_per_dose': REGULATORY_LIMIT_NG_PER_DOSE,
        })
    summary_df = pd.DataFrame(summary_rows)
    return summary_df, contam_df



def main():
    print('Building baseline distribution from natural viral genomes...')
    baseline_df, baseline_names = build_baseline()
    # Use ONLY viral sequences for Z-score baseline (not Borrelia)
    viral_baseline = baseline_df[baseline_df['baseline_class'] == 'natural_viral']
    print(f'  Viral baseline: {len(viral_baseline)} natural viral genomes')
    print(f'  Borrelia ref (not used for Z): '
          f'{len(baseline_df) - len(viral_baseline)}')
    if len(viral_baseline) < 8:
        print('  WARNING: small viral baseline — Z-scores may be unstable')

    baseline_df_for_scoring = viral_baseline

    # Test set: vaccine vectors + positive controls + held-out naturals
    test_paths = [
        # Engineered (positives)
        ('Pfizer_2P_vaccine',  BASE / 'validation_data/2P_Vaccine.fasta'),
        ('HexaPro_6P',         BASE / 'validation_data/HexaPro.fasta'),
        ('Pfizer_vector_OR134577',
         BASE / 'data/all_sequences/pfizer_bnt162b2_vector_OR134577.fasta'),
        ('Moderna_vector_OR134578',
         BASE / 'data/all_sequences/moderna_mrna1273_vector_OR134578.fasta'),
        # Synthetic constructs from test_blind (true positives)
        ('synthetic_Vector_035_pUC',
         BASE / 'test_blind/synthetic_constructs/synthetic_Vector_035_pUC.fasta'),
        ('synthetic_GoldenGate_001_AarI',
         BASE / 'test_blind/synthetic_constructs/synthetic_GoldenGate_001_AarI.fasta'),
        ('synthetic_PromoterFusion_016_SV40',
         BASE / 'test_blind/synthetic_constructs/synthetic_PromoterFusion_016_SV40.fasta'),
        # Held-out: Wuhan-Hu-1 (it's in baseline but let's confirm Z=0)
        ('Wuhan_Hu1_full',
         BASE / 'data/sequences/NC_045512_Wuhan-Hu-1_reference.fasta'),
    ]

    test_df, feature_cols, bm, bs = score_test_set(
        baseline_df_for_scoring, test_paths)

    # Empirical p-value: rank test composite vs baseline composites
    baseline_composites = score_baseline(
        baseline_df_for_scoring, feature_cols, bm, bs)
    test_df['empirical_p'] = test_df['composite_z_sum'].apply(
        lambda s: (np.sum(baseline_composites >= s) + 1) /
                  (len(baseline_composites) + 1))
    test_df['empirical_p'] = test_df['empirical_p'].round(4)

    # Classification
    def classify(row):
        z = row['composite_z_sum']
        n = row['n_features_anomalous']
        p = row['empirical_p']
        if z >= 15 and n >= 4 and p <= 0.05:
            return 'HIGH_CONFIDENCE_SYNTHETIC'
        if z >= 8 and n >= 3:
            return 'LIKELY_SYNTHETIC'
        if z >= 4 and n >= 2:
            return 'POSSIBLE_ANOMALY'
        return 'NATURAL'
    test_df['classification'] = test_df.apply(classify, axis=1)

    # Report
    print('\n' + '=' * 90)
    print('CALIBRATED COMPOSITE SCORE')
    print('=' * 90)
    show = ['label', 'composite_z_sum', 'n_features_anomalous',
            'empirical_p', 'classification']
    print(test_df[show].sort_values(
        'composite_z_sum', ascending=False).to_string(index=False))

    print('\nPer-feature Z-scores (only engineered direction counts):')
    z_cols = [c for c in test_df.columns if c.endswith('_z')]
    print(test_df[['label'] + z_cols].sort_values(
        'label').to_string(index=False))

    print(f'\nBaseline distribution of composite (natural):')
    print(f'  n = {len(baseline_composites)}')
    print(f'  mean = {baseline_composites.mean():.2f}')
    print(f'  std  = {baseline_composites.std():.2f}')
    print(f'  max  = {baseline_composites.max():.2f}')
    print(f'  95th pct = {np.percentile(baseline_composites, 95):.2f}')

    baseline_df.to_csv(BASE / 'calibrated_baseline_features.csv', index=False)
    test_df.to_csv(BASE / 'calibrated_composite_scores.csv', index=False)
    print(f'\nFiles:')
    print(f'  calibrated_baseline_features.csv ({len(baseline_df)} rows)')
    print(f'  calibrated_composite_scores.csv  ({len(test_df)} rows)')

    # ---- Combined view: scorer Z + contamination prior ----
    contam_df = load_contamination_measurements()
    if contam_df is not None:
        summary_df, contam_annotated = combine_scorer_with_contamination(
            test_df, contam_df)
        contam_annotated.to_csv(
            BASE / 'combined_contamination_by_lab.csv', index=False)
        summary_df.to_csv(
            BASE / 'combined_per_brand_summary.csv', index=False)
        print('\n' + '=' * 90)
        print('COMBINED VIEW  (sequence Z-score + vial contamination prior)')
        print('=' * 90)
        print(summary_df.to_string(index=False))
        print(f'\nFiles:')
        print(f'  combined_contamination_by_lab.csv ({len(contam_annotated)} rows)')
        print(f'  combined_per_brand_summary.csv    ({len(summary_df)} rows)')
    else:
        print(f'\n[combined view skipped: {CONTAMINATION_CSV.name} not found]')


if __name__ == '__main__':
    main()
