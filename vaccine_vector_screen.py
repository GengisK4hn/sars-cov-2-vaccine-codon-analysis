#!/usr/bin/env python3
"""
vaccine_vector_screen.py
========================
Scan every vaccine FASTA on disk for known vector / backbone components:
  - SV40 promoter / enhancer / origin
  - SV40 polyA signal
  - CMV / EF1a / CAG promoters
  - pUC / ColE1 origin
  - KanR / NeoR / AmpR (antibiotic selection markers)
  - BGH polyA
  - T7 / SP6 promoter

For each vaccine file, report which vector components are present, their
coordinates, and % identity. Exact substring search + fuzzy 30-mer matching.

This is the "Step C" of the forensic plan: characterize what's in the
vaccine references beyond just the spike gene.
"""

from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
import json

BASE = Path('/home/dad/supracode-tool')


# Known vector components (canonical sequences from NCBI / Addgene / standard
# plasmid maps). These are short enough for substring matching.
VECTOR_COMPONENTS = {
    # SV40 promoter region (enhancer + early promoter)
    'SV40_enhancer_72bp': 'GTTCCGGTGGCCTGTAGAAATCCCTATCAGTGGGTGTCGGTAACGAAAGCAGTGCCAGCTGACCTTGAGG',
    'SV40_promoter_core': 'GGTAGGAGTAGGGCGGAGCGGGCGGGTGGCGGCGGTGGCGGCGGTGGCGGCGGTGGCGGCGG',
    # SV40 origin of replication
    'SV40_ori': 'GAGGCCCTGGAATCGGTGGCGGTGCCTGTAGAAATCCCTATCAGTGGGTGTCGGTAACGAAAG',
    # SV40 polyadenylation signal
    'SV40_polyA': 'AATTTGTGTTTTATTTGTTAACTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTGAATTTGT',
    # pCMV-BD backbone variants (the actual backbone family used in
    # BNT162b2/mRNA1273, not canonical SV40). Extracted from OR134577.1
    # at positions verified by comprehensive_vector_analysis.json to
    # share 100% identity with pCMV-BD AF151088.1.
    # See VECTOR_CORRECTED_VERDICT_20260712.txt for chain of evidence.
    'pCMV-BD_SV40_early_polyA': 'ATCGCCCTTCCCAACAGTTGCGCAGCCTGAATGGCGAATGGAGATCCAATTTTTAAGTGTATAATGTGTTAAACTACTGATTCTAATTGTTTGTGTATTTTAGATTCACAGTCCCAAGGCTCATTTCAGGCCCCTCAGTCCTCACAGTCTGTTCATGATCATAATCAGCCATACCACATTTGTAGAGGTTTTACTTGCTTTAAAAAACCTCCCACACCTCCCCCTGAACCTGAAACATAAAATGAATGCAATTGTTGTTGTTAACTTGTTTATTGCAGCTTATAATGGTTACAAATAAAGCAATAGCATCACAAATTTCACAAATAAAGCATTTTTTTCACTGCATTCTAGTTGTGGTTTGTCCAAACTCATCAATGTATCTTAA',
    'pCMV-BD_f1_origin': 'AACGCGTAAATTGTAAGCGTTAATATTTTGTTAAAATTCGCGTTAAATTTTTGTTAAATCAGCTCATTTTTTAACCAATAGGCCGAAATCGGCAAAATCCCTTATAAATCAAAAGAATAGACCGAGATAGGGTTGAGTGTTGTTCCAGTTTGGAACAAGAGTCCACTATTAAAGAACGTGGACTCCAACGTCAAAGGGCGAAAAACCGTCTATCAGGGCGATGGCCCACTACGTGAACCATCACCCTAATCAAGTTTTTTGGGGTCGAGGTGCCGTAAAGCACTAAATCGGAACCCTAAAGGGAGCCCCCGATTTAGAGCTTGACGGGGAAAGCCGGCGAACGTGGCGAGAAAGGAAGGGAAGAAAGCGAAAGGAGCGGGCGCTAGGGCGCTGGCAAGTGTAGCGGTCACGCTGCGCGTAACCACCACACCCGCCGCGCTTAATGCGCCGCTACAGGGCGCGT',
    'pCMV-BD_SV40_early_promoter_ori': 'CCTGAGGCGGAAAGAACCAGCTGTGGAATGTGTGTCAGTTAGGGTGTGGAAAGTCCCCAGGCTCCCCAGCAGGCAGAAGTATGCAAAGCATGCATCTCAATTAGTCAGCAACCAGGTGTGGAAAGTCCCCAGGCTCCCCAGCAGGCAGAAGTATGCAAAGCATGCATCTCAATTAGTCAGCAACCATAGTCCCGCCCCTAACTCCGCCCATCCCGCCCCTAACTCCGCCCAGTTCCGCCCATTCTCCGCCCCATGGCTGACTAATTTTTTTTATTTATGCAGAGGCCGAGGCCGCCTCGGCCTCTGAGCTATTCCAGAAGTAGTGAGGAGGCTTTTTTGGAGGCCTAGGCTTTTGCAAAGATCGATCAAGAGACAGGATGAGGATCGTTTCGC',
    # BGH polyA (bovine growth hormone)
    'BGH_polyA': 'GCCTCGACTGTGCCTTCTAGTTGCCAGCCATCTGTTGTTTGCCCCTCCCCCGTGCCTTCCTTGACCCTGG',
    # CMV immediate early promoter (canonical)
    'CMV_promoter': 'CCCCATTGACGTCAATGGGAGTTTGTTTTGGCACCAAAATCAACGGGACTTTCCAAAATGTCGTAACAACT',
    # EF1a promoter core
    'EF1a_promoter': 'AGCGGCCGCGGCCCGGGCGGTGGCGCGTGCCTGCTGGGGGCCCCTCGGGCCCCGGCGGGCCCCGCGGGGC',
    # Kanamycin resistance (CDS片段)
    'KanR_NeoR_marker': 'ATGATTGAACAAGATGGATTGCACGCAGGTTCTCCGGCCGCTTGGGTGGAGAGGCTATTCGGCTATGACTG',
    # Beta-lactamase / AmpR (CDS片段)
    'AmpR_blactamase': 'ATGAGTATTCAACATTTCCGTGTCGCCCTTATTCCCTTTTTTGCGGCATTTTGCCTTCCTGTTTTTGCTCA',
    # T7 promoter
    'T7_promoter': 'TAATACGACTCACTATAGGG',
    # SP6 promoter
    'SP6_promoter': 'ATTTAGGTGACACTATAG',
    # pUC origin region
    'pUC_ori': 'GCGGTGTTGGGCGGAGTTGAGTTGACAGCGTTACGCGTTTAAGCGGTGCCGGTTTACGGCCGGCGAATTC',
}


def rc(s):
    return str(Seq(s).reverse_complement())


def find_matches(query, seq, min_len=20):
    """Return list of (pos, strand, matched_len, identity) for query in seq.
    Uses sliding window of length len(query), step 5, allows mismatches."""
    hits = []
    L = len(query)
    seq = seq.upper()
    query = query.upper()
    # Forward exact
    p = 0
    while True:
        i = seq.find(query, p)
        if i == -1:
            break
        hits.append((i, '+', L, 100.0))
        p = i + 1
    # Reverse complement exact
    qrc = rc(query)
    p = 0
    while True:
        i = seq.find(qrc, p)
        if i == -1:
            break
        hits.append((i, '-', L, 100.0))
        p = i + 1
    # Fuzzy: 30-mers from query
    if L >= 30:
        for start in range(0, L - 30 + 1, 10):
            kmer = query[start:start + 30]
            for i in range(0, len(seq) - 30 + 1, 5):
                sub = seq[i:i + 30]
                matches = sum(1 for a, b in zip(kmer, sub) if a == b)
                if matches >= 27:  # >=90% identity
                    hits.append((i, '+', 30, 100.0 * matches / 30))
            kmer_rc = rc(kmer)
            for i in range(0, len(seq) - 30 + 1, 5):
                sub = seq[i:i + 30]
                matches = sum(1 for a, b in zip(kmer_rc, sub) if a == b)
                if matches >= 27:
                    hits.append((i, '-', 30, 100.0 * matches / 30))
    # Deduplicate (keep best per position)
    seen = {}
    for pos, strand, ml, ident in hits:
        key = (pos, strand)
        if key not in seen or ident > seen[key][2]:
            seen[key] = (pos, strand, ident)
    return list(seen.values())


def scan_vaccine(path):
    """Return list of hit dicts for a vaccine file."""
    out = []
    try:
        records = list(SeqIO.parse(str(path), 'fasta'))
    except Exception:
        return out
    for rec in records:
        seq = str(rec.seq).upper()
        if len(seq) < 50:
            continue
        for comp_name, comp_seq in VECTOR_COMPONENTS.items():
            hits = find_matches(comp_seq, seq)
            for pos, strand, ident in hits:
                ctx_l = seq[max(0, pos - 30):pos]
                ctx_r = seq[pos + len(comp_seq):pos + len(comp_seq) + 30] \
                    if strand == '+' else seq[pos + len(rc(comp_seq)):pos + len(rc(comp_seq)) + 30]
                out.append({
                    'file': str(path.relative_to(BASE)),
                    'record': rec.id[:80],
                    'record_len': len(seq),
                    'component': comp_name,
                    'pos': pos,
                    'strand': strand,
                    'identity': round(ident, 1),
                    'context_l': ctx_l,
                    'context_r': ctx_r,
                })
    return out


def main():
    # Find all vaccine FASTAs
    vaccine_patterns = [
        'pfizer', 'bnt', 'moderna', 'mrna', 'janssen', 'jNJ', 'astrazeneca',
        'oxford', '2p_vaccine', 'hexapro', '6p', 'bnt162', 'mrna1273',
        'pfizer_bnt', 'moderna_mrna', 'vaccine_vector', '2p_spike',
    ]

    all_fastas = list(BASE.rglob('*.fasta')) + list(BASE.rglob('*.fa')) + \
                 list(BASE.rglob('*.fna'))
    vaccine_files = []
    for p in all_fastas:
        name = p.name.lower()
        if any(v in name for v in vaccine_patterns) and p.is_file():
            vaccine_files.append(p)
    vaccine_files = sorted(set(vaccine_files))

    print(f'Found {len(vaccine_files)} vaccine-related FASTA files')
    for p in vaccine_files:
        print(f'  {p.relative_to(BASE)}')

    all_hits = []
    for p in vaccine_files:
        hits = scan_vaccine(p)
        if hits:
            all_hits.extend(hits)
            print(f'\n{p.relative_to(BASE)}: {len(hits)} component hits')
            comps = sorted(set(h['component'] for h in hits))
            for c in comps:
                print(f'  - {c}')

    df = pd.DataFrame(all_hits)
    if df.empty:
        print('\nNO vector component hits found.')
        return

    # Summary by file x component
    print('\n' + '=' * 70)
    print('SUMMARY: vector components per vaccine file')
    print('=' * 70)
    pivot = df.groupby(['file', 'component']).size().unstack(fill_value=0)
    print(pivot.to_string())

    df.to_csv(BASE / 'vaccine_vector_screen.csv', index=False)
    print(f'\nFile: vaccine_vector_screen.csv ({len(df)} rows)')

    # Critical components present?
    print('\nCritical findings:')
    critical = ['SV40_enhancer_72bp', 'SV40_promoter_core', 'SV40_ori',
                'SV40_polyA', 'CMV_promoter', 'EF1a_promoter',
                'KanR_NeoR_marker', 'AmpR_blactamase', 'BGH_polyA',
                # pCMV-BD backbone variants (added 2026-07-13 to fix
                # the false-negative documented in VECTOR_CORRECTED_VERDICT_20260712.txt)
                'pCMV-BD_SV40_early_polyA',
                'pCMV-BD_f1_origin',
                'pCMV-BD_SV40_early_promoter_ori']
    for comp in critical:
        hits = df[df['component'] == comp]
        if len(hits):
            files = sorted(set(hits['file']))
            print(f'  {comp}: {len(hits)} hits in {len(files)} files')
            for f in files[:3]:
                print(f'      {f}')


if __name__ == '__main__':
    main()
