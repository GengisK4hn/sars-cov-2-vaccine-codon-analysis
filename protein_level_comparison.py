#!/usr/bin/env python3
"""
Protein-level comparison: Pfizer 2P spike vs Wuhan spike.
Runs prion-domain analysis (PLAAC + local Q/N) and NLS scanner
on BOTH sequences to establish whether codon optimization or
the 2P mutations change protein-level features.

Key principle: the 2P vaccine spike and Wuhan spike encode
proteins that are 99.8% identical (only K986P, V987P differ).
So any DIFFERENCE in prion/NLS signal between them must come
from those two mutations -- NOT from codon usage. This is the
null-model control the user's findings need.
"""
import sys
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq

# Make repo modules importable
REPO = Path("/home/dad/supracode-tool")
sys.path.insert(0, str(REPO))

from spike_prion_analyzer import (
    _extract_spike_protein, locate_prd_and_mutations, run_full_plaac
)
from scan_protein_nls import translate_dna, scan_protein_for_nls

PFIZER_FASTA = REPO / "validation_data/2P_Vaccine.fasta"
WUHAN_FASTA  = REPO / "data/all_sequences/sars_cov_2_wuhan_spike.fasta"

def load_protein(fasta_path, label):
    """Load a FASTA and extract the spike protein."""
    rec = next(SeqIO.parse(fasta_path, "fasta"))
    dna = str(rec.seq).upper()
    prot = _extract_spike_protein(dna)
    # Truncate at first stop
    stop = prot.find("*")
    if stop > 0:
        prot = prot[:stop]
    print(f"{label}: {len(dna)} nt DNA -> {len(prot)} aa protein")
    return dna, prot

def prion_comparison(proteins, dnas):
    """Run PLAAC + local Q/N prion analysis on both proteins."""
    print("\n" + "=" * 72)
    print("PRION-DOMAIN ANALYSIS (PLAAC HMM + local Q/N scan)")
    print("=" * 72)
    # Build variant dict for PLAAC
    variant_dict = {label: prot for label, prot in proteins.items()}
    try:
        plaac_df = run_full_plaac(variant_dict)
        print("\n--- PLAAC HMM results ---")
        cols = [c for c in ["variant","max_LLR","COREscore","PrD_start","PrD_end","has_PrD"] if c in plaac_df.columns]
        print(plaac_df[cols].to_string(index=False))
    except Exception as e:
        print(f"PLAAC run failed: {e}")
        plaac_df = None
    # Local Q/N analysis
    print("\n--- Local Q/N + mutation analysis (Tetz PrD 473-510) ---")
    local_rows = []
    for label, prot in proteins.items():
        # locate_prd_and_mutations expects DNA (it extracts the spike
        # protein internally via _extract_spike_protein). Pass the DNA.
        r = locate_prd_and_mutations(Seq(dnas[label]), label)
        local_rows.append(r)
        print(f"\n{label}:")
        print(f"  input_background  : {r['input_background']}")
        print(f"  PrD sequence     : {r['prd_sequence']}")
        print(f"  Q+N content      : {r['qn_content_prd_percent']}%")
        print(f"  mutations in PrD : {r['mutations_in_prd']}")
    return plaac_df, local_rows

def nls_comparison(proteins, dnas):
    """Run NLS scan on the CORRECT frame (frame +1) for both.
    Also scan all 6 frames to show how many spurious motifs appear
    in non-translated frames -- this is the null-model point."""
    print("\n" + "=" * 72)
    print("NLS SCAN (nuclear localization signals)")
    print("=" * 72)
    for label in proteins:
        prot = proteins[label]
        dna = dnas[label]
        print(f"\n--- {label} ---")
        # Correct-frame scan (frame +1, i.e. the actual translated protein)
        result_correct = scan_protein_for_nls(prot, "+", 1)
        print(f"  Correct frame (+1, {len(prot)} aa): {result_correct['total_count']} NLS motifs")
        for nls in result_correct["nls_found"]:
            print(f"    {nls['name']}: {nls['count']}x  pattern={nls['motif']}  pos={nls['positions'][:3]}")
        # All-6-frame scan (shows the null-model problem)
        all_prots = translate_dna(dna)
        total_all = 0
        for fdir, fnum, pseq in all_prots:
            r = scan_protein_for_nls(pseq, fdir, fnum)
            total_all += r["total_count"]
        print(f"  All 6 frames combined: {total_all} NLS motifs (incl. non-translated frames)")
        print(f"  -> Spurious (non-translated) motifs: {total_all - result_correct['total_count']}")

def aa_diff(proteins):
    """Show the exact aa differences between the two proteins."""
    print("\n" + "=" * 72)
    print("AA-LEVEL DIFFERENCES (vaccine vs Wuhan)")
    print("=" * 72)
    p = proteins.get("Pfizer_2P")
    w = proteins.get("Wuhan")
    if not p or not w:
        print("missing protein")
        return
    n = min(len(p), len(w))
    diffs = [(i+1, w[i], p[i]) for i in range(n) if w[i] != p[i]]
    print(f"Lengths: Pfizer={len(p)} aa, Wuhan={len(w)} aa")
    print(f"Differences: {len(diffs)} positions")
    for pos, wa, pa in diffs[:20]:
        flag = "  <-- 2P mutation" if pos in (986, 987) else ""
        print(f"  pos {pos}: Wuhan={wa} -> Pfizer={pa}{flag}")
    if len(diffs) > 20:
        print(f"  ... and {len(diffs)-20} more")

def main():
    print("=" * 72)
    print("PROTEIN-LEVEL FORENSIC COMPARISON")
    print("Pfizer BNT162b2 2P spike  vs  Wuhan-Hu-1 spike (natural control)")
    print("=" * 72)
    p_dna, p_prot = load_protein(PFIZER_FASTA, "Pfizer_2P")
    w_dna, w_prot = load_protein(WUHAN_FASTA,  "Wuhan")
    proteins = {"Pfizer_2P": p_prot, "Wuhan": w_prot}
    dnas     = {"Pfizer_2P": p_dna, "Wuhan": w_dna}
    aa_diff(proteins)
    prion_comparison(proteins, dnas)
    nls_comparison(proteins, dnas)
    print("\n" + "=" * 72)
    print("DONE")
    print("=" * 72)

if __name__ == "__main__":
    main()
