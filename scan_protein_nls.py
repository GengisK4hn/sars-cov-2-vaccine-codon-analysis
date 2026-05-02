#!/usr/bin/env python3
"""
Scan TRANSLATED protein sequences for Nuclear Localization Signals (NLS)
"""

import sys
import re
from Bio import SeqIO
from Bio.Seq import Seq

# Known NLS motifs in protein sequences
PROTEIN_NLS_PATTERNS = [
    # Classical monopartite NLS
    ("PKKKRKV", "SV40_large_T_antigen_NLS"),
    ("KKKRK", "SV40_type_short"),
    
    # Bipartite NLS (2 basic clusters separated by 10-12 residues)
    (r"KR.{10,12}KKK", "bipartite_NLS"),
    
    # Other NLS variants
    ("KKKR", "basic_NLS_1"),
    ("RKKR", "reverse_NLS"),
    ("KKRR", "basic_NLS_2"),
    ("RRRR", "polyarginine_NLS"),
    ("KKKK", "polylysine_NLS"),
    ("KRRR", "three_arginine"),
    ("RRKR", "arginine_lysine_cluster"),
]

def translate_dna(dna_seq):
    """Translate DNA to protein (all 6 frames)"""
    proteins = []
    
    for frame in range(3):
        seq = dna_seq[frame:]
        # Translate in forward direction
        try:
            protein = str(Seq(seq).translate(to_stop=False))
            proteins.append(('frame+', frame+1, protein))
        except:
            pass
    
    # Also translate reverse complement
    rev_comp = str(Seq(dna_seq).reverse_complement())
    for frame in range(3):
        seq = rev_comp[frame:]
        try:
            protein = str(Seq(seq).translate(to_stop=False))
            proteins.append(('frame-', frame+1, protein))
        except:
            pass
    
    return proteins

def scan_protein_for_nls(protein_seq, frame_dir, frame_num):
    """Scan protein sequence for NLS motifs"""
    
    results = {
        'frame': f"{frame_dir}{frame_num}",
        'protein_length': len(protein_seq),
        'nls_found': [],
        'total_count': 0
    }
    
    for pattern, name in PROTEIN_NLS_PATTERNS:
        matches = list(re.finditer(pattern, protein_seq))
        if matches:
            results['nls_found'].append({
                'motif': pattern,
                'name': name,
                'count': len(matches),
                'positions': [(m.start(), m.group()) for m in matches[:5]]
            })
            results['total_count'] += len(matches)
    
    return results

def main():
    if len(sys.argv) < 2:
        print("Usage: scan_protein_nls.py <dna_fasta_file> [output_file]")
        sys.exit(1)
    
    fasta_file = sys.argv[1]
    output_file = sys.argv[2] if len(sys.argv) > 2 else None
    
    print("=" * 80)
    print("PROTEIN-LEVEL NLS SCAN (6-frame translation)")
    print("=" * 80)
    print(f"Input: {fasta_file}")
    
    records = list(SeqIO.parse(fasta_file, "fasta"))
    
    if not records:
        print("Error: No sequences found")
        sys.exit(1)
    
    all_results = []
    
    for record in records:
        print(f"\nAnalyzing {record.id}...")
        dna_seq = str(record.seq)
        
        # Translate all 6 frames
        proteins = translate_dna(dna_seq)
        
        print(f"  Translated {len(proteins)} reading frames")
        
        record_results = {
            'sequence_id': record.id,
            'dna_length': len(dna_seq),
            'frames': []
        }
        
        for frame_dir, frame_num, protein_seq in proteins:
            result = scan_protein_for_nls(protein_seq, frame_dir, frame_num)
            if result['total_count'] > 0:
                record_results['frames'].append(result)
                print(f"  Frame {result['frame']}: {result['total_count']} NLS motifs found")
        
        all_results.append(record_results)
    
    # Generate report
    report = []
    report.append("=" * 80)
    report.append("PROTEIN-LEVEL NLS SCAN RESULTS")
    report.append("=" * 80)
    
    for record_results in all_results:
        report.append(f"\nSequence: {record_results['sequence_id']}")
        report.append(f"DNA Length: {record_results['dna_length']} bp")
        
        if not record_results['frames']:
            report.append("\n✓ No NLS motifs detected in any reading frame")
            report.append("  This suggests the translated spike protein does NOT contain")
            report.append("  classical nuclear localization signals.")
        else:
            report.append(f"\n⚠ NLS MOTIFS DETECTED in {len(record_results['frames'])} reading frame(s):")
            
            for frame in record_results['frames']:
                report.append(f"\n  Frame {frame['frame']} (Protein length: {frame['protein_length']} aa):")
                report.append(f"  Total NLS motifs: {frame['total_count']}")
                
                for nls in frame['nls_found']:
                    report.append(f"    • {nls['name']}: {nls['count']} occurrences")
                    report.append(f"      Pattern: {nls['motif']}")
                    report.append(f"      Sample positions: {nls['positions']}")
    
    report.append("\n" + "=" * 80)
    report.append("INTERPRETATION")
    report.append("=" * 80)
    report.append("If NLS motifs are present in translated proteins:")
    report.append("  • Vaccine-derived spike protein could enter nucleus")
    report.append("  • Increased risk of genomic integration")
    report.append("  • Potential interference with nuclear processes")
    report.append("\nIf NO NLS motifs are detected:")
    report.append("  • Spike protein likely remains cytoplasmic/membrane-bound")
    report.append("  • Lower nuclear entry risk")
    report.append("  • Consistent with natural SARS-CoV-2 spike localization")
    report.append("=" * 80)
    
    report_text = "\n".join(report)
    print("\n" + report_text)
    
    if output_file:
        with open(output_file, 'w') as f:
            f.write(report_text)
        print(f"\nReport saved to: {output_file}")

if __name__ == "__main__":
    main()
