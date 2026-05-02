#!/usr/bin/env python3
"""
CODON OPTIMIZATION VERIFIER
============================

Verifies codon optimization patterns by comparing:
1. Target sequence vs. natural SARS-CoV-2 codon usage
2. Target sequence vs. human codon usage
3. Identifies which codons are optimized
4. Calculates RSCU (Relative Synonymous Codon Usage)

Answers: Is this sequence optimized for human expression?
"""

import sys
import os
import json
import logging
from pathlib import Path
from typing import Dict, List, Tuple
from collections import defaultdict
import numpy as np
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class CodonOptimizationVerifier:
    """Verify codon optimization patterns."""

    # Standard genetic code
    CODON_TABLE = {
        'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
        'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
        'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
        'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
        'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
        'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'CGA': 'R', 'CGG': 'R', 'CGA': 'R', 'CGG': 'R',
        'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
        'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
        'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
    }

    # Human preferred codons (high CAI)
    HUMAN_PREFERRED = {
        'F': 'TTT', 'L': 'CTG', 'I': 'ATC', 'M': 'ATG',
        'V': 'GTG', 'S': 'AGC', 'P': 'CCC', 'T': 'ACC',
        'A': 'GCT', 'Y': 'TAC', 'H': 'CAC', 'Q': 'CAG',
        'N': 'AAC', 'K': 'AAG', 'D': 'GAC', 'E': 'GAG',
        'C': 'TGC', 'W': 'TGG', 'R': 'CGC', 'G': 'GGC',
        'S': 'AGC', 'R': 'CGC', 'L': 'CTG'
    }

    def __init__(self):
        self.results = {}

    def build_codon_usage(self, sequences: Dict[str, str]) -> Dict[str, Dict[str, int]]:
        """Build codon usage table from reference sequences."""
        codon_counts = defaultdict(lambda: defaultdict(int))

        for seq_id, seq in sequences.items():
            seq_upper = seq.upper()
            # Convert T to U if RNA
            seq_upper = seq_upper.replace('T', 'U')

            for i in range(0, len(seq_upper) - 2, 3):
                codon = seq_upper[i:i+3]
                if len(codon) == 3 and all(c in 'AUCG' for c in codon):
                    if codon in self.CODON_TABLE:
                        aa = self.CODON_TABLE[codon]
                        if aa != '*':  # Skip stop codons
                            codon_counts[aa][codon] += 1

        return dict(codon_counts)

    def calculate_rscu(self, target_seq: str, reference_seqs: Dict[str, str]) -> Dict[str, float]:
        """
        Calculate Relative Synonymous Codon Usage (RSCU).

        CORRECT FORMULA:
        RSCU_ij = (observed count of codon j) / (expected count if all synonymous codons used equally)
        Expected = (total count of amino acid i) / (number of synonymous codons)

        RSCU > 1: codon used more often than expected (preferred)
        RSCU < 1: codon used less often than expected (avoided)
        RSCU = 1: codon used as expected (neutral)
        """
        # Count codons in target
        target_codon_counts = defaultdict(int)
        target_aa_counts = defaultdict(int)
        target_upper = target_seq.upper().replace('T', 'U')

        for i in range(0, len(target_upper) - 2, 3):
            codon = target_upper[i:i+3]
            if len(codon) == 3 and all(c in 'AUCG' for c in codon):
                if codon in self.CODON_TABLE:
                    aa = self.CODON_TABLE[codon]
                    if aa != '*':
                        target_codon_counts[codon] += 1
                        target_aa_counts[aa] += 1

        # Build reference codon usage (for comparison)
        ref_codon_usage = self.build_codon_usage(reference_seqs)

        # Calculate RSCU for each codon in target
        rscu_values = {}
        comparison_to_natural = {}

        for aa, codons in self._get_synonymous_codons().items():
            if aa not in target_aa_counts:
                continue

            total_aa_count = target_aa_counts[aa]
            num_synonymous = len(codons)

            if num_synonymous == 0 or total_aa_count == 0:
                continue

            # Expected count if all synonymous codons used equally
            expected_count = total_aa_count / num_synonymous

            for codon in codons:
                if codon in target_codon_counts:
                    observed_count = target_codon_counts[codon]
                    # CORRECT RSCU CALCULATION
                    rscu_values[codon] = observed_count / expected_count if expected_count > 0 else 0

                    # Also compare to natural reference frequency
                    if aa in ref_codon_usage:
                        ref_total = sum(ref_codon_usage[aa].values())
                        if ref_total > 0:
                            ref_freq = ref_codon_usage[aa].get(codon, 0) / ref_total
                            target_freq = observed_count / total_aa_count
                            # Ratio of target frequency to natural frequency
                            if ref_freq > 0:
                                comparison_to_natural[codon] = target_freq / ref_freq

        return rscu_values, comparison_to_natural

    def _get_synonymous_codons(self) -> Dict[str, List[str]]:
        """Get synonymous codons for each amino acid."""
        aa_to_codons = defaultdict(list)

        for codon, aa in self.CODON_TABLE.items():
            if aa != '*':
                aa_to_codons[aa].append(codon)

        return dict(aa_to_codons)

    def detect_optimization(self, target_seq: str, natural_seqs: Dict[str, str],
                          human_pref: bool = True) -> Dict:
        """
        Detect codon optimization patterns.

        Returns:
            Analysis of optimized codons and overall assessment.
        """
        logger.info("Detecting codon optimization...")

        # Calculate RSCU (within target) and comparison to natural
        rscu_values, comparison_to_natural = self.calculate_rscu(target_seq, natural_seqs)

        # Find optimized codons based on RSCU (>1.5 = strongly preferred)
        optimized_codons = []
        for codon, rscu in rscu_values.items():
            # Check if it matches human preferred codons
            aa = self.CODON_TABLE.get(codon.replace('U', 'T'), '')
            is_human_preferred = self.HUMAN_PREFERRED.get(aa, '') == codon.replace('U', 'T')

            # Get comparison to natural if available
            vs_natural = comparison_to_natural.get(codon, None)

            if rscu > 1.5:
                optimized_codons.append({
                    'codon': codon,
                    'rscu': rscu,
                    'amino_acid': aa,
                    'level': 'STRONGLY_OPTIMIZED',
                    'human_preferred': is_human_preferred,
                    'vs_natural_ratio': vs_natural
                })
            elif rscu > 1.2:
                optimized_codons.append({
                    'codon': codon,
                    'rscu': rscu,
                    'amino_acid': aa,
                    'level': 'MODERATELY_OPTIMIZED',
                    'human_preferred': is_human_preferred,
                    'vs_natural_ratio': vs_natural
                })

        # Calculate optimization score based on RSCU
        if rscu_values:
            avg_rscu = np.mean(list(rscu_values.values()))
            max_rscu = max(rscu_values.values())
            optimized_count = len([c for c in rscu_values.values() if c > 1.5])

            # Assessment based on RSCU values
            if avg_rscu > 1.3:
                assessment = 'HIGHLY_OPTIMIZED'
            elif avg_rscu > 1.15:
                assessment = 'OPTIMIZED'
            elif avg_rscu > 1.05:
                assessment = 'PARTIALLY_OPTIMIZED'
            else:
                assessment = 'NOT_OPTIMIZED'
        else:
            avg_rscu = 0
            max_rscu = 0
            optimized_count = 0
            assessment = 'UNKNOWN'

        # Calculate deviation from natural sequences
        if comparison_to_natural:
            avg_vs_natural = np.mean(list(comparison_to_natural.values()))
            max_vs_natural = max(comparison_to_natural.values())
            # Count codons with >2x deviation from natural
            high_deviation_count = len([c for c in comparison_to_natural.values() if c > 2.0 or c < 0.5])
        else:
            avg_vs_natural = 0
            max_vs_natural = 0
            high_deviation_count = 0

        return {
            'assessment': assessment,
            'avg_rscu': avg_rscu,
            'max_rscu': max_rscu,
            'optimized_codon_count': optimized_count,
            'optimized_codons': optimized_codons[:20],  # Top 20
            'vs_natural': {
                'avg_ratio': avg_vs_natural,
                'max_ratio': max_vs_natural,
                'high_deviation_count': high_deviation_count
            }
        }

    def compare_to_natural(self, target_seq: str, natural_seqs: Dict[str, str]) -> Dict:
        """Compare target to natural SARS-CoV-2 codon usage."""
        # Build natural codon usage
        natural_codon_usage = self.build_codon_usage(natural_seqs)

        # Analyze target
        target_upper = target_seq.upper().replace('T', 'U')
        target_codon_counts = defaultdict(int)

        for i in range(0, len(target_upper) - 2, 3):
            codon = target_upper[i:i+3]
            if len(codon) == 3 and all(c in 'AUCG' for c in codon):
                if codon in self.CODON_TABLE:
                    aa = self.CODON_TABLE[codon]
                    if aa != '*':
                        target_codon_counts[codon] += 1

        # Compare
        comparisons = []
        for aa in self._get_synonymous_codons().keys():
            if aa not in natural_codon_usage:
                continue

            nat_counts = natural_codon_usage[aa]
            nat_total = sum(nat_counts.values())

            if nat_total == 0:
                continue

            # Find most used codon in natural
            nat_most_common = max(nat_counts.items(), key=lambda x: x[1])

            # Find most used in target
            aa_codons = self._get_synonymous_codons()[aa]
            target_counts = {c: target_codon_counts.get(c.replace('T', 'U'), 0)
                           for c in aa_codons}
            target_total = sum(target_counts.values())

            if target_total == 0:
                continue

            target_most_common = max(target_counts.items(), key=lambda x: x[1])

            comparisons.append({
                'amino_acid': aa,
                'natural_most_common': nat_most_common[0],
                'natural_frequency': nat_most_common[1] / nat_total,
                'target_most_common': target_most_common[0],
                'target_frequency': target_most_common[1] / target_total,
                'uses_human_preferred': target_most_common[0] == self.HUMAN_PREFERRED.get(aa, ''),
                'differs_from_natural': target_most_common[0] != nat_most_common[0]
            })

        return comparisons

    def generate_report(self, optimization: Dict, comparisons: List[Dict],
                       seq_id: str) -> str:
        """Generate codon optimization report."""
        report = []
        report.append("=" * 80)
        report.append("CODON OPTIMIZATION VERIFICATION (CORRECTED RSCU)")
        report.append("=" * 80)
        report.append(f"Sequence: {seq_id}")
        report.append("")

        # Overall assessment
        report.append("-" * 80)
        report.append("OPTIMIZATION ASSESSMENT")
        report.append("-" * 80)
        report.append(f"Overall: {optimization['assessment']}")
        report.append(f"Average RSCU: {optimization['avg_rscu']:.4f}")
        report.append(f"Max RSCU: {optimization['max_rscu']:.4f}")
        report.append(f"Optimized Codons: {optimization['optimized_codon_count']}")
        report.append("")

        # Add comparison to natural if available
        if 'vs_natural' in optimization:
            vs_nat = optimization['vs_natural']
            report.append("-" * 80)
            report.append("DEVIATION FROM NATURAL SARS-CoV-2")
            report.append("-" * 80)
            report.append(f"Average frequency ratio: {vs_nat['avg_ratio']:.4f}")
            report.append(f"Max frequency ratio: {vs_nat['max_ratio']:.4f}")
            report.append(f"High deviation codons (>2x): {vs_nat['high_deviation_count']}")
            report.append("")

        # Optimized codons
        if optimization['optimized_codons']:
            report.append("-" * 80)
            report.append("OPTIMIZED CODONS (RSCU > 1.2)")
            report.append("-" * 80)
            report.append("RSCU > 1.0 means codon used more than expected by chance")
            report.append("")

            for i, codon in enumerate(optimization['optimized_codons'][:15], 1):
                human_pref = " ✓" if codon.get('human_preferred') else ""
                vs_nat = codon.get('vs_natural_ratio')
                vs_nat_str = f" (vs natural: {vs_nat:.2f}x)" if vs_nat else ""

                report.append(f"{i}. {codon['codon']} ({codon['amino_acid']}){human_pref}")
                report.append(f"   RSCU: {codon['rscu']:.4f}{vs_nat_str}")
                report.append(f"   Level: {codon['level']}")
            report.append("")

        # Comparisons
        if comparisons:
            report.append("-" * 80)
            report.append("COMPARISON TO NATURAL SARS-CoV-2")
            report.append("-" * 80)
            report.append("Amino acids with codon usage changes:")
            report.append("")

            diff_count = 0
            for comp in comparisons:
                if comp['differs_from_natural']:
                    diff_count += 1
                    report.append(f"{comp['amino_acid']}:")
                    report.append(f"  Natural: {comp['natural_most_common']} ({comp['natural_frequency']:.2%})")
                    report.append(f"  Target:  {comp['target_most_common']} ({comp['target_frequency']:.2%})")
                    report.append(f"  Human preferred: {comp['uses_human_preferred']}")
                    report.append("")

            report.append(f"Total amino acids with different codon preference: {diff_count}/20")
            report.append("")

        # Interpretation
        report.append("-" * 80)
        report.append("INTERPRETATION")
        report.append("-" * 80)

        if optimization['assessment'] in ['HIGHLY_OPTIMIZED', 'OPTIMIZED']:
            report.append("⚠️  CODON OPTIMIZATION DETECTED")
            report.append("")
            report.append("Evidence:")
            report.append(f"  • Average RSCU: {optimization['avg_rscu']:.4f} (>1.0 indicates bias)")
            if 'vs_natural' in optimization:
                vs_nat = optimization['vs_natural']
                if vs_nat['avg_ratio'] > 1.5 or vs_nat['avg_ratio'] < 0.67:
                    report.append(f"  • {vs_nat['high_deviation_count']} codons deviate >2x from natural")
            if diff_count > 10:
                report.append(f"  • {diff_count} amino acids use different codons than natural virus")
            report.append("")
            report.append("RECOMMENDATION: Sequence appears optimized for human expression")
            report.append("  • May increase translation efficiency")
            report.append("  • May alter immune response")
            report.append("  • Indicates laboratory engineering")
        elif optimization['assessment'] == 'PARTIALLY_OPTIMIZED':
            report.append("⚠️  PARTIAL CODON OPTIMIZATION")
            report.append("")
            report.append("Some codon optimization detected but not extensive.")
            report.append("RECOMMENDATION: Monitor expression levels")
        else:
            report.append("✓ NO SIGNIFICANT OPTIMIZATION")
            report.append("")
            report.append("Codon usage appears similar to natural sequences.")
            report.append("RECOMMENDATION: Standard expression expected")

        report.append("")
        report.append("=" * 80)

        return "\n".join(report)


def main():
    """Main workflow."""
    import argparse

    parser = argparse.ArgumentParser(
        description='Verify codon optimization patterns'
    )
    parser.add_argument('target', help='Target sequence FASTA')
    parser.add_argument('--natural', nargs='+', help='Natural SARS-CoV-2 FASTA files')
    parser.add_argument('--output', help='Output report file')

    args = parser.parse_args()

    verifier = CodonOptimizationVerifier()

    # Load target
    logger.info(f"Loading target from {args.target}")
    target_records = list(SeqIO.parse(args.target, "fasta"))
    target_seq = str(target_records[0].seq)
    target_id = target_records[0].id

    # Load natural references
    natural_seqs = {}
    if args.natural:
        for nat_file in args.natural:
            for record in SeqIO.parse(nat_file, "fasta"):
                natural_seqs[record.id] = str(record.seq)
    else:
        # Use default
        data_dir = Path('data/sequences')
        if data_dir.exists():
            for fasta_file in data_dir.glob('*.fasta'):
                try:
                    for record in SeqIO.parse(fasta_file, "fasta"):
                        if 'Wuhan' in record.id or 'MT020880' in record.id:
                            natural_seqs[record.id] = str(record.seq)
                except:
                    pass

    if not natural_seqs:
        logger.error("No natural reference sequences found")
        return 1

    logger.info(f"Loaded {len(natural_seqs)} natural reference sequences")

    # Detect optimization
    optimization = verifier.detect_optimization(target_seq, natural_seqs)

    # Compare to natural
    comparisons = verifier.compare_to_natural(target_seq, natural_seqs)

    # Generate report
    report = verifier.generate_report(optimization, comparisons, target_id)
    print(report)

    if args.output:
        with open(args.output, 'w') as f:
            f.write(report)
        logger.info(f"Report saved to {args.output}")

    return 0


if __name__ == '__main__':
    sys.exit(main())
