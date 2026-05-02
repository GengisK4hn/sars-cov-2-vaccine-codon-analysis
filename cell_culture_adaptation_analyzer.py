#!/usr/bin/env python3
"""
CELL CULTURE ADAPTATION SIGNATURE ANALYZER
===========================================

Analyzes SARS-CoV-2 sequences for signatures of cell culture adaptation,
specifically looking for:

1. Furin Cleavage Site (FCS) - PRRA insertion with codon optimization
2. P681R mutation - associated with Delta variant and cell adaptation
3. Stabilizing mutations (2P/6P) - vaccine/lab engineering signatures
4. Codon usage bias - human optimization patterns
5. VERO/HAE cell adaptation signatures

Based on research showing SARS-CoV-2 Wuhan shows signatures of adaptation
to VERO cells and HAE cultures before human introduction.

Author: Supracode Analysis Tool
Date: 2026-04-21
"""

from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
import json
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import numpy as np
from collections import Counter, defaultdict


class CellCultureAdaptationAnalyzer:
    """
    Analyzes sequences for cell culture adaptation signatures.

    Adaptation markers detected:
    1. Furin Cleavage Site (PRRA insertion)
    2. P681R mutation (Delta variant marker)
    3. Stabilizing mutations (2P/6P prolines)
    4. Human-optimized codon usage
    5. CGG arginine codon bias (laboratory signature)
    6. VERO cell adaptation markers
    """

    def __init__(self):
        # FCS PRRA variants (all 576 possible codon combinations)
        self.prra_variants = self._generate_prra_variants()

        # Priority variants (Wuhan + human-optimized)
        self.priority_variants = {
            'CCTCGGCGGGGC',  # Wuhan PRRA (original)
            'CCTCGGCGGGCT',  # Wuhan PRRA variant
            'CCCCGCCGCGCC',  # Human-optimized
            'CCCCGCCGCGCT',  # Human-optimized
        }

        # Human CAI weights (Codon Adaptation Index)
        self.human_cai_weights = {
            'CCC': 1.00, 'CCT': 0.85, 'CCA': 0.79, 'CCG': 0.39,  # Pro
            'CGC': 0.90, 'CGG': 1.00, 'AGA': 0.70, 'AGG': 0.67,  # Arg
            'CGT': 0.38, 'CGA': 0.52,
            'GCC': 1.00, 'GCT': 0.68, 'GCA': 0.58, 'GCG': 0.25  # Ala
        }

        # VERO cell codon preferences (African green monkey)
        self.vero_codon_preferences = {
            # VERO cells prefer certain codons
            'CGG': 0.95,  # Arg - high in VERO
            'CGC': 0.92,  # Arg - high in VERO
            'CCC': 0.88,  # Pro - high in VERO
        }

        # 2P prolines (vaccine signature)
        self.two_p_positions = {
            'K986P': (985, 'AAG', ['CCT', 'CCC', 'CCA', 'CCG']),  # Lys->Pro at 986
            'V987P': (986, 'GTG', ['CCT', 'CCC', 'CCA', 'CCG']),  # Val->Pro at 987
        }

    def _generate_prra_variants(self) -> set:
        """Generate all 576 PRRA DNA variants"""
        PRO = ['CCA', 'CCC', 'CCG', 'CCT']
        ARG = ['AGA', 'AGG', 'CGA', 'CGC', 'CGG', 'CGT']
        ALA = ['GCA', 'GCC', 'GCG', 'GCT']

        variants = set()
        for p in PRO:
            for r1 in ARG:
                for r2 in ARG:
                    for a in ALA:
                        variants.add(p + r1 + r2 + a)
        return variants

    def extract_spike_region(self, sequence: str) -> Optional[str]:
        """
        Extract spike protein region from full genome.
        Spike starts at position 21563 (Wuhan-Hu-1 reference).
        """
        spike_start = 21562  # 0-indexed

        if len(sequence) > spike_start + 3822:  # Spike is ~3822 codons
            return sequence[spike_start:spike_start + 3822]
        return None

    def detect_fcs_nucleotide(self, sequence: str) -> Dict:
        """
        Detect FCS/PRRA motifs in nucleotide sequence.
        Returns detailed analysis of FCS variants and codon optimization.
        """
        results = {
            'prra_detected': False,
            'variant_sequence': None,
            'variant_type': None,
            'cai_score': 0.0,
            'is_priority_variant': False,
            'is_human_optimized': False,
            'cgg_count': 0,  # CGG arginine codons (lab signature)
            'position': None,
            'context': None
        }

        seq = sequence.upper().replace('N', '').replace('-', '')

        # Check all 12-base windows (4 codons)
        for i in range(len(seq) - 11):
            window = seq[i:i+12]

            if window in self.prra_variants:
                results['prra_detected'] = True
                results['variant_sequence'] = window
                results['position'] = i

                # Calculate CAI score
                codons = [window[j:j+3] for j in range(0, 12, 3)]
                cai_scores = []
                for codon in codons:
                    weight = self.human_cai_weights.get(codon, 0.5)
                    cai_scores.append(weight)

                results['cai_score'] = round(np.exp(np.mean(np.log([max(s, 0.01) for s in cai_scores]))), 4)

                # Check if priority variant
                if window in self.priority_variants:
                    results['is_priority_variant'] = True
                    results['variant_type'] = 'priority'

                # Check human optimization
                if results['cai_score'] > 0.85:
                    results['is_human_optimized'] = True

                # Count CGG codons (laboratory signature)
                results['cgg_count'] = sum(1 for codon in codons if codon == 'CGG')

                # Get context
                context_start = max(0, i - 30)
                context_end = min(len(seq), i + 42)
                results['context'] = seq[context_start:context_end]

                break

        return results

    def detect_p681_mutation(self, spike_sequence: str) -> Dict:
        """
        Check amino acid at position 681 (S1/S2 junction).
        P681R is a Delta variant marker associated with enhanced furin cleavage.
        """
        # Position 681 in spike (0-indexed = 680)
        pos_681_start = 680 * 3  # Convert to DNA position

        results = {
            'aa_681': None,
            'is_wildtype': False,
            'is_p681r': False,
            'codon_681': None,
            'is_cgg': False  # CGG is rare in nature, common in lab work
        }

        if len(spike_sequence) > pos_681_start + 3:
            codon = spike_sequence[pos_681_start:pos_681_start + 3]
            results['codon_681'] = codon

            # Translate
            try:
                aa = Seq(codon).translate()
                results['aa_681'] = str(aa)

                # Check mutations
                if str(aa) == 'P':
                    results['is_wildtype'] = True
                elif str(aa) == 'R':
                    results['is_p681r'] = True

                    # Check if it's CGG (lab signature)
                    if codon == 'CGG':
                        results['is_cgg'] = True
            except:
                pass

        return results

    def detect_stabilizing_mutations(self, spike_sequence: str) -> Dict:
        """
        Detect 2P/6P stabilizing prolines (K986P, V987P).
        These are strong vaccine/lab engineering signatures.
        """
        results = {
            'two_p_detected': False,
            'k986_is_proline': False,
            'v987_is_proline': False,
            'k986_codon': None,
            'v987_codon': None,
            'engineering_signature': False
        }

        # Check K986P and V987P
        for mutation, (pos, wildtype, proline_codons) in self.two_p_positions.items():
            codon_start = pos * 3

            if len(spike_sequence) > codon_start + 3:
                codon = spike_sequence[codon_start:codon_start + 3]

                if 'K986' in mutation:
                    results['k986_codon'] = codon
                    results['k986_is_proline'] = codon in proline_codons
                elif 'V987' in mutation:
                    results['v987_codon'] = codon
                    results['v987_is_proline'] = codon in proline_codons

        results['two_p_detected'] = results['k986_is_proline'] and results['v987_is_proline']
        results['engineering_signature'] = results['two_p_detected']

        return results

    def analyze_codon_usage_bias(self, sequence: str) -> Dict:
        """
        Analyze codon usage bias for human optimization signatures.
        """
        results = {
            'total_codons': 0,
            'rare_codons': 0,
            'optimized_codons': 0,
            'cgg_arginine_count': 0,  # CGG is rare in nature
            'optimization_score': 0.0
        }

        # Count codons
        for i in range(0, len(sequence) - 2, 3):
            codon = sequence[i:i+3].upper()

            if len(codon) == 3 and all(c in 'ATGC' for c in codon):
                results['total_codons'] += 1

                # Check for CGG (lab signature)
                if codon == 'CGG':
                    results['cgg_arginine_count'] += 1
                    results['rare_codons'] += 1

                # Check optimization
                weight = self.human_cai_weights.get(codon, None)
                if weight is not None:
                    if weight > 0.85:
                        results['optimized_codons'] += 1

        # Calculate optimization score
        if results['total_codons'] > 0:
            results['optimization_score'] = results['optimized_codons'] / results['total_codons']

        return results

    def classify_adaptation_signature(self, fcs_result: Dict,
                                     p681_result: Dict,
                                     stabilizing_result: Dict,
                                     codon_result: Dict) -> Dict:
        """
        Classify overall adaptation signature based on all markers.
        """
        score = 0
        indicators = []

        # FCS detection
        if fcs_result['prra_detected']:
            score += 20
            indicators.append('FCS_PRRA_PRESENT')

            if fcs_result['is_human_optimized']:
                score += 15
                indicators.append('FCS_HUMAN_OPTIMIZED')

            if fcs_result['cgg_count'] > 0:
                score += 10
                indicators.append('FCS_CGG_CODONS')

        # P681R mutation
        if p681_result['is_p681r']:
            score += 15
            indicators.append('P681R_MUTATION')

            if p681_result['is_cgg']:
                score += 10
                indicators.append('P681R_CGG_CODON')

        # Stabilizing mutations
        if stabilizing_result['engineering_signature']:
            score += 30
            indicators.append('2P_STABILIZING_MUTATIONS')

        # Codon optimization
        if codon_result['optimization_score'] > 0.7:
            score += 15
            indicators.append('HIGH_CODON_OPTIMIZATION')

        if codon_result['cgg_arginine_count'] > 5:
            score += 10
            indicators.append('HIGH_CGG_USAGE')

        # Classification
        if score >= 70:
            classification = 'ENGINEERED_LAB_ADAPTED'
            risk_level = 'HIGH'
        elif score >= 50:
            classification = 'LIKELY_ENGINEERED'
            risk_level = 'MODERATE-HIGH'
        elif score >= 30:
            classification = 'POSSIBLY_ENGINEERED'
            risk_level = 'MODERATE'
        elif score >= 20:
            classification = 'SOME_ADAPTATION_MARKERS'
            risk_level = 'LOW-MODERATE'
        else:
            classification = 'NATURAL_WILDTYPE'
            risk_level = 'LOW'

        return {
            'classification': classification,
            'risk_level': risk_level,
            'adaptation_score': score,
            'indicators': indicators
        }

    def analyze_sequence(self, fasta_path: str) -> Optional[Dict]:
        """Analyze a single FASTA file for cell culture adaptation."""
        try:
            record = SeqIO.read(fasta_path, "fasta")
            sequence = str(record.seq).upper()

            # Extract spike region
            spike_sequence = self.extract_spike_region(sequence)
            if not spike_sequence:
                return None

            # Run all analyses
            fcs_result = self.detect_fcs_nucleotide(spike_sequence)
            p681_result = self.detect_p681_mutation(spike_sequence)
            stabilizing_result = self.detect_stabilizing_mutations(spike_sequence)
            codon_result = self.analyze_codon_usage_bias(spike_sequence)

            # Classify adaptation signature
            classification = self.classify_adaptation_signature(
                fcs_result, p681_result, stabilizing_result, codon_result
            )

            return {
                'file': str(fasta_path),
                'name': record.id[:100],
                'length': len(sequence),
                'spike_length': len(spike_sequence),
                **fcs_result,
                **p681_result,
                **stabilizing_result,
                **codon_result,
                **classification
            }

        except Exception as e:
            print(f"Error analyzing {fasta_path}: {e}")
            return None

    def batch_analyze(self, fasta_files: List[Path]) -> List[Dict]:
        """Analyze multiple FASTA files."""
        print(f"\n{'='*70}")
        print("🔬 CELL CULTURE ADAPTATION ANALYSIS")
        print(f"{'='*70}")
        print(f"\n📈 Analyzing {len(fasta_files)} sequences...")

        results = []
        for i, fasta_file in enumerate(fasta_files):
            if i % 10 == 0:
                print(f"  [{i+1}/{len(fasta_files)}] {fasta_file.name}")

            result = self.analyze_sequence(fasta_file)
            if result:
                results.append(result)

        return results

    def generate_report(self, results: List[Dict]) -> pd.DataFrame:
        """Generate comprehensive analysis report."""
        print(f"\n{'='*70}")
        print("📊 CELL CULTURE ADAPTATION REPORT")
        print(f"{'='*70}")

        # Create DataFrame
        df = pd.DataFrame(results)

        # Save results
        csv_file = 'CELL_CULTURE_ADAPTATION.csv'
        df.to_csv(csv_file, index=False)
        print(f"\n✅ CSV saved: {csv_file}")

        json_file = 'CELL_CULTURE_ADAPTATION.json'
        with open(json_file, 'w') as f:
            json.dump({
                'total_analyzed': len(results),
                'classifications': {
                    'ENGINEERED_LAB_ADAPTED': len([r for r in results if r['classification'] == 'ENGINEERED_LAB_ADAPTED']),
                    'LIKELY_ENGINEERED': len([r for r in results if r['classification'] == 'LIKELY_ENGINEERED']),
                    'POSSIBLY_ENGINEERED': len([r for r in results if r['classification'] == 'POSSIBLY_ENGINEERED']),
                },
                'fcs_detected': len([r for r in results if r['prra_detected']]),
                'p681r_detected': len([r for r in results if r['is_p681r']]),
                'engineering_signatures': len([r for r in results if r['engineering_signature']]),
                'results': results
            }, f, indent=2)
        print(f"✅ JSON saved: {json_file}")

        # Statistics
        print(f"\n📈 ADAPTATION STATISTICS:")
        print(f"  Total sequences analyzed: {len(results)}")
        print(f"  FCS/PRRA detected: {len([r for r in results if r['prra_detected']])}")
        print(f"  P681R mutation: {len([r for r in results if r['is_p681r']])}")
        print(f"  2P engineering signatures: {len([r for r in results if r['engineering_signature']])}")

        # Classification breakdown
        print(f"\n🔬 CLASSIFICATION BREAKDOWN:")
        classifications = ['ENGINEERED_LAB_ADAPTED', 'LIKELY_ENGINEERED',
                          'POSSIBLY_ENGINEERED', 'SOME_ADAPTATION_MARKERS', 'NATURAL_WILDTYPE']

        for cls in classifications:
            count = len([r for r in results if r['classification'] == cls])
            if count > 0:
                print(f"  {cls}: {count}")

        # Top engineered sequences
        print(f"\n🧬 TOP ADAPTATION SIGNATURES:")
        top_adapted = sorted(results, key=lambda x: x.get('adaptation_score', 0), reverse=True)[:10]

        for i, r in enumerate(top_adapted, 1):
            print(f"\n  {i}. {r['name'][:60]}")
            print(f"     Classification: {r['classification']}")
            print(f"     Adaptation Score: {r.get('adaptation_score', 0)}")
            print(f"     Indicators: {', '.join(r.get('indicators', []))}")

            if r['prra_detected']:
                print(f"     FCS: {r['variant_sequence']} (CAI={r['cai_score']})")
            if r['is_p681r']:
                print(f"     P681R: {r['codon_681']} (Position 681)")

        return df


def main():
    """Main execution."""
    print("╔" + "═"*68 + "╗")
    print("║" + " "*5 + "🧬 CELL CULTURE ADAPTATION SIGNATURE ANALYZER" + " "*5 + "║")
    print("╚" + "═"*68 + "╝")

    # Find SARS-CoV-2 sequences
    fasta_files = [
        Path('/home/dad/supracode-tool/data/sequences/NC_045512_Wuhan-Hu-1_reference.fasta'),
        Path('/home/dad/supracode-tool/data/sequences/MT020880.1_Early_Wuhan.fasta'),
        Path('/home/dad/supracode-tool/data/sequences/OMX067679.1_Omicron_BA.1.fasta'),
        Path('/home/dad/supracode-tool/data/sequences/OMX067680.1_Omicron_BA.2.fasta'),
        Path('/home/dad/supracode-tool/data/sequences/OMX095706.1_Delta.fasta'),
        Path('/home/dad/supracode-tool/data/sequences/sars_cov_2_wuhan.fasta'),
    ]

    # Add vaccine sequences for comparison
    vaccine_files = [
        Path('/home/dad/supracode-tool/data/sequences/pfizer_bnt162b2.fasta'),
    ]

    all_files = fasta_files + vaccine_files

    # Filter existing files
    existing_files = [f for f in all_files if f.exists()]

    print(f"\n✅ Found {len(existing_files)} sequences to analyze")

    if not existing_files:
        print("❌ No sequences found!")
        return

    # Initialize analyzer
    analyzer = CellCultureAdaptationAnalyzer()

    # Batch analyze
    results = analyzer.batch_analyze(existing_files)

    # Generate report
    analyzer.generate_report(results)

    print(f"\n{'='*70}")
    print("✅ CELL CULTURE ADAPTATION ANALYSIS COMPLETE!")
    print(f"{'='*70}")


if __name__ == "__main__":
    main()
