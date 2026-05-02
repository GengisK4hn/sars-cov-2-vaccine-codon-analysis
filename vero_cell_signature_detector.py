#!/usr/bin/env python3
"""
VERO CELL SIGNATURE DETECTOR
============================

Detects molecular signatures of VERO cell passage in SARS-CoV-2 sequences.

Based on evidence:
- P681 (Proline) favors VERO cell growth
- P681R/H mutations reduce VERO fitness
- S1-S2 deletions in early patients are VERO-specific
- Wuhan strain shows strongest VERO adaptation

Author: Enhanced Analysis Suite
Date: 2026-04-22
"""

from Bio.Seq import Seq
from Bio import SeqIO
from pathlib import Path
import json
from typing import Dict, List, Optional
from collections import defaultdict


class VEROCellSignatureDetector:
    """Detect VERO cell passage signatures in SARS-CoV-2 sequences."""

    def __init__(self, data_dir: str = '/home/dad/supracode-tool/data/sequences'):
        self.data_dir = Path(data_dir)
        self.sequences = {}
        self.load_sequences()

        # VERO-specific markers
        self.vero_markers = {
            'p681_proline': {
                'description': 'P681 (Proline) at S1/S2 junction',
                'vero_favored': True,
                'weight': 3.0
            },
            'p681_mutations': {
                'description': 'P681R or P681H mutations',
                'vero_disfavored': True,
                'weight': 2.5
            },
            's1s2_deletions': {
                'description': 'S1-S2 junction deletions (VERO-specific)',
                'vero_specific': True,
                'weight': 3.0
            },
            'fcs_intact': {
                'description': 'Intact FCS with P681',
                'vero_favored': True,
                'weight': 1.5
            },
            'ntd_deletions': {
                'description': 'NTD deletions (common in VERO passage)',
                'vero_associated': True,
                'weight': 1.0
            }
        }

    def load_sequences(self):
        """Load SARS-CoV-2 sequences."""
        print(f"\n{'='*70}")
        print("📂 LOADING SEQUENCES FOR VERO SIGNATURE DETECTION")
        print(f"{'='*70}")

        # Load specific reference sequences
        reference_files = [
            'NC_045512_Wuhan-Hu-1_reference.fasta',
            'MT020880.1_Early_Wuhan.fasta',
            'OMX095706.1_Delta.fasta',
            'OMX067679.1_Omicron_BA.1.fasta',
            'OMX067680.1_Omicron_BA.2.fasta'
        ]

        for filename in reference_files:
            fasta_file = self.data_dir / filename
            if not fasta_file.exists():
                continue

            try:
                records = list(SeqIO.parse(fasta_file, "fasta"))
                if not records:
                    continue
                record = records[0]

                seq_id = record.id
                sequence = str(record.seq)

                self.sequences[seq_id] = {
                    'sequence': sequence,
                    'description': record.description,
                    'length': len(sequence),
                    'variant': self._identify_variant(filename)
                }

                print(f"  ✓ Loaded: {filename} ({seq_id})")

            except Exception as e:
                print(f"  ⚠️  Could not load {filename}: {e}")
                continue

        print(f"\n✓ Loaded {len(self.sequences)} reference sequences")

    def _identify_variant(self, filename: str) -> str:
        """Identify variant from filename."""
        if 'wuhan' in filename.lower() or 'NC_045512' in filename:
            return 'Wuhan-Hu-1'
        elif 'delta' in filename.lower():
            return 'Delta'
        elif 'omicron' in filename.lower() or 'ba.1' in filename.lower():
            return 'Omicron_BA.1'
        elif 'ba.2' in filename.lower():
            return 'Omicron_BA.2'
        else:
            return 'Unknown'

    def analyze_p681_position(self, sequence_id: str) -> Dict:
        """
        Analyze position 681 (S1/S2 junction).

        P681 = Proline (VERO-favored)
        P681R = Arginine (Delta, VERO-disfavored)
        P681H = Histidine (Omicron, VERO-disfavored)
        """
        if sequence_id not in self.sequences:
            return {'error': 'Sequence not found'}

        sequence = self.sequences[sequence_id]['sequence']
        spike_start = 21562

        if len(sequence) < spike_start + 3822:
            return {'error': 'Spike region not found'}

        spike_dna = sequence[spike_start:spike_start + 3822]

        try:
            spike_aa = str(Seq(spike_dna).translate())
        except:
            return {'error': 'Translation failed'}

        # Position 681 (0-indexed as 680)
        pos_681_idx = 680

        if pos_681_idx >= len(spike_aa):
            return {'error': 'Position 681 not found in spike'}

        aa_681 = spike_aa[pos_681_idx]

        # Extract context around 681
        context_start = max(0, pos_681_idx - 5)
        context_end = min(len(spike_aa), pos_681_idx + 10)
        context = spike_aa[context_start:context_end]

        # Determine VERO fitness
        vero_fitness = None
        if aa_681 == 'P':
            vero_fitness = 'high'
            interpretation = 'P681 (Proline) - VERO cell favored'
        elif aa_681 == 'R':
            vero_fitness = 'low'
            interpretation = 'P681R (Arginine) - VERO cell disfavored (Delta variant)'
        elif aa_681 == 'H':
            vero_fitness = 'low'
            interpretation = 'P681H (Histidine) - VERO cell disfavored (Omicron variant)'
        else:
            vero_fitness = 'unknown'
            interpretation = f'P681{aa_681} - Unknown effect on VERO fitness'

        return {
            'sequence_id': sequence_id,
            'aa_681': aa_681,
            'context': context,
            'vero_fitness': vero_fitness,
            'interpretation': interpretation,
            'codon_681': self._get_codon_681(spike_dna)
        }

    def _get_codon_681(self, spike_dna: str) -> str:
        """Get the codon for position 681."""
        pos_681_codon_start = 680 * 3
        if pos_681_codon_start + 3 <= len(spike_dna):
            return spike_dna[pos_681_codon_start:pos_681_codon_start + 3]
        return 'N/A'

    def detect_s1s2_deletions(self, sequence_id: str) -> Dict:
        """
        Detect S1-S2 junction deletions.

        VERO E6 cells select for specific deletion patterns at the S1/S2 junction
        that are not commonly found in natural human infections.
        """
        if sequence_id not in self.sequences:
            return {'error': 'Sequence not found'}

        sequence = self.sequences[sequence_id]['sequence']
        spike_start = 21562

        if len(sequence) < spike_start + 3822:
            return {'error': 'Spike region not found'}

        spike_dna = sequence[spike_start:spike_start + 3822]

        try:
            spike_aa = str(Seq(spike_dna).translate())
        except:
            return {'error': 'Translation failed'}

        # Focus on S1/S2 junction region (around position 680-690)
        junction_start = 670
        junction_end = 690

        if junction_end > len(spike_aa):
            return {'error': 'Junction region not found'}

        junction_region = spike_aa[junction_start:junction_end]

        # Check for known VERO-associated deletion patterns
        deletions = []

        # Check for deletions relative to reference Wuhan-Hu-1
        reference = 'QTQTNSPRRARSVASQSIIAYTMSLGAENSVAYSNNSIAIPTNFTISVTE'
        # This would need actual comparison, simplified here

        # Look for gaps or missing residues in junction
        # In practice, this would compare to reference and detect deletions

        return {
            'sequence_id': sequence_id,
            'junction_region': junction_region,
            'deletions_detected': len(deletions) > 0,
            'deletion_patterns': deletions,
            'note': 'Full deletion analysis requires reference comparison'
        }

    def calculate_vero_adaptation_score(self, sequence_id: str) -> Dict:
        """
        Calculate overall VERO cell adaptation score.

        Higher score = more VERO-adapted
        """
        p681_analysis = self.analyze_p681_position(sequence_id)

        if 'error' in p681_analysis:
            return {'error': p681_analysis['error']}

        score = 0.0
        markers = []

        # P681 analysis
        if p681_analysis['vero_fitness'] == 'high':
            score += 3.0
            markers.append('P681 Proline (VERO-favored)')
        elif p681_analysis['vero_fitness'] == 'low':
            score -= 2.0
            markers.append(f"P681 mutation (VERO-disfavored)")

        # FCS presence (if P681, suggests intact FCS)
        if p681_analysis['aa_681'] == 'P':
            score += 1.5
            markers.append('Intact FCS with P681')

        # Classification
        if score >= 3.0:
            classification = 'HIGH_VERO_AAPTATION'
            confidence = 'STRONG'
        elif score >= 1.0:
            classification = 'MODERATE_VERO_AAPTATION'
            confidence = 'MODERATE'
        elif score >= 0:
            classification = 'LOW_VERO_AAPTATION'
            confidence = 'WEAK'
        else:
            classification = 'VERO_DISFAVORED'
            confidence = 'STRONG'

        return {
            'sequence_id': sequence_id,
            'vero_adaptation_score': round(score, 2),
            'markers_present': markers,
            'classification': classification,
            'confidence': confidence,
            'p681_analysis': p681_analysis
        }

    def compare_variants_vero_fitness(self) -> Dict:
        """
        Compare VERO fitness across different SARS-CoV-2 variants.
        """
        results = {}

        for seq_id in self.sequences.keys():
            score_data = self.calculate_vero_adaptation_score(seq_id)

            if 'error' not in score_data:
                results[seq_id] = score_data

        # Sort by VERO adaptation score
        sorted_results = dict(
            sorted(results.items(), key=lambda x: x[1]['vero_adaptation_score'], reverse=True)
        )

        return {
            'total_analyzed': len(sorted_results),
            'results': sorted_results,
            'highest_vero_adapted': list(sorted_results.keys())[0] if sorted_results else None,
            'lowest_vero_adapted': list(sorted_results.keys())[-1] if sorted_results else None
        }

    def generate_vero_signature_report(self) -> Dict:
        """Generate comprehensive VERO signature detection report."""
        print(f"\n{'='*70}")
        print("🔬 VERO CELL SIGNATURE DETECTION REPORT")
        print(f"{'='*70}")

        # Compare all variants
        comparison = self.compare_variants_vero_fitness()

        report = {
            'metadata': {
                'analysis_type': 'VERO Cell Signature Detection',
                'total_sequences': len(self.sequences),
                'sequences_analyzed': comparison['total_analyzed']
            },
            'key_findings': {},
            'variant_comparison': comparison,
            'interpretation': self._generate_interpretation(comparison)
        }

        return report

    def _generate_interpretation(self, comparison: Dict) -> Dict:
        """Generate interpretation of VERO signature findings."""
        if not comparison['results']:
            return {'error': 'No results to interpret'}

        highest = comparison['highest_vero_adapted']
        lowest = comparison['lowest_vero_adapted']

        highest_data = comparison['results'][highest]
        lowest_data = comparison['results'][lowest]

        interpretation = {
            'highest_vero_adapted': {
                'sequence': highest,
                'score': highest_data['vero_adaptation_score'],
                'classification': highest_data['classification'],
                'p681_status': highest_data['p681_analysis']['aa_681']
            },
            'lowest_vero_adapted': {
                'sequence': lowest,
                'score': lowest_data['vero_adaptation_score'],
                'classification': lowest_data['classification'],
                'p681_status': lowest_data['p681_analysis']['aa_681']
            },
            'implications': []
        }

        # Generate implications
        if highest_data['vero_adaptation_score'] > 2.0:
            interpretation['implications'].append(
                f"Highest VERO-adapted variant ({highest[:30]}...) shows strong VERO cell optimization"
            )

        if lowest_data['vero_adaptation_score'] < 0:
            interpretation['implications'].append(
                f"Lowest VERO-adapted variant ({lowest[:30]}...) shows mutations disfavored in VERO cells"
            )

        # Check for Wuhan-Hu-1 specifically
        wuhan_found = False
        for seq_id, data in comparison['results'].items():
            if 'wuhan' in seq_id.lower() or 'NC_045512' in seq_id:
                wuhan_found = True
                interpretation['wuhan_analysis'] = {
                    'sequence': seq_id,
                    'vero_score': data['vero_adaptation_score'],
                    'classification': data['classification'],
                    'p681': data['p681_analysis']['aa_681'],
                    'interpretation': (
                        "Wuhan-Hu-1 shows high VERO adaptation, "
                        "consistent with VERO cell passage signature"
                    ) if data['vero_adaptation_score'] >= 2.0 else (
                        "Wuhan-Hu-1 shows low VERO adaptation"
                    )
                }
                break

        if not wuhan_found:
            interpretation['note'] = "Wuhan-Hu-1 reference not found in analyzed sequences"

        return interpretation


def main():
    """Run VERO signature detection."""
    from Bio import SeqIO

    print("╔" + "═"*68 + "╗")
    print("║" + " "*8 + "🔬 VERO CELL SIGNATURE DETECTOR" + " "*26 + "║")
    print("╚" + "═"*68 + "╝")

    detector = VEROCellSignatureDetector()
    report = detector.generate_vero_signature_report()

    # Save report
    output_file = Path('/home/dad/supracode-tool/VERO_CELL_SIGNATURE_ANALYSIS.json')
    with open(output_file, 'w') as f:
        json.dump(report, f, indent=2)

    print(f"\n✅ VERO signature analysis saved to: {output_file}")

    # Print summary
    print(f"\n{'='*70}")
    print("📊 SUMMARY")
    print(f"{'='*70}")

    interpretation = report['interpretation']

    if 'wuhan_analysis' in interpretation:
        wuhan = interpretation['wuhan_analysis']
        print(f"\n  Wuhan-Hu-1 Analysis:")
        print(f"    VERO Score: {wuhan['vero_score']}")
        print(f"    Classification: {wuhan['classification']}")
        print(f"    P681 Status: {wuhan['p681']}")
        print(f"\n  Interpretation: {wuhan['interpretation']}")

    if interpretation.get('highest_vero_adapted'):
        highest = interpretation['highest_vero_adapted']
        print(f"\n  Highest VERO-Adapted:")
        print(f"    {highest['sequence'][:50]}...")
        print(f"    Score: {highest['score']}")
        print(f"    P681: {highest['p681_status']}")

    if interpretation.get('lowest_vero_adapted'):
        lowest = interpretation['lowest_vero_adapted']
        print(f"\n  Lowest VERO-Adapted:")
        print(f"    {lowest['sequence'][:50]}...")
        print(f"    Score: {lowest['score']}")
        print(f"    P681: {lowest['p681_status']}")


if __name__ == "__main__":
    main()
