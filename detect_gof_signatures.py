#!/usr/bin/env python3
"""
Enhanced GOF Signature Detection
Incorporating HKU4r-HZAU-2020 findings
"""

import re
from Bio import SeqIO
from collections import defaultdict

class GOFSignatureDetector:
    def __init__(self):
        # Polybasic FCS motifs (from HKU4r findings + SARS-CoV-2)
        self.fcs_patterns = {
            "PRRAR": "SARS-CoV-2_FCS",
            "KQQR": "HKU4r_FCS_like",
            "RRAR": "FCS_variant",
            "RX[RK]R": "General_polybasic",
        }
        
        # Restriction sites for infectious clone assembly
        self.restriction_sites = {
            "BsaI": "GGTCTC",
            "BsmBI": "CGTCTC",
            "BbsI": "GAAGAC",
            "SapI": "GCTCTTC",
        }
        
        # Codon patterns associated with lab engineering
        self.suspicious_codons = {
            "CGGCGG": "Double_CGG_arginine",
            "CGG": "Rare_mammalian_codon",
        }
        
    def scan_fcs(self, seq):
        """Scan for furin cleavage site motifs"""
        results = []
        for pattern, name in self.fcs_patterns.items():
            if "X" in pattern:
                # Regex pattern
                matches = list(re.finditer(pattern, seq))
                if matches:
                    results.append({
                        'type': 'FCS_motif',
                        'pattern': name,
                        'positions': [m.start() for m in matches]
                    })
            else:
                # Literal pattern
                if pattern in seq:
                    positions = [m.start() for m in re.finditer(f'(?={pattern})', seq)]
                    results.append({
                        'type': 'FCS_motif',
                        'pattern': name,
                        'positions': positions
                    })
        return results
    
    def scan_restriction_sites(self, seq):
        """Scan for infectious clone assembly sites"""
        results = []
        for site_name, site_seq in self.restriction_sites.items():
            if site_seq in seq:
                count = seq.count(site_seq)
                results.append({
                    'type': 'restriction_site',
                    'site': site_name,
                    'count': count
                })
        return results
    
    def scan_codon_anomalies(self, seq):
        """Scan for suspicious codon patterns"""
        results = []
        for codon, name in self.suspicious_codons.items():
            if codon in seq:
                count = seq.count(codon)
                results.append({
                    'type': 'codon_anomaly',
                    'pattern': name,
                    'count': count
                })
        return results
    
    def detect_chimeric_rbd(self, seq_record, references):
        """Detect RBD swapping (chimera signature)"""
        # This would require BLAST against reference database
        # Placeholder for implementation
        pass
    
    def comprehensive_scan(self, fasta_file):
        """Run comprehensive GOF signature scan"""
        records = list(SeqIO.parse(fasta_file, "fasta"))
        
        results = {
            'file': fasta_file,
            'sequences': []
        }
        
        for record in records:
            seq = str(record.seq)
            seq_result = {
                'id': record.id,
                'length': len(seq),
                'signatures': []
            }
            
            # Scan for all signatures
            seq_result['signatures'].extend(self.scan_fcs(seq))
            seq_result['signatures'].extend(self.scan_restriction_sites(seq))
            seq_result['signatures'].extend(self.scan_codon_anomalies(seq))
            
            results['sequences'].append(seq_result)
        
        return results
    
    def generate_report(self, results):
        """Generate human-readable report"""
        report = []
        report.append("=" * 80)
        report.append("GOF SIGNATURE DETECTION REPORT")
        report.append("Incorporating HKU4r-HZAU-2020 findings")
        report.append("=" * 80)
        
        for seq in results['sequences']:
            report.append(f"\nSequence: {seq['id']}")
            report.append(f"Length: {seq['length']} bp")
            
            if not seq['signatures']:
                report.append("✓ No GOF signatures detected")
            else:
                report.append(f"⚠ {len(seq['signatures'])} signatures detected:")
                
                for sig in seq['signatures']:
                    if sig['type'] == 'FCS_motif':
                        report.append(f"  • FCS Motif: {sig['pattern']} at positions {sig['positions']}")
                    elif sig['type'] == 'restriction_site':
                        report.append(f"  • Restriction Site: {sig['site']} ({sig['count']} occurrences)")
                    elif sig['type'] == 'codon_anomaly':
                        report.append(f"  • Codon Anomaly: {sig['pattern']} ({sig['count']} occurrences)")
        
        return "\n".join(report)

def main():
    import sys
    
    if len(sys.argv) < 2:
        print("Usage: detect_gof_signatures.py <fasta_file>")
        sys.exit(1)
    
    detector = GOFSignatureDetector()
    results = detector.comprehensive_scan(sys.argv[1])
    
    report = detector.generate_report(results)
    print(report)

if __name__ == "__main__":
    main()
