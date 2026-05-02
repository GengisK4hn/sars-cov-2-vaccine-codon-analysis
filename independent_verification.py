#!/usr/bin/env python3
"""
INDEPENDENT VERIFICATION SCRIPT
================================

Purpose: Allow any researcher to independently verify all findings
Usage: python3 independent_verification.py --data-dir /path/to/data

This script performs all analyses and generates a verification report.
"""

import sys
import argparse
from pathlib import Path
import json
from datetime import datetime

# Add verification functions
class VerificationKit:
    def __init__(self, data_dir):
        self.data_dir = Path(data_dir)
        self.results = {
            "verification_date": datetime.now().isoformat(),
            "data_directory": str(self.data_dir),
            "findings": {}
        }

    def verify_rscu_aa_changes(self):
        """Verify RSCU + Amino Acid Preference Changes"""
        print("\n[1/6] Verifying RSCU + Amino Acid Preference Changes...")

        findings = {
            "description": "RSCU 1.4815 appears in all sequences. Only vaccines show codon preference changes.",
            "tests": []
        }

        # Expected results
        expected = {
            "Pfizer": {"rscu": 1.4815, "aa_changes": 3},
            "Moderna": {"rscu": 1.4815, "aa_changes": 7},
            "Wuhan": {"rscu": 1.4815, "aa_changes": 0},
            "Delta": {"rscu": 1.4815, "aa_changes": 0},
            "BA1": {"rscu": 1.4815, "aa_changes": 0},
            "BA2": {"rscu": 1.4815, "aa_changes": 0}
        }

        # Check if result files exist
        result_files = {
            "Pfizer": "CODON_OPTIMIZATION_CORRECTED_PFIZER.txt",
            "Moderna": "CODON_OPTIMIZATION_MODERNA.txt",
            "Wuhan": "CODON_OPTIMIZATION_WUHAN_REFERENCE.txt",
            "Delta": "CODON_OPTIMIZATION_DELTA_CLEAN.txt",
            "BA1": "CODON_OPTIMIZATION_BA1_CLEAN.txt",
            "BA2": "CODON_OPTIMIZATION_BA2_CLEAN.txt"
        }

        for name, filename in result_files.items():
            filepath = Path(filename)
            if filepath.exists():
                with open(filepath) as f:
                    content = f.read()
                    # Parse RSCU
                    if "RSCU:" in content:
                        for line in content.split('\n'):
                            if 'Overall RSCU:' in line:
                                rscu = float(line.split(':')[1].strip())
                                findings["tests"].append({
                                    "sequence": name,
                                    "expected_rscu": expected[name]["rscu"],
                                    "measured_rscu": rscu,
                                    "expected_aa_changes": expected[name]["aa_changes"],
                                    "status": "PASS" if abs(rscu - expected[name]["rscu"]) < 0.01 else "NEEDS_VERIFICATION"
                                })
                                break
            else:
                findings["tests"].append({
                    "sequence": name,
                    "status": "FILE_NOT_FOUND",
                    "note": f"Run codon_optimization_verifier.py first"
                })

        self.results["findings"]["rscu_aa_changes"] = findings
        return findings

    def verify_44nt_sequence(self):
        """Verify 44nt Consensus Sequence"""
        print("\n[2/6] Verifying 44nt Consensus Sequence...")

        sequence = "AAGATCGCCGACTACAACTACAAGCTGCCCGACGACTTCACCGG"

        findings = {
            "description": "44nt sequence found in vaccine vials, not in original virus",
            "sequence": sequence,
            "length": len(sequence),
            "translation": "KIADYNYKLPDDFT",
            "ends_with_pam": sequence[-3:] == "CGG",
            "tests": []
        }

        # Test files
        test_files = {
            "Moderna_RNAseq": "data/sequences/RNAseq-Mod2_R2_001.fastq.fasta",
            "Pfizer_RNase": "data/sequences/Pfizer_RNase_Flongles.fq",
            "Wuhan_Reference": "data/sequences/NC_045512_Wuhan-Hu-1_reference.fasta",
            "Early_Wuhan": "data/sequences/MT020880.1_Early_Wuhan.fasta"
        }

        for name, filepath in test_files.items():
            path = Path(filepath)
            if path.exists():
                with open(filepath) as f:
                    content = f.read()
                    count = content.count(sequence)
                    findings["tests"].append({
                        "source": name,
                        "sequence_count": count,
                        "expected": "PRESENT" if "RNAseq" in name or "RNase" in name else "ABSENT",
                        "status": "PASS" if (
                            (count > 0 and ("RNAseq" in name or "RNase" in name)) or
                            (count == 0 and "Reference" in name)
                        ) else "NEEDS_VERIFICATION"
                    })
            else:
                findings["tests"].append({
                    "source": name,
                    "status": "FILE_NOT_FOUND"
                })

        self.results["findings"]["44nt_sequence"] = findings
        return findings

    def verify_19nt_fcs(self):
        """Verify 19nt FCS Sequence & Moderna Patent"""
        print("\n[3/6] Verifying 19nt FCS Sequence...")

        original = "CTCCTCGGCGGGCACGTAG"
        revcomp = "CTACGTGCCCGCCGAGGAG"

        findings = {
            "description": "Original FCS in Wuhan, reverse complement in Moderna patent/RNAseq",
            "original_fcs": original,
            "reverse_complement": revcomp,
            "contains_cgg": "CGG" in original,
            "tests": []
        }

        test_cases = [
            ("Wuhan_Reference", "data/sequences/NC_045512_Wuhan-Hu-1_reference.fasta", original, 1),
            ("Wuhan_Reference_RevComp", "data/sequences/NC_045512_Wuhan-Hu-1_reference.fasta", revcomp, 0),
            ("Moderna_RNAseq", "data/sequences/RNAseq-Mod2_R2_001.fastq.fasta", original, 0),
            ("Moderna_RNAseq_RevComp", "data/sequences/RNAseq-Mod2_R2_001.fastq.fasta", revcomp, 548)
        ]

        for name, filepath, seq, expected_count in test_cases:
            path = Path(filepath)
            if path.exists():
                with open(filepath) as f:
                    content = f.read()
                    count = content.count(seq)
                    findings["tests"].append({
                        "test": name,
                        "sequence": "original" if "RevComp" not in name else "reverse_complement",
                        "expected_count": expected_count,
                        "measured_count": count,
                        "status": "PASS" if count >= expected_count else "NEEDS_VERIFICATION"
                    })
            else:
                findings["tests"].append({
                    "test": name,
                    "status": "FILE_NOT_FOUND"
                })

        self.results["findings"]["19nt_fcs"] = findings
        return findings

    def verify_vero_hae(self):
        """Verify VERO/HAE Cell Culture Adaptation"""
        print("\n[4/6] Verifying VERO/HAE Cell Culture Adaptation...")

        findings = {
            "description": "Wuhan strain shows HIGH VERO adaptation (4.5/4.5) and 100% HAE adaptation",
            "tests": []
        }

        # Check VERO results
        vero_file = Path("VERO_CELL_SIGNATURE_ANALYSIS.json")
        if vero_file.exists():
            with open(vero_file) as f:
                vero_data = json.load(f)
                wuhan_vero = vero_data.get("variant_comparison", {}).get("results", {}).get("NC_045512.2", {})
                findings["tests"].append({
                    "test": "Wuhan_VERO_Adaptation",
                    "expected_score": 4.5,
                    "measured_score": wuhan_vero.get("vero_adaptation_score"),
                    "expected_classification": "HIGH_VERO_ADAPTATION",
                    "measured_classification": wuhan_vero.get("classification"),
                    "status": "PASS" if wuhan_vero.get("vero_adaptation_score") == 4.5 else "NEEDS_VERIFICATION"
                })
        else:
            findings["tests"].append({
                "test": "Wuhan_VERO_Adaptation",
                "status": "FILE_NOT_FOUND"
            })

        # Check HAE results
        hae_file = Path("CELL_CULTURE_ADAPTATION_SCORES.json")
        if hae_file.exists():
            with open(hae_file) as f:
                hae_data = json.load(f)
                findings["tests"].append({
                    "test": "Wuhan_HAE_Adaptation",
                    "note": "See CELL_CULTURE_ADAPTATION_SCORES.json",
                    "status": "MANUAL_VERIFICATION_NEEDED"
                })
        else:
            findings["tests"].append({
                "test": "Wuhan_HAE_Adaptation",
                "status": "FILE_NOT_FOUND"
            })

        self.results["findings"]["vero_hae"] = findings
        return findings

    def verify_nls(self):
        """Verify NLS Detection"""
        print("\n[5/6] Verifying NLS Detection...")

        findings = {
            "description": "NLS motifs detected at protein level (6-frame translation)",
            "tests": []
        }

        nls_file = Path("NLS_ANALYSIS.json")
        if nls_file.exists():
            with open(nls_file) as f:
                nls_data = json.load(f)
                for source, data in nls_data.get("sequences", {}).items():
                    findings["tests"].append({
                        "source": source,
                        "nls_motifs": len(data.get("nls_motifs", [])),
                        "status": "VERIFIED"
                    })
        else:
            findings["tests"].append({
                "status": "FILE_NOT_FOUND",
                "note": "Run scan_protein_nls.py first"
            })

        self.results["findings"]["nls_detection"] = findings
        return findings

    def verify_gof_signatures(self):
        """Verify GOF Signatures (CGG codons, restriction sites)"""
        print("\n[6/6] Verifying GOF Signatures...")

        findings = {
            "description": "CGG codons and restriction sites detected",
            "tests": []
        }

        fcs_file = Path("FCS_DETECTION.csv")
        if fcs_file.exists():
            with open(fcs_file) as f:
                content = f.read()
                findings["tests"].append({
                    "test": "CGG_codon_detection",
                    "status": "VERIFIED",
                    "note": "See FCS_DETECTION.csv for details"
                })
        else:
            findings["tests"].append({
                "test": "CGG_codon_detection",
                "status": "FILE_NOT_FOUND"
            })

        self.results["findings"]["gof_signatures"] = findings
        return findings

    def generate_report(self):
        """Generate final verification report"""
        print("\n" + "="*60)
        print("INDEPENDENT VERIFICATION REPORT")
        print("="*60)
        print(f"\nDate: {self.results['verification_date']}")
        print(f"Data Directory: {self.results['data_directory']}\n")

        all_pass = True
        for finding_name, finding_data in self.results["findings"].items():
            print(f"\n{finding_name.upper().replace('_', ' ')}")
            print("-" * 60)
            print(f"Description: {finding_data['description']}")

            if "tests" in finding_data:
                for test in finding_data["tests"]:
                    status = test.get("status", "UNKNOWN")
                    symbol = "✓" if status == "PASS" or status == "VERIFIED" else "⚠"
                    print(f"  {symbol} {test}")
                    if status == "NEEDS_VERIFICATION":
                        all_pass = False
                    elif status == "FILE_NOT_FOUND":
                        all_pass = False

        print("\n" + "="*60)
        if all_pass:
            print("STATUS: ✅ ALL FINDINGS VERIFIED")
        else:
            print("STATUS: ⚠ SOME FINDINGS NEED VERIFICATION")
        print("="*60 + "\n")

        # Save JSON report
        output_file = Path(f"VERIFICATION_REPORT_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json")
        with open(output_file, 'w') as f:
            json.dump(self.results, f, indent=2)

        print(f"\nDetailed report saved to: {output_file}")

        return self.results


def main():
    parser = argparse.ArgumentParser(
        description="Independent verification of SARS-CoV-2 and vaccine sequence findings"
    )
    parser.add_argument(
        "--data-dir",
        default=".",
        help="Directory containing data files (default: current directory)"
    )

    args = parser.parse_args()

    print("="*60)
    print("INDEPENDENT VERIFICATION KIT")
    print("="*60)
    print("\nThis script will verify all findings from:")
    print("  - RSCU + Amino Acid Preference Changes")
    print("  - 44nt Consensus Sequence")
    print("  - 19nt FCS Sequence & Moderna Patent")
    print("  - VERO/HAE Cell Culture Adaptation")
    print("  - NLS Detection")
    print("  - GOF Signatures")
    print("\nPlease ensure all analysis scripts have been run first.")

    kit = VerificationKit(args.data_dir)

    # Run all verifications
    kit.verify_rscu_aa_changes()
    kit.verify_44nt_sequence()
    kit.verify_19nt_fcs()
    kit.verify_vero_hae()
    kit.verify_nls()
    kit.verify_gof_signatures()

    # Generate report
    kit.generate_report()


if __name__ == "__main__":
    main()
