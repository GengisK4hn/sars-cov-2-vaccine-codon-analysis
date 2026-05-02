# Quick Start Guide - Reproduce the Findings

This guide will walk you through reproducing the key findings in **under 10 minutes**.

## Prerequisites

```bash
# Install Python dependencies
pip install -r requirements.txt
```

## Step 1: Download Reference Sequences (2 minutes)

```bash
# Download all 7 required sequences from NCBI
bash download_all_sequences.sh
```

**Expected output:**
```
========================================
SARS-CoV-2 Sequence Download Script
========================================

--- Vaccine References ---
Downloading Pfizer (OR134577.1)...
  ✓ Downloaded Pfizer (4.3K)
Downloading Moderna (OR134578.1)...
  ✓ Downloaded Moderna (3.8K)

--- Natural Variants ---
Downloading Wuhan Reference (NC_045512.2)...
  ✓ Downloaded Wuhan Reference (9.8K)
...

✅ All sequences downloaded successfully!
```

## Step 2: Verify the 44nt Consensus Sequence (30 seconds)

```bash
# Count in Moderna RNAseq data
grep -c "AAGATCGCCGACTACAACTACAAGCTGCCCGACGACTTCACCGG" \
  data/sequences/RNAseq-Mod2_R2_001.fastq.fasta
```

**Expected output:**
```
156086
```

```bash
# Verify ABSENCE in natural variants
grep -c "AAGATCGCCGACTACAACTACAAGCTGCCCGACGACTTCACCGG" \
  data/sequences/NC_045512_Wuhan-Hu-1_reference.fasta
```

**Expected output:**
```
0
```

## Step 3: Verify the 19nt FCS Reverse Complement (30 seconds)

```bash
# Count reverse complement in Moderna RNAseq
grep -c "CTACGTGCCCGCCGAGGAG" \
  data/sequences/RNAseq-Mod2_R2_001.fastq.fasta
```

**Expected output:**
```
548
```

```bash
# Verify ABSENCE in Wuhan
grep -c "CTACGTGCCCGCCGAGGAG" \
  data/sequences/NC_045512_Wuhan-Hu-1_reference.fasta
```

**Expected output:**
```
0
```

## Step 4: Run Codon Optimization Analysis (5 minutes)

```bash
# Analyze Pfizer
python codon_optimization_verifier.py \
  data/vaccine_references/pfizer_bnt162b2_vector_OR134577.fasta \
  --natural data/sequences/MT020880.1_Early_Wuhan.fasta data/sequences/SARS-CoV-2.fasta \
  --output CODON_OPTIMIZATION_PFIZER.txt
```

**Expected output (last few lines):**
```
Overall: HIGHLY_OPTIMIZED
Average RSCU: 1.4815

Amino acids with codon usage changes:
R: Natural AGA (54.11%) → Target CGA (27.68%)
T: Natural ACA (56.26%) → Target ACC (43.94%)
A: Natural GCA (57.26%) → Target GCC (38.58%)

Total amino acids with different codon preference: 3/20
```

```bash
# Analyze Moderna
python codon_optimization_verifier.py \
  data/vaccine_references/moderna_mrna1273_vector_OR134578.fasta \
  --natural data/sequences/MT020880.1_Early_Wuhan.fasta data/sequences/SARS-CoV-2.fasta \
  --output CODON_OPTIMIZATION_MODERNA.txt
```

**Expected output (last line):**
```
Total amino acids with different codon preference: 7/20
```

```bash
# Analyze Delta variant
python codon_optimization_verifier.py \
  data/sequences/OMX095706.1_Delta.fasta \
  --natural data/sequences/MT020880.1_Early_Wuhan.fasta \
  --output CODON_OPTIMIZATION_DELTA.txt
```

**Expected output (last line):**
```
Total amino acids with different codon preference: 0/20
```

## Step 5: Verify All Findings at Once (1 minute)

```bash
# Run the comprehensive verification script
python independent_verification.py --data-dir data
```

**Expected output:**
```
============================================================
INDEPENDENT VERIFICATION REPORT
============================================================

[1/6] Verifying RSCU + Amino Acid Preference Changes...
[2/6] Verifying 44nt Consensus Sequence...
[3/6] Verifying 19nt FCS Sequence...
[4/6] Verifying VERO/HAE Cell Culture Adaptation...
[5/6] Verifying NLS Detection...
[6/6] Verifying GOF Signatures...

============================================================
STATUS: ✅ ALL FINDINGS VERIFIED
============================================================
```

## Summary of Verified Findings

| Finding | Pfizer | Moderna | All Variants | Status |
|---------|--------|--------|--------------|--------|
| **AA Changes** | 3/20 (15%) | 7/20 (35%) | 0/20 (0%) | ✅ |
| **44nt Sequence** | 3 reads | 156,086 reads | 0 reads | ✅ |
| **19nt FCS Revcomp** | 0 reads | 548 reads | 0 reads | ✅ |

## Troubleshooting

**Problem:** `python: command not found`
```bash
# Use python3 instead
python3 codon_optimization_verifier.py ...
```

**Problem:** `ModuleNotFoundError: No module named 'Bio'`
```bash
# Install BioPython
pip install biopython
```

**Problem:** Download script fails
```bash
# Check internet connection
# Or download sequences manually from NCBI:
# https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2
```

## Expected Time to Reproduce

- **Fast verification (grep commands):** 1 minute
- **Full codon analysis:** 10 minutes
- **Complete verification:** 15 minutes

## Need Help?

- Open an issue on GitHub
- Check `VERIFIED_FINDINGS_SUMMARY.md` for detailed methodology
- Review analysis scripts for implementation details

---

**Total time to reproduce all findings: ~15 minutes**
**Difficulty level: Easy (basic command line skills required)**
